"""Run parent planning, the RF/chemistry worker, and retained B/C planning.

This module is a real-execution orchestration path. Tests inject a backend and
CPU stage; the command-line path always uses ResidentGraph and rf_worker.
One process owns an output directory. Stored errors and incomplete searches are
retained rather than automatically retried.
"""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import signal
import subprocess
import time
import traceback

from .graph_core import NORMAL, POLICIES, QueryStore, ResidentGraph
from .graph_core import evaluate_policy, read, save, summary
from .cpu_costs import summarize_cpu_work, actual_work_components
from .qualification import settings, evaluator, report_qualification

ROOT = Path(__file__).resolve().parents[1]
ARMS = {"B_joint6plus6": "task_contribution",
        "C_joint6plus6": "permuted_contribution",
        "B_fixedtask6plus6": "task_contribution",
        "C_fixedtask6plus6": "permuted_contribution"}


def versions():
    return {name: importlib.metadata.version(name)
            for name in ("syntheseus", "torch", "rdkit")}


class RoleAwareGraph:
    """Preserve the distinct, existing parent and candidate stock semantics.

    Actual parent code (v2_graph2edits_probe.py) accepts a stock target and removes
    its standard InChIKey from inventory. Formal512 candidate code rejects stock
    targets before search. ResidentGraph implements both via target_role. Role is
    written into the real request record before backend execution.
    """
    def __init__(self, backend, context):
        self.backend, self.context = backend, context

    @property
    def load_seconds(self):
        return self.backend.load_seconds

    def search(self, request):
        request["target_role"] = self.context["target_role"]
        return self.backend.search(request)

    def audit(self, graph, smiles):
        return self.backend.audit(graph, smiles)


def graph_config(config, installed_versions):
    qualification, policies = settings(config)
    expected = {"parent_calls": 12, "parent_soft_seconds": 45.,
                "candidate_calls": [6, 6], "candidate_shared_soft_seconds": 45.,
                "candidate_score_slots_per_arm": 8,
                "candidate_forward_attempts_per_arm": 16,
                "policies": list(policies)}
    if any(config.get(key) != value for key, value in expected.items()):
        raise ValueError("Method configuration differs from the retained B/C procedure")
    compatible = {"stock": config["stock_keys"], "model_dir": config["graph_model_dir"],
                  "num_top_results": 20, "max_expansion_depth": 12,
                  "reaction_model_calls_limit": 12, "time_limit_s": 45,
                  "target_exclusion": "remove_target_standard_inchikey_from_inventory_per_search",
                  "cache": "reset_between_every_target_including_duplicates",
                  "versions": installed_versions}
    return {"scope": config["scope"], "package_version": config["package_version"],
            "stock": config["stock_keys"], "model_dir": config["graph_model_dir"],
            "device": config["device"], "versions": installed_versions,
            "compatible12_settings": compatible, "policies": list(policies),
            **({"qualification": qualification} if "qualification" in config else {}),
            "explicit_reuse_study": config.get("reuse_study"),
            "shared_soft_seconds": 45.,
            "time_note": "Resident search time, including reset, consumes the logical allowance; "
                         "model load, CPU stage, audit, graph save and retrieval are separate.",
            "current_study_reuse_selection": "First path for exact molecule/call cap/soft time, "
                "including negative/error/incomplete queries; adopted seed is retained.",
            "parent_target_role": "Stock parent allowed; target standard InChIKey excluded",
            "candidate_target_role": "Stock candidate rejected; target standard InChIKey excluded"}


def run_rf_candidates(output, config):
    """Use the configured RF interpreter, never the Graph environment's sklearn."""
    command = [config["rf_python"], "-m", "lfr_method.rf_worker", "candidates",
               "--output", str(output)]
    env = dict(os.environ, PYTHONPATH=str(ROOT) + os.pathsep + os.environ.get("PYTHONPATH", ""),
               OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", MKL_NUM_THREADS="2", PYTHONUTF8="1")
    with (output / "candidate_worker.log").open("a", encoding="utf-8") as stream:
        code = subprocess.call(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
                               stdout=stream, stderr=subprocess.STDOUT)
    if code:
        raise RuntimeError(f"RF candidate worker exited {code}; retained log: {output / 'candidate_worker.log'}")


def parent_row(source, query):
    """Preserve reported work while refusing unaudited or over-budget positives."""
    effective = dict(query)
    if int(query.get("model_calls", 0)) > 12:
        effective["status"] = "reaction_call_budget_limitation"
    elif query.get("status") == "solved" and not query.get("verified_solved"):
        effective["status"] = "route_audit_limitation"
    return dict(source, requested_search_seed=source["search_seed"],
                search_seed=int(query.get("adopted_seed", source["search_seed"])),
                graph_verified=effective.get("status") == "solved" and bool(query.get("verified_solved")),
                graph_result=effective, reported_parent_result=query,
                graph_result_path=query.get("result_path"), graph_audit=query.get("audit"),
                prior_graph_reused=bool(query.get("history_path")),
                query_seed_adopted=int(query.get("adopted_seed", source["search_seed"])) != source["search_seed"])


def stage_index(stage, sources):
    entries = stage["sources"]
    result = {entry["source"]["sample_id"]: entry for entry in entries}
    if len(result) != len(entries) or set(result) != {s["sample_id"] for s in sources}:
        raise ValueError("Candidate stage does not cover exactly the selected raw source IDs")
    for source in sources:
        if result[source["sample_id"]]["source"] != source:
            raise ValueError("Candidate stage source identity differs from the parent stage")
    return result


def run(output, *, backend_factory=None, candidate_runner=run_rf_candidates,
        version_provider=versions, auditor_factory=None):
    output = Path(output)
    config = read(output / "method_config.json")
    qualification, policies = settings(config)
    evaluate = evaluator(qualification)
    selected = read(output / "selected_inputs.json")
    scored = read(output / "sources_scored.json")
    if scored["supplied_inputs"] != selected["rows"] or scored["rf_model_dir"] != config["rf_model_dir"]:
        raise ValueError("RF source scoring belongs to different supplied inputs or model directory")
    raw = scored["rows"]
    if len(raw) != len(selected["rows"]) or len({s["sample_id"] for s in raw}) != len(raw):
        raise ValueError("Normalized sources lost raw slots or repeat an identifier")
    context = {"target_role": "parent"}
    graph_output = output / "graph_measurements"
    installed = version_provider()
    factory = backend_factory or (lambda: ResidentGraph(config["stock_keys"],
                                                       config["graph_model_dir"], config["device"]))
    # With no explicit reuse study, scanning our own excluded directory yields
    # an empty history. Current-directory query continuation remains enabled.
    history = Path(config["reuse_study"]) if config.get("reuse_study") else graph_output
    if config.get("reuse_study"):
        from .exact_study_store import ExactStudyStore
        store = ExactStudyStore(graph_output, graph_config(config, installed), history,
                                lambda: RoleAwareGraph(factory(), context), auditor_factory)
    else:
        store = QueryStore(graph_output, graph_config(config, installed), history,
                           lambda: RoleAwareGraph(factory(), context))
    state_path = output / "method_execution.json"
    state = {"pid": os.getpid(), "status": "running", "stage": "parents",
             "started_unix": time.time(), "raw_n": len(raw)}
    save(state_path, state)
    try:
        parents = []
        for source in raw:
            path = output / "parent_outcomes" / (source["sample_id"] + ".json")
            if path.exists():
                entry = read(path)
                if entry["source_scored"] != source:
                    raise ValueError("Saved parent belongs to a different normalized source")
                row = entry["source_with_parent"]
            else:
                if source["valid"]:
                    query = store(source["canonical_smiles"], source["search_seed"], 12, 45.)
                    row = parent_row(source, query)
                else:
                    row = dict(source, graph_verified=False, graph_result={"status": "invalid"},
                               graph_result_path=None, graph_audit=None)
                save(path, {"source_scored": source, "source_with_parent": row})
            parents.append(row)
            save(output / "method_progress.json", {"stage": "parents", "last_source": source["sample_id"],
                 "parent_raw_rows": len(parents), "new_graph_calls_this_invocation": store.new_calls_this_invocation})
        source_report = {"scope": config["scope"], "rows": parents,
                         "raw_n": len(parents), "training": 0, "sampling": 0}
        parent_path = output / "sources_with_parents.json"
        if parent_path.exists() and read(parent_path) != source_report:
            raise ValueError("Completed parent source identities differ")
        save(parent_path, source_report)
        state.update(stage="cpu_candidates")
        save(state_path, state)
        candidate_runner(output, config)
        candidate_report = read(output / "candidate_stage.json")
        candidates = stage_index(candidate_report, parents)
        context["target_role"] = "candidate"
        state.update(stage="candidate_planning")
        save(state_path, state)
        rows = []
        for source in parents:
            for policy in policies:
                path = output / "policy_outcomes" / policy / (source["sample_id"] + ".json")
                previous = read(path) if path.exists() else None
                entry = candidates[source["sample_id"]]
                arm = entry.get("arms", {}).get(ARMS[policy])
                row = evaluate(source, policy, arm, store, saved=previous,
                                      persist=lambda result, p=path: save(p, result))
                report_qualification(row, qualification)
                save(path, row)
                rows.append(row)
            save(output / "method_progress.json", {"stage": "candidate_planning",
                 "last_source": source["sample_id"], "policy_rows": len(rows),
                 "new_graph_calls_this_invocation": store.new_calls_this_invocation})
        records = [record for _, record in store.cache.values()]
        new = [r for r in records if r.get("cache_mode") == "new"]
        reused = [r for r in records if r.get("cache_mode") == "explicit_study_exact"]
        actual = {"training": 0, "sampling": 0,
                  "RF_source_predictions": scored["actual_new_source_predictions"],
                  "RF_candidate_predictions": candidate_report["actual_new_candidate_predictions"],
                  "forward_reactions": candidate_report["actual_new_forward_attempts"],
                  "Graph_queries": len(new), "Graph_model_calls": sum(r.get("model_calls", 0) for r in new)}
        actual["RF_predictions"] = actual["RF_source_predictions"] + actual["RF_candidate_predictions"]
        # Recompute from cumulative completed records, including old candidate
        # reports without the new summary. Do not clear earlier spend on resume
        # or assume missing legacy root/attribution counters are measured zero.
        cpu_work = summarize_cpu_work(candidate_report["sources"])
        actual.update(actual_work_components(cpu_work))
        limited = any(not r["scientific_endpoint_available"] for r in rows)
        report = {"scope": config["scope"], "mode": "real_method_execution",
                  "selected_raw_n": len(parents), "rows": rows,
                  "qualification": qualification,
                  "policies": summary(rows, policies), "actual_new_work": actual,
                  "actual_new_work_scope": "Cumulative observed new work in this output directory; not replay work",
                  "actual_new_work_field_scopes": {
                      "forward_reactions": "Legacy candidate-formation-only field; excludes source replay. "
                          "Use forward_reactions_total when its component observations are complete.",
                      "source_replay_forward_reactions": "Direct native source-recovery RunReactants invocations.",
                      "candidate_forward_reactions": "New native candidate-formation RunReactants invocations; excludes cached attempts.",
                      "forward_reactions_total": "Source replay plus candidate formation; null if either observation is incomplete.",
                      "retro_template_calls": "Native reverse-template RunReactants invocations for root discovery.",
                      "source_RF_path_attributions": "Physical source RF path decompositions, shared by both arms; not RF predictions."},
                  "actual_new_cpu_work": cpu_work,
                  "new_graph_calls_this_invocation": store.new_calls_this_invocation,
                  "new_search_seconds_this_invocation": store.actual_search_seconds_this_invocation,
                  "stored_incomplete_cost_queries": sum(not r.get("cost_observation_complete", False) for r in records),
                  "new_graph_work_cost_observation_complete": all(r.get("cost_observation_complete", False) for r in new),
                  "historical12_queries": sum(bool(r.get("history_path")) and r["request"]["calls_limit"] == 12
                                              for r in records),
                  "exact_study_reused_queries": len(reused),
                  "old_graph_work": {"model_calls": sum(r.get("model_calls", 0) for r in reused),
                                     "search_seconds": sum(r.get("search_seconds", 0.) for r in reused),
                                     "cost_observation_complete": all(r.get("cost_observation_complete", False) for r in reused),
                                     "interpretation": "Observed old work; incomplete-cost totals are lower bounds, not measured zero effort"},
                  "new_saved_graph_audit_seconds_this_invocation": getattr(store, "new_reuse_audit_seconds_this_invocation", 0.),
                  "new_saved_graph_read_seconds_this_invocation": getattr(store, "new_reuse_graph_read_seconds_this_invocation", 0.),
                  "saved_graph_auditor_setup_seconds_this_invocation": getattr(getattr(store, "auditor", None), "setup_seconds", 0.),
                  "parent_verified_routes": sum(s["graph_verified"] for s in parents),
                  "parent_unknown_raw_slots": sum(s["valid"] and s["graph_result"]["status"] not in NORMAL for s in parents),
                  "versions": installed, "completed_unix": time.time(),
                  "route_evidence": "Each new solved search is saved and reloaded for native-route, stock-leaf and target-exclusion audit; each explicitly imported solved query is separately re-audited from its original saved graph in this package without loading an inference model.",
                  "interpretation": "A selected-source software round trip, not a new efficacy cohort. Logical shared45/6+6 budgets and actual directory costs differ."}
        save(output / "method_result.json", report)
        state.update(status="complete_with_execution_limits" if limited else "complete",
                     ended_unix=time.time(), new_graph_calls_this_invocation=store.new_calls_this_invocation)
        save(state_path, state)
        return int(limited)
    except BaseException:
        state.update(status="interrupted_or_failed", error=traceback.format_exc(), ended_unix=time.time(),
                     new_graph_calls_this_invocation=store.new_calls_this_invocation)
        save(state_path, state)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    def interrupted(signum, frame):
        raise KeyboardInterrupt("Termination signal " + str(signum))
    signal.signal(signal.SIGTERM, interrupted)
    code = run(args.output)
    result = read(args.output / "method_result.json")
    print(json.dumps({"mode": result["mode"], "policies": result["policies"],
                      "actual_new_work": result["actual_new_work"]}), flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
