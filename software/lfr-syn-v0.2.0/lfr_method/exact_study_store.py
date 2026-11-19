"""Reuse exact Graph measurements from one explicitly supplied study.

The key is readable molecule/call-cap/requested-soft-time data, not a digest.
All saved states participate in first-path selection; unknowns are not skipped
in favour of later completed searches. Solved witnesses are audited afresh in
the current package without constructing an inference model.
"""
import copy
import json
import math
from pathlib import Path
import pickle
import time
import traceback

from .graph_core import NORMAL, QueryStore, ResidentGraph, read, save

PROJECT_DEV = Path("/home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914/development")
PROJECT_USER = Path("/home/chu")


def index_key(request):
    return json.dumps(list(QueryStore.chemical_key(request)), ensure_ascii=False,
                      separators=(",", ":"), allow_nan=False)


class SavedGraphAuditor:
    """The retained native audit, with stock/RDKit only and no model weights."""
    def __init__(self, stock):
        begin = time.monotonic()
        from rdkit import Chem
        from lfr_syn.shared_stock import InchiKeyInventory
        self.Chem = Chem
        self.stock = InchiKeyInventory.load_from_file(stock)
        self.setup_seconds = time.monotonic() - begin

    audit = ResidentGraph.audit


class ExactStudyStore(QueryStore):
    def scan_history(self, root):
        # Replace the generic terminal-only old12 scanner for explicit reuse.
        return {}

    def __init__(self, output, config, study, backend_factory, auditor_factory=None):
        self.study = Path(study).resolve()
        if not self.study.is_dir():
            raise ValueError("Explicit reuse study is missing: " + str(self.study))
        self.auditor_factory = auditor_factory or (lambda: SavedGraphAuditor(config["stock"]))
        self.auditor = None
        self.new_reuse_audit_seconds_this_invocation = 0.
        self.new_reuse_graph_read_seconds_this_invocation = 0.
        super().__init__(output, config, self.study, backend_factory)
        path = self.output / "explicit_study_exact_index.json"
        if path.exists():
            retained = read(path)
            if (retained["study"] != str(self.study)
                    or retained["settings"] != config["compatible12_settings"]):
                raise ValueError("Retained exact study index belongs to different settings/study")
            self.exact_index = retained["index"]
        else:
            self.exact_index, diagnostics = self.scan_explicit_study()
            save(path, {"study": str(self.study), "settings": config["compatible12_settings"],
                        "index": self.exact_index, "diagnostics": diagnostics,
                        "selection": "Lexicographically first physical query path for exact molecule/call cap/soft time, "
                                     "including unknown/error/incomplete states; equivalent included aliases resolve to origin.",
                        "created_unix": time.time()})

    def in_study(self, path):
        return Path(path).resolve().is_relative_to(self.study)

    def allowed_declared_path(self, path):
        """Only precise links to resources authorized by the existing project.

        This does not enumerate PROJECT_DEV or a trace directory. The root task
        explicitly authorized following an in-study record's physical/history
        origin into those already used project resources.
        """
        path = Path(path).resolve()
        if self.in_study(path) or path.is_relative_to(PROJECT_DEV.resolve()):
            return True
        if path.is_relative_to(PROJECT_USER.resolve()):
            relative = path.relative_to(PROJECT_USER.resolve())
            return bool(relative.parts and relative.parts[0].startswith("trace_lfr_"))
        return False

    @staticmethod
    def validate_request(request, path):
        if (not isinstance(request, dict) or type(request.get("calls_limit")) is not int
                or request["calls_limit"] not in (6, 12)
                or type(request.get("soft_seconds")) not in (int, float)
                or not math.isfinite(request["soft_seconds"]) or not 0 < request["soft_seconds"] <= 45.
                or type(request.get("seed")) is not int
                or not isinstance(request.get("smiles"), str) or not request["smiles"]):
            raise ValueError("Malformed exact saved Graph request: " + str(path))

    def normalized_request(self, record, current, path):
        if "request" in record:
            request = record["request"]
        else:
            # Actual legacy parent format. Prefix/mixed-controller records are
            # deliberately inapplicable; a cap12 measurement is never cap6.
            if "old_model_calls" in record or "new_model_calls" in record:
                return None
            request = {"smiles": record.get("canonical_smiles"), "seed": record.get("seed"),
                       "calls_limit": current.get("reaction_model_calls_limit"),
                       "soft_seconds": current.get("time_limit_s")}
        self.validate_request(request, path)
        return dict(request)

    def compatible_config(self, path, modern, cache):
        cfg = path.parent.parent / "config.json" if modern else path.parent / "config.json"
        if not self.allowed_declared_path(cfg) or not cfg.is_file():
            raise ValueError("Saved Graph request has no authorized corresponding configuration: " + str(path))
        if cfg not in cache:
            cache[cfg] = read(cfg)
        current = cache[cfg]
        wanted = self.config["compatible12_settings"]
        if modern:
            compatible = (current.get("compatible12_settings") == wanted
                          and current.get("shared_soft_seconds") == 45.)
        else:
            compatible = all(current.get(k) == v for k, v in wanted.items())
        return current if compatible else None, cfg

    @staticmethod
    def spent(record):
        return (record.get("model_calls", 0),
                record.get("search_seconds", record.get("total_seconds", 0.)),
                bool(record.get("cost_observation_complete", record.get("status") in NORMAL)))

    @staticmethod
    def semantic_audit(audit):
        if not isinstance(audit, dict):
            return audit
        # A re-audit's duration is new work; it is not the identity of the route.
        return {key: value for key, value in audit.items() if key not in ("seconds",)}

    def equivalent_alias(self, alias, origin):
        a, b = alias["record"], origin["record"]
        if (index_key(alias["request"]) != index_key(origin["request"])
                or alias["request"]["seed"] != origin["request"]["seed"]
                or self.spent(a) != self.spent(b)
                or a.get("status") != b.get("status")
                or a.get("solved", False) != b.get("solved", False)
                or bool(a.get("audit_error")) != bool(b.get("audit_error"))):
            return False
        if "request" in b:
            if (a.get("verified_solved") != b.get("verified_solved")
                    or a.get("audit_complete") != b.get("audit_complete")
                    or self.semantic_audit(a.get("audit")) != self.semantic_audit(b.get("audit"))):
                return False
        # Either the declared import includes the exact origin record, or all
        # modern request/work/endpoint/audit fields above establish equivalence.
        return a.get("original_result") == b or "request" in b

    @staticmethod
    def origin_name(record):
        if record.get("cache_mode") not in ("joint512_exact", "historical12", "explicit_study_exact"):
            return None
        for key in ("physical_origin_path", "study_query_path", "old512_query_path", "history_path"):
            value = record.get(key)
            if isinstance(value, str) and value:
                return value
        return None

    def scan_explicit_study(self):
        cache, compatible, diagnostics = {}, {}, {"query_files_seen": 0, "compatible_query_files": 0,
            "equivalent_import_aliases_skipped": 0, "incompatible_configs": [],
            "declared_external_origin_files_read": 0, "declared_origin_limitations": []}
        for path in sorted(self.study.rglob("*.result.json")):
            if self.output.resolve() in path.resolve().parents:
                continue
            if not self.in_study(path):
                raise ValueError("Explicit study contains an external query symlink: " + str(path))
            diagnostics["query_files_seen"] += 1
            record = read(path)
            modern = "request" in record
            current, cfg = self.compatible_config(path, modern, cache)
            if current is None:
                if str(cfg) not in diagnostics["incompatible_configs"]:
                    diagnostics["incompatible_configs"].append(str(cfg))
                continue
            request = self.normalized_request(record, current, path)
            if request is None:
                continue
            compatible[path.resolve()] = {"record": record, "request": request, "modern": modern}
            diagnostics["compatible_query_files"] += 1
        initial_paths = list(compatible)
        aliases, resolved = {}, {}
        def physical(path, seen=None):
            if path in resolved:
                return resolved[path]
            seen = set() if seen is None else set(seen)
            if path in seen:
                raise ValueError("Cyclic saved Graph import aliases")
            seen.add(path)
            entry = compatible[path]
            name = self.origin_name(entry["record"])
            origin = Path(name).resolve() if name else None
            if origin is None:
                resolved[path] = path
                return path
            try:
                if not origin.name.endswith(".result.json") or not self.allowed_declared_path(origin):
                    raise ValueError("Declared origin lies outside authorized project resources")
                if origin not in compatible:
                    # Read this specific declared file and its corresponding
                    # config only. Never rglob its directory or discover peers.
                    record = read(origin)
                    modern = "request" in record
                    current, cfg = self.compatible_config(origin, modern, cache)
                    if current is None:
                        raise ValueError("Declared origin has incompatible settings")
                    request = self.normalized_request(record, current, origin)
                    if request is None:
                        raise ValueError("Declared origin is a prefix/mixed measurement")
                    compatible[origin] = {"record": record, "request": request, "modern": modern}
                    if not self.in_study(origin):
                        diagnostics["declared_external_origin_files_read"] += 1
                if not self.equivalent_alias(entry, compatible[origin]):
                    raise ValueError("Declared alias/origin query, seed, work or state differs")
                target = physical(origin, seen)
                aliases[path] = origin
                diagnostics["equivalent_import_aliases_skipped"] += 1
                resolved[path] = target
                return target
            except (OSError, ValueError, TypeError, KeyError) as exc:
                # Keep the selected alias; import will retain an unknown if its
                # exact origin cannot be established, never try a better record.
                diagnostics["declared_origin_limitations"].append({"path": str(path), "origin": str(origin), "error": str(exc)})
                resolved[path] = path
                return path
        grouped = {}
        for path in initial_paths:
            origin = physical(path)
            entry = compatible[origin]
            grouped.setdefault(index_key(entry["request"]), set()).add(origin)
        index = {}
        for key, paths in grouped.items():
            chosen = min(paths, key=str)
            entry = compatible[chosen]
            index[key] = {"path": str(chosen), "request_at_index": entry["request"],
                          "status_at_index": entry["record"].get("status"),
                          "modern": entry["modern"], "saved_physical_query_count": len(paths)}
        diagnostics["exact_settings_keys"] = len(index)
        return index, diagnostics

    def selected_graph_path(self, selected_path, old):
        sibling = selected_path.with_name(selected_path.name.replace(".result.json", ".graph.pkl"))
        if sibling.is_file() and self.allowed_declared_path(sibling):
            return sibling
        raise ValueError("The physical query's corresponding saved Graph witness is unavailable")

    def verify_declared_origin(self, selected_path, old, old_request):
        path, entry, seen = selected_path.resolve(), {"record": old, "request": old_request}, set()
        while True:
            if path in seen:
                raise ValueError("Cyclic selected Graph import aliases")
            seen.add(path)
            name = self.origin_name(entry["record"])
            if not name:
                return path, entry["record"]
            origin = Path(name).resolve()
            if not origin.name.endswith(".result.json") or not self.allowed_declared_path(origin):
                raise ValueError("Selected import alias origin lies outside authorized project resources")
            physical = read(origin)
            modern = "request" in physical
            current, cfg = self.compatible_config(origin, modern, {})
            if current is None:
                raise ValueError("Selected import alias origin has incompatible settings")
            request = self.normalized_request(physical, current, origin)
            other = {"record": physical, "request": request}
            if request is None or not self.equivalent_alias(entry, other):
                raise ValueError("Selected import alias differs from its physical query; no prefix or status substitution")
            path, entry = origin, other

    def ensure_auditor(self):
        if self.auditor is None:
            self.auditor = self.auditor_factory()
        return self.auditor

    def import_selected(self, selected, requested):
        begin = time.monotonic()
        old_path = Path(selected["path"])
        row = self.next_row
        self.next_row += 1
        path = self.queries / f"query-{row:05d}.result.json"
        record = {"row": row, "request": dict(requested), "requested_seed": requested["seed"],
                  "adopted_seed": requested["seed"], "current_study_seed_adopted": False,
                  "status": "study_read_error", "model_calls": 0, "search_seconds": 0.,
                  "verified_solved": False, "cost_observation_complete": False,
                  "cache_mode": "explicit_study_exact", "history_path": str(old_path),
                  "study_query_path": str(old_path), "physical_origin_path": str(old_path),
                  "result_path": str(path), "actual_new_calls": 0,
                  "selected_study_record": selected, "audit_seconds": 0.,
                  "audit_graph_read_seconds": 0., "new_work": {"model_calls": 0}}
        try:
            if not self.allowed_declared_path(old_path):
                raise ValueError("Preselected query is outside authorized project resources")
            old = read(old_path)
            current, cfg = self.compatible_config(old_path, selected["modern"], {})
            if current is None:
                raise ValueError("Preselected saved Graph configuration changed")
            old_request = self.normalized_request(old, current, old_path)
            self.validate_request(old_request, old_path)
            if (old_request != selected["request_at_index"] or index_key(old_request) != index_key(requested)
                    or old.get("status") != selected["status_at_index"]):
                raise ValueError("Preselected saved Graph identity/status changed")
            seed = old_request["seed"]
            calls, elapsed, complete = self.spent(old)
            record.update(request=dict(requested, seed=seed), adopted_seed=seed,
                          study_seed_adopted=seed != requested["seed"], original_result=copy.deepcopy(old),
                          status=old.get("status"), solved=old.get("solved", False),
                          model_calls=calls, search_seconds=elapsed, cost_observation_complete=complete,
                          first_solution_calls=old.get("first_solution_calls"), graph_nodes=old.get("graph_nodes"),
                          original_graph_path=old.get("graph_path"), original_result_path=old.get("result_path"),
                          original_audit=copy.deepcopy(old.get("audit")))
            if (type(calls) is not int or calls < 0 or type(elapsed) not in (int, float)
                    or elapsed < 0 or not math.isfinite(elapsed)):
                raise ValueError("Preselected saved Graph cost is invalid")
            record["old_work"] = {"model_calls": calls, "search_seconds": elapsed,
                                  "cost_observation_complete": complete, "query_path": str(old_path)}
            if old.get("status") in NORMAL:
                physical_path, physical_record = self.verify_declared_origin(old_path, old, old_request)
                record["physical_origin_path"] = str(physical_path)
                if (old.get("solved") is not (old["status"] == "solved")
                        or not complete or calls > requested["calls_limit"]):
                    raise ValueError("Preselected normal query has invalid endpoint/cost observation")
                if selected["modern"]:
                    audit = old.get("audit")
                    if (not old.get("audit_complete") or old.get("verified_solved") is not old["solved"]
                            or not isinstance(audit, dict) or audit.get("verified_saved_graph") is not old["solved"]
                            or old.get("audit_error") or audit.get("audit_error")
                            or (old["solved"] and audit.get("target") != requested["smiles"])):
                        raise ValueError("Preselected query's saved endpoint/audit is unavailable or inconsistent")
                else:
                    record["historical_audit"] = self.historical_audit(old_path, old)
                graph_path = self.selected_graph_path(physical_path, physical_record)
                record["graph_path"] = str(graph_path)
                if old["solved"]:
                    reading = time.monotonic()
                    try:
                        with graph_path.open("rb") as stream:
                            graph = pickle.load(stream)
                    finally:
                        record["audit_graph_read_seconds"] = time.monotonic() - reading
                        self.new_reuse_graph_read_seconds_this_invocation += record["audit_graph_read_seconds"]
                    auditor = self.ensure_auditor()
                    auditing = time.monotonic()
                    try:
                        record["audit"] = auditor.audit(graph, requested["smiles"])
                    finally:
                        record["audit_seconds"] = time.monotonic() - auditing
                        self.new_reuse_audit_seconds_this_invocation += record["audit_seconds"]
                    if record["audit"].get("verified_saved_graph") is not True:
                        raise ValueError("Current package's saved-graph audit did not verify the route")
                    record["verified_solved"] = True
                    record["study_saved_graph_reaudited"] = True
                else:
                    record["audit"] = {"reported_solved": False, "verified_saved_graph": False, "seconds": 0.}
                record["audit_complete"] = True
            else:
                if old.get("status") == "running":
                    record.update(status="saved_incomplete_query", cost_observation_complete=False,
                                  cost_note="Unobserved interruption: counters are lower bounds, not measured zero work")
                elif old.get("status") == "search_complete":
                    record["status"] = "saved_study_audit_incomplete"
        except BaseException as exc:
            record.update(status="study_reuse_limitation", verified_solved=False,
                          error=traceback.format_exc(), interrupted=not isinstance(exc, Exception))
            # Preserve the raw record; invalid counters are explicitly unknown,
            # not available as measured logical zero cost or a retry trigger.
            if type(record["model_calls"]) is not int or record["model_calls"] < 0:
                record["model_calls"] = 0
                record["cost_observation_complete"] = False
            if (type(record["search_seconds"]) not in (int, float) or record["search_seconds"] < 0
                    or not math.isfinite(record["search_seconds"])):
                record["search_seconds"] = 0.
                record["cost_observation_complete"] = False
        record["new_work"].update(graph_read_seconds=record["audit_graph_read_seconds"],
                                    graph_audit_seconds=record["audit_seconds"])
        save(path, record)
        self.cache[self.key(record["request"])] = (path, record)
        if record.get("interrupted"):
            raise KeyboardInterrupt("Interrupted saved-graph import retained at " + str(path))
        return dict(record, actual_new_calls=0, retrieval_seconds=time.monotonic() - begin)

    def __call__(self, smiles, seed, calls, seconds):
        request = {"smiles": smiles, "seed": int(seed), "calls_limit": int(calls), "soft_seconds": float(seconds)}
        if self.first_current_query(request):
            return super().__call__(smiles, seed, calls, seconds)
        selected = self.exact_index.get(index_key(request))
        if selected is not None:
            return self.import_selected(selected, request)
        return super().__call__(smiles, seed, calls, seconds)
