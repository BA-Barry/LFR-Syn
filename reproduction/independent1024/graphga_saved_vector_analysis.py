"""Analyze fixed GraphGA Original/B/C rows; no scientific model is imported.

Input rows retain terminal draw positions, including duplicate molecules. The
uncertainty unit is the complete optimizer run, never an individual molecule.
"""
import argparse
from collections import Counter, defaultdict
import csv
import itertools
import json
import math
from pathlib import Path
import platform
import statistics
import time

POLICIES = ("Original", "B_joint6plus6", "C_joint6plus6")
ENDPOINTS = ("external_joint", "joint", "solved")
TASK_NAMES = ("gsk3b_activity", "jnk3_activity")
MODEL_REFERENCE = {
    "repository": "https://github.com/Tibogoss/CombiMOTS",
    "published_commit": "3f644b727982a95b4f87115ccb8d6b4c918a02ba",
    "software": "Chemprop1.5.2", "tasks": list(TASK_NAMES), "models": 10,
}


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def finite_probability(value):
    return (isinstance(value, (float, int)) and not isinstance(value, bool)
            and math.isfinite(value) and 0 <= value <= 1)


def endpoint_bounds(row, endpoint):
    value = row[endpoint]
    if isinstance(value, bool):
        return [int(value), int(value)]
    if value is not None:
        raise ValueError("Endpoint must be an explicit boolean or null: " + endpoint)
    lower = row.get("endpoint_lower_" + endpoint, 0)
    upper = row.get("endpoint_upper_" + endpoint, 1)
    if lower not in (0, 1) or upper not in (0, 1) or lower > upper:
        raise ValueError("Invalid endpoint bounds: " + endpoint)
    return [int(lower), int(upper)]


def check_membership(rows, design):
    """Ordinary ID/position checks prevent lost duplicate draw positions."""
    if tuple(design["policies"]) != POLICIES:
        raise ValueError("Expected prospective Original/B/C policies")
    specs = {s["cohort_id"]: s for s in design["runs"]}
    if len(specs) != 8 or Counter(s["task"] for s in specs.values()) != {"gsk3b": 4, "jnk3": 4}:
        raise ValueError("Expected eight optimizer runs, four per fixed task")
    n = design["source_positions_per_run"]
    if n != 64 or design["raw_source_positions"] != 512:
        raise ValueError("This planned analysis uses 64 terminal draws per run")
    grouped = defaultdict(list)
    identities = {}
    seen = set()
    for row in rows:
        key = (row["policy"], row["sample_id"])
        if key in seen or row["policy"] not in POLICIES:
            raise ValueError("Duplicate or unexpected policy/sample identifier")
        seen.add(key)
        spec = specs.get(row["cohort_id"])
        if spec is None or (row["task"], row["replicate"]) != (spec["task"], spec["replicate"]):
            raise ValueError("Row does not belong to a specified optimizer run")
        if type(row["position"]) is not int or row["position"] not in range(n):
            raise ValueError("Invalid raw terminal draw position")
        identity = tuple(row[k] for k in ("cohort_id", "task", "replicate", "position",
                                         "source_smiles", "source_valid", "source_score", "source_qed"))
        previous = identities.setdefault(row["sample_id"], identity)
        if previous != identity:
            raise ValueError("Policy source identity mismatch")
        grouped[(row["policy"], row["cohort_id"])].append(row)
        if type(row["source_valid"]) is not bool or type(row["scientific_endpoint_available"]) is not bool:
            raise ValueError("Explicit source-valid and route-execution flags required")
        if not finite_probability(row["source_score"]) or not finite_probability(row["final_score"]):
            raise ValueError("Invalid saved RF probability")
        for name in ("source_qed", "final_qed"):
            if row[name] is not None and not finite_probability(row[name]):
                raise ValueError("Invalid saved QED")
        if row["source_valid"] != bool(row["source_smiles"]):
            raise ValueError("Valid source must have its canonical structure")
        if row["final_smiles"] is not None and (not isinstance(row["final_smiles"], str) or not row["final_smiles"]):
            raise ValueError("Final structure must be a canonical string or null")
        for endpoint in ("joint", "solved"):
            endpoint_bounds(row, endpoint)
            if row["scientific_endpoint_available"] != isinstance(row[endpoint], bool):
                raise ValueError("Route/RF endpoint availability mismatch")
        if row["scientific_endpoint_available"]:
            if row["solved"] and not row["final_smiles"]:
                raise ValueError("Audited route endpoint lacks a final structure")
            expected = bool(row["solved"] and row["source_valid"]
                            and row["final_score"] >= design["method"]["source_RF_cutoffs"][row["task"]])
            if row["joint"] != expected:
                raise ValueError("Saved guiding-RF joint endpoint mismatch")
        if row["replaced"] and (not row["scientific_endpoint_available"] or not row["solved"]
                or row["final_smiles"] == row["source_smiles"]
                or row["final_score"] < row["source_score"]
                or row["final_qed"] is None or row["source_qed"] is None
                or row["final_qed"] < row["source_qed"] - .05):
            raise ValueError("Replacement violates saved source-relative qualification")
    if len(rows) != 1536 or len(grouped) != 24:
        raise ValueError("Need the complete 512 terminal positions for each of three policies")
    for policy in POLICIES:
        for cohort in specs:
            selected = grouped[(policy, cohort)]
            if len(selected) != n or {r["position"] for r in selected} != set(range(n)):
                raise ValueError("Incomplete or duplicated terminal draw position")
    # A different sample ID at a matching position must not evade pairing checks.
    for cohort in specs:
        reference = {(r["position"], r["sample_id"]) for r in grouped[("Original", cohort)]}
        if any({(r["position"], r["sample_id"]) for r in grouped[(p, cohort)]} != reference for p in POLICIES[1:]):
            raise ValueError("Policies do not have the same raw draw identifiers")
    return {p: [r for r in rows if r["policy"] == p] for p in POLICIES}


def vector_records(external):
    """Accept saved 10x2 vectors or the existing predictions.json format.

    Fold vectors, rather than rounded/means-only scores, are authoritative here.
    Model compatibility during cache acquisition remains the retained reader's
    responsibility; a supplied model reference, when present, must agree.
    """
    tasks = external["tasks"]
    if tasks not in (list(TASK_NAMES), ["gsk3b", "jnk3"]):
        raise ValueError("Wrong or reordered two-task vector columns")
    if external.get("model_reference") is not None and external["model_reference"] != MODEL_REFERENCE:
        raise ValueError("External scores use a different published model reference")
    if external.get("checkpoint_count", external.get("models", 10)) != 10:
        raise ValueError("Expected all ten checkpoint predictions")
    records = external.get("predictions", external.get("vectors"))
    if not isinstance(records, dict):
        raise ValueError("External JSON needs predictions or vectors by canonical structure")
    parsed = {}
    for smiles, record in records.items():
        if not isinstance(smiles, str) or not smiles:
            raise ValueError("Invalid external structure key")
        folds = record.get("folds") if isinstance(record, dict) else record
        if (not isinstance(folds, list) or len(folds) != 10
                or any(not isinstance(v, list) or len(v) != 2
                       or not all(map(finite_probability, v)) for v in folds)):
            raise ValueError("External score needs ten finite two-task fold vectors: " + smiles)
        means = [statistics.mean(v[i] for v in folds) for i in range(2)]
        if isinstance(record, dict) and "mean" in record:
            if (not isinstance(record["mean"], list) or len(record["mean"]) != 2
                    or not all(map(finite_probability, record["mean"]))
                    or any(abs(record["mean"][i] - means[i]) > 1e-7 for i in range(2))):
                raise ValueError("Saved external means conflict with fold vectors")
        parsed[smiles] = {"mean": means, "folds": folds,
                          "fold_provenance": record.get("fold_provenance", []) if isinstance(record, dict) else []}
    return parsed


def enrich(rows, external):
    vectors = vector_records(external)
    missing = set()
    enriched = []
    for original in rows:
        row = dict(original)
        index = ("gsk3b", "jnk3").index(row["task"])
        for prefix, smiles in (("source_", row["source_smiles"]), ("", row["final_smiles"])):
            record = vectors.get(smiles)
            row[prefix + "external_probability"] = record["mean"][index] if record else (0. if not smiles else None)
            if smiles and record is None:
                missing.add(smiles)
        row["external_score_available"] = row["final_smiles"] in vectors or not row["final_smiles"]
        probability = row["external_probability"]
        for name, profile_key in (("final_scaffold", "descriptors"), ("source_scaffold", "source_descriptors")):
            if isinstance(row.get(profile_key), dict):
                row[name] = row[profile_key].get("murcko_scaffold")
        if row["scientific_endpoint_available"] and row["solved"] is False:
            row["external_joint"] = False
        elif row["scientific_endpoint_available"] and probability is not None:
            # Guiding RF qualification is deliberately absent from this endpoint.
            row["external_joint"] = bool(row["solved"] and probability >= .5)
        else:
            row["external_joint"] = None
            row["endpoint_lower_external_joint"] = 0
            row["endpoint_upper_external_joint"] = endpoint_bounds(row, "solved")[1]
        enriched.append(row)
    return enriched, {"supplied_complete_vectors": len(vectors), "missing_required_structures": sorted(missing),
                      "arithmetic": "Unrounded arithmetic mean of the ten saved task-column probabilities.",
                      "scientific_model_calls": 0}


def percentile(values, probability):
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    low = int(position)
    fraction = position - low
    return ordered[low] * (1 - fraction) + ordered[min(low + 1, len(ordered) - 1)] * fraction


def exact_run_bootstrap(effects):
    """All 4**4 resamples per fixed task, equally weighted: 65,536."""
    tasks = sorted(effects)
    if len(tasks) != 2 or any(len(effects[t]) != 4 for t in tasks):
        raise ValueError("Four complete optimizer-run effects per fixed task required")
    means = [[statistics.mean(sample) for sample in itertools.product(effects[t], repeat=4)] for t in tasks]
    distribution = [(a + b) / 2 for a in means[0] for b in means[1]]
    return {"effect_pp": statistics.mean(statistics.mean(effects[t]) for t in tasks) * 100,
            "interval_95_pp": [percentile(distribution, .025) * 100, percentile(distribution, .975) * 100],
            "resample_combinations": len(distribution),
            "method": "Empirical percentile bootstrap of paired complete optimizer runs within two fixed tasks; equal task weighting; all ordered resamples enumerated.",
            "scope": "Descriptive uncertainty conditional on these predictors/assets and four runs per task; not a randomization test, power guarantee or biological validation."}


def summarize(rows):
    result = {"raw_n": len(rows), "route_execution_unknown_slots": sum(not r["scientific_endpoint_available"] for r in rows),
              "missing_external_score_slots": sum(not r["external_score_available"] for r in rows),
              "distinct_source_task_structures": len({(r["task"], r["source_smiles"]) for r in rows if r["source_smiles"]}),
              "distinct_final_task_structures_observed": len({(r["task"], r["final_smiles"]) for r in rows if r["final_smiles"]}),
              "replaced_slots": sum(r["replaced"] for r in rows)}
    for endpoint in ENDPOINTS:
        bounds = [endpoint_bounds(r, endpoint) for r in rows]
        successful = {(r["task"], r["final_smiles"]) for r in rows if r[endpoint] is True}
        unknown_upper = sum(r[endpoint] is None and endpoint_bounds(r, endpoint)[1] for r in rows)
        result[endpoint] = {"known_successes": sum(r[endpoint] is True for r in rows),
                            "raw_count_bounds": [sum(b[i] for b in bounds) for i in (0, 1)],
                            "raw_rate_bounds": [sum(b[i] for b in bounds) / len(rows) for i in (0, 1)],
                            "unknown_slots": sum(r[endpoint] is None for r in rows),
                            "distinct_task_structure_known_successes": len(successful),
                            "distinct_task_structure_count_bounds": [len(successful), len(successful) + unknown_upper]}
    # Distinct output structures across tasks are also descriptive, without task pairing.
    result["distinct_count_note"] = "Task/canonical-structure identities; duplicates remain in raw yields. Unknown final structures are not inferred from operational fallback identities. Upper bounds allow one additional structure per potentially successful unknown slot."
    for name in ("source_qed", "final_qed"):
        values = [r[name] for r in rows if r[name] is not None]
        result[name + "_observed"] = {"n": len(values), "mean": statistics.mean(values) if values else None,
                                      "median": statistics.median(values) if values else None}
    if any("final_scaffold" in r for r in rows):
        result["nonempty_final_scaffolds_observed"] = len({r["final_scaffold"] for r in rows if r.get("final_scaffold")})
        result["scaffold_changed_replacements"] = sum(r.get("source_scaffold") != r.get("final_scaffold") for r in rows if r["replaced"] and "source_scaffold" in r and "final_scaffold" in r)
        for endpoint in ENDPOINTS:
            result[endpoint]["nonempty_known_success_scaffolds"] = len({r["final_scaffold"] for r in rows if r[endpoint] is True and r.get("final_scaffold")})
    profiles = [r["descriptors"] for r in rows if isinstance(r.get("descriptors"), dict)]
    if profiles:
        complete = [p for p in profiles if all(k in p for k in ("molecular_weight", "pains", "brenk", "components"))]
        result["observed_profile_diagnostics"] = {"profile_slots": len(profiles), "complete_filter_profile_slots": len(complete),
            "pains_flagged_slots": sum(bool(p["pains"]) for p in complete),
            "brenk_flagged_slots": sum(bool(p["brenk"]) for p in complete),
            "mw_above500_slots": sum(p["molecular_weight"] > 500 for p in complete),
            "multiple_component_slots": sum(p["components"] > 1 for p in complete),
            "note": "Observed final/fallback structure diagnostics; no row is removed from the primary endpoint and unknown replacements may have other descriptors."}
    labeled = [r for r in rows if r.get("public_overlap_available")]
    if labeled:
        result["public_table_overlap_diagnostics"] = {"observed_slots_with_lookup": len(labeled),
            "task_labeled_overlap_slots": sum(bool(r.get("public_task_labeled_overlap")) for r in labeled),
            "note": "Descriptive matches to the available public table; not complete training membership or independence evidence."}
    return result


def contrast(first, second, endpoint, policies, design):
    indexes = {p: {r["sample_id"]: r for r in policies[p]} for p in (first, second)}
    grouped = defaultdict(list)
    for row in policies[first]:
        grouped[row["cohort_id"]].append(row)
    runs, effects = [], defaultdict(list)
    lower = upper = wins = losses = unknown = 0
    for spec in design["runs"]:
        selected = sorted(grouped[spec["cohort_id"]], key=lambda r: r["position"])
        lo = hi = rw = rl = ru = 0
        for a in selected:
            b = indexes[second][a["sample_id"]]
            al, au = endpoint_bounds(a, endpoint)
            bl, bu = endpoint_bounds(b, endpoint)
            lo += al - bu
            hi += au - bl
            if a[endpoint] is None or b[endpoint] is None:
                ru += 1
            else:
                rw += int(a[endpoint] and not b[endpoint])
                rl += int(b[endpoint] and not a[endpoint])
        lower += lo; upper += hi; wins += rw; losses += rl; unknown += ru
        runs.append({"cohort_id": spec["cohort_id"], "task": spec["task"], "replicate": spec["replicate"],
                     "raw_n": len(selected), "net_count_bounds": [lo, hi], "unknown_pairs": ru,
                     "known_wins": rw, "known_losses": rl})
        effects[spec["task"]].append(lo / len(selected))
    result = {"first": first, "second": second, "endpoint": endpoint, "raw_n": len(policies[first]),
              "net_count_bounds": [lower, upper], "effect_pp_bounds": [lower * 100 / 512, upper * 100 / 512],
              "known_wins": wins, "known_losses": losses, "unknown_pairs": unknown,
              "runs": runs, "run_level_bootstrap": None}
    if unknown == 0:
        result["run_level_bootstrap"] = exact_run_bootstrap(effects)
        result.update(positive_runs=sum(r["net_count_bounds"][0] > 0 for r in runs),
                      negative_runs=sum(r["net_count_bounds"][0] < 0 for r in runs),
                      tied_runs=sum(r["net_count_bounds"][0] == 0 for r in runs))
    else:
        result["interval_scope"] = "Technical bounds only; no complete-endpoint bootstrap is reported."
    return result


def analyze(rows, design, external):
    check_membership(rows, design)
    enriched, external_audit = enrich(rows, external)
    policies = {p: [r for r in enriched if r["policy"] == p] for p in POLICIES}
    result = {"experiment": design["experiment"], "scope": "Fixed completed GraphGA output analysis; no generation, inference, chemistry, route search or witness re-audit.",
              "raw_source_positions": 512, "independent_optimizer_runs": 8, "policy_rows": len(rows),
              "primary_contrast": "B minus C: ten-checkpoint mean task probability >= .5 AND audited route; guiding RF is absent.",
              "policies": {p: summarize(selected) for p, selected in policies.items()},
              "by_optimizer_run": {s["cohort_id"]: {p: summarize([r for r in policies[p] if r["cohort_id"] == s["cohort_id"]]) for p in POLICIES} for s in design["runs"]},
              "primary": contrast("B_joint6plus6", "C_joint6plus6", "external_joint", policies, design),
              "secondary_B_minus_C": {e: contrast("B_joint6plus6", "C_joint6plus6", e, policies, design) for e in ("joint", "solved")},
              "descriptive_B_minus_Original": {e: contrast("B_joint6plus6", "Original", e, policies, design) for e in ENDPOINTS},
              "external_vector_arithmetic": external_audit,
              "limitations": ["Ten Chemprop checkpoints share the broad ExCAPE source with the optimized RF; incomplete original training composition prevents an independent biological or reconstructed out-of-training claim.",
                              "Terminal positions are uniform draws with replacement from each optimized complete terminal population. They are not independent runs or raw mutation proposals.",
                              "Distinct structures and scaffolds are descriptive yield checks; raw positions remain the prespecified paired denominator.",
                              "Route evidence is inherited from saved method results; this analysis does not reload or re-audit native graphs."]}
    original = policies["Original"]
    fixed = sum(r["source_valid"] and r["source_score"] >= design["method"]["source_RF_cutoffs"][r["task"]] for r in original)
    result["fixed_source_RF_ceiling"] = {"raw_qualified_positions": fixed,
        "distinct_qualified_task_structures": len({(r["task"], r["source_smiles"]) for r in original if r["source_valid"] and r["source_score"] >= design["method"]["source_RF_cutoffs"][r["task"]]}),
        "B_joint_minus_raw_ceiling_bounds": [v - fixed for v in result["policies"]["B_joint6plus6"]["joint"]["raw_count_bounds"]],
        "interpretation": "Maximum raw RF-qualified yield reachable by retaining only these source positions, even if all their routes were solved. It is not an external-score ceiling."}
    return result, enriched


def load_rows(path):
    path = Path(path)
    if path.suffix == ".jsonl":
        return [json.loads(line) for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    data = read(path)
    return data["rows"] if isinstance(data, dict) else data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("design", "rows", "external", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    result, rows = analyze(load_rows(args.rows), read(args.design), read(args.external))
    result["runtime"] = {"python": platform.python_version(), "elapsed_seconds": time.perf_counter() - started,
                         "scientific_model_calls": 0}
    result["input_paths"] = {name: str(getattr(args, name).resolve()) for name in ("design", "rows", "external")}
    args.output.mkdir(parents=True, exist_ok=True)
    save(args.output / "analysis.json", result)
    with (args.output / "enriched_rows.jsonl").open("w", encoding="utf-8") as stream:
        for row in rows:
            stream.write(json.dumps(row, allow_nan=False) + "\n")
    with (args.output / "policy_yields.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("policy", "endpoint", "raw_n", "known_successes", "raw_lower", "raw_upper", "distinct_known_successes", "unknown_slots"))
        writer.writeheader()
        for policy, data in result["policies"].items():
            for endpoint in ENDPOINTS:
                counts = data[endpoint]
                writer.writerow({"policy": policy, "endpoint": endpoint, "raw_n": data["raw_n"],
                                 "known_successes": counts["known_successes"], "raw_lower": counts["raw_count_bounds"][0],
                                 "raw_upper": counts["raw_count_bounds"][1], "distinct_known_successes": counts["distinct_task_structure_known_successes"],
                                 "unknown_slots": counts["unknown_slots"]})
    print(json.dumps({"output": str(args.output), "primary": result["primary"], "runtime": result["runtime"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
