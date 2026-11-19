"""Summarize saved CPU counters without chemistry or inference imports.

Counts are cumulative for this output directory, including completed records
loaded on continuation. Missing legacy fields are unknown, not measured zero.
"""

ARMS = ("task_contribution", "permuted_contribution")
SKIPPED_FORMATION = {
    "retained_invalid", "retained_joint_success", "parent_execution_limitation"
}
SOURCE_COUNTERS = {
    "source_replay_forward_reactions": "source_replay_forward_calls",
    "retro_template_calls": "retro_template_calls",
    "source_RF_path_attributions": "source_RF_path_attributions",
}


def _count(value):
    return value if type(value) is int and value >= 0 else None


def _component(observations):
    missing = [location for location, value in observations if value is None]
    observed = sum(value for _, value in observations if value is not None)
    return {"observed_count": observed, "count": None if missing else observed,
            "cost_observation_complete": not missing, "missing_counter_fields": missing}


def _skipped_formation(record):
    if record.get("status") not in SKIPPED_FORMATION:
        return False
    costs = record.get("costs", {})
    if any(costs.get(field) not in (None, 0) for field in SOURCE_COUNTERS.values()):
        return False
    # A retained label with actual formation data is not evidence of zero work.
    for arm in record.get("arms", {}).values():
        if any(arm.get(field) for field in
               ("planned_candidates", "scored_candidates", "forward_attempts")):
            return False
        if arm.get("costs", {}).get("new_forward_attempts") not in (None, 0):
            return False
    return True


def summarize_cpu_work(sources):
    """Read existing source/arm counters; never reconstruct native operations.

    Root discovery has direct native RunReactants counters; candidate formation
    has separate new/cached counters. A terminal retained-source branch skips
    formation, so absent root counters there denote known zero. For a formed
    source, an absent counter is unknown even if its scientific outcome is known.
    """
    observations = {name: [] for name in SOURCE_COUNTERS}
    observations["candidate_forward_reactions"] = []
    for index, record in enumerate(sources):
        sid = record.get("source", {}).get("sample_id", f"source-index-{index}")
        skipped = _skipped_formation(record)
        costs = record.get("costs", {})
        for name, field in SOURCE_COUNTERS.items():
            value = _count(costs.get(field))
            if field not in costs and skipped:
                value = 0
            observations[name].append((f"{sid}/costs/{field}", value))
        for arm in ARMS:
            arm_costs = record.get("arms", {}).get(arm, {}).get("costs", {})
            field = "new_forward_attempts"
            value = _count(arm_costs.get(field))
            if field not in arm_costs and skipped:
                value = 0
            observations["candidate_forward_reactions"].append(
                (f"{sid}/arms/{arm}/costs/{field}", value))
    components = {name: _component(rows) for name, rows in observations.items()}
    forward = _component(observations["source_replay_forward_reactions"] +
                         observations["candidate_forward_reactions"])
    return {"scope": "Cumulative recorded new physical CPU work in this output directory; "
                     "completed source records retain their original costs on continuation.",
            "source_records": len(sources), "components": components,
            "forward_reactions_total": forward,
            "cost_observation_complete": all(
                part["cost_observation_complete"] for part in components.values()),
            "interpretation": "count is null if a required legacy counter is absent or invalid; "
                              "observed_count is only the sum of available nonnegative integer "
                              "counters, not a complete total. Candidate-only legacy aggregates "
                              "are retained separately and do not fill missing source counters."}


def actual_work_components(cpu_work):
    """Project complete component counts while retaining unknowns as null."""
    fields = {name: part["count"] for name, part in cpu_work["components"].items()}
    fields["forward_reactions_total"] = cpu_work["forward_reactions_total"]["count"]
    return fields
