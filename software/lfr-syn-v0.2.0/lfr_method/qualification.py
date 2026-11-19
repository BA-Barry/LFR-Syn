"""Named qualification modes; missing historical configuration remains strict."""
POLICIES_BY_QUALIFICATION = {
    "fixed_task": ("B_fixedtask6plus6", "C_fixedtask6plus6"),
    "source_preserving": ("B_joint6plus6", "C_joint6plus6"),
}


def settings(config):
    mode = config.get("qualification", "source_preserving")
    if mode not in POLICIES_BY_QUALIFICATION:
        raise ValueError("Unknown candidate qualification mode")
    policies = POLICIES_BY_QUALIFICATION[mode]
    if "policies" in config and tuple(config["policies"]) != policies:
        raise ValueError("Policy names conflict with qualification mode")
    return mode, policies


def evaluator(mode):
    if mode == "fixed_task":
        from .fixedtask_kernel import evaluate_policy
    elif mode == "source_preserving":
        from .graph_core import evaluate_policy
    else:
        raise ValueError("Unknown candidate qualification mode")
    return evaluate_policy


def report_qualification(row, mode):
    row.update(qualification=mode,
               guiding_RF_nondegradation_guarantee=mode == "source_preserving",
               guiding_RF_change=(row["final_score"] - row["source_score"]
                                  if row["scientific_endpoint_available"] else None))
    return row
