"""Experiment-local eligibility change; all policy execution code is retained.

No legacy module globals or candidate records are mutated. The original policy
function's code object runs with a private global namespace containing only two
replacements: the candidate eligibility function and the two revised arm IDs.
"""
import math
from types import FunctionType

from lfr_method import graph_core as legacy

POLICIES = ("B_fixedtask6plus6", "C_fixedtask6plus6")
ARMS = {POLICIES[0]: "task_contribution", POLICIES[1]: "permuted_contribution"}
LEGACY_POLICIES = {POLICIES[0]: "B_joint6plus6", POLICIES[1]: "C_joint6plus6"}
ELIGIBILITY = "candidate_RF >= fixed task cutoff; candidate_QED >= source_QED - 0.05"


def eligibility_records(source, arm):
    """Requalify saved measurements without consulting the strict legacy flag."""
    scored = arm.get("scored_candidates", [])
    if len(scored) > 8:
        raise ValueError("More than eight scored candidates in a policy source")
    seen, result = set(), []
    floor = legacy.THRESHOLDS[source["task"]]
    for candidate in scored:
        if candidate["smiles"] in seen:
            raise ValueError("Duplicate scored candidate identity")
        seen.add(candidate["smiles"])
        score, qed = float(candidate["task_score"]), float(candidate["qed"])
        if not math.isfinite(score) or not math.isfinite(qed):
            raise ValueError("Non-finite candidate score/QED")
        eligible = (score >= floor and qed >= float(source["qed"]) - .05
                    and candidate["smiles"] != source["canonical_smiles"])
        # task_pass stays the saved legacy flag. fixed_task_pass is the revised
        # scientific decision, so a false legacy flag never vetoes this arm.
        result.append(dict(candidate, legacy_task_pass=candidate.get("task_pass"),
                           fixed_task_pass=eligible, fixed_task_cutoff=floor))
    return result


def candidate_pool(source, arm):
    pool = [candidate for candidate in eligibility_records(source, arm)
            if candidate["fixed_task_pass"]]
    return sorted(pool, key=lambda c: (-c["task_score"], -c["qed"], c["smiles"]))[:2]


_local_globals = dict(legacy.evaluate_policy.__globals__,
                      POLICIES=POLICIES, candidate_pool=candidate_pool)
evaluate_policy = FunctionType(legacy.evaluate_policy.__code__, _local_globals,
                               "evaluate_policy", legacy.evaluate_policy.__defaults__,
                               legacy.evaluate_policy.__closure__)
evaluate_policy.__kwdefaults__ = legacy.evaluate_policy.__kwdefaults__
evaluate_policy.__doc__ = legacy.evaluate_policy.__doc__
