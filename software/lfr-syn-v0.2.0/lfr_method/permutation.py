"""Core extracted from retained joint_validation512_candidates.py; see code_origins.json."""
import random


def permutation(values, source, root_index, group_index):
    seed = (10202000 + source['task_index'] * 1000000 + source['replicate'] * 10000
            + source['position'] * 100 + root_index * 10 + group_index)
    result = list(values)
    random.Random(seed).shuffle(result)
    return result, seed
