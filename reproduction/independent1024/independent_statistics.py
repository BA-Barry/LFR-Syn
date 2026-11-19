"""Prespecified paired optimizer-run statistics; standard library only.

Positions and checkpoints are never resampled as independent observations.
The two fixed tasks have equal weight. Technical unknowns retain identification
bounds and do not enter a conventional complete-data test or interval.
"""
import math
import random
import statistics

TASKS = ('gsk3b', 'jnk3')
BOOTSTRAP_REPLICATES = 100000
BOOTSTRAP_SEED = 2026101801


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _beta_fraction(a, b, x):
    tiny = 1e-300
    c = 1.
    d = 1. - (a + b) * x / (a + 1.)
    d = 1. / (d if abs(d) >= tiny else tiny)
    h = d
    for m in range(1, 501):
        aa = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
        d = 1. + aa * d
        c = 1. + aa / c
        d = 1. / (d if abs(d) >= tiny else tiny)
        c = c if abs(c) >= tiny else tiny
        h *= d * c
        aa = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1. + aa * d
        c = 1. + aa / c
        d = 1. / (d if abs(d) >= tiny else tiny)
        c = c if abs(c) >= tiny else tiny
        change = d * c
        h *= change
        if abs(change - 1.) < 3e-14:
            return h
    raise ArithmeticError('Incomplete-beta continued fraction did not converge')


def regularized_beta(x, a, b):
    require(0 <= x <= 1 and a > 0 and b > 0, 'Invalid beta parameters')
    if x == 0 or x == 1:
        return float(x)
    scale = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                     + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1) / (a + b + 2):
        return scale * _beta_fraction(a, b, x) / a
    return 1. - scale * _beta_fraction(b, a, 1. - x) / b


def student_tail(t, df):
    """One tail at |t|, computed without cancellation for large |t|."""
    require(df > 0 and math.isfinite(df) and math.isfinite(t), 'Invalid Student t input')
    square = t * t
    if square <= df:
        # df/(df+t*t) can round to one near zero, losing the entire t term.
        # Its complement has no subtraction and this branch's tail is not tiny.
        return .5 * (1. - regularized_beta(square / (df + square), .5, df / 2.))
    return .5 * regularized_beta(df / (df + square), df / 2., .5)


def student_quantile(probability, df):
    require(.5 < probability < 1 and math.isfinite(probability)
            and df > 0 and math.isfinite(df), 'Upper Student quantile required')
    lo, hi = 0., 1.
    while 1. - student_tail(hi, df) < probability:
        hi *= 2
    for _ in range(90):
        mid = (lo + hi) / 2
        if 1. - student_tail(mid, df) < probability:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def percentile(ordered, q):
    require(bool(ordered) and 0 <= q <= 1, 'Nonempty percentile input required')
    position = (len(ordered) - 1) * q
    lo, hi = math.floor(position), math.ceil(position)
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (position - lo)


def _task_summary(values):
    mean = statistics.fmean(values)
    variance = statistics.variance(values)
    result = {'optimizer_runs': len(values), 'effects_pp': [100 * x for x in values],
              'mean_effect_pp': 100 * mean, 'run_variance': variance,
              'interval_95_pp': None, 'p_two_sided': None}
    if variance == 0:
        result['interval_status'] = 'not_estimable_zero_observed_run_variance'
        return result
    se = math.sqrt(variance / len(values))
    critical = student_quantile(.975, len(values) - 1)
    result.update(interval_status='complete_data_paired_run_t',
                  interval_95_pp=[100 * (mean - critical * se), 100 * (mean + critical * se)],
                  p_two_sided=min(1., 2 * student_tail(mean / se, len(values) - 1)))
    return result


def welch_task_equal(effects):
    require(set(effects) == set(TASKS) and all(len(effects[t]) == 8 for t in TASKS),
            'Eight complete paired optimizer-run effects per fixed task required')
    require(all(math.isfinite(x) for t in TASKS for x in effects[t]), 'Nonfinite run effect')
    tasks = {task: _task_summary(effects[task]) for task in TASKS}
    mean = statistics.fmean(statistics.fmean(effects[t]) for t in TASKS)
    variance_terms = [statistics.variance(effects[t]) / len(effects[t]) for t in TASKS]
    sum_variance = sum(variance_terms)
    result = {'task_equal_mean_effect_pp': 100 * mean, 'task_summaries': tasks,
              'interval_95_pp': None, 'p_two_sided': None, 'alpha': .05,
              'method': 'Task-equal paired optimizer-run mean; Welch-Satterthwaite t interval/test',
              'standard_error_pp': 50 * math.sqrt(sum_variance),
              'scope': 'Independent optimizer runs conditional on the two fixed tasks and retained predictors/assets; no biological validation claim.'}
    if sum_variance == 0:
        result.update(interval_status='not_estimable_zero_observed_run_variance', degrees_of_freedom=None)
        return result
    df = sum_variance ** 2 / sum(term ** 2 / 7 for term in variance_terms)
    se = math.sqrt(sum_variance) / 2
    critical = student_quantile(.975, df)
    result.update(interval_status='complete_data_welch', degrees_of_freedom=df,
                  interval_95_pp=[100 * (mean - critical * se), 100 * (mean + critical * se)],
                  p_two_sided=min(1., 2 * student_tail(mean / se, df)),
                  t_statistic=mean / se)
    return result


def stratified_run_bootstrap(effects, *, replicates=BOOTSTRAP_REPLICATES, seed=BOOTSTRAP_SEED):
    require(set(effects) == set(TASKS) and all(len(effects[t]) == 8 for t in TASKS),
            'Bootstrap requires eight paired optimizer-run effects per task')
    require(type(replicates) is int and replicates > 0, 'Positive bootstrap replicate count required')
    rng = random.Random(seed)
    first, second = (tuple(effects[t]) for t in TASKS)
    values = []
    for _ in range(replicates):
        mean_first = sum(first[rng.randrange(8)] for _ in range(8)) / 8
        mean_second = sum(second[rng.randrange(8)] for _ in range(8)) / 8
        values.append((mean_first + mean_second) / 2)
    values.sort()
    return {'method': 'Percentile sensitivity bootstrap of paired optimizer-run effects within each fixed task, equal task weights',
            'replicates': replicates, 'seed': seed,
            'interval_95_pp': [100 * percentile(values, .025), 100 * percentile(values, .975)],
            'task_equal_mean_effect_pp': 50 * (statistics.fmean(first) + statistics.fmean(second)),
            'resampling_unit': 'complete paired optimizer run; never raw positions, structures or checkpoints'}


def infer_from_runs(runs, *, bootstrap=False, replicates=BOOTSTRAP_REPLICATES, seed=BOOTSTRAP_SEED):
    grouped = {task: [r for r in runs if r['task'] == task] for task in TASKS}
    require(all(len(grouped[t]) == 8 for t in TASKS) and len(runs) == 16, 'All sixteen independent runs required')
    lower = statistics.fmean(statistics.fmean(r['effect_bounds'][0] for r in grouped[t]) for t in TASKS)
    upper = statistics.fmean(statistics.fmean(r['effect_bounds'][1] for r in grouped[t]) for t in TASKS)
    unknown = sum(r['unknown_pairs'] for r in runs)
    result = {'task_equal_identification_bounds_pp': [100 * lower, 100 * upper],
              'technical_unknown_pairs': unknown, 'welch': None, 'bootstrap_sensitivity': None,
              'task_effect_bounds_pp': {t: [100 * statistics.fmean(r['effect_bounds'][i] for r in grouped[t]) for i in (0, 1)] for t in TASKS}}
    if unknown:
        result['interval_status'] = 'technical_unknowns_identification_bounds_only'
        return result
    effects = {t: [r['effect_bounds'][0] for r in grouped[t]] for t in TASKS}
    result['welch'] = welch_task_equal(effects)
    result['interval_status'] = result['welch']['interval_status']
    if bootstrap and result['welch']['degrees_of_freedom'] is not None:
        result['bootstrap_sensitivity'] = stratified_run_bootstrap(effects, replicates=replicates, seed=seed)
    return result
