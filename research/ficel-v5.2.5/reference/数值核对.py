"""V5.2.5 的有限数学/实现核对；不是科研性能实验。

仅处理本文件内的手工四格例子，不读取科研数据，不接服务器，不改旧程序。
参考集合语义：同一次观察以已有 observation_ref 为引用；重读同一引用不加证据。
这里的引用去重只防止重复读取，并不把共享控制转换成独立多动作似然。
"""

from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import json

MAX_TERMS = 4096
MAX_BITS = 65536
MAX_OPS = 2_000_000
BISECTION_STEPS = 24


class NumericalUnresolved(Exception):
    pass


class Work:
    def __init__(self, max_ops=MAX_OPS, max_terms=MAX_TERMS, max_bits=MAX_BITS):
        self.ops = 0
        self.max_ops, self.max_terms, self.max_bits = max_ops, max_terms, max_bits

    def check(self, value):
        self.ops += 1
        if self.ops > self.max_ops:
            raise NumericalUnresolved("rational_operation_limit")
        value = F(value)
        if max(abs(value.numerator).bit_length(), value.denominator.bit_length()) > self.max_bits:
            raise NumericalUnresolved("integer_bit_limit")
        return value

    def add(self, a, b):
        return self.check(a + b)

    def mul(self, a, b):
        return self.check(a * b)

    def div(self, a, b):
        return self.check(a / b)

    def power(self, a, exponent):
        answer = F(1)
        for _ in range(exponent):
            answer = self.mul(answer, a)
        return answer

    def fact(self, n):
        answer = F(1)
        for j in range(2, n + 1):
            answer = self.mul(answer, j)
        return answer

    def choose(self, n, k):
        # 正因子逐项形成，不在位数检查前生成一个任意大的组合数。
        k = min(k, n - k)
        answer = F(1)
        for j in range(1, k + 1):
            answer = self.div(self.mul(answer, n - k + j), j)
        return answer


def canonical_observations(records):
    """普通引用集合；同引用内容冲突是输入错误，不被静默合并。"""
    unique = {}
    for r in records:
        ref, source, counts = r
        counts = tuple(counts)
        if len(counts) != 4 or any(type(n) is not int or n < 0 for n in counts):
            raise ValueError("invalid four-cell counts")
        payload = (source, counts)
        if ref in unique and unique[ref] != payload:
            raise ValueError("same observation reference has incompatible content")
        unique[ref] = payload
    return tuple((ref, *unique[ref]) for ref in sorted(unique))


def source_counts(records):
    by_source = {}
    for _, source, counts in canonical_observations(records):
        if source not in by_source:
            by_source[source] = [0] * 4
        for k, n in enumerate(counts):
            by_source[source][k] += n
    return tuple(tuple(by_source[s]) for s in sorted(by_source))


def posterior(records, kappa=F(1), work=None):
    """固定 α=(1,1,1,1) 的原 Dirichlet 层级模型的精确有限混合。"""
    w = work or Work()
    kappa = F(kappa)
    if kappa <= 0:
        raise ValueError("kappa must be positive")
    counts = source_counts(records)
    if all(sum(n) == 1 for n in counts):
        # K_i=1 的真正共轭特例。κ 完全抵消，不能称被数据识别。
        totals = tuple(sum(n[k] for n in counts) for k in range(4))
        alpha = tuple(1 + n for n in totals)
        z = w.check(6)
        for n in totals:
            z = w.mul(z, w.fact(n))
        z = w.div(z, w.fact(sum(totals) + 3))
        return {alpha: F(1)}, z, w

    polynomials = []
    for k in range(4):
        poly = [F(1)]
        for n in counts:
            for j in range(n[k]):
                new = [F(0)] * (len(poly) + 1)
                for exponent, coefficient in enumerate(poly):
                    new[exponent] = w.add(new[exponent], w.mul(coefficient, j))
                    new[exponent + 1] = w.add(new[exponent + 1], w.mul(coefficient, kappa))
                poly = new
        polynomials.append([(e, v) for e, v in enumerate(poly) if v])

    term_count = 1
    for p in polynomials:
        term_count *= len(p)
        if term_count > w.max_terms:
            raise NumericalUnresolved("mixture_term_limit")

    terms = {(): F(1)}
    for p in polynomials:
        terms = {old + (e,): w.mul(c, v) for old, c in terms.items() for e, v in p}
    unnormalised = {}
    normaliser = F(0)
    for exponents, coefficient in terms.items():
        value = w.mul(coefficient, 6)
        for exponent in exponents:
            value = w.mul(value, w.fact(exponent))
        value = w.div(value, w.fact(sum(exponents) + 3))
        alpha = tuple(1 + e for e in exponents)
        unnormalised[alpha] = value
        normaliser = w.add(normaliser, value)
    mixture = {a: w.div(v, normaliser) for a, v in unnormalised.items()}

    # 模型证据包括多项式似然的 θ 无关常数；不能只比较归一化混合权重。
    prefactor = F(1)
    for n in counts:
        total = sum(n)
        mult = w.fact(total)
        for cell in n:
            mult = w.div(mult, w.fact(cell))
        rising = F(1)
        for j in range(total):
            rising = w.mul(rising, kappa + j)
        prefactor = w.mul(prefactor, w.div(mult, rising))
    return mixture, w.mul(normaliser, prefactor), w


def mean_delta(mixture, work):
    answer = F(0)
    for a, weight in mixture.items():
        answer = work.add(answer, work.mul(weight, F(a[1] - a[2], sum(a))))
    return answer


def beta_survival_integer(t, a, b, work):
    """1-I_t(a,b)，整数参数时的非负有限二项式和。"""
    if t <= 0:
        return F(1)
    if t >= 1:
        return F(0)
    n = a + b - 1
    answer = F(0)
    for h in range(a):
        term = work.mul(work.choose(n, h), work.power(t, h))
        term = work.mul(term, work.power(1 - t, n - h))
        answer = work.add(answer, term)
    return answer


def positive_delta_survival(t, a, b, c, work):
    answer = F(0)
    for j in range(a):
        coefficient = work.div(work.choose(b + j - 1, j), work.power(F(2), b + j))
        term = work.mul(coefficient, beta_survival_integer(t, a - j, b + j + c, work))
        answer = work.add(answer, term)
    return answer


def delta_cdf(t, mixture, work):
    t = F(t)
    if t <= -1:
        return F(0)
    if t >= 1:
        return F(1)
    answer = F(0)
    for alpha, weight in mixture.items():
        a, b, c = alpha[1], alpha[2], alpha[0] + alpha[3]
        if t >= 0:
            value = work.add(1, -positive_delta_survival(t, a, b, c, work))
        else:
            value = positive_delta_survival(-t, b, a, c, work)
        answer = work.add(answer, work.mul(weight, value))
    return answer


def lcb_interval(mixture, tail, work):
    tail = F(tail)
    if not 0 < tail < 1:
        raise ValueError("tail must lie inside (0,1)")
    lo, hi = F(-1), F(1)
    for _ in range(BISECTION_STEPS):
        middle = (lo + hi) / 2
        cdf = delta_cdf(middle, mixture, work)
        if cdf == tail:
            return middle, middle
        if cdf < tail:
            lo = middle
        else:
            hi = middle
    return lo, hi


def exact_score_sign(mixture, tail, declared_relative_cost, work):
    """连续严格递增 CDF 下 LCB-cost 的精确符号；费用为声明有理数。"""
    x = delta_cdf(declared_relative_cost, mixture, work)
    return "POSITIVE" if x < tail else "NEGATIVE" if x > tail else "ZERO"


def interval_gate(intervals):
    """仅用于核对原有 0 回退和数值未决分支；不实现科研 Gate。"""
    if all(hi <= 0 for _, lo, hi in intervals):
        return "A0"
    for name, lo, hi in intervals:
        if lo > 0 and all(name == other or lo > other_hi for other, _, other_hi in intervals):
            return "APPLY:" + name
    return "NUM_UNRESOLVED"


def fraction_text(x):
    return str(F(x))


def main():
    tests = []
    def check(name, passed, details=None):
        tests.append({"name": name, "passed": bool(passed), "details": details})
        if not passed:
            raise AssertionError(name)

    one = [("one-source-block", "source-A", (0, 2, 0, 0))]
    mix, evidence, work = posterior(one)
    check("hierarchical_exact_mixture", mix == {(1, 2, 1, 1): F(5, 7), (1, 3, 1, 1): F(2, 7)},
          {str(k): str(v) for k, v in mix.items()})
    check("hierarchical_mean_delta_5_over_21", mean_delta(mix, work) == F(5, 21))
    check("hierarchical_evidence_7_over_40", evidence == F(7, 40))

    two = [("source-A-block", "source-A", (0, 1, 0, 0)),
           ("source-B-block", "source-B", (0, 1, 0, 0))]
    results = [posterior(two, k) for k in (F(1, 4), F(1), F(4))]
    check("Ki1_kappa_cancels_exactly", all(r[0] == results[0][0] and r[1] == results[0][1] for r in results),
          {"posterior": {str(a): str(v) for a, v in results[0][0].items()}, "evidence": str(results[0][1])})
    check("independent_sources_are_not_within_source_repeats", mean_delta(results[0][0], results[0][2]) == F(1, 3))

    prior, _, prior_work = posterior([])
    grid = [F(-1), F(-3, 4), F(-1, 2), F(-1, 4), F(0), F(1, 4), F(1, 2), F(3, 4), F(1)]
    values = [delta_cdf(t, prior, prior_work) for t in grid]
    check("cdf_boundaries_and_known_value", values[0] == 0 and values[-1] == 1 and values[4] == F(1, 2)
          and delta_cdf(F(1, 2), prior, prior_work) == F(15, 16),
          {str(t): str(v) for t, v in zip(grid, values)})
    check("cdf_symmetry", all(delta_cdf(-t, prior, prior_work) == 1 - delta_cdf(t, prior, prior_work) for t in grid))
    check("cdf_monotone", all(a <= b for a, b in zip(values, values[1:])))
    reflected = {tuple((a[0], a[2], a[1], a[3])): weight for a, weight in mix.items()}
    check("asymmetric_mixture_reflection", all(delta_cdf(-t, mix, work) == 1 - delta_cdf(t, reflected, work) for t in grid))

    lo, hi = lcb_interval(prior, F(1, 20), prior_work)
    check("lcb_exact_cdf_bracket", delta_cdf(lo, prior, prior_work) <= F(1, 20) <= delta_cdf(hi, prior, prior_work)
          and hi - lo <= F(1, 2 ** 23),
          {"lo": str(lo), "hi": str(hi), "width": str(hi - lo)})
    check("cold_start_negative_lcb", hi < 0 and exact_score_sign(prior, F(1, 20), F(0), prior_work) == "NEGATIVE")
    check("exact_zero_score_and_a0", exact_score_sign(prior, F(1, 16), F(-1, 2), prior_work) == "ZERO"
          and interval_gate([("zero", F(0), F(0))]) == "A0")

    positive_records = [("positive-" + str(i), "positive-source-" + str(i), (0, 1, 0, 0)) for i in range(8)]
    pos, _, pos_work = posterior(positive_records)
    plo, phi = lcb_interval(pos, F(1, 20), pos_work)
    check("positive_lcb_and_separated_selection", plo > 0
          and exact_score_sign(pos, F(1, 20), F(0), pos_work) == "POSITIVE"
          and interval_gate([("positive", plo, phi), ("prior", lo, hi)]) == "APPLY:positive",
          {"lo": str(plo), "hi": str(phi)})
    check("cross_zero_or_overlapping_intervals_unresolved",
          interval_gate([("cross-zero", F(-1, 100), F(1, 100))]) == "NUM_UNRESOLVED"
          and interval_gate([("a", F(1, 10), F(1, 5)), ("b", F(3, 20), F(1, 4))]) == "NUM_UNRESOLVED")

    one_step = []
    for k in range(4):
        counts = tuple(int(j == k) for j in range(4))
        _, extended_z, _ = posterior(one + [("future-" + str(k), "future-source", counts)])
        one_step.append(extended_z / evidence)
    check("evidence_ratio_one_step_predictive_sums_one", sum(one_step) == 1, [str(p) for p in one_step])
    two_step = []
    for n0 in range(3):
        for n1 in range(3 - n0):
            for n2 in range(3 - n0 - n1):
                counts = (n0, n1, n2, 2 - n0 - n1 - n2)
                _, extended_z, _ = posterior(one + [("future-pair", "future-source", counts)])
                two_step.append(extended_z / evidence)
    check("evidence_ratio_two_repeat_count_predictive_sums_one", len(two_step) == 10 and sum(two_step) == 1,
          [str(p) for p in two_step])

    reread = one + one
    dup_mix, dup_evidence, _ = posterior(reread)
    check("rereading_reference_set_does_not_add_information", dup_mix == mix and dup_evidence == evidence
          and len(canonical_observations(reread)) == 1)
    conflict_handled = False
    try:
        canonical_observations(one + [("one-source-block", "source-A", (0, 1, 0, 0))])
    except ValueError:
        conflict_handled = True
    check("incompatible_same_reference_not_silently_deduplicated", conflict_handled)

    limits = []
    for kind, instance in [("mixture_term_limit", Work(max_terms=1)),
                           ("integer_bit_limit", Work(max_bits=1)),
                           ("rational_operation_limit", Work(max_ops=0))]:
        try:
            posterior(one, work=instance)
        except NumericalUnresolved as exc:
            limits.append(str(exc))
    check("all_resource_limits_return_numerical_unresolved", set(limits) == {
          "mixture_term_limit", "integer_bit_limit", "rational_operation_limit"}, limits)

    report = {
        "kind": "deterministic_mathematical_and_implementation_checks_not_performance_experiment",
        "data": "handwritten_toy_counts_only_no_research_data_no_server",
        "four_cell_order": ["11", "10", "01", "00"],
        "reference_semantics": "set of existing observation_ref values; no duplicated likelihood from rereading; incompatible content raises ordinary input error; shared controls do not become independent joint observations",
        "limits": {"max_nonzero_mixture_terms": MAX_TERMS, "max_integer_bits": MAX_BITS,
                   "max_rational_operations_per_query": MAX_OPS, "lcb_bisection_steps": BISECTION_STEPS},
        "all_passed": all(t["passed"] for t in tests),
        "test_count": len(tests), "tests": tests,
        "scope": "conditional population-mean theta posterior, mean LCB, marginal paired prediction and count-block prediction; not q_new full-effect CDF, not shared-control multi-action joint IG, not effectiveness evidence",
    }
    output = Path(__file__).with_name("数值核对结果.json")
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"all_passed": report["all_passed"], "test_count": len(tests), "output": str(output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
