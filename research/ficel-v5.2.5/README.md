# FICEL-Syn V5.2.5

Last research revision: 7 November 2026. Public deposit: 19 November 2026.

- [Complete editable six-algorithm scheme](FICEL-Syn_V5.2.5_六算法完整方案_计算时序与接口补充.docx)
- [Full readable text](修订全文.txt): text extraction linearizes some native Word equations; use Word for full notation.
- [Supplement summary](本次补充说明.md)
- [Existing implementation audit](旧实现接口核对.md)
- [Reference numerical checks](reference/数值核对.py) and [recorded result](reference/数值核对结果.json)

Run from the repository root using Python 3.10+:

```sh
python research/ficel-v5.2.5/reference/数值核对.py
```

This writes the numerical result JSON beside the script. The 19 checks use handcrafted counts only; they do not measure molecular efficacy. The numerical reference preserves within-source repeats, the hierarchical four-cell posterior, mean-LCB interpretation and bounded computation. It is not a completed six-algorithm runtime or a large-data inference guarantee.

The finite instance fixes the cross-source state during one development batch, updates CEMA at the batch boundary, uses one declared statistical layer for the course, then allows one finite G0→G1 TL phase. Negative/zero evidence changes CEMA and selection; TL receives actual course-distribution changes, not an automatic extra failure gradient. Candidate-formation and chemical identity interfaces still need the specific connections documented in the audit.

Historical evidence is retained without relabeling it as new FICEL efficacy. No server or new chemistry/model experiment was invoked by the theory revision or this public deposit.
