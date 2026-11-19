# LFR-Syn

Research manuscripts, method code and saved-study reproduction for reaction-constrained molecular replacement. Earlier completed writing is the starting point of this repository; later research will be recorded as an explicit continuation.

## Read the completed writing

| Version | Main text | Supplement | Scope |
| --- | --- | --- | --- |
| 2 October 2026 | [PDF](manuscripts/2026-10-02/main.pdf) · [Word](manuscripts/2026-10-02/main.docx) · [source](manuscripts/2026-10-02/main.md) | [source](manuscripts/2026-10-02/supplement.md) | Task-guided reaction-constrained replacement on fixed generator outputs. The dated local file contains its later retained corrections; this is a historical working manuscript. |
| 6 October 2026 | [PDF](manuscripts/2026-10-06/main.pdf) · [Word](manuscripts/2026-10-06/main.docx) · [source](manuscripts/2026-10-06/main.md) | [PDF](manuscripts/2026-10-06/supplement.pdf) · [Word](manuscripts/2026-10-06/supplement.docx) · [source](manuscripts/2026-10-06/supplement.md) | Unified English working manuscript, independent 16-run / 1,024-position five-arm study, and tradeoffs. |

The manuscripts are retained research drafts, not journal publications or a claim of submission readiness. Author and affiliation information was deferred. Historical decisions describe their own dates and do not automatically define the latest research claim.

## Existing evidence and executable material

- [Historical method v0.2.0](software/lfr-syn-v0.2.0/README.md): source-preserving and fixed-task policies, real-execution adapters, tests and saved replay examples. New model execution needs external assets and the documented environment.
- [Portable independent-study analysis](reproduction/independent1024/README.md): supplied saved rows, checkpoint predictions, descriptors and retained audit reports. Standard-library-only reconstruction does not rerun chemistry or native-route audits.
- [Data and implementation scope](DATA_AND_CODE.md).

The 6 October study reported external predicted-task qualification plus an audited route for 883/1,024 fixed-task B outputs versus 202/1,024 source-preserving B outputs. Fixed-task B/C was 883/880, which establishes neither universal B superiority nor equivalence. These are computational endpoints under declared predictors, stocks, planners and budgets; they do not establish measured activity or laboratory synthesis. See the manuscript for run-level inference, losses and diversity.

## Quick checks without models or GPUs

Python 3.10+:

```sh
cd software/lfr-syn-v0.2.0
python run_method.py replay --sources examples/fixed_task/source_saved_independent.json --cached-records examples/fixed_task/candidate_records.json --cached-outcomes examples/fixed_task/policy_outcomes.json --output ../../validation-output/fixed-task-replay
python -m unittest discover -v
```

From the repository root:

```sh
python -I -S -B reproduction/independent1024/reproduce_independent1024.py --data reproduction/independent1024/data --reference reproduction/independent1024/data/reference --output validation-output/independent1024
python -I -S -B reproduction/independent1024/test_saved_inputs.py -v
```

The complete saved analysis includes 100,000 bootstrap replicates and is CPU-only. It reconstructs saved-record statistics; it creates no new efficacy observations.

## Public-deposit notes

First public deposit: 19 November 2026, authorized by the project owner. Text inside dated source documents saying “private” or “not yet uploaded” describes the earlier local state. Original scientific values and negative results are retained. Connection credentials, machine-login records, model weights, complete stock/template distributions and third-party corpora are excluded. Some original Word/PDF documents and provenance strings retain historical machine paths; those are not required for the portable saved analysis.

No project-wide reuse license has been selected. Public visibility is not a new MIT/Apache license, and no additional rights over external assets are asserted.
