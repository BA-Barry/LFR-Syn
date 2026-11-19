# Local method packaging — fixed-task qualification adopted (2026-10-06)

The formally adopted candidate qualification is now the default of
`run_method.py execute`: `--qualification fixed_task`. Candidate RF must reach
the fixed task floor (GSK3B 0.10 or JNK3 0.05), and candidate QED must remain at
least source QED minus 0.05. The complete retained scored pool is requalified
by the exact `fixedtask_kernel.py` used in the completed 16-run, 1,024-position
five-arm validation. Source exclusion, ordering, eight scored candidates per
arm, sixteen forward attempts, two candidate queries of at most six calls and
one shared 45-second allowance are unchanged. Outputs are named
`B_fixedtask6plus6` and `C_fixedtask6plus6`, record the qualification mode and
RF change, and make **no source-RF nondegradation guarantee**. Qualification
does not establish measured activity or a globally best B policy.

Use `--qualification source_preserving` for the historical strict rule
`candidate_RF >= max(source_RF, task_floor)` and the retained
`B_joint6plus6`/`C_joint6plus6` names. Only that mode guarantees candidate RF
nondegradation. Existing output configurations without a mode remain strict;
changing the mode requires a separate output directory. Historical bundled
replay examples below explicitly select strict mode, and replay never searches
for an unmeasured fixed-policy candidate. Neither the new CLI integration nor
its injected tests repeat model experiments. The original chemistry, score,
reaction, search, audit and strict policy kernels remain unchanged.

```powershell
python run_method.py execute --help
python -m unittest test_qualification_modes -v
```

Portability update, 4 October 2026: caller-relative inputs, outputs, assets and
explicit interpreter paths are now resolved before worker subprocesses change
directory; bare interpreter names retain normal PATH lookup. The original
failure was reproduced in a temporary copy, and 53 path/worker/reuse tests passed
with Python 3.12.14, without scientific model calls. The separate
`submission_preparation_20261004/software_setup/portable_fresh_quickstart.md`
contains a proposed relocated two-environment installation recipe based on the
observed study versions. On 5 October 2026, the two ordinary isolated virtual
environments were actually installed and checked, and the existing one-source
candidate fixture completed a real no-history execution in those new
environments. See
`submission_preparation_20261004/software_setup/clean_fresh_one_source_software_validation_20261005.md`
and its private receipt archive. This installation/execution evidence is
distinct from the earlier round trips reported below and is not an efficacy
cohort. The RF environment reports PyTDC metadata as null; this does not relax
the historical RF-cache compatibility checks.

This B/C method package is a private preparation artifact. It has not been
published. The previously missing Graph orchestration is implemented. An actual
one-source model round trip and a separate two-source explicit-history round
trip completed on A6000 on 2026-10-03/04. These selected inputs check software
branches; they are not an additional scientific effect cohort.

The fresh branch run made 17 RF predictions, 13 candidate-formation forward
reaction calls (excluding source replay) and four
Graph queries comprising 25 reaction-model calls. Its parent and both B
candidates remained unsolved; C's first candidate solved and its saved graph
passed native-route, stock-leaf and target-exclusion audit. The explicit-history
run covered the retained no-qualified-candidate input and the same candidate
branch: five old Graph measurements were reused, with zero new RF predictions,
forward reaction calls or Graph model calls. Its solved C witness was loaded
from the original physical file and audited afresh. Neither run trained or
sampled a generator. Full outputs are in the project directories
`method_roundtrip_branch_20261003` and
`method_roundtrip_reuse_two_branches_20261003`.

The first explicit-history invocation stopped before scientific execution
because distribution metadata `2023.9.6` was compared literally with RDKit's
runtime version `2023.09.6`. Numeric release comparison now accepts these
equivalent spellings and separately verifies the loaded runtime. A different
release remains incompatible. The original failure record is retained;
seventeen targeted reuse tests passed locally and in the actual RF environment.

## What is operational

`run_method.py replay` recomputes retained B/C policy decisions from real saved
candidate scores and query records using the retained policy evaluator. It
requires only Python's standard library and makes no model calls:

```powershell
python run_method.py replay --qualification source_preserving --output local_replay_first
python run_method.py replay --qualification source_preserving --sources examples/source_candidate_branch.json --output local_replay_candidate
python run_method.py replay --sources examples/fixed_task/source_saved_independent.json --cached-records examples/fixed_task/candidate_records.json --cached-outcomes examples/fixed_task/policy_outcomes.json --output local_replay_fixed_task
```

Both examples were run locally with Python 3.14 on 2026-10-03. Their scientific
endpoint fields matched the saved outcomes exactly (zero differences). This
accepts the original route-audit fields; it does not independently deserialize
or audit the route graphs.

The fixed-task example is the first saved independent-study source with both
fixed-policy replacement branches, `graphga-qualconfirm-gsk3b-rep0-r000`.
It uses complete measured candidates and query outcomes; default fixed-task
replay matched both original saved policy outcomes with zero scientific-field
differences and zero new scientific work. It is a software fixture from the
completed study, adding no efficacy cohort or independent statistical unit.
Its source paths and selection purpose are in
`examples/fixed_task/example_provenance.json`.

The default example is `gsk3b-rep0-r000`, the first raw source, retaining its real
task index 0, replicate 0, position 0 and search seed 10302000. Its saved B/C
policies had no qualified candidate. The explicitly selected branch example
`gsk3b-rep0-r005` is the first source in that run with qualified candidates for
both arms, and retains position 5 and seed 10302005. It was chosen for software
branch coverage, not for a new effect estimate. No extra molecules were sampled.
Default selection is the first one supplied row, including invalids/duplicates.

## Real-execution worker and local validation

`run_method.py execute` now connects RF source scoring, 12-call/45-second parent
searches, the separate RF/chemistry interpreter, and B/C candidate planning.
The command-line worker calls the real retained `ResidentGraph`; there is no
replay substitute on this path. Every reported solved search is saved as a graph
pickle, reloaded, and audited using native route extraction, stock membership,
target identity and target-leaf exclusion. Query errors, audit errors, call
overruns and unobserved interruptions stay unknown rather than becoming failed
scientific outcomes. Resuming the same output never silently repeats a saved
incomplete or failed query.

Local control-flow tests use injected backends and CPU stages and load no
third-party models. Run them from this package directory:

```powershell
python -m unittest lfr_method.test_graph_worker lfr_method.test_exact_study_store -v
```

Thirty Graph worker/reuse tests passed with Python 3.14 on 2026-10-03. They cover the full stage
sequence and interpreter separation, shared remaining time, current-output
query reuse, real-request cost accounting, saved-witness auditing, invalid-slot
retention, parent retention, audit/parent errors, parent call overrun,
interruption/resume, explicit first-path history selection and source identity.
Exact reuse tests also cover cap12/cap6 incompatibility, exact second-branch
time, unknown/error/incomplete retention, cross-directory declared origins
without directory scans, alias equality and new saved-witness re-auditing. Both
included software examples use their real source/candidate identities and saved
logical costs in injected tests: five saved Graph measurements are reused with
zero new model calls, and the C solved witness is re-audited by a synthetic
auditor. These tests use synthetic graph objects and therefore do not verify
the actual archived scientific witnesses.

The actual retained search method is also exercised against fake installed APIs
to verify per-request reset and parent/candidate stock semantics. This is
software validation, not GPU execution or evidence of an efficacy increase.

The retained search and saved-graph audit implementation is present in
`lfr_method/graph_core.py`: RetroStarSearch, Graph2Edits, 20 expansions, depth 12,
target InChIKey exclusion, reset per target, saved-witness extraction and stock
checks. `lfr_method/graph_worker.py` supplies parent/candidate scheduling, the
role adapter, stage persistence and cumulative actual-work reporting. Its role
adapter was checked against the actual parent `v2_graph2edits_probe.py` and
formal candidate `joint_validation512_graph.py`: a stock parent is allowed but
its key is removed from inventory; a stock candidate is rejected before search.
Both roles remove the target key and reset the reaction-model cache between
new chemical searches. The parent stage and candidate stage have separate
logical budgets.

## Minimal actual model round trip

### CPU cost reporting

The legacy `actual_new_work.forward_reactions` and candidate-stage
`actual_new_forward_attempts` fields count **candidate formation only**. Their
explicit scope labels now prevent reading them as full chemistry totals. The
new `source_replay_forward_reactions`, `candidate_forward_reactions`, and
`forward_reactions_total` fields separately report source replay, newly formed
candidates, and their sum. `retro_template_calls` reports reverse-template
native calls, and `source_RF_path_attributions` counts actual path
decompositions shared by both arms; it is not an additional RF prediction.

`actual_new_cpu_work` retains the component counters, their observation status,
and the sum of available counters. Counts are cumulative in the output
directory: continuing from a completed candidate stage retains its original
source/arm costs rather than reporting them as new zero-cost measurements. A
formed legacy source with a missing counter has an unknown (`null`) component
and no complete total; `observed_count` records only the available subtotal.
Only a terminal branch that actually skipped formation has known zero root,
retro and attribution work without those fields.

For the actual fresh isolated fixture on 5 October, saved counters give 121
source-recovery `RunReactants` calls plus 13 new candidate calls, hence 134
native forward calls, together with 1,071 retro calls and one shared physical
RF path decomposition. The local bookkeeping repair reproduced these counts
from the existing saved records without rerunning prediction, reactions or
searches. The already executed source package, original results, relocated
snapshot and receipt archives were preserved. This reporting repair has not
been deployed to a scientific worker or used for an additional model run.

```powershell
python -m unittest lfr_method.test_cpu_costs lfr_method.test_graph_worker lfr_method.test_exact_study_store test_cli_paths test_rf_reuse -v
```

On Windows Python 3.12.14, 65 tests completed with 64 passes and one expected
Unix-venv-symlink skip. These cover complete-record continuation, missing
legacy counters, measured cached-zero work, actual counter scopes and the
existing workflow/reuse safety semantics without real model inference.

Copy the package to the indicated private directory on A6000, preserving its
subdirectories. The following Linux command uses the existing project assets
and one existing, already measured source. It creates an independent output
directory, performs no training or generation, and deliberately omits
`--reuse-study` so the Graph searches on this first run are actual new calls:

```bash
cd /home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914/development/method_reproduction_20261003
/home/chu/lfr-envs/pytdc-1.1.15-sklearn-1.2.2/bin/python run_method.py execute \
  --sources examples/source_candidate_branch.json --limit 1 \
  --output /home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914/development/method_roundtrip_branch_20261003 \
  --rf-python /home/chu/lfr-envs/pytdc-1.1.15-sklearn-1.2.2/bin/python \
  --graph-python /home/chu/lfr-envs/syntheseus-formal-py311-0.7.2-graph2edits-0.2.0-torch2.2.2-cu121/bin/python \
  --rf-model-dir /home/chu/lfr-assets/tdc-oracle-interface-20260908/oracle \
  --forward-templates /home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914/development/external_baselines/SynTwins/data/templates.csv \
  --retro-templates /home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914/development/external_baselines/SynTwins/data/retro_templates.csv \
  --stock-smiles /home/chu/trace_lfr_g0_mapping_20260721/raw_external/paroutes_v2_7341155/stock_n1.txt \
  --stock-keys /home/chu/trace_lfr_g0_stock_identity_20260721/runs/n1_aiz_20260721/stock_inchikeys.txt \
  --graph-model-dir /home/chu/trace_lfr_syntheseus_smoke_20260721/model-cache/Graph2Edits_backward \
  --device cuda:0
```

The source `gsk3b-rep0-r005` is already included because it has qualified
candidates in both arms and an archived C replacement route. The selection is
for software branch coverage, not for a positive effect estimate. One source
exercises source RF/QED, reaction-constrained candidate formation, eight logical
score slots, attribution/permutation, both arms' 6+6 scheduling, and saved-route
audit if the measured search solves. Saved source and candidate scientific
fields can be compared with the retained records; timing, physical cache costs
and search outcomes must be reported as actually measured. A saved success is
not copied to manufacture a successful round trip.

Inspect `method_execution.json`, `method_result.json`, both interpreter logs,
`sources_with_parents.json`, `candidate_stage.json`, per-policy outcomes and
`graph_measurements/queries/*.result.json` with their graph witnesses. The result
labels costs as cumulative actual new work in this independent directory.
Zero or nonzero successes are both usable observations; unexplained execution
limits leave the real round trip incomplete. Repeating the same command is a
continuation/cache check, not an independent experiment. Additional raw512
generation or repeating the entire cohort is unnecessary for this software
check. The existing first-raw example can separately test real retention if a
further execution branch is needed.

## Current reuse boundaries

No historical Graph lookup is enabled unless `--reuse-study` is explicitly
supplied. `ExactStudyStore` then reads the selected study's parent and policy
query records at identical molecule, call cap and requested soft time. Six-call
candidate queries and 12-call parent queries are separately applicable; there
is no 12-call prefix converted into a six-call result, including through an
import alias. The exact measured first-branch duration determines the second
branch's remaining allowance, with no approximate matching. Selection is the
lexicographically first physical query path, including unsolved, error and
incomplete states, and its seed is adopted. A later solved or completed record
does not replace an earlier unknown. Current-output continuation likewise
keeps the first exact request and does not silently repeat its search.

Declared import aliases are folded back only when their request, seed, spent
work, status and audit evidence are equivalent. The root task authorized
following the explicit study's precise physical/history links into the existing
project's `BASE/development` or `/home/chu/trace_lfr_*` resources. Only those
declared files and corresponding configurations are read; their directories
are never scanned for alternative or more favourable queries. Every imported
solved query is reloaded from the physical query's original graph pickle and
audited in the current package. `SavedGraphAuditor` loads stock/RDKit and uses
the retained native audit, without constructing Graph2Edits or loading its
weights. An incompatible, incomplete or failed selected witness stays unknown.
Original search work, zero new inference calls, and new graph-read/audit/setup
costs are reported separately. Unobserved expenditure is labelled incomplete,
with observed totals treated as lower bounds.

To exercise both existing software examples through explicit reuse, use the
full `execute` command above with a new output directory and add:

```bash
--reuse-study /home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914/development/joint_validation512_v1
```

Use `examples/source_candidate_branch.json` for one output and
`examples/source_first_raw.json` for another. This preserves the software
examples' identities and allows their existing measurements to be reused.
Missing exact observations are actually measured, not inferred from replay.
The independent no-reuse command above remains the small actual-search check.

CPU reuse is supported specifically for the original `joint_validation512_v1`
layout, original RF/template/stock asset paths, RDKit 2023.09.6 and scikit-learn
1.2.2. `rf_reuse.py` reads 456 internal ledger files: 17 source/training/new-score
ledgers and 439 completed source records. It uses all observed planned candidate
scores, including low scores, and does not follow external `score_provenance`
paths or directories named elsewhere in those ledgers. Source scoring is
missing-only. Candidate preparation reuses saved root discoveries, compatible
forward reactions and all available RF observations while recomputing
attribution, ranking and enrollment. Logical score/forward budgets remain the
same when physical measurements are cached. Seventeen dedicated RF reuse tests
passed locally and in the actual retained RF environment. The explicit-history
run then measured zero new RF/forward/Graph calls. All 16 candidate RF values,
planned identities and qualified outputs in its selected candidate branch
matched the deliberately fresh branch run:

```powershell
python -m unittest test_rf_reuse -v
```

The old ledgers do not record general SMARTS/model content versions. Matching
the original paths and environments cannot detect assets replaced at those
same paths, and this package does not claim content identity verification or
support relocated/arbitrary studies. A deliberately fresh bounded source round
trip and the reuse checks are reproduction checks, not new efficacy cohorts.
Retained historical data are read only. End-to-end reproducibility is assessed
after the actual model and witness checks, without repeating raw512 generation
or the entire scientific cohort.

## Source and scientific semantics

`code_origins.json` records original function names, line locations and explicit
packaging changes. `sources/` contains retained project scripts acquired through
read-only access to their A6000 paths on 2026-10-03. These are retained current
project sources, not newly asserted frozen execution snapshots. The full new
packaging version is 0.1.0 and is distinct from the original 512 controller.

The core uses the original root ordering, two one-missing-role and two
two-missing-role roots with unused-capacity filling; forward replay must recover
the source. B/C permutation groups and the New512 seed formula retain the
source's actual task/replicate/position. Four precursor ranks at most four roots
give at most 16 attempts; the first eight unique QED-qualified products consume
logical scoring slots, including cached scores. Fixed-task qualification is RF
at least the task cutoff, with QED at least source QED−0.05. Explicit historical
source-preserving mode uses the greater of source score and task cutoff.
Candidate planning is six calls then six calls with one shared 45-second soft
allowance; the second search gets remaining time and cannot inherit unused
calls. Parent 12-call/45-second work is separate. Search failures retain the
source; execution/audit errors remain scientifically unknown.

Only project source and existing project-generated example records are copied.
No third-party model weights, template tables or stock files are bundled. Real
execution is intended to use explicit paths for the RF model directory,
forward and reverse template CSVs, SMILES stock, InChIKey stock and Graph2Edits
model directory, plus separate RF and Graph interpreters. Their provenance,
permissions and exact training composition are not resolved by this package.

The original science used RF/RDKit in the PyTDC 1.1.15, scikit-learn 1.2.2,
RDKit 2023.9.6 environment, and planning in Syntheseus 0.7.2, Graph2Edits 0.2.0,
Torch 2.2.2+cu121, RDKit 2025.9.2. This local replay needs none of these packages.
