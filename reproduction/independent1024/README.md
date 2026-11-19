# Private saved-data analysis: independent 16-run five-arm qualification study

This separate package reconstructs the completed 1,024-position study using
Python's standard library. It contains all 16 new optimizer runs (eight per
fixed task), the same 64 raw positions per run and all five policy rows per
position. Earlier development runs remain excluded from the primary inference.
The fixed-task recipe was formally adopted on 6 October 2026. Its primary is
fixed-task B minus source-preserving B for predicted external qualification
and an audited route. This is computational evidence; the package makes no
claim of biological confirmation or a globally best B policy.

Run from any working directory with Python 3.10 or later:

```sh
python -I -S -B /path/to/analysis_reproduction_independent1024/reproduce_independent1024.py \
  --data /path/to/analysis_reproduction_independent1024/data \
  --reference /path/to/analysis_reproduction_independent1024/data/reference \
  --output /path/to/chosen/results

python -I -S -B /path/to/analysis_reproduction_independent1024/test_saved_inputs.py -v
```

`-I -S` disables environment/user-site imports and site packages. The entry
point explicitly loads only the delivered project modules beside itself. It
never opens original absolute provenance paths. The output may be outside the
package; no A6000 layout, GPU, internet access or model installation is needed.

| Input | Scope and authoritative origin |
| --- | --- |
| `data/design.json` | The original declared 16-run five-arm design from terminal `external_fixed_outputs_v1/fixed_outputs.json` |
| `data/sources_with_parents.json` | All 1,024 unchanged source objects, including raw positions and original identity/path metadata |
| `data/policy_rows.json` | All 5,120 fixed policy rows; documented field projection, with saved attempt status used for diagnostics |
| `data/predictions.json` | Full ten-by-two checkpoint vectors and provenance for 1,118 structures; means and checkpoint dispersion are recomputed |
| `data/profiles.json` and `data/structure_descriptors.json` | Already independently audited saved QED, molecular weight, alert/scaffold profiles and source/final similarities |
| `data/candidate_stage.json` | The unchanged source objects and saved arm/native-status/planned/scored candidates consumed by the original mechanism analysis |
| `data/reference/analysis.json` | Original full analysis, first read after reconstruction, used only for comparison |
| `data/data_provenance.json` | Every original data/code path and retained-field projection; original paths are descriptive strings |
| `data/completed_audits/` | The four actual completed independent audit reports and aggregate receipt; included as retained evidence |
| `data/acquisition_metadata.json` | Original execution/reuse/runtime metadata, preserved separately and not presented as newly measured costs |

The original project `evaluate_independent.py`, `independent_statistics.py`
and helper modules are copied unchanged. The entry point uses the new study's
`analyze` and 16-run statistics; it never invokes the old eight-run analysis,
new scientific `execute`, predictor acquisition or chemistry classes. The
sole primary uses the prespecified 100,000 within-task paired-run bootstrap,
seed 2026101801, and the already independently checked Welch calculation.
Ten checkpoint predictions, 1,024 positions, molecule duplicates and scaffold
unions do not become independent experimental runs.

`results/analysis.json` reconstructs the primary, every prespecified secondary,
all cohort/task summaries, RF/external/QED/mass changes, output concentration,
alert sensitivities, raw-input duplicates and candidate-pool mechanisms. It
reconstructs the external endpoint from the unrounded intended-task mean at
0.5 and the saved audited-route status; a guiding-RF conjunct is absent. Raw
positions and any technical-unknown bounds remain in place. Successful unique
structure/scaffold counts are descriptive unions, not additional run-level
probability estimates. `policy_yields.csv` is a convenience display.

`results/comparison.json` retains every nonzero floating difference. Counts,
memberships and nonnumeric fields are exact; floating differences above
1e-10 are reported as failures. `execution_receipt.json` records the actual
Python/caller/output paths, imported scientific modules and zero scientific
model/reaction work. The actual Python 3.14.4 run reconstructed **13,757 scalar
fields**, with zero comparison failures and 16 final-bit floating differences
(maximum 1.28×10⁻¹⁵). The eight malformed-input tests protect raw membership,
RF qualification, task-column identity, full fold vectors and mechanism source
identity. Relocated private extraction is recorded separately in the parent
folder; retained results are inspection examples, not a substitute for running.

This is saved-data analysis reproduction. It consumes saved native structures,
scaffolds, profiles, similarities and route outcomes; it does not recompute
canonical identities, descriptors, template chemistry or deserialize/audit
native route graphs. Those independent scientific checks already completed
in the retained audit reports. No RF, Chemprop, Graph2Edits, reaction templates,
full stock, generator corpus or third-party weights are included or executed.
Public hosting, redistribution permission and the final archive citation remain
separate release decisions. Nothing here authorizes a public deposit.
