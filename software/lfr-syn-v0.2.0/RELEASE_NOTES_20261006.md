# Qualification integration, version 0.2.0

After formal adoption on 6 October, the command-line default is
`--qualification fixed_task`: GSK3B RF >=0.10 or JNK3 RF >=0.05,
candidate QED >=source QED minus0.05, and source structure exclusion.
The exact project `fixedtask_kernel.py` used by the completed independent
16-run study is included unchanged. It requalifies the full saved/scored pool
without interpreting the historical strict `task_pass` flag as a veto.

The original strict policy evaluator, candidate score/chemistry kernels,
retrieval/root ranking, source attribution, reaction settings, inventory and
native route audit implementation are unchanged. Both modes retain eight
scored candidates per arm, sixteen forward attempts, score/QED/SMILES ordering,
two candidate branches of at most six model calls, one shared 45-second logical
allowance and the same audited-route requirement. This change integrates the
already validated recipe; it adds no algorithm branch or efficacy experiment.

`--qualification source_preserving` retains the historical
`max(source_RF, task_floor)` rule and `B_joint6plus6`/`C_joint6plus6` names.
Missing historical configuration mode remains strict. Explicit fixed mode uses
`B_fixedtask6plus6`/`C_fixedtask6plus6`. Mode is part of the ordinary output
configuration, so another mode needs a separate output directory. RF change
is reported for known final outcomes; only strict mode promises guiding-RF
nondegradation. Neither fixed policy is presented as globally best, and neither
score is measured biological activity.

The bundled historical replay examples use explicit strict mode. A separate
fixed-mode fixture uses the first measured source with both completed fixed
replacement branches in the independent study. Its default fixed-mode replay
matches both original saved outcomes exactly, adding zero scientific work and
zero independent efficacy units. Unmeasured candidates are never searched by
the replay command.

The final relevant ordinary suite ran 74 tests: 73 passed, with one expected
Windows skip for Unix interpreter-symlink semantics; zero failures/errors.
Eight qualification tests and the new default fixed-mode real-record replay
test cover the changed option, while the retained worker/reuse/path/CPU-cost
tests cover its integration. Backends, predictor and chemistry doubles are
used for execution tests. No model weights or native reactions were executed.
The fresh native installation/execution evidence from 5 October remains the
original strict-mode software proof; it has not been relabeled as new fixed
execution or repeated.

Receipts and the private relocated source archive are under
`submission_preparation_20261006/software_validation_20261006/`.
The independent study's complete saved-data analysis has a separate package,
`submission_preparation_20261006/analysis_reproduction_independent1024/`.
No public upload, project-license grant or third-party weight redistribution
is performed by this release preparation.
