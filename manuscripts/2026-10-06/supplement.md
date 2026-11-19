# Supplementary methods and reproduction notes for LFR-Syn

This Supplement specifies the adopted fixed-task qualification recipe and preserves the complete evidence used to develop and test it. Old256, New512, HIV384, the earlier eight-run GraphGA512 comparison and the prospective sixteen-run GraphGA1024 qualification validation are separate populations. Historical source-relative qualification is identified explicitly; its previously selected outputs and primary contrasts have not been rewritten as fixed-qualification experiments. Negative Mix, null HIV LR, concentration and chemical-alert sensitivities, unsuccessful allocation/precursor revisions and eight historical native execution limitations remain reported. S23 defines the prospective qualification test; S24 records development alternatives, tradeoffs, audits and complete cost scopes. Computational route and predictor results do not establish biological activity or laboratory synthesis.

## S1. Evaluation populations and the original two evidence batches

The scientific question is whether bounded, reaction-constrained replacement can increase the proportion of original generator output positions whose final structure satisfies both a fixed task predictor and an audited computational route criterion. Every position has exactly one final output. Successful replacement changes that structure; unsuccessful intervention preserves its source. The evaluation denominator includes invalid structures, duplicate chemistry, stock members and low task scores. Generating additional samples to fill failed positions is outside the recorded design.

Table S1 distinguishes the completed development batch from the new validation. Their results must not be pooled. Existing chemistry measurements can be reused with provenance, but reuse does not turn an old measurement or a different search seed into a new independent training replicate.

**Table S1. Original development and validation populations.** These historical populations remain separate from the subsequent GraphGA qualification validation in S23.

| Feature | Old256: `independent_training256_v1` | New512: `joint_validation512_v1` |
|---|---|---|
| Source training runs | Four: two per fixed task | Eight: four per fixed task |
| Terminal raw draws | 64 per run; 256 total | 64 per run; 512 total |
| Main replacement trigger | Valid source with unresolved initial route | Valid source with failed initial joint endpoint; the route-only policy is a specific comparator |
| Main candidate-route allocation | First qualified candidate, up to 12 calls / 45 soft seconds | Seven recorded policies, including first-candidate 12 calls and two-candidate 6+6 calls with **one shared 45-second soft allowance** |
| Main mechanism/composition contrast | Contribution versus permuted contribution | Mix versus contribution under the same 6+6 allocation |
| Current evidence status | Completed generation, planning, quality/external checks and historical-pool diagnostics | Eight source runs, seven route policies, external/structure and saved-output sensitivity analyses complete; eight native endpoints remain bounded unknowns |
| Relation to old budget diagnostic | Diagnostic used saved first-candidate prefixes and measured missing second candidates, with a separate 45-second limit per candidate | Fresh 6-call queries and the actual remaining shared time; no promotion of old trajectory prefixes to new prospective evidence |

The Old256 saved-record evidence has 256 positions, 247 valid structures and nine invalid structures. Its completed joint counts are 21 for original outputs, 93 for Morgan, 126 for contribution and 97 for permuted contribution; the operational SynTwins-compatible count is 70 with a four-slot execution-uncertainty upper bound of 74. The external fixed-output threshold-0.5 counts for contribution, Morgan and permutation are 84, 55 and 56. These are Old256 results, retained as development evidence rather than inserted into New512 estimates. Relevant records are listed in S12.

**Implementation evidence:** [Old256 main working manuscript](G:/逆合成结构化反馈/LFR-Syn_Bioinformatics_fixed_slot_draft_20261001.md), [New512 design](G:/逆合成结构化反馈/joint_validation512_design.json), [New512 generation](G:/逆合成结构化反馈/joint_validation512_run.py), [New512 preparation](G:/逆合成结构化反馈/joint_validation512_prepare.py).

The same-source extensions are described in S14–S17 and S19; HIV and the earlier strict-qualification GraphGA study are specified in S18 and S21. The private software evidence in S20 is not an effect cohort. S22 reanalyzes the completed New512 allocation comparison. The sixteen new GraphGA runs in S23 are the inference units for the prospectively specified qualification contrast. The eight earlier GraphGA runs, including their subsequent fixed-qualification development outputs, are excluded from that new primary sample. Exact compatible measurements can be reused without turning their acquisition into a new source run.

## S2. Complete source training, terminal sampling and source identities

All eight New512 source runs start from the same existing REINVENT prior and optimize only their corresponding fixed GSK3B or JNK3 task predictor. REINVENT4 supplies the staged-learning and sampling framework (Loeffler et al., 2024). The configuration uses 128 steps, batch size 16, DAP, sigma 128 and learning rate 0.0001. Prior and initial agent are both `/home/chu/ai4s-tyk2/external/REINVENT4/priors/reinvent.prior`. There is no route or replacement feedback in source training.

The historical trainer retains two ExternalProcess endpoint labels with weights 0.8 and 0.2. Under `--task-only`, `v2_score_task_only.py` returns the same raw task probability under both labels. Their arithmetic mean is therefore the task score, rather than a mixture with a synthesis feedback reward. Training feedback traces record the SMILES and raw oracle values. Nominal steps times batch size is not substituted for observed oracle-row expenditure.

Only `final.chkpt` supplies the 64 terminal draws. Sampling sets `unique_molecules=false` and `randomize_smiles=true`. No checkpoint ranking, successful-run selection, replacement of invalid draws or duplicate removal determines enrollment. Table S2 contains the saved seeds; physical GPU preferences can change when only one card is available, and actual allocation is retained in job metadata.

**Table S2. New512 training, sampling and preferred GPU assignments.** Actual allocation is retained separately from design preferences.

| Task/run | Training seed | Sampling seed | Preferred physical GPU in design | Raw positions |
|---|---:|---:|---:|---|
| GSK3B rep0 | 10102001 | 10103001 | 0 | 0–63 |
| JNK3 rep0 | 10102101 | 10103101 | 1 | 0–63 |
| GSK3B rep1 | 10102002 | 10103002 | 1 | 0–63 |
| JNK3 rep1 | 10102102 | 10103102 | 0 | 0–63 |
| GSK3B rep2 | 10102003 | 10103003 | 0 | 0–63 |
| JNK3 rep2 | 10102103 | 10103103 | 1 | 0–63 |
| GSK3B rep3 | 10102004 | 10103004 | 1 | 0–63 |
| JNK3 rep3 | 10102104 | 10103104 | 0 | 0–63 |

Source preparation parses each raw SMILES with the RF environment's RDKit, requires one connected component and saves isomeric canonical SMILES. Invalid positions have no canonical structure, task score zero and no positive route or joint endpoint. No salt stripping, tautomer enumeration or resampling is introduced here. Identity records include `sample_id`, task, task index, replicate, cohort, raw position, raw/canonical SMILES, validity, score, QED, training/sampling seeds, requested/adopted search seed and the saved parent graph/audit references.

A complete training run is the stochastic source replicate. Its 64 positions share a checkpoint. Duplicate molecules retain separate positions, while their measurements may be shared. Source preparation reports within-run/task duplicates and exact source overlap with saved history. The completed `generation_complete.json` covers all eight runs and 512 positions; partial startup directories are not substituted for these scheduled runs.

## S3. Fixed RF tasks and source-conditioned contribution weights

The fixed PyTDC GSK3B and JNK3 random-forest assets are `gsk3b_current.pkl` and `jnk3_current.pkl` in the existing oracle directory. The task resources identify ExCAPE-DB as their broad activity-data origin (Huang et al., 2021; Sun et al., 2017). Inputs are RDKit radius-2, 2,048-bit Morgan fingerprints, converted to float32 arrays; chirality is not enabled in the fingerprint calls. The score is `predict_proba` for class 1. This explicit implementation is used instead of inferring fingerprint parameters from an ECFP label in general documentation (Rogers and Hahn, 2010).

The task thresholds are 0.10 for GSK3B and 0.05 for JNK3. They define computational qualification and do not represent measured potency. Two candidate qualification modes are retained:

```text
Historical source-relative (strict) mode:
    T_task(z) >= max(T_task(x), task_threshold)
Adopted fixed-task mode:
    T_task(z) >= task_threshold
Both modes:
    QED(z) >= QED(x) - 0.05
```

Old256, New512 and the extensions in S14–S22 used the strict mode unless explicitly identified otherwise. The adopted recipe uses the fixed mode after the prospective comparison in S23. It permits a lower guiding RF score than the source; the score loss is measured rather than presented as nondegradation. Candidate pool formation, source contribution calculation, RF ordering, reaction constraints, stock exclusions and route budgets are unchanged between the paired qualification arms.

Contribution retrieval uses a probability-form adaptation of random-forest decision-path feature contributions (Palczewska et al., 2014). For each tree, the algorithm follows the source's actual branch path. At a split on bit b, it assigns the difference between the visited child's class-1 probability and its parent's class-1 probability to b. It averages these increments across trees. The code checks that the mean root-node probability plus all bit contributions reproduces the mean leaf probability, and that this prediction matches the saved source task score.

This directly decomposes the sklearn forest's mean probability. It does not compute SHAP values, a fragment-masking explanation, or a biological causal effect. The original classification-vote interpretation of the cited method has leaf-purity qualifications; the present implementation uses node probabilities.

Let H(x) be the set of distinct unfolded radius-2 environment identifiers in source x, and let phi_b(x) be the folded bit contribution. For each environment e:

```text
b = e modulo 2048
n_b = number of distinct identifiers in H(x) mapping to b
w_e = max(0, phi_b(x)) / n_b
```

The division is over distinct source environment identifiers, not their occurrence frequencies. It apportions a folded bit's positive contribution across the observed unfolded identifiers and retains collision ambiguity. No atom-mapped pharmacophore interpretation is implied.

For each missing precursor role j, define E_j as the intersection between H(x) and that original precursor's unfolded environments; E_j is empty at preserved stock roles. Contribution retention for a compatible inventory precursor z is the sum of w_e over e in E_j also present in z. Inventory candidates are ordered by retention descending, precursor Morgan Tanimoto descending and canonical SMILES ascending. Morgan retrieval uses the same role matches, with Tanimoto and canonical tie-breaking.

For the permutation control, E_0/E_1 define three disjoint groups: environments in both missing roles, in role 0 only, or in role 1 only. The algorithm shuffles weight values within each group's sorted identifiers, preserving group membership and the weight multiset. New512 seeds are:

```text
10202000 + task_index*1000000 + replicate*10000
         + position*100 + root_index*10 + group_index
```

There is one recorded permutation per group. Some permutations need not alter ranks, particularly for empty groups or equal weights. Such positions remain in the main population.

**Implementation evidence:** `local_contributions`, `fingerprints` and `score_plan` in [common candidate code](G:/逆合成结构化反馈/v2_common_slot_retrieval_candidates.py); `generate_plan` in [Old256 candidate code reused by New512](G:/逆合成结构化反馈/independent_training256_candidates.py); `permutation` in [New512 candidate code](G:/逆合成结构化反馈/joint_validation512_candidates.py).

## S4. Reaction roots, inventory retrieval and score-blind Mix enrollment

The same project inventory and SynTwins forward/reverse templates are used by the matched internal policies. The inventory is derived from the retained PaRoutes stock file (Genheden and Bjerrum, 2022). Its canonical connected SMILES inventory contains 13,432 entries in the native initialization record. Planner stock membership uses standard InChIKey representation; this representation, rather than the source's RF score, determines whether a precursor role is missing.

Root discovery considers reverse-template applications whose corresponding forward reaction has two reactant roles. Both precursor structures must be valid, and applying the forward reaction to the ordered precursor pair must recover the canonical source among its connected products. Roots are deduplicated by reaction reference and ordered precursor pair. The procedure retains up to two distinct-reference roots with one missing-stock role and two with two missing roles; remaining capacity is filled from the ordered extra roots, to at most four. The final selected root order is missing-role count, reaction reference and ordered precursor SMILES.

A source disconnection with zero missing-stock roles is recorded among replay roots but is not silently added to the replacement root pool. This matters for an already route-solved, low-task source: joint-failure intervention does not guarantee that the current missing-role operation can supply a root. A no-root or no-qualified-candidate result is a measured coverage limit of the specified operation, not proof that no chemical improvement exists.

Stocked roles preserve their original precursor. At missing roles, an inventory molecule must match the forward template's reactant substructure query. One missing role receives up to four ranked alternatives. Two missing roles receive two alternatives each and their four combinations. Enumeration uses rank as the outer loop and root as the inner loop. Forward products are sanitized, require one connected component and are inspected in canonical SMILES order. The source, invalid products, duplicate products and products in the planner stock are excluded. The first eight unique QED-qualified products are enrolled before candidate task-score lookup. At most four roots times four combinations are attempted; cache hits consume the logical combination allowance.

Mix operates on precursor-combination proposals before candidate RF values are read:

```text
1. Create the unchanged B-contribution and A-Morgan proposal streams.
2. Alternate B, A, B, A, ...; start with B.
3. Deduplicate (reaction reference, ordered reactants).
   A duplicate advances its stream but makes no extra forward attempt.
4. Apply the proposal and the shared source/stock/validity/QED exclusions.
5. Enroll canonical-unique products, stopping at 8 enrolled products
   or 16 executed logical combination proposals.
6. If one stream ends, consume the other within these same ceilings.
7. Save the complete enrollment plan, then obtain missing RF values.
```

“Score-blind” refers to **candidate** task-score enrollment and alternation. B still uses source predictor contributions, and the source trigger uses its recorded task qualification. Mix does not form two separately scored eight-product pools, switch streams according to candidate RF success, or take a retrospective union of successful outputs.

After enrollment, matched arms within each study apply its explicitly specified RF mode and the common QED tolerance. Eligible candidates are ordered by RF descending, QED descending and canonical SMILES ascending. The route stage takes at most the first two; a one-candidate policy takes only the first. A forward-template relation proposes a candidate, while a separate saved Graph search establishes its reported route qualification; the reported Graph route need not be the candidate's original proposal root.

**Implementation evidence:** [root discovery](G:/逆合成结构化反馈/v2_unified_raw8_candidates.py), [internal plans and RF lookup](G:/逆合成结构化反馈/independent_training256_candidates.py), [Mix proposal streams and enrollment](G:/逆合成结构化反馈/joint_validation512_candidates.py).

## S5. Seven New512 policies and shared planning allowance

For New512, let r(x) denote the original audited route indicator and T(x) the task score. Joint failure means valid(x) and [not r(x) or T(x) below the task threshold]. Invalid sources and initially guiding-RF/route-successful sources are retained. The policy table gives the original strict-qualification identifiers and settings; it does not imply that these seven policies were rerun with the subsequently adopted fixed mode. The five prospective qualification arms are defined separately in S23.

| Identifier | Trigger | Retrieval | Candidate-route calls |
|---|---|---|---|
| `B_route12` | Original route unresolved | Contribution | First candidate: at most 12 |
| `B_joint12` | Original joint failed | Contribution | First candidate: at most 12 |
| `B_joint6plus6` | Original joint failed | Contribution | First at most 6; second at most 6 if needed |
| `Mix_joint6plus6` | Original joint failed | Score-blind B/A proposal alternation | Same 6+6 |
| `C_joint6plus6` | Original joint failed | Permuted contribution | Same 6+6 |
| `A_joint6plus6` | Original joint failed | Morgan | Same 6+6 |
| `SynTwins_joint6plus6` | Original joint failed | Disclosed native depth3/k10 adaptation | Same 6+6 after native enumeration |

Every policy/source has one 45-second soft route-search allowance after the resident Graph model is available. For 6+6, the second receives 45 minus the recorded first-search seconds and is skipped if no time remains. Unspent calls from the first candidate do not transfer to the second. A policy with one qualifying candidate has at most six calls under 6+6. A verified first-candidate solution stops the policy. A search/audit execution limitation yields an unknown policy endpoint rather than continuing until a favorable result is found.

Time is checked between expansions, and an individual model call may overrun. Actual search duration and overrun are recorded; 45 seconds is not a hard process deadline. Query setup, model loading, graph serialization, audit and cache retrieval are recorded separately. A reused query's historical duration and calls still consume the logical policy allowance; an old record without a pure search timer is explicitly labeled as including other overhead. Equal ceilings do not establish equal observed wallclock.

Policy order rotates cyclically with the complete cohort index and is saved in worker configuration. Eight cohorts and seven policies imply a near-balanced cyclic schedule, not an exactly balanced factorial order. The controller's two workers process the same policies across four complete cohorts each; they do not select successful sources.

## S6. Graph2Edits search, target exclusion, auditing and reuse

The planner uses Syntheseus 0.7.2, syntheseus-graph2edits 0.2.0 and the retained Graph2Edits backward checkpoint (Zhong et al., 2023; Maziarz et al., 2025). Multi-step search is Syntheseus `RetroStarSearch` with `ConstantNodeEvaluator(0.0)` and `ReactionModelLogProbCost` (Chen et al., 2020). No learned Retro* molecular-value checkpoint is loaded.

Search requests 20 one-step predictions, sets maximum expansion-depth parameter 12, stops on the first solution, prevents repeated molecules in a tree and resets reaction-model/search cache between genuinely executed targets. Initial parents have at most 12 calls and 45 soft seconds. Each target's standard InChIKey is removed from inventory before search. Replacement candidates already in stock are excluded before scoring and rejected if they nevertheless reach the candidate planner.

A claimed positive requires a saved graph and route audit. The audit checks graph/synthesis-graph validity, route connectivity and acyclicity, canonical target identity, nonempty stock-supported starting molecules, exclusion of the target as a starting molecule and alternating graph-path depth at most 12. That last number is the implemented graph-depth condition, not a declaration of 12 laboratory reaction steps. A graph with a reported solved status but failed audit is an execution-limited endpoint.

The global source-parent plan groups all raw positions by canonical molecule across tasks. It selects the first lexicographic exact-compatible historical terminal parent record, irrespective of solved versus budget-unsolved status, and adopts its recorded seed. If none exists, one query supplies every duplicate source position in that group. Requested parent seeds are `10302000 + task_index*1000000 + replicate*10000 + position`. Sources retain requested and adopted seeds and their audit references. Parent preparation does not complete `sources.json` by converting interrupted searches to ordinary failures.

Candidate-route history reuse requires matching stock/model paths, package versions, top-results/depth settings, target exclusion, reset behavior and the 12-call/45-second configuration, with a corresponding saved graph. Historical first-path selection does not rank results by success. Missing audit can be completed from the preselected saved graph; audit disagreement remains unavailable rather than selecting another favorable historical path.

For current candidate queries, the key is exact canonical molecule, call cap and exact requested soft seconds; the selected stored query contributes its recorded seed. The path rule includes current unsolved, error and incomplete records. A six-call request cannot be supplied by truncating an old twelve-call trajectory, and different remaining second-branch times are different request settings. Search seeds are measurement provenance, not extra training replicates.

In the historical New512 execution, candidate query stores are **worker-local**. The history scanner excludes `policy_graph` directories, so the two GPU workers do not guarantee global cross-worker deduplication of identical new candidate queries. Parent grouping, the single internal RF writer, the native coordinator and external fold caches have broader respective sharing scopes. Actual per-worker totals must be retained rather than described as globally unique Graph work.

**Implementation evidence:** [source parent plan and audit mapping](G:/逆合成结构化反馈/joint_validation512_prepare.py), [policy evaluation, ResidentGraph and QueryStore](G:/逆合成结构化反馈/joint_validation512_graph.py).

## S7. SynTwins-compatible native comparator

The candidate generator uses SynTwins official commit `53d6bbcfc540971460278ad85e31d312b749fa55`, depth 3, k=10, 1,024-bit native fingerprints and `get_twins(..., fast=True, return_as_list=False)` (Chen et al., 2026). Each newly enumerated molecule has a fresh tree and a requested enumeration allowance of 300 seconds. This is a native execution allowance, not a limit on the number of RDKit template calls. Initialization and seed/tree construction are separate from the timed enumeration; observed elapsed duration is retained.

The original checkout is left intact. The copied `syntree_compatible.py` changes exactly two previously disclosed operations:

1. A negative `filter_k` uses the full ranked building-block list instead of Python's negative end slice.
2. Reaction-role access skips role indices outside the available stock KD-tree/building-block lists.

The implementation compares the relevant original source/templates to the recorded Git version, compares retained preprocessing templates to that checkout and compares the preprocessing building-block sequence to the canonical original inventory. These are existing asset-identity checks, not a new experiment. It reuses the retained preprocessing under `official_syntwins_small_comparison_v1/data` and records initialization costs.

Native outputs are traversed in the returned similarity-ranked insertion order. After canonical uniqueness, source/stock exclusions and the same QED tolerance, only the first eight products receive RF lookup. Cache hits consume the eight-score allowance. The common RF qualification/order then chooses up to two for Graph 6+6/45 evaluation.

Requested new seeds are `10402000 + task_index*1000000 + replicate*10000 + position` and are saved before enumeration. Compatible complete historical `.native.json` records can be reused by exact source and demonstrated version/template/stock/depth/k/bit/fast/300-second settings; the first lexicographic file is selected without choosing outcomes. Recorded execution-limited historical observations remain unknown. If compatibility cannot be demonstrated, that directory is explicitly listed as not reused. Old final “best candidate” selections are not imported as the new policy result.

Within the native coordinator, identical source chemistry adopts the already recorded representative seed/output. With `--workers 2`, workers enumerate only; RF lookup and aggregate/per-source writes remain under one coordinator. Per-source plans, scores and completion files support restart. SIGINT/SIGTERM or an observed timeout preserves a terminal unknown and observed template-call/time data. Diagnostic partial outputs are not treated as a completed candidate pool; a hard interruption with incomplete cost observation is not assigned measured zero expenditure.

The completed New512 native candidate stage retains all 512 source rows, with eight execution-limited rows and 172 sources having at least one RF/QED-qualified candidate. The 172 is candidate-qualification coverage, not an audited route or joint-success count. The eight limitations remain scientific unknowns; finishing enumeration does not resolve them into chemical failures. These counts come from the saved [native candidate stage](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/joint_validation512_v1/syntwins_depth3_k10/candidate_stage.json). Its downstream route stage is now complete: known RF-and-audited-route joint successes are 153/512 with a conservative upper bound of 161/512. All eight limitations are recorded NativeTimeout observations; they remain unknown in the final comparison.

This comparator substitutes the project inventory and common downstream RF/Graph policy. It is not a reproduction of the original full inventory, GA/REINVENT multi-property optimization or its oracle budget. Native enumeration is outside the common candidate-route 45 seconds and differs from the internal 16-combination allowance.

**Implementation evidence:** [New512 native stage](G:/逆合成结构化反馈/joint_validation512_native.py), [retained Old256 native driver](G:/逆合成结构化反馈/independent_training256_native.py), [native preprocessing record](G:/逆合成结构化反馈/official_syntwins_small_comparison_v1/initialization.json).

## S8. Endpoints, training-run comparisons and unknown bounds

For a known final output y at raw position i:

```text
route_i = audited_route(y)
joint_i = audited_route(y) AND T_task(y) >= task_threshold
yield = sum(endpoint_i) / 512
distinct joint outputs = unique (task, canonical final SMILES) among known successes
```

Invalid positions have known zero endpoints. Route-only and joint endpoints, task/run effects, source strata and distinct known outputs are reported together. The source strata are defined before policy outcome: invalid; high/low original task score crossed with original route status. They do not select the primary population.

The New512 primary contrast is `Mix_joint6plus6 - B_joint6plus6`. Recorded secondary contrasts are B−C under 6+6, joint-trigger versus route-trigger B12, B6+6 versus B12, and B/Mix versus Morgan/native. For each contrast, the same source positions are paired inside each of eight complete runs. Wins/losses describe the slot-level mechanism but are not independent training replicates.

For complete endpoints, each run contributes its paired joint-count difference divided by 64. The empirical percentile bootstrap resamples four runs with replacement within each of the two fixed tasks, averages within task and weights the two tasks equally. Enumerating all 4^4 ordered resamples per task gives 65,536 combined resamples. The 2.5th and 97.5th percentiles use the saved analysis's linear interpolation. This is a conditional empirical bootstrap, not an exact randomization test, a power guarantee or cross-target population inference. There are only four runs per task, with shared fixed predictor/model/stock and chemical measurement reuse.

Execution-unavailable rows retain `solved=None`, `joint=None` and separate operational fallback fields. Their RF/route lower bounds are the known source fallback endpoints, and upper bounds permit success. For policy p versus q, the per-position difference bounds are [L_p−U_q, U_p−L_q]; summation and division use all 512 positions. If any contrast endpoint is unavailable, the code reports its bounds and unknown-pair count and does not bootstrap a falsely completed endpoint vector. Bounds are identification limits for unobserved execution, not confidence intervals.

An internal report can contain all seven policy labels with `pending_candidate_stage` for intervened native sources. Completed internal contrasts can be summarized with that status clearly retained. The final external procedure requires 3,584 policy-position records and refuses pending candidate-stage entries; genuine native execution limitations remain allowed as unknowns.

The saved [final seven-policy analysis](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/joint_validation512_v1/analysis.json) covers all 512 positions and all eight scheduled runs for all seven policies. No candidate-stage placeholders remain. The six internal policies have complete scientific endpoints; native retains eight execution-unavailable endpoints. The earlier `analysis_internal_policies_20261002.json` remains a phase snapshot, not the final native result. Table S3 reports final audited outcomes, rather than RF/QED candidate-qualification coverage.

**Table S3. Complete New512 fixed-output outcomes.** Bounds concern execution-unavailable native decisions, not sampling uncertainty. The external endpoint uses the released ten-model mean at least 0.5 and the same audited route, as defined in S9.

| Output policy | RF-and-audited-route joint successes / 512 | Audited-route successes / 512 | External-and-audited-route joint successes / 512 |
|---|---:|---:|---:|
| Original source outputs | 52 | 133 | 23 |
| `B_route12` | 254 | 335 | 147 |
| `B_joint12` | 282 | 335 | 155 |
| `B_joint6plus6` | 311 | 361 | 176 |
| `Mix_joint6plus6` | 296 | 347 | 156 |
| `C_joint6plus6` | 249 | 303 | 112 |
| `A_joint6plus6` | 219 | 276 | 89 |
| `SynTwins_joint6plus6` | 153–161 | 218–225 | 58–66 |

The prespecified primary Mix−B contrast is −15 RF joint successes (18 wins, 33 losses; three positive and five negative runs). Its conditional run-bootstrap interval is −6.445 to +0.586 percentage points around −2.930 percentage points. The fixed-output external threshold-0.5 difference is −20 (18 wins, 38 losses; one positive and seven negative runs), with an interval of −6.445 to −0.977 percentage points around −3.906 percentage points. These observations do not support superior joint yield from this Mix construction. The primary designation is retained, and favorable recorded secondary results do not convert the unsupported Mix-superiority hypothesis into a positive primary result.

The recorded component contrasts are directionally consistent across the optimized RF endpoint and the external architecture. Contribution B−permutation C adds 62 RF joint successes and 64 external joint successes; both differences are positive in all eight runs. Joint-trigger B12 versus route-trigger B12 adds 28 RF joint successes (seven positive runs, one tie) and eight external successes (five positive runs, three ties). B6+6 versus B12 adds 29 RF and 21 external successes, both positive in all eight runs. B−Morgan adds 92 RF and 87 external successes, both positive in all eight runs. B−native RF bounds are +150 to +158, and external bounds are +110 to +118; native contrasts have no complete-endpoint bootstrap interval. These are recorded secondary computational comparisons and do not establish measured biological potency or laboratory route feasibility.

The completed [fixed-output RF and historical-pool analysis](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/joint_validation512_v1/fixed_output_sensitivity_20261002.json) changes endpoint thresholds without recomputing triggers, candidate ranks, chosen branches or final structures. At common RF thresholds 0.10, 0.25, 0.50 and 0.75, B6+6 has 251, 127, 48 and nine joint successes; Mix has 234, 111, 36 and four; C has 201, 82, 26 and three; Morgan has 166, 63, 20 and three. These are exploratory sensitivities, not separate confirmatory claims. The original task-specific threshold remains primary.

At the primary task-specific threshold, B has 310 known unique (task,canonical) joint outputs and only one duplicate joint-slot occupancy. All 310 are outside the specified historical solved pool of 117 structures from `reused_raw96_v1` and `syntwins_depth3_raw96_v1`; this establishes additions relative to that recorded pool only. Mix has 295 known unique task-structure joint outputs, 294 absent by task identity and 293 absent by global canonical identity in that pool. These checks execute no new model inference and do not establish novelty against the entire chemical literature or original model-training set.

**Implementation evidence:** [New512 analysis](G:/逆合成结构化反馈/joint_validation512_analysis.py), [external unknown handling](G:/逆合成结构化反馈/joint_validation512_external.py).

## S9. Fixed-output external property and structure sensitivity

The property sensitivity uses the released CombiMOTS two-task classification checkpoints at repository commit `3f644b727982a95b4f87115ccb8d6b4c918a02ba`: `models/gsk3b_jnk3/fold_0` through `fold_9/model_0/model.pt` (Southiratn et al., 2025). The property architecture follows the directed-message-passing framework described by Yang et al. (2019). This study uses the released property models and does not run the CombiMOTS molecular-search algorithm.

Task-column order is exactly `gsk3b_activity, jnk3_activity`. Chemprop 1.5.2 loads models on CPU, with two Torch/BLAS threads and batches of 32 structures. The external procedure scores all fixed policy outputs and also sources for paired source change. Identical canonical chemistry shares both task outputs for each checkpoint. No external score changes an enrolled pool, candidate rank, chosen output or existing Graph/RF endpoint.

Historical raw96, Old256 and budget-diagnostic complete/per-fold caches are merged with same-model/task-order provenance. Ten complete fold vectors are required to recreate a mean; a historical mean alone is insufficient. New512 predicts only missing fold–structure pairs and saves each batch. The historical union currently supplies 710 complete ten-fold structures (701 from the raw96/Old256 union plus nine budget-diagnostic additions); this is cache coverage, not the number of new validation molecules or an outcome. `prepare` reports actual required coverage without loading property models, `run` fills missing vectors and `analyze` requires a complete cache and loads no property models.

The completed New512 external run requires 1,225 unique canonical structures including sources and all seven fixed outputs. Its logical coverage is 12,250 checkpoint–structure pairs; 140 pairs from 14 complete historical structures are reused, and 12,110 missing pairs are newly predicted. All ten checkpoints are loaded, each new fold ledger contains 1,211 structures, and no historical fold prediction is recomputed. The invocation's recorded 62.718 seconds includes its external/descriptor workflow; the enclosing process stage records 70.095 seconds including process startup. These timings are different observations, not additive measurements of property inference alone.

For known endpoints, external joint qualification is audited route AND the relevant mean external probability at least 0.5; 0.3 and 0.7 are exploratory thresholds. It replaces the optimized RF threshold for this sensitivity definition rather than adding a triple RF/external/route requirement. Invalid outputs retain zero in full-denominator means. Fold population SD is reported as model dispersion. Paired policy contrasts follow the eight-run method above; continuous probability and route-weighted probability differences are descriptive per-run/overall quantities with execution bounds, without an additional reported bootstrap interval.

Unknown historical strict-mode replacements can have a worse external score than their known source even when the guiding RF rule enforces nondegradation; fixed-mode replacements can also have lower guiding RF scores. Consequently an unknown row's scientific external probability has bounds [0,1], and external route-joint has bounds [0, route upper]. The observed fallback external score is labeled operational and is not used as a scientific success lower bound. If the external endpoint is known, its value remains a point. Unknown quality-filtered outcomes are likewise not inferred from fallback descriptors.

The checkpoint arguments reference `data/excape+aixfuse_curated.csv`, fold-specific `train_data_size` values of 77,581–77,583 and a ten-fold configuration. In fold order the values are 77581, 77581, 77581, 77581, 77582, 77583, 77583, 77582, 77581 and 77581. The published `GSK3B_JNK3.csv` contains 95,159 rows, but the complete original curated training rows/order have not been established. Exact canonical overlap with the public table is a descriptive provenance check; absence from that table does not prove exclusion from the original training set or its scaffolds. Averaging all checkpoints is not reconstructed out-of-fold prediction. The RF and external models share broad ExCAPE origin. This is cross-architecture sensitivity, not independent biological validation; probability thresholds are uncalibrated and do not establish potency.

Structural sensitivity uses RDKit QED and molecular weight, PAINS_A/B/C and BRENK catalogs, Murcko scaffolds with `includeChirality=False`, source/output Morgan radius-2/2048 Tanimoto, source-relative QED/MW changes and unique (task, canonical) outputs (Bickerton et al., 2012; Baell and Holloway, 2010; Brenk et al., 2008; Bemis and Murcko, 1996). Empty acyclic scaffolds are counted separately. Fixed-output filters include no PAINS, no PAINS/BRENK, MW≤500 and the source QED tolerance. The script reports original RF-joint quality sensitivity and external threshold/quality sensitivity. Alerts are heuristic descriptors, not safety validation. Operational fallback descriptors remain flagged when the policy's scientific endpoint is unknown.

The final B6+6 operational outputs contain 490 unique structures among 491 valid positions, 471 nonempty unique Murcko scaffolds, one empty scaffold and a largest scaffold occupancy of five. Among 259 replacements, the median source Morgan similarity is 0.252, median QED change is +0.110 and median molecular-weight change is −77.749 Da; 258 replacements change Murcko scaffold. These distributions describe the implemented, often scaffold-changing replacement and do not imply local medicinal-chemistry equivalence. Saved Graph audit witnesses for B's 259 replacements contain 221 one-reaction routes and 38 routes with two to five reaction records. Mix has 207 one-reaction routes among 244 replacements; native has 44 among 101. Reaction-record counts are not necessarily the longest route depth. The result describes improved final output identities under this inventory/planner, not newly discovered routes to the unchanged original molecules. The RF joint count remains 306 without PAINS alerts, 214 without PAINS/BRENK alerts and 287 with MW≤500, against the unfiltered 311. At external threshold 0.5, corresponding counts are 173, 132 and 167, against 176. Source-QED tolerance preserves all 311 RF and 176 external successes because it is part of replacement eligibility. No filter is used to reselect the final outputs.

The Mix−B external difference remains negative under no-PAINS (−19), no-PAINS/BRENK (−12), MW≤500 (−17) and source-QED-tolerance (−20) filters. There is one exact public-table overlap among B's final output slots; the corresponding operational count is one for Mix, three for C, two for Morgan and two for native. These counts refer to the released public `GSK3B_JNK3.csv`, whose original curated-training coverage remains unresolved; neither the largely distinct outputs nor the structural filters supply independent experimental activity evidence.

The completed [public-label exclusion sensitivity](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/joint_validation512_v1/public_label_overlap_sensitivity_20261002.json) retains all 512 raw positions and the already selected outputs, but does not count a joint success when the final canonical structure overlaps the available public labeled index. B6+6's only overlapping output, `jnk3-rep1-r029`, is a non-joint source fallback; its RF count remains 311 and external threshold-0.5 count remains 176. Morgan remains 219/89, Mix remains 296/156, and native bounds remain 153–161/58–66. C decreases from 249/112 to 248/111 because its successful replacement at `gsk3b-rep0-r012` overlaps the public table. No model is rerun, no fallback or replacement is reselected, and this exploratory calculation does not establish complete exclusion from the original model-training data: the available public index is not the complete curated checkpoint-training dataset.

**Implementation evidence:** [New512 external/structure script](G:/逆合成结构化反馈/joint_validation512_external.py), [released checkpoint inspection](G:/逆合成结构化反馈/external_property_raw96_v1/model_inspection.json), [download provenance](G:/逆合成结构化反馈/external_property_raw96_v1/download_metadata.json), [Old256 external script](G:/逆合成结构化反馈/external_property_independent256.py), [budget external cache procedure](G:/逆合成结构化反馈/budget_replay256_external.py).

## S10. Dependency and asset provenance

Tables S4–S5 report observed package metadata and retained paths; they do not prescribe replacing assets with newer releases. The environment metadata was read without importing/running property or route models. The RF virtual environment inherits system-site packages, so inherited package versions were resolved through its interpreter metadata rather than inferred from the virtual-environment directory.

**Table S4. Recorded execution environments.**

| Role | Python and packages | Existing interpreter/location | Scope and remaining provenance |
|---|---|---|---|
| Generator | Python 3.11.15; reinvent 4.7.15; Torch 2.9.1+cu128; RDKit 2026.3.3; NumPy 1.26.4 | `/home/chu/ai4s-tyk2/.venv_reinvent/bin/python`; sampling executable in the same environment | Fixed task reward is delegated to the RF interpreter; generator sklearn 1.9.0 is not used to deserialize the fixed RF assets |
| RF, internal/native chemistry, external property | Python 3.11.15; PyTDC 1.1.15; sklearn 1.2.2; RDKit 2023.9.6; NumPy 1.26.4; Chemprop 1.5.2; Torch 2.4.1+cu121 | `/home/chu/lfr-envs/pytdc-1.1.15-sklearn-1.2.2/bin/python` | Source/candidate identities and structural descriptors use this environment; external inference is CPU-only |
| Graph route execution/audit | Python 3.11.15; Syntheseus 0.7.2; syntheseus-graph2edits 0.2.0; Torch 2.2.2+cu121; RDKit 2025.9.2; NumPy 1.26.4 | `/home/chu/lfr-envs/syntheseus-formal-py311-0.7.2-graph2edits-0.2.0-torch2.2.2-cu121/bin/python` | Saved parent configuration separately records Syntheseus/Torch/RDKit versions |
| Process scheduling | Existing A6000 host and foreground WSL/SSH controller; two GPU workers allowed | `joint_validation512_execution.py` | Hardware assignment/concurrency and per-invocation durations belong in execution records; no isolated throughput benchmark is asserted |

The differing RDKit versions are disclosed. Saved canonical identities remain authoritative for pairing; downstream checks compare them rather than silently rewriting selected structures. A version/identity disagreement is a recorded execution/provenance problem.

**Table S5. Retained assets and their roles.** Paths are current project locations; public release-relative paths and permissions remain to be finalized.

| Asset | Version/source and location | Recorded role or limitation |
|---|---|---|
| REINVENT prior | `/home/chu/ai4s-tyk2/external/REINVENT4/priors/reinvent.prior`; installed generator version above | Public-source correspondence to REINVENT prior v1, Zenodo 15641297, is documented with Apache 2.0 notices; the retained prior's exact original acquisition and binary equivalence remain unverified (S20) |
| Fixed task forests | `/home/chu/lfr-assets/tdc-oracle-interface-20260908/oracle/{gsk3b,jnk3}_current.pkl` | Existing PyTDC-compatible forests; tested full-size filenames/IDs are documented in Harvard Dataverse V85 (CC0); original detailed training composition and fresh current-endpoint recovery remain unresolved (S20) |
| Canonical stock input | `/home/chu/trace_lfr_g0_mapping_20260721/raw_external/paroutes_v2_7341155/stock_n1.txt` | PaRoutes-derived; 13,432 connected canonical SMILES; source access is Zenodo 7341155 with CC BY 4.0 (S20) |
| Planner stock representation | `/home/chu/trace_lfr_g0_stock_identity_20260721/runs/n1_aiz_20260721/stock_inchikeys.txt` | 13,415 standard InChIKeys; target removed per query; this representation count differs from canonical SMILES count |
| Native templates/implementation | `development/external_baselines/SynTwins` at commit `53d6bbcfc540971460278ad85e31d312b749fa55`; `data/templates.csv` and `data/retro_templates.csv` | Retained preprocessing compared to exact version; 1,163 reverse templates; applicable redistribution grant for the root code/templates remains unknown (S20) |
| Native preprocessing | `development/official_syntwins_small_comparison_v1/data` | Ordered building blocks and both template tables; historical first initialization is reported separately from reuse |
| Graph2Edits weights | `/home/chu/trace_lfr_syntheseus_smoke_20260721/model-cache/Graph2Edits_backward`; retained manifest lists `epoch_123.pt` | Checkpoint access is Figshare 25047314/file 44194301 with a deposit MIT declaration; the held-out-role label does not certify exclusion from every generator/property-training source (S20) |
| External ensemble | CombiMOTS commit `3f644b727982a95b4f87115ccb8d6b4c918a02ba`; `development/external_property_raw96_v1/published/models/gsk3b_jnk3/fold_{0..9}/model_0/model.pt` | Released ten-model ensemble; complete original curated training table/order and weight-specific redistribution grant remain unresolved; use exact upstream access (S20) |

The Graph weight manifest is retained at [recovered asset record](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/assets/graph2edits_backward.asset-manifest.json). Native and external Git versions are also documented in the existing [bibliography/provenance record](G:/逆合成结构化反馈/bibliography_20261002.json). Public papers provide component credit; the saved project asset records establish the tested implementation.

## S11. Computational work, shared caches and restart boundaries

Table S6 specifies what must be counted. The source-generation and initial parent costs are common preparation. They are not charged only to the successful method. The candidate-route ceiling excludes native enumeration, template/root work, RF decomposition/prediction, initialization, audits and serialization.

**Table S6. Historical New512 allowance and accounting definitions.** Logical policy exposure and newly executed physical work are distinct quantities.

| Stage | Logical per-policy/source allowance | Newly executed work and reused work to report |
|---|---|---|
| Source training/sampling | Eight complete scheduled runs, 64 raw terminal draws each | Observed training/sampling times, oracle rows/batches, initialization and controller walltime |
| Initial source routes | 12 calls / 45 soft seconds for each valid source measurement | Raw-slot logical requests/calls; globally grouped unique queries; compatible history reused; fresh calls/time; audit/setup/save work |
| Shared internal roots/weights | At most four selected replay roots; root discovery itself is not included in the 16 candidate combinations | Reverse-template applications, source-replay forward applications, source RF decompositions, matching/ranking time and reused complete roots |
| Internal candidate formation | At most 16 combinations and 8 RF score slots per source/policy | Logical cache hits still consume allowances; actual new forward applications/errors/time and unique task predictions are reported separately |
| Native candidate formation | Depth3/k10/1024/fast; 300-second enumeration request and 8 scored outputs | Template application counts, completed/limited enumeration, initialization/worker durations, full-output reuse, requested/adopted seeds and incomplete cost observations |
| Candidate routes | 12 or 6+6 calls; one shared 45-second allowance | Each policy's logical calls/search duration; each worker's new calls/duration, historical/current query reuse, audits/model loading/setup/graph-read/save costs |
| External fixed-output check | Ten fold vectors per required canonical structure | Missing-only checkpoint–structure predictions; historical/current cache reuse; model loads and per-batch times |
| Saved-output analyses | Fixed final structures and all positions | Deterministic postprocessing/descriptors; no new training, RF or Graph calls |

Actual Graph work is counted from each worker's unique stored `queries/query-*.result.json` ledger. Reused policy outcomes can retain a former query's `actual_new_calls`; summing that field across policy outcomes would charge the same invocation repeatedly. Policy logical work and trigger-filtered carried shared observations are exported separately. Native raw records distinguish new enumeration from historical and current-study reuse by their saved provenance, even when their logical time/call costs are identical.

The completed internal phase retains `internal_policy_results_20261002.json` paired with `internal_execution_20261002.json` in each worker directory. These preserve the first resident Graph model load before a later native invocation writes `policy_results.json` and `execution.json` again. The cost reader deduplicates captures of one invocation by output directory, PID, start time and cohort membership; separate native/restart invocation loads remain separate. Missing or interrupted load observations remain unavailable or lower bounds. The two first internal loads are recorded as 1.8035 and 1.7934 seconds; neither is an isolated end-to-end policy latency. [Cost accounting implementation](G:/逆合成结构化反馈/joint_validation512_costs.py) accepts the paired snapshots explicitly or discovers saved internal report snapshots in the worker directories.

The completed cost summary is based on actual raw query/native/fold ledgers and keeps unknown instrumentation separate. Newly executed source planning consists of 491 queries and 4,987 Graph calls; all 133 reported solved source routes have saved verified witnesses. Candidate planning stores 1,353 unique worker-local queries: 1,342 new queries and 11 compatible historical twelve-call queries. New candidate work uses 3,918 Graph calls and 7,632.439 search seconds. All 983 solved candidate query records have `verified_solved=True` and `verified_saved_graph=True`; the other 370 are budget-unsolved, with no saved audit disagreement. Source and candidate Graph calls sum to 8,905 newly executed calls. The sum of all policy logical candidate calls is 5,571, including reused measurements; it is not an additional bill to add to the 3,918 actual calls.

| Completed shared stage | Actual observed work | Recorded time or observation boundary |
|---|---|---|
| Eight source training/sampling runs | 15,993 training RF rows; 512 raw terminal positions | 2,909.437 training-process seconds and 46.072 sampling-process seconds; parallel generation controller walltime 1,481.737 seconds |
| Source RF | 490 new task predictions and one cached source prediction | 0.746 scoring seconds; 4.311 model-load seconds |
| Original source Graph queries | 491 new queries; 4,987 calls | 9,135.930 search seconds; 19.000 observed resident-model-load seconds |
| Shared internal root formation | 470,169 reverse-template applications; at least 22,599 source-replay forward applications | Source-replay counts are missing for ten reused root records; the observed sum is a lower bound, not an exact total |
| Internal candidate forward/RF work | 10,779 new forward applications; 6,646 new RF predictions, including one disclosed repeat | 6.675 forward seconds; 15.237 RF scoring seconds; 608.120 enclosing process-stage seconds |
| Native candidate formation | 439 new enumerations including eight limitations; 14,973,625 reverse and 7,496,772 forward template applications; 3,138 new RF predictions | 10,544.890 summed enumeration seconds; 12.517 RF scoring seconds; 8,807.248 enclosing process-stage seconds with two enumeration workers |
| All policy candidate Graph queries | 1,342 new and 11 historical queries; 3,918 new calls | 7,632.439 new search seconds; 23.527 observed resident-model-load seconds is a lower bound because an interrupted invocation lacks a load observation |
| Fixed-output external evaluation | 12,110 new and 140 historical checkpoint–structure pairs for 1,225 structures | 62.718 workflow seconds; 70.095 enclosing process-stage seconds |

Policy logical candidate calls are 898 for B_route12, 1,013 for B_joint12, 812 for B6+6, 782 for Mix, 722 for C, 405 for Morgan and 939 for native. B_route12 acts on 358 valid route-unresolved sources, while joint-trigger policies act on 439 joint-failed sources. Root/attribution work carried into a policy comparison is filtered by its actual trigger. The common per-source ceiling permits different realized use; neither these logical values nor shared actual timings establish an isolated equal-total-cost advantage.

The original controller terminated with a recorded BrokenPipeError after completed internal/native formation stages. Its completed scientific records were retained. A separate recorded completion controller resumed the same native policy directories, completed both native route workers, analysis and external evaluation; neither eight source training runs nor six completed internal policies were regenerated. Both `pipeline_execution.json`/archived invocations and `completion_execution.json`/archived completion invocations must be read and deduplicated by PID, process start and start time. Their ended controller intervals total 10,265.688 seconds; this excludes the gap between invocations and is not uninterrupted end-to-end walltime. The nine stage observations include the brief interrupted native launch, the two completed recovery native workers, analysis and external execution. The 439 `.native_request.json` files are saved requests, with zero pending corresponding terminal native results; eight scientific unknowns remain terminal. They must not be labeled 439 pending jobs.

After the two accounting corrections and the cross-stage RF audit were incorporated, the final cost CLI was deployed and rerun on the existing A6000 host. It exited 0 with all 512 positions and all seven policy reports complete, zero input-read issues, zero pending native requests, and both controller/nine stage observations included. The refreshed `costs_20261002.json` and its two CSV exports have been copied back to the recovered evidence directory. This was a read-only accounting refresh; the scientific effects, actual 6,646/3,138 RF-stage costs and native unknown bounds were preserved.

RF caches use task plus canonical structure; external caches store both ordered task outputs per canonical structure and fold. Internal forward-cache identities include reaction reference and ordered reactants, with source-exclusion provenance when an old cache stores filtered products. Only complete unified root records are eligible for root-list reuse; earlier incomplete one-role lists are not substituted for the two-role operation.

Current internal candidate lookup explicitly includes historical RF records, Old256 training history, New512 raw-source scores and the current candidate ledger. Source preparation and native scoring explicitly also ingest all eight New512 training feedback traces. The internal candidate loader does not explicitly bridge every New512 training-only score into its lookup. The completed [canonical overlap audit](G:/逆合成结构化反馈/joint_validation512_training_score_overlap_20261002.json) compared all 15,993 training score rows from the eight traces with the 6,646 newly executed candidate RF predictions using the source preparation's connected isomeric canonical-SMILES rule. It confirmed **one repeated JNK3 prediction** (1/6,646, 0.0150%), zero GSK3B repeats and no score conflict. The repeated structure `COc1ccc(Cc2ccccc2)cc1` had score 0.0 in both records and failed candidate task qualification. This already executed prediction remains included in the actual 6,646-prediction cost; no saving or individual prediction time is inferred retrospectively. The present batch's candidate enrollment and completed endpoints are not regenerated. The missing history bridge is a future cache correction, and exhaustive avoidance of repeated RF work is not claimed for this execution.

The completed [internal/native RF stage-overlap audit](G:/逆合成结构化反馈/joint_validation512_RF_stage_overlap_20261002.json) also identifies 315 repeated task-structure predictions between the parallel internal and native stages: 132 GSK3B and 183 JNK3, with zero score conflicts. The internal stage executes 6,646 new prediction rows and native executes 3,138, for 9,784 actual rows and 9,469 unique task-structure identities. Native reads its RF history once at startup while internal candidate scoring is still writing later records, so it does not ingest those later cache additions. These already executed repeated predictions remain charged. The one internal/training repeat above is a different overlap boundary, and the two audits must not be conflated into one inferred saving or a new scientific endpoint. Global zero-repeat execution is not claimed. For future parallel execution, the ordinary RF cache must be refreshed or coordinated before prediction; completed outcomes are preserved without retrospective reruns.

Per-source plans are saved before candidate RF lookup; completed source records are consumed on restart. External fold batches save at 32-structure boundaries. Native and Graph observed interruptions preserve unknown records. Unobserved process interruption can leave only a lower bound on spent work, which is not presented as measured zero. Reproduction preserves output directories and ordinary identity/status comparisons; it does not discard unfavorable records or rerun old measurements to manufacture new replication.

Summing stage durations under overlapping CPU/GPU execution is not the controller's walltime and not isolated end-to-end runtime of each policy. Shared cache hits remove new invocation work but do not remove historical acquisition costs. Different candidate coverage can cause policies to use different amounts of their common ceiling. Native and internal formation budgets are different. These boundaries prevent an equal-total-cost efficiency claim from being inferred solely from slot yield or call caps.

## S12. Reproduction from saved sources and existing evidence

The entry point for reproducing the policy comparison is the completed, unchanged `sources.json` and recorded `design.json`. It does not require rerunning eight source trainings or the Old256 studies. Restore the retained stock, template/preprocessing, RF/checkpoint assets and historical measurement files at their recorded locations or provide the corresponding CLI overrides. The code directory must include the new stage scripts and their unchanged common/Old256 helpers; `PYTHONPATH` also includes the project's source package and historical helper location.

The following Bash commands describe reproduction or interruption recovery from saved sources on the existing host; they were not executed to prepare this supplement. All actual scientific stages in this batch have completed. For inspection and manuscript reanalysis, use saved records and analysis-only commands. Stage commands retain their existing output directories/cohort assignments and should not be launched concurrently with an active controller.

```bash
LFR_BASE=/home/chu/lfr-syn/artifacts/lfr_syn_v2_20260914
LFR_CODE="$LFR_BASE/code"
LFR_STUDY="$LFR_BASE/development/joint_validation512_v1"
LFR_RF_PY=/home/chu/lfr-envs/pytdc-1.1.15-sklearn-1.2.2/bin/python
LFR_GRAPH_PY=/home/chu/lfr-envs/syntheseus-formal-py311-0.7.2-graph2edits-0.2.0-torch2.2.2-cu121/bin/python
export PYTHONPATH="$LFR_CODE:/mnt/c/Users/Chu:/home/chu/lfr-syn/src"

# Complete/resume internal candidates from the saved 512 source positions.
"$LFR_RF_PY" "$LFR_CODE/joint_validation512_candidates.py" \
  --output "$LFR_STUDY" --sources "$LFR_STUDY/sources.json"

# Complete/resume the disclosed native candidate stage; RF stays single-writer.
"$LFR_RF_PY" "$LFR_CODE/joint_validation512_native.py" candidates \
  --workers 2 --sources "$LFR_STUDY/sources.json" \
  --design "$LFR_STUDY/design.json" \
  --output "$LFR_STUDY/syntwins_depth3_k10"

# Preserve the controller's worker directories/cohort coverage.
CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 \
"$LFR_GRAPH_PY" "$LFR_CODE/joint_validation512_graph.py" \
  --study "$LFR_STUDY" --output "$LFR_STUDY/policy_graph/gpu0" --device cuda:0 \
  --native "$LFR_STUDY/syntwins_depth3_k10/candidate_stage.json" \
  --cohort gsk3b-rep0 --cohort jnk3-rep1 --cohort gsk3b-rep2 --cohort jnk3-rep3

CUDA_VISIBLE_DEVICES=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 \
"$LFR_GRAPH_PY" "$LFR_CODE/joint_validation512_graph.py" \
  --study "$LFR_STUDY" --output "$LFR_STUDY/policy_graph/gpu1" --device cuda:0 \
  --native "$LFR_STUDY/syntwins_depth3_k10/candidate_stage.json" \
  --cohort jnk3-rep0 --cohort gsk3b-rep1 --cohort jnk3-rep2 --cohort gsk3b-rep3

"$LFR_GRAPH_PY" "$LFR_CODE/joint_validation512_analysis.py" --study "$LFR_STUDY" \
  --results "$LFR_STUDY/policy_graph/gpu0/policy_results.json" \
            "$LFR_STUDY/policy_graph/gpu1/policy_results.json"

# Read saved cost ledgers only; no model, chemistry or search execution.
# Retain the existing canonical training-score overlap audit in LFR_STUDY
# or its parent directory so it is carried into the cost report.
"$LFR_RF_PY" "$LFR_CODE/joint_validation512_costs.py" --study "$LFR_STUDY" \
  --graph-output "$LFR_STUDY/policy_graph/gpu0" \
  --graph-output "$LFR_STUDY/policy_graph/gpu1" \
  --graph-snapshots "$LFR_STUDY/policy_graph/gpu0/internal_policy_results_20261002.json" \
                    "$LFR_STUDY/policy_graph/gpu0/internal_execution_20261002.json" \
                    "$LFR_STUDY/policy_graph/gpu1/internal_policy_results_20261002.json" \
                    "$LFR_STUDY/policy_graph/gpu1/internal_execution_20261002.json" \
  --output "$LFR_STUDY/costs_20261002.json"

# First inspect missing external fold coverage; no property model is loaded.
"$LFR_RF_PY" "$LFR_CODE/joint_validation512_external.py" --mode prepare \
  --study "$LFR_STUDY" \
  --results "$LFR_STUDY/policy_graph/gpu0/policy_results.json" \
            "$LFR_STUDY/policy_graph/gpu1/policy_results.json"

# Fill missing external vectors only after final fixed outputs are complete.
"$LFR_RF_PY" "$LFR_CODE/joint_validation512_external.py" --mode run \
  --study "$LFR_STUDY" \
  --results "$LFR_STUDY/policy_graph/gpu0/policy_results.json" \
            "$LFR_STUDY/policy_graph/gpu1/policy_results.json"

# Future repeated analysis can use --mode analyze instead of run:
# complete caches are required, and no property model is loaded.
```

The commands above are sequential for clarity. The existing controller can overlap internal Graph work with native CPU enumeration and later supply native candidates to the same outcome directories. For one physical GPU, process the two recorded worker groups sequentially with the available physical card while preserving logical `cuda:0`, directory and cohort assignments, and record the changed schedule. Creating a new all-cohort worker directory would lose worker-local current-query reuse.

The Graph script may return status 1 after writing complete policy results with execution limitations; its `execution.json` distinguishes that state from an interrupted process. The controller checks actual stage completion before downstream work. A native placeholder is replaced when its real candidate stage is supplied; a genuinely executed limitation is preserved. Final external analysis accepts all 512×7 rows only after pending candidate-stage placeholders are resolved.

**Table S7. Minimum saved files and interpretation.**

| File relative to New512 study | Purpose |
|---|---|
| `design.json`, `generation_complete.json` and each `{task}-rep{replicate}/oracle-only/{summary.json,sampling_summary.json,samples.csv,feedback.jsonl}` | Source design, complete training/sampling membership, raw draw order and task traces |
| `raw_scores.json`, `parent_graph_plan.json`, `parent_graph_records.json`, `sources.json` | Score reuse, canonical source grouping, saved Graph/audit provenance and all final source identities |
| `candidate_stage.json`, `source_records/*.plan.json`, `source_records/*.complete.json`, `new_scores.json`, `forward_cache.json` | Before-score enrollment, selected roots/weights/permutations/Mix proposals and recorded RF/forward work |
| `syntwins_depth3_k10/{design.json,assets.json,native_history.json,candidate_stage.json,row-*.native.json,row-*.scored.json}` | Native settings/compatibility/reuse, output order, limitations, adopted seeds and RF scores |
| `policy_graph/gpu{0,1}/{config.json,queries/*,outcomes/*,policy_results.json,execution.json}` | Exact request budgets, saved graph/audit witnesses, all policy outcomes and completion/cost status |
| `policy_graph/gpu{0,1}/{internal_policy_results_20261002.json,internal_execution_20261002.json}` | Paired completed internal-phase snapshots; preserve scientific rows and first resident model-load costs before native continuation |
| `analysis_internal_policies_20261002.json` | Retained internal-phase snapshot; its native placeholders are historical phase status rather than final scientific outcomes |
| `analysis.json`, `fixed_output_sensitivity_20261002.json` | Completed seven-policy eight-run paired contrasts/strata/bounds and fixed-output RF/historical-pool sensitivities |
| `pipeline_execution.json`, `pipeline_invocation-*.json`, `completion_execution.json`, `completion_invocation-*.json` | Original and recovery process intervals, observed stages and resume provenance; ordinary PID/start identity removes duplicate archived/current observations |
| `joint_validation512_training_score_overlap_20261002.json`, `costs_20261002.{json,csv}`, `costs_20261002_policy_sources.csv` | Canonical training/candidate RF overlap audit, actual/shared versus policy logical cost accounting and per-source provenance |
| `external_fixed_outputs_v1/{prediction_plan.json,new_fold_*.json,predictions.json,structure_descriptors.json,analysis.json,slot_outcomes.csv}` | Missing-only external fold work, fixed-output probabilities/descriptors and full-slot exports |

Pure-mock tests exist for candidate, Graph and analysis control flow, native interruption/parallel restart, external fixed-output/cache/unknown semantics, and cost accounting for archived/recovered invocation coverage and terminal native requests. These establish implementation behavior and restart accounting; they are not molecular efficacy evidence. The supplement was prepared without model inference, resource benchmarking or additional experiments.

Completed Old256 evidence can be re-read from the recovered snapshot without new calls: [primary evidence](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/independent_training256_v1/analysis.json), [historical-pool analysis](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/historical_pool_diagnostic256_v1/analysis.json), [reuse sensitivity](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/independent_training256_v1/reuse_sensitivity_20261001.json), [budget diagnostic](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/budget_replay256_v1/analysis.json) and [budget fixed-output external sensitivity](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/budget_replay256_v1/external_sensitivity.json).

## S13. Reproducibility scope and remaining access limitations

Each completed computational-effect cohort retains its own source population, qualification mode and contrast labels. The negative New512 Mix primary, null HIV LR primary, earlier GraphGA task heterogeneity/concentration/alert findings and historical native unknown bounds remain reported. The New512 allocation component uses existing matched measurements (S22). The adopted fixed-qualification result comes from the sixteen new optimizer runs in S23; earlier qualification development and unsuccessful alternatives remain exploratory (S24).

Private code, environment recipes, saved-analysis inputs, selected fresh/history executions and relocated replays are available locally. On the existing A6000 host, separate clean RF/chemistry and Graph environments were installed, and an existing one-source software fixture completed native execution and independent saved-output checks (S20). That software fixture is excluded from efficacy samples. Installation on the existing host is not a new-machine deployment, and the complete scientific cohorts have not been regenerated from scratch in a second independent installation. Saved-record replay and native-witness auditing are also distinct from rerunning predictor inference.

The tested RF, Graph2Edits and stock deposit licenses are documented in S10/S20. Root SynTwins/template and CombiMOTS weight redistribution grants remain unresolved; the retained REINVENT prior's original acquisition and deposited-binary correspondence remain unverified. Original kinase predictor/prior training compositions, historical HIV table mapping and complete training-data exclusion cannot be reconstructed from available exact public-table checks. Stable public code/data deposition and its final availability statement remain pending. Private archives do not establish public access or permission to redistribute upstream assets.

Historical New512 retains 315 internal/native repeated RF scores and one training repeat, with all actual expenditure charged. The prospective GraphGA1024 implementation uses cumulative missing-only observation records, exact compatible origins and shared complete saved plans; its measured physical costs and remaining reuse scopes are stated in S24. No globally duplicate-free execution or uncached per-arm runtime advantage is inferred from a common route-call ceiling. Fixed qualification trades guiding RF retention and output diversity against qualified-route yield. The substantial shared qualification gain does not establish a general B-over-permutation advantage or cross-task biological transfer.

## S14. Exact TreeSHAP in an adaptive matched extension of New512

### S14.1. Question, matched population and unchanged workflow

After completion of the original seven-policy comparison, an adaptive extension tested whether a different RF attribution definition could improve the same bounded replacement workflow. The extension compares `TreeSHAP_joint6plus6` with the already completed `B_joint6plus6`; it is exploratory and uses the same source molecules rather than an independently generated validation cohort. The original prespecified Mix−B primary contrast and its negative result remain unchanged (S8). The extension does not designate an additional confirmatory primary success.

All 512 original raw positions from the eight completed training runs were retained, including 21 invalid positions and 52 original RF-qualified/audited-route successes. Attribution and candidate intervention were applied to the 439 valid original joint-failed positions. There was no new source training, source sampling or retrospective selection of successful positions. Original source identity, slot order and run/task membership were retained. The same fixed RF assets, class-1 probability, binary radius-2/2048 Morgan inputs, source-relative RF qualification, QED tolerance, inventory, forward reactions, Graph2Edits planner, stock and saved-route audit standard were used.

The changed component was the source feature attribution used to rerank the complete reaction-compatible inventory. The extension reused each source's original selected exact-replay roots; it did not introduce extra disconnection opportunities or restrict the alternative ranking to B's selected products. Positive contributions were mapped to present unfolded source environments, divided over folded-bit collisions and propagated through the same shared source/precursor/candidate environments. Absent or negative environments did not gain a new role. Retained-weight, Morgan-similarity and canonical-SMILES tie rules were unchanged. Enrollment preceded task-score lookup, with at most 16 forward attempts and eight QED-qualified unique scoring positions per source. RF-qualified candidates were ranked by RF score, QED and canonical SMILES; the first two were attempted under the same six-plus-six call ceilings and shared 45-second soft planning budget. Replacement required a verified saved route; otherwise the original output was retained. External predictions were added only to the fixed final outputs and never influenced candidate choice.

Historical chemistry was reused only under compatible reaction/reactant/source-exclusion, task/structure or exact Graph request identities. An eligible saved Graph record was selected by the recorded lexical rule, independently of its success, and its saved seed was adopted explicitly. A completed 12-call request was not truncated into a fictitious six-call observation. Original B outcomes were read as the reference and were not regenerated. These constraints make the attribution comparison matched within the tested cohort; they do not make an adaptive same-source analysis independent of earlier results.

### S14.2. Numerical failure and exact path-dependent repair

The initial implementation used `shap==0.46.0`, `TreeExplainer(data=None, feature_perturbation="tree_path_dependent", model_output="raw")` and non-approximate values. Tree-path-dependent SHAP uses the model's stored path covers without a separate background sample (Lundberg et al., 2018; [official TreeExplainer documentation](https://shap.readthedocs.io/en/latest/generated/shap.TreeExplainer.html)). On the unchanged first source, `gsk3b-rep0-r000`, both the original RF and SHAP's prediction representation returned class probabilities approximately `[0.67, 0.33]`; the source's saved score was 0.33. However, the initial SHAP reconstruction was `[0.794612546941983, 0.3306130965957135]` with expected values `[0.9374907935439072, 0.06250920645609306]`. The additivity check failed before any scientific candidate plan was produced. A diagnostic-only extraction of the failed vector did not enter candidate ranking or route evaluation. Tree-cover and depth inspection did not establish an asset/prediction mismatch: the first forest contained 100 trees, maximum depth 78 and 385,384 nodes, with conserved positive covers. The failed invocation is retained as an implementation failure, not a molecular result.

The replacement implementation evaluates the same path-dependent Shapley game by a positive leaf-polynomial integral, without changing the RF, coalition definition, inputs or `1e-8` additivity tolerance. Quadrature reformulation of path-dependent TreeSHAP has been described by Wettenstein, Mitchell and Yu (2026); that preprint reports a practical fixed eight-point formulation. The present implementation uses path-specific order to integrate each finite polynomial exactly in real arithmetic and credits that quadrature principle as prior work. It does not claim a new attribution algorithm or use a fixed eight-point approximation.

For a leaf $\ell$, repeated splits on feature $i$ are merged into one coalition variable. Let $z_{\ell i}$ be the product of child/parent weighted-cover fractions at its splits, $h_{\ell i}(x)$ the product of the corresponding branch-consistency indicators for the source, $U_\ell$ the unique split features and $v_\ell$ the class-1 leaf probability. Its feature contribution is

$$
\phi_i^{(\ell)}(x)=v_\ell\bigl(h_{\ell i}(x)-z_{\ell i}\bigr)
\int_0^1\prod_{j\in U_\ell\setminus\{i\}}
\bigl[(1-t)z_{\ell j}+t h_{\ell j}(x)\bigr]\,dt.
$$

The integrand has degree at most $d-1$, where $d=|U_\ell|$. Gauss-Legendre quadrature with $\lceil d/2\rceil$ points integrates that polynomial exactly in real arithmetic; constant leaves contribute only to the baseline. Positive prefix/suffix products and `long double` path products/accumulation avoid subtraction-based path unwinding. Contributions and leaf baselines are summed within trees and averaged across the unchanged forest. Ordinary finite-precision output error remains measurable; algebraic exactness is not a claim of bitwise arithmetic exactness. The baseline is checked against the forest's root class-1 probability, and each reconstruction is checked against both the current RF prediction and the saved source score. The implementation uses no attribution rescaling, interventional background, RF refit or approximate single-order attribution.

Three distinct validations were retained. An independent local calculation covered 80 small-tree cases, including both classes, repeated path features and constant trees: exhaustive rational coalition values and the rational leaf formula agreed exactly; quadrature's maximum feature, reconstruction and baseline errors were respectively `7.772e-16`, `8.882e-16` and `4.441e-16`. On three shallow reference forests, all 15 full-vector comparisons with SHAP 0.46.0 had maximum absolute feature error `1.665e-16`. The unchanged first real source reconstructed 0.33000000000000007 from baseline 0.06250920645609306, with error `5.551e-17`. Finally, all 439 intervention sources passed saved-score/current-prediction/additivity checks; their largest reconstruction error was `1.665e-16`. The 80-case proof check and shallow/first-source implementation tests are numerical evidence, not additional candidate or joint-yield observations.

### S14.3. Completed endpoints and conditional paired contrasts

Both compared arms have complete scientific endpoints at all 512 positions, with zero execution-unknown or pending positions. The external endpoint follows S9: audited route AND the relevant ten-checkpoint mean probability at least the stated threshold. It replaces the RF threshold for this fixed-output sensitivity; it is not a triple RF/external/route requirement. Canonical task–structure uniqueness is distinguished from raw-slot success counts.

**Table S8. Fixed-output outcomes in the matched attribution extension.** Each arm retains denominator 512. B values are the completed original reference; source values are the unchanged original outputs.

| Endpoint | Original sources | Saved B6+6 | TreeSHAP6+6 |
|---|---:|---:|---:|
| RF-qualified AND audited-route slots | 52 | 311 | 302 |
| Audited-route slots | 133 | 361 | 356 |
| Unique RF-joint task–structures | 52 | 310 | 300 |
| Replaced slots | 0 | 259 | 250 |
| External mean ≥0.5 AND audited-route slots | 23 | 176 | 175 |
| Unique external-0.5-joint task–structures | 23 | 175 | 173 |
| External mean ≥0.3 AND audited-route slots, exploratory | — | 197 | 205 |
| External mean ≥0.7 AND audited-route slots, exploratory | — | 150 | 150 |

TreeSHAP increases RF joint yield by 250 slots relative to the original sources, while the matched B reference increases it by 259. For RF joint yield, TreeSHAP−B has 27 winning and 36 losing source positions, net −9/512 (−1.758 percentage points). For audited routes the corresponding values are 22 wins, 27 losses and net −5/512 (−0.977 percentage points). Among GSK3B positions RF joint counts are 165 versus 164 and route counts are 185 versus 185; among JNK3 positions RF joint counts are 137 versus 147 and route counts are 171 versus 176. Thus the RF disadvantage is concentrated in the tested JNK3 runs rather than uniformly negative across both tasks.

Paired uncertainty was computed by the same conditional method as S8: enumerate all 65,536 ordered resamples of four complete runs within each of the two fixed tasks, maintaining equal task weight. These empirical percentile intervals describe only these predictors, assets and four runs per task; they are not molecule-independent inference, multiplicity-adjusted confirmation or a guarantee of biological generality. Wins and losses refer to matched raw positions, whereas uncertainty resamples complete paired runs.

**Table S9. Exploratory TreeSHAP−B paired contrasts.** Effects and 95% conditional intervals are percentage points; run directions are positive/negative/tied across eight complete runs.

| Endpoint | Slot wins/losses | Net slots | Effect | Conditional interval | Run directions |
|---|---:|---:|---:|---:|---:|
| RF joint | 27/36 | −9 | −1.758 | [−2.734, −0.781] | 2/5/1 |
| Audited route | 22/27 | −5 | −0.977 | [−2.344, +0.391] | 3/4/1 |
| External-0.5 joint | 27/28 | −1 | −0.195 | [−2.148, +1.953] | 3/4/1 |
| External-0.3 joint, exploratory | 35/27 | +8 | +1.563 | [+0.391, +2.734] | 6/1/1 |
| External-0.7 joint, exploratory | 24/24 | 0 | 0.000 | [−1.367, +1.758] | 2/5/1 |

In replicate order 0–3, GSK3B RF net differences are `[+1, 0, +1, −1]` and JNK3 RF differences are `[−3, −4, −1, −2]`. External-0.5 differences are `[−2, 0, +2, −1]` for GSK3B and `[−1, +3, −3, +1]` for JNK3. The external-0.5 interval includes zero; this does not establish equivalence. The threshold-0.3 result is a disclosed exploratory sensitivity and does not replace either the original RF endpoint or the threshold-0.5 comparison. Overall mean external probabilities are 0.403366 for TreeSHAP and 0.396637 for B, a descriptive +0.006729 difference; this continuous summary does not establish superior thresholded joint yield.

At external threshold 0.5, fixed-output no-PAINS, no-PAINS/BRENK and MW≤500 counts are respectively 170, 125 and 165 for TreeSHAP, compared with 173, 132 and 167 for B. QED-tolerance filtering preserves all 175 and 176 successes. These are descriptor sensitivities of unchanged selected outputs, without reselection or additional activity validation. The available 95,159-row public labeled table, representing 95,157 unique valid structures, overlaps two TreeSHAP output positions and one B position. They are `jnk3-rep1-r029` in both arms and additionally `jnk3-rep2-r047` in TreeSHAP; all are non-RF-joint source fallbacks below the external-0.5 threshold. Excluding their success contribution while retaining all 512 positions leaves counts at 302/175 and 311/176. This supports robustness to these exact public matches, while the full curated checkpoint-training table and scaffold membership remain unavailable (S9); it does not establish original training-set exclusion.

### S14.4. Complementarity is an observed diagnostic, not an achieved combined policy

RF joint success is shared at 275 positions, TreeSHAP-only at 27, B-only at 36 and absent in both at 174. The outcome-informed union is 338/512. External-0.5 success is shared at 148 positions, TreeSHAP-only at 27, B-only at 28 and absent in both at 309, giving union 203/512. These unions are diagnostic potential upper bounds for selecting between the two already executed outputs with knowledge of their outcomes. No combined policy achieved 338 or 203 under the shared score/search budget, and those counts must not be presented as a new method's result.

The discordances are mechanistically informative. Among the 27 TreeSHAP-only RF successes, B has no RF/QED-qualified candidate at 15 positions and exhausted candidate planning budgets at 12. Among the 36 B-only RF successes, TreeSHAP has no qualified candidate at 22 and budget-unresolved candidates at 14. For external-0.5-only successes, TreeSHAP's 27 positions correspond to B statuses of six no-qualified, seven budget-unresolved and 14 replacement-solved outputs below the external threshold; B's 28 positions correspond to TreeSHAP statuses of six no-qualified, six budget-unresolved and 16 replacement-solved outputs below threshold. Qualification, bounded planning and predictor disagreement therefore all contribute to discordance.

The ordered first-two qualified candidate lists differ at 221 positions and agree at 291; 178 positions have any shared candidate within the two lists. TreeSHAP has at least one qualified candidate at 283 positions versus 286 for B; both have at least two at 223. The final outputs agree at 350 positions and differ at 162. Successful task–structure sets share 176 identities, with 124 identities unique to TreeSHAP and 134 unique to B. Distinct chemical outputs and exclusive per-slot successes are different quantities. These observations motivate testing attribution selection or candidate allocation using information available before search; they do not validate outcome-based selection or establish a completed attribution combination. The current data do not justify replacing B on the basis of superior overall RF joint yield.

### S14.5. Actual new work, reuse and recovery accounting

The successful candidate stage computed 439 source attributions in 107.995 recorded seconds and made 439 separate source RF prediction/additivity checks. It used the existing exact source roots with zero new retrospective-template discovery or source replay-forward calls and no source-model training. Full inventory matching took 331.139 seconds. Among 4,715 logical candidate forward attempts, 1,628 were physically new and 3,087 were reused; new forward reaction execution took 1.098 seconds. Among 2,350 logical candidate score positions, 813 required new task predictions (361 GSK3B and 452 JNK3) and 1,537 used saved scores. These logical score positions are not 2,350 additional physical model evaluations. The 439 source reconstruction checks are additional physical predictions, accounted separately from candidate scores.

Graph evaluation used 323 distinct stored requests: 164 new requests with 456 physical model calls and 943.682 seconds of new search, plus 159 exactly compatible original-New512 requests carrying 330 logical saved calls and 648.661 saved search seconds. There were zero unknown reused requests and zero historical 12-call requests substituted for six-call requests. Five compatible same-structure requests also adopted current-stage saved results. Accordingly TreeSHAP's full logical policy cost is 808 calls and 1,635.283 search seconds, distinct from its actual new Graph work; model resident loading added 10.256 seconds. The corresponding already-completed B reference costs are not newly charged to the extension.

The fixed-output external stage required 871 canonical structures including sources, or 8,710 logical checkpoint–structure vectors. It reused 7,710 historical vectors and newly predicted 1,000 vectors for 100 missing structures across ten checkpoints; no historical fold prediction was recomputed. All 4,910 B final-valid-slot fold vectors matched the completed original reference. The recorded new inference batches total 2.762 seconds, excluding checkpoint loading and descriptor work; the external workflow invocation totals 35.571 seconds. Descriptors were computed for all 871 required structures, so reuse of molecular inference does not imply zero repeated metadata work.

The first failed candidate controller stage lasted 45.467 seconds and yielded no scientific candidates. Independent mathematical and numerical diagnostics, including three shallow forest fits, 15 comparison rows and one stable first-real-source attribution/prediction, are retained separately from the completed 439-source candidate stage. After the numerical repair, the successful candidate and Graph stages lasted 645.621 and 969.294 seconds. A subsequent analysis invocation failed after 0.571 seconds; recovery resumed analysis and external postprocessing, with completed stages of 6.099 and 39.393 seconds, and did not repeat candidates or Graph. The completed candidate/Graph/analysis/external stage windows sum to approximately 1,660.406 seconds (27.673 minutes). This sum excludes the failed stages, diagnostic/repair work and gaps between invocations; it is not a continuous end-to-end run or total allocated GPU time. The latest controller record is terminal `complete`, and all scientific rows are available.

**Evidence and implementation:** [matched design](G:/逆合成结构化反馈/attribution_extension512_design.json); [candidate summary](G:/逆合成结构化反馈/attribution_extension512_v1/attribution512_candidate_summary.json); [matched RF/route analysis](G:/逆合成结构化反馈/attribution_extension512_v1/attribution512_analysis.json); [fixed-output external analysis](G:/逆合成结构化反馈/attribution_extension512_v1/attribution512_external_analysis.json); [complementarity](G:/逆合成结构化反馈/attribution_extension512_v1/complementarity_summary.json); [initial numerical failure](G:/逆合成结构化反馈/attribution_extension512_v1/attribution512_TreeSHAP_first_source_failure_diagnostic.json); [stable RF implementation checks](G:/逆合成结构化反馈/attribution_extension512_v1/attribution512_TreeSHAP_stable_implementation_validation.json); [independent small-tree mathematical checks](G:/逆合成结构化反馈/attribution_extension512_v1/gauss_math_check.json); [controller histories](G:/逆合成结构化反馈/attribution_extension512_v1/attribution512_pipeline_execution.json); [positive quadrature implementation](G:/逆合成结构化反馈/tree_shap_positive_quadrature.py); [native leaf integration](G:/逆合成结构化反馈/tree_shap_positive_quadrature.cpp); [candidate](G:/逆合成结构化反馈/attribution_extension512_candidates.py), [Graph](G:/逆合成结构化反馈/attribution_extension512_graph.py), [analysis](G:/逆合成结构化反馈/attribution_extension512_analysis.py) and [external](G:/逆合成结构化反馈/attribution_extension512_external.py) scripts. These local records establish the tested extension; external component references below provide attribution credit rather than project result evidence.

## S15. Source-only attribution diagnostics and saved-output routing replay

### S15.1. Feature boundary, paired role concentration and multiple exploration

The completed attribution extension was followed by deterministic local analyses of its saved source inputs and single-arm outcomes. These analyses made zero new model, chemical or Graph executions. They retained all 512 original positions and all eight complete 64-position training runs. The 439 intervention-eligible positions entered feature bins; the 73 retained positions did not. Source-only here means information available after the original source RF/QED, exact root discovery and relevant inventory ranking, but before candidate forward generation, candidate RF scoring or candidate Graph search. The original exact source-replay chemistry had already occurred during root discovery. The analyses therefore do not describe information available before any chemical computation.

Although saved plans are named `before_oracle_lookup`, those files already contain forward products and QED. The feature extractor read only the source/selected-root/role/environment and precursor-ranking fields within the stated time boundary. It excluded candidate products, candidate QED, qualified-candidate counts, candidate scores, final structures, route outcomes, external probabilities and measured execution costs from predictor inputs. Source identity, selected roots, compatible inventory counts and source environment sets were matched across arms for 2,160 selected-root/missing-role records. The full source dictionaries and all saved endpoint labels were also reconciled against the paired records; none of these identity checks supplies independent predictive validation.

For arm a and missing role k, positive source-present environment weights are normalized within the identical saved role set: `p(a,k,i) = w(a,k,i) / sum_i w(a,k,i)`. The role concentration is `HHI(a,k) = sum_i p(a,k,i)^2`. The routing feature first computes `HHI(TreeSHAP,k) − HHI(B,k)` on each identical role, then averages these differences over roles with nonzero positive totals in both arms. It is not the difference of arm means calculated over different role subsets. It describes source environments covered by the saved roles rather than full-source global concentration or candidate-product attribution.

Of the 2,160 role records, 2,124 have nonzero positive totals in both arms. Twenty-four have zero B total but nonzero TreeSHAP total, two have the reverse, and ten have zero total in both arms. Thus 26 have discordant zero/nonzero status, with another ten jointly zero; all 36 are missing for the paired HHI difference and are not imputed as zero. A source with no usable paired role remains missing. In the present records, 409 positions have usable routing inputs and 103 default to B: 73 retained positions and 30 intervention positions without selected roots. Missing inputs, no roots and nonintervened sources default to B for every replay rule.

The diagnostic report contains 62 bins across 13 source or precursor-ranking features. Natural task thresholds, root integers and the HHI zero point were used where applicable; continuous quartile boundaries were derived solely from nonmissing intervention inputs before outcome labels were joined. Feature identification and interpretation nevertheless followed inspection of the complete 512 outcomes, so the exercise remains post hoc and subject to multiple exploration and selective emphasis. It used no significance filtering or multiplicity-adjusted claim. Natural and quartile bins overlap and must not be added as separate successes. All bins, exact boundaries and task/run summaries are retained in the diagnostic JSON.

For illustration, HHI difference >0 covers 109 positions with RF TreeSHAP-only/B-only counts 9/7 (net +2) and external-0.5-only counts 13/10 (net +3). HHI difference ≤0 covers 300 positions with corresponding counts 18/29 (net −11) and 14/18 (net −4). The lowest input quartile, ≤−0.0786312605, covers 103 positions and has RF net −11 and external net −6; the highest quartile, >0.000762526994, covers 102 positions and has RF net +4 and external net +5. The latter RF difference is positive/negative/tied in 5/2/1 runs, and the external difference in 4/1/3 runs. An original exact-replay-root count >56 covers 104 positions with RF net 0 and external net +6 (external run directions 7/1/0). These associations suggest candidate hypotheses about role concentration and available root opportunities; they do not establish concentration as a causal mechanism or identify a validated switching rule.

### S15.2. Three limited rules and the RF-only fitting objective

The next local replay restricted attention to one feature, the paired role HHI difference, and three explicit rules. Every rule selects one already saved single-arm output per original position; it does not union successes or generate a new candidate pool. Rule 1 selects TreeSHAP if the usable difference is strictly >0 and B otherwise. Rule 2 uses the full-input upper-quartile cutoff 0.000762527. Its exact unrounded input cutoff is 0.0007625269941309001; the rounded value makes the same 512 choices. This full-data cutoff and the selected HHI feature are post hoc and cannot serve as an independent validation rule.

Rule 3 leaves out one complete training run at a time. In each of the eight folds, the remaining seven runs (448 raw positions) supply training inputs and RF joint labels only. For this replay, a restricted candidate family was fixed: always B, the natural zero cutoff, and cutoffs at the 0th, 25th, 50th, 75th and 100th percentiles of the training folds' usable HHI input distribution. Each fold evaluates seven candidate rules; it does not scan further features, lambda grids or complex selectors. The objective averages the RF net-rate difference within each task across its available complete runs, then averages the two task means with equal task weights. This prevents the task with four remaining runs from outweighing the task with three. Candidate objectives are compared exactly using rational arithmetic. Equal RF objective favors fewer TreeSHAP selections, then always B and then the larger cutoff. External outcome labels are unavailable to the fitting projection and do not fit or break ties between thresholds. The chosen rule is applied to the 64 positions of the held-out run before their labels are used for replay summaries.

**Table S10. Saved-output routing replays on the unchanged 512 positions.** Reference B has 311 RF joint and 176 external-0.5 joint successes. These are output replays rather than achieved results of a newly executed routing policy; no independent-repeat interval is supplied.

| Rule | TreeSHAP selections | Replayed RF joint | RF wins/losses/net versus B | Replayed external-0.5 joint | External wins/losses/net versus B |
|---|---:|---:|---:|---:|---:|
| Natural HHI difference >0 | 109 | 313 | 9/7/+2 | 179 | 13/10/+3 |
| Full-data upper-quartile cutoff, post hoc | 102 | 315 | 8/4/+4 | 181 | 12/7/+5 |
| Leave-one-run-out RF-fitted cutoff | 169 | 311 | 10/10/0 | 179 | 16/13/+3 |

The leave-one-run-out replay's RF net is zero, with 3/3/2 positive/negative/tied held-out runs; its external net is +3, with directions 4/3/1. GSK3B contributes RF +3 and external +3, whereas JNK3 contributes RF −3 and external 0. The selected cutoffs range from −0.06983849798 to +0.005508448510, indicating unstable cutoff choice across the small number of source runs. All threshold candidates, RF training objectives, selections and held-out results are retained. The natural-zero and full-data upper-quartile apparent RF improvements are +2 and +4 positions, respectively; they do not remain a positive aggregate RF gain when the cutoff is fitted without the held-out run's RF labels.

**Table S11. Leave-one-complete-run-out threshold choices and saved held-out outcomes.** Training objective values are percentage points under the equal-task RF objective. Cutoffs are dimensionless HHI differences.

| Held-out run | Chosen cutoff | Training RF objective (pp) | Held-out TreeSHAP selections | Held-out RF net versus B | Held-out external net versus B |
|---|---:|---:|---:|---:|---:|
| gsk3b-rep0 | -0.02879520806 | 0.7812 | 31 | +1 | +0 |
| gsk3b-rep1 | 0.0003202429121 | 0.5208 | 16 | +1 | +2 |
| gsk3b-rep2 | -0.02809390081 | 0.5208 | 30 | +2 | +2 |
| gsk3b-rep3 | -0.02850385098 | 1.3021 | 33 | -1 | -1 |
| jnk3-rep0 | -0.06983849798 | 0.7812 | 28 | -2 | -1 |
| jnk3-rep1 | 0.00550844851 | 0.8464 | 3 | +0 | +1 |
| jnk3-rep2 | -0.02657723434 | 0.5208 | 22 | +0 | +1 |
| jnk3-rep3 | 0.00382314845 | 1.1068 | 6 | -1 | -1 |

### S15.3. Evidence boundary and implications for further experiments

Leaving out a run controls only the RF-label threshold-fitting step. HHI itself was selected after observing all 512 outcomes, including those later designated held out; this feature-selection information remains shared. The analysis is therefore not a fully independent test or a complete nested assessment of the discovery procedure. Training sets also overlap between folds, so the eight fitted-rule results are not eight independent replications of one fixed prespecified policy. The original 65,536 independent-run bootstrap is not applied to these replays, and no apparent replay benefit is reported with that interval. An external label excluded from threshold fitting still provides a sensitivity on the same sources and broadly shared prediction provenance, not independent biological evidence.

No new routing implementation or its budgeted execution was measured. The replays read outcomes from both completed single-arm workflows, whose individual route budgets and seeds were already recorded. They do not establish that a future selector will obtain the same outputs, scheduling, wall time or total work under its own execution. Computing both source attributions and required role or inventory information entails additional preselection work; if both rankings are computed, that sorting cost also belongs in the accounting. A prospective implementation could generate only the selected candidate pool after computing the source feature, but the resulting compute savings and net efficiency have not been measured. Zero new inference in this local analysis is not zero deployment overhead.

These results did not support adopting the fitted attribution switch. The original Mix−B primary comparison and the historical B outputs remain intact, and the leave-one-run-out routing RF aggregate remains zero. No new selector is included in the adopted fixed-qualification recipe. The later prospective test changes qualification within the original B/C formation operation (S23), rather than treating these post hoc replay gains as achieved combined-policy evidence. These observations establish neither reliable switching nor general attribution superiority, laboratory synthesis success or biological potency.

**Evidence and reproduction:** [source feature diagnostic JSON](G:/逆合成结构化反馈/attribution_extension512_v1/source_presearch_diagnostics.json), [diagnostic reading report](G:/逆合成结构化反馈/attribution_extension512_v1/source_presearch_diagnostics.md), [source feature CSV](G:/逆合成结构化反馈/attribution_extension512_v1/source_presearch_features.csv), [routing replay JSON](G:/逆合成结构化反馈/attribution_extension512_v1/source_routing_replay.json), [routing reading report](G:/逆合成结构化反馈/attribution_extension512_v1/source_routing_replay.md) and [standard-library routing script](G:/逆合成结构化反馈/attribution_extension512_v1/source_routing_replay.py). The JSON records the complete rule families and source identities; CSV `outcome_` columns are labels, never deployment inputs. All outputs are derived from saved records, and these analyses do not repeat the molecular experiments.

## S16. Completed second-planner sensitivity on 128 fixed saved outputs

### S16.1. Same-source population and fixed structures

This post hoc analysis asks whether the comparative pattern of the completed original, `B_joint6plus6` and `C_joint6plus6` outputs remains visible when their final structures are evaluated by a second retrosynthesis planner. It uses raw positions 0–15 from each of the eight completed New512 training runs: four GSK3B runs and four JNK3 runs, with 16 positions per run and a total denominator of 128. This positional selection applies no route-success or property-score filter. Its definition follows observation of the completed New512 experiment, however, so it is not a new independent source cohort or a prespecified confirmation of planner generality.

The original source identities and each arm's final canonical structure, RF score, QED, ten external fold probabilities and saved Graph audit remain unchanged. All six invalid source positions are retained as known zeros in all three arms. Each arm therefore has 122 valid logical query references; the 366 references share 232 distinct valid canonical structures. Query identifiers follow lexicographic canonical-SMILES order, starting at zero. Each query receives seed `10203000 + query_id`; the assigned identifiers are shuffled for execution with seed 2026100201. A shared canonical structure receives one AiZ search, with its result referenced at every corresponding source-position/arm occurrence. The source positions remain the outcome denominator; neither the 232 structures nor successful queries replace it.

The analysis performs whole-target searches on these saved final structures. It does not regenerate proposals, recompute contribution or permutation weights, rescore candidates, or rerun the 6+6 replacement decision under AiZ. Thus the arm names identify the original output-producing policies; they do not describe an AiZ implementation of those policies. The full 512-position Graph analysis and its Mix−B primary contrast are unchanged.

### S16.2. Actual AiZ configuration, budgets and route audit

The recorded execution uses AiZynthFinder 4.4.1, an MCTS-based retrosynthesis system (Genheden et al., 2020), with the saved USPTO expansion model, template library and filter model. The actual protocol records RDKit 2023.9.6 and ONNX Runtime 1.25.1. Both ONNX models use `CPUExecutionProvider`; the scientific searches in this sensitivity did not use GPU inference. The n1 stock contains 13,415 standard InChIKeys. That representation count is distinct from the canonical connected-SMILES inventory count reported in S4.

Every genuinely executed target starts with a reset query cache, a stock exclusion for its standard InChIKey and an independent model allowance. The recorded configuration allows at most 64 search iterations and 45 soft seconds per target, with `max_transforms=12`, cycle pruning and `return_first=false`. The MCTS exploration constant is 1.4, the default prior is 0.5, priors are enabled and the search reward is the state score. The expansion allowance is 32 attempted prediction calls and 64 attempted input rows; the filter allowance is 256 attempted calls and 256 attempted rows. The meter checks the next call before execution and retains attempted and failed work if prediction raises an exception. A denied call is recorded separately rather than added to attempted calls. Timing is checked between iterations, so 45 seconds is a soft search allowance rather than a hard process timeout. The two planners' depth limits have different units: AiZ counts reaction transforms, whereas the saved Graph audit counts alternating And/Or topology edges. Their numeric limits of 12 do not establish equal search depth or cost.

The query records preserve the complete saved search tree, reported solved-node paths and route-audit evidence. A positive endpoint requires a matching saved solved node, canonical target identity, a valid occurrence-based structural replay and supported template replay for every recorded step. The structural checks follow actual parent/child molecule occurrences, preserve repeated precursor occurrences, check route connectivity and reaction depth, and require all terminal occurrences to be stock-supported while excluding the target from the starting material. Template replay uses the saved reaction action without an additional learned-model prediction. A planner-reported solution with missing or contradictory saved evidence remains `unknown_audit`; it is not counted from the reported Boolean alone.

`resource_exhausted` denotes a normal pre-over-limit model-budget stop with complete saved evidence and an audit containing no verified route. It is a known negative within this allowance, not proof that the target has no synthesis route. An ordinary iteration, soft-time or search-exhaustion stop with the same complete negative audit is `unsolved_within_budget`. Missing outputs, OOM or other execution exceptions remain `unknown_execution`; incomplete audit of a reported solution remains `unknown_audit`. These rules retain uncertainty rather than treating technical failures as chemical negatives.

All 232 scientific queries finished: 131 have `verified_solved`, 100 have `resource_exhausted` and one has `unsolved_within_budget`. There are no missing, execution-unknown or audit-unknown queries. These are counts of distinct physical queries, not success counts for an arm's 128 source positions. A subsequent independent identity/inventory check passed for all 232 queries, covering 1,481 saved routes, 8,240 saved steps, 8,026 leaf-occurrence references and 13,475 distinct serialized molecules. It independently recomputed canonical identities, parent transform relationships and actual standard-InChIKey stock membership, with zero new model calls and zero new route searches. This verifies the saved identity/inventory evidence; it does not establish laboratory reaction feasibility or measured synthesis success.

### S16.3. Endpoints, paired counts and conditional uncertainty

The route endpoint is the audited AiZ whole-target result. `rf_task_threshold_and_aiz` combines that route with the saved RF score at the original task thresholds, 0.10 for GSK3B and 0.05 for JNK3. The separate external endpoint combines the audited route with the arithmetic mean of the saved ten fold probabilities at least 0.5; it does not add an RF qualification requirement. No RF or external-model prediction is rerun. External cutoffs 0.3 and 0.7 remain exploratory outputs in the analysis JSON.

**Table S12. Completed second-planner outcomes and saved Graph reference on the same 128 positions.** All entries are success counts with denominator 128, retaining six invalid source positions as zeros. Graph rows are unchanged saved reference results, not new Graph searches. No union of the two planners' successful routes is used.

| Planner/result source | Saved output arm | Audited route | RF task threshold and route | External mean ≥0.5 and route |
|---|---|---:|---:|---:|
| AiZ, newly executed whole-target queries | original | 43 | 21 | 9 |
| AiZ, newly executed whole-target queries | B_joint6plus6 | 80 | 68 | 42 |
| AiZ, newly executed whole-target queries | C_joint6plus6 | 65 | 49 | 20 |
| Graph, saved reference | original | 28 | 10 | 6 |
| Graph, saved reference | B_joint6plus6 | 88 | 76 | 48 |
| Graph, saved reference | C_joint6plus6 | 73 | 58 | 26 |

Against C, B has 23 RF-joint wins and four losses, a net gain of 19 positions or 14.84375 percentage points. The conditional 95% empirical percentile interval is [9.375, 20.3125] percentage points. For the external-0.5 joint endpoint, B has 24 wins and two losses, a net gain of 22 positions or 17.1875 percentage points, with conditional interval [10.15625, 21.875]. RF joint counts for GSK3B are 34 versus 27 and for JNK3 are 34 versus 22. External joint counts are 27 versus 14 and 15 versus six, respectively. RF net differences are strictly positive in all eight runs; external net differences are positive in seven and tied in one, GSK3B rep2. The route-only B−C contrast has 20 wins and five losses, net +15/128 or 11.71875 percentage points, with conditional interval [5.46875, 17.96875].

Against original outputs, B has 47 RF-joint wins and zero losses, net +36.71875 percentage points with conditional interval [28.125, 42.96875]. The external-0.5 joint contrast has 33 wins and zero losses, net +25.78125 percentage points with interval [21.09375, 30.46875]. The route-only contrast has 38 wins and one loss, net +28.90625 percentage points with interval [20.3125, 36.71875]. Complete slot expressions, paired wins/losses, task strata and all eight run strata are retained in the analysis JSON.

Intervals are supplied only because all expected queries have terminal, known endpoints. The analysis resamples the four complete training-run differences within each fixed task with replacement, gives the two task means equal weight and enumerates all `4^4 × 4^4 = 65,536` combinations, using linear empirical percentile interpolation. A replicate here contributes its 16 selected positions, not a new 16-position training run. The intervals describe conditional variability among the existing runs under this saved-output sensitivity. They are not independent-source validation or multiplicity-adjusted confirmation. Had any query remained unknown, count and paired-difference bounds would have retained one shared binary route variable per canonical query; unknown references would not have been treated as independent arm outcomes, and these intervals would not have been supplied.

### S16.4. Reuse lookup and observed computational work

Before new search, the recorded historical lookup scanned 8,713 files in the existing remote artifact roots. It found 125 matching Graph target records and zero exact-compatible AiZ query records. A matching canonical target alone did not establish compatibility across planner, protocol and seed. The completed progress record therefore reports 232 new searches and zero reused AiZ searches. Reuse of the original sources, final outputs and RF/external measurements is distinct from reuse of a historical AiZ route query; the former occurred and the latter did not.

Each distinct canonical query is counted once in physical work. Total recorded target-search time is 458.843 seconds. Expansion inference accounts for 7,241 attempted and successful calls, 8,440 input rows and 35.637 seconds; filter inference accounts for 28,434 attempted and successful calls, 28,434 rows and 18.412 seconds. No prediction call failed. There are 210 expansion and 19 filter budget denials, which do not enter attempted-call totals. Model timings are components of the search time and are not added again to it. The progress record separately reports 12.032 seconds for resident model loading and 639.301 seconds for the controller including history lookup and loading; the enclosing foreground execution was 640.754 seconds. The subsequent independent inventory check took 28.040 seconds and is separate from that search-controller timer.

For logical arm references, saved target-search time sums are 273.433 seconds for original, 244.637 for B and 241.464 for C. Their expansion call/row counts are 3,773/4,542, 3,761/4,409 and 3,777/4,469; filter calls/rows are 15,689/15,689, 15,522/15,522 and 15,184/15,184. These sums include a shared query at each logical reference and must not be added to obtain physical work. The analyzer keeps new-search, compatible-historical and missing-provenance work separate; all recorded query work in this run is in the new-search category. Missing costs would remain unobserved rather than zero. This analysis adds zero RF and zero external-model inference and does not reexecute the source training or candidate-generation stages.

### S16.5. Interpretation and reproduction boundary

The observed B−C advantage persists for both fixed tasks and the RF and external joint endpoints under the recorded AiZ configuration. This is positive evidence that the saved contribution-selected outputs' comparative pattern is not confined to their original Graph route labels. The scope remains a post hoc 128-position subset of the same eight New512 source runs, using fixed output structures originally selected by the Graph workflows. It does not show that an AiZ-guided 6+6 replacement policy would select these same outputs, establish an optimal policy, prove planner-independent behavior or establish equal cost between AiZ and Graph. The different route counts and depth units are reported directly; successful routes are not combined retrospectively. Computed routes and predicted task scores remain computational endpoints rather than measured synthesis or activity.

**Evidence and reproduction:** [fixed-input preparation and source/arm identities](G:/逆合成结构化反馈/second_planner128_v1/inputs.json), [actual executed protocol](G:/逆合成结构化反馈/second_planner128_v1/second_planner128_actual_protocol.json), [historical query lookup](G:/逆合成结构化反馈/second_planner128_v1/second_planner128_history_lookup.json), [completed controller progress](G:/逆合成结构化反馈/second_planner128_v1/second_planner128_progress.json), [complete analysis JSON](G:/逆合成结构化反馈/second_planner128_v1/second_planner128_analysis.json), [analysis reading report](G:/逆合成结构化反馈/second_planner128_v1/second_planner128_analysis.md), [independent identity/inventory check](G:/逆合成结构化反馈/second_planner128_v1/independent_stock_identity_check.json) and [standard-library analysis code](G:/逆合成结构化反馈/second_planner128_analyze.py). The preparation-state label in `inputs.json` is a preparation snapshot; completion is established by the terminal query records, analysis and progress. The run's `executed_code` snapshot preserves the runner/auditor/meter code that produced the scientific queries. Subsequent local strengthening of runner exception classification and resume checks is not substituted as the code executed in this completed run.


A separately executed saved-tree checker, without importing the original route auditor, found all 232 queries and 1,481 reported verified routes consistent. A later incremental check covered 19,869 trace nodes and 95,059 molecule references, including strict state/stock fields and reconstruction of the resident cumulative model counters; all queries remained consistent. It took 1.222 seconds, with zero tree rereads, route reaudit, model calls or new queries. Its first check had incorrectly compared cumulative absolute budget-denial counters with per-query increments; the initial report and corrected source are retained, and no scientific endpoint was changed. These checks and their actual source versions are separately saved in [evidence-check documentation](G:/逆合成结构化反馈/second_planner128_v1/evidence_checks/README.md).

The complete scientific archive has been transferred locally (143,165,518 bytes), fully decompressed for integrity inspection, and safely extracted. It contains 714 members, including all 232 result, trace and tree triplets; the [archive inspection record](G:/逆合成结构化反馈/second_planner128_v1/archive_inventory_20261003.json) supplies the observed counts. Later independent check reports and code are separate additions rather than silently treated as members of that original archive. The local working and archived actual input JSON have zero semantic differences; their byte lengths differ only through serialization.

## S17. Task-free precursor-complexity control on the unchanged New512 sources

### S17.1. Question and matched intervention

The complexity control tests whether B's advantage can be reproduced by preferentially using simpler stocked precursors without source-task information. It is an adaptive, same-source extension of New512, not a new generator cohort. The complete 512-position denominator and the original negative Mix−B primary comparison are retained.

S ranks the full reaction-compatible precursor inventory by the SA method of Ertl and Schuffenhauer (2009), implemented in RDKit Contrib `SA_Score.calculateScore`, ascending, followed by molecular weight ascending and canonical SMILES ascending. The ranking uses precursor properties before product formation or candidate RF lookup. It does not use source contributions, candidate RF values, product complexity, or future route outcomes. Inventory descriptor observations are shared across source roles rather than recomputed per source.

All other components match B: joint-failure trigger, the same complete saved replay roots and inventory, preserved stocked roles, up to four alternatives for one missing role or two per role for two missing roles, rank-then-root proposal order, sixteen forward combinations, eight QED-qualified scoring slots, source-relative RF preservation, and the existing 6+6/shared-45-second route policy. A solved candidate replaces the source only after the same native route audit. B's completed outcomes are read as the reference; they are not regenerated.

### S17.2. Completed outcomes and interpretation

**Table S13. Complexity-control outcomes on 512 unchanged raw positions.** All endpoints are known. External qualification uses the unchanged ten-model mean at least 0.5 in place of RF qualification.

| Endpoint | Saved B6+6 | S6+6 | B minus S |
|---|---:|---:|---:|
| RF qualification and audited route | 311 | 83 | +228 |
| Audited route | 361 | 156 | +205 |
| External qualification and audited route | 176 | 35 | +141 |
| Distinct RF-joint task/structures | 310 | 73 | +237 |
| Distinct external-joint task/structures | 175 | 29 | +146 |
| Replaced raw positions | 259 | 31 | +228 |

B−S RF joint yield is +44.531 percentage points, with the conditional run-bootstrap interval [39.258, 49.414]; the external difference is +27.539 points [24.609, 30.469]. Every one of the eight completed source runs has a positive RF and external difference. Uncertainty uses complete paired runs within each fixed task, as in S8, rather than independently resampling molecules.

S has an RF/QED-qualified candidate on 31 source positions, compared with 286 under B. Its lower final yield is therefore accompanied by lower task-preserving candidate coverage under this specific precursor-complexity ranking. This supports the usefulness of task information relative to the tested task-free alternative. It does not exclude all complexity-related explanations, establish superiority to every synthetic-accessibility algorithm, or establish a total-cost advantage: S performs much less candidate route work because it seldom supplies qualified candidates.

The saved candidate stage adds 1,446 forward applications and 992 RF predictions, with no new training, terminal sampling, root discovery or source attribution. S's logical candidate route work is 35 model calls and 56.208 summed recorded search seconds; B's already measured reference carries 812 calls and 1,619.591 seconds. These are separate policy work views, not a claim that B is cheaper.

**Evidence:** [prospective extension design](G:/逆合成结构化反馈/complexity_control512_v1/design.json), [completed RF/route analysis](G:/逆合成结构化反馈/complexity_control512_v1/archive_contents/complexity_control512_v1/analysis.json), [retained candidate stage](G:/逆合成结构化反馈/complexity_control512_v1/archive_contents/complexity_control512_v1/candidate_stage.json), and [completed fixed-output external analysis](G:/逆合成结构化反馈/complexity_control512_v1/external_fixed_outputs_v1/analysis.json).

## S18. New HIV source runs with role-separated guiding and evaluation predictors

### S18.1. Assay definition, molecular inputs and four data roles

This extension asks whether the contribution mechanism transfers to a non-kinase prediction task when the guiding and evaluation models have disclosed, nonoverlapping training roles. It uses the retained MoleculeNet HIV table (Wu et al., 2018) from the official DeepChem loader's object URL. CA and CM are positive, and CI is negative. The endpoint concerns predicted cellular HIV-replication screen classification; it is not a protein-affinity assay or a measurement of generated-compound activity.

The retained table has 41,127 valid unique canonical rows, comprising 1,443 positives and 39,684 negatives. Atom maps are cleared, but no salt stripping, tautomer normalization or largest-fragment relabeling is applied. Full labeled multicomponent molecules remain in predictor development; a label is never assigned to an extracted fragment. Generated and replacement structures still require one connected component. Predictor assessment reports the full test distribution and its single-component stratum.

Whole achiral Murcko scaffolds define initial groups, with all empty scaffolds in one group. Groups sharing identical full guiding-RF or evaluation-LR fingerprint vectors are merged, using the actual vectors rather than a new digest. There are 19,089 initial scaffold groups and 18,950 merged groups. One seeded allocation orders groups by descending size, with a shuffled tie order, and balances row/positive/negative counts toward 35%/35%/15%/15%. Split seed 2026100301 is fixed without performance-based seed selection. No retained scaffold group or identical full predictor vector crosses roles.

**Table S14. Actual HIV predictor-development roles.**

| Role | Rows | Positive | Negative | Merged groups |
|---|---:|---:|---:|---:|
| Guiding RF training | 14,394 | 505 | 13,889 | 7,986 |
| Evaluation LR training | 14,395 | 505 | 13,890 | 7,987 |
| Threshold calibration | 6,169 | 216 | 5,953 | 1,487 |
| Untouched test | 6,169 | 217 | 5,952 | 1,490 |

The RF uses achiral radius-2/2,048-bit binary Morgan inputs, 200 trees, maximum depth 20, minimum leaf size 2, square-root feature sampling, bootstrap sampling, balanced-subsample class weights, seed 2026100302 and eight CPU jobs. Its raw class-one probability supplies source optimization, contributions, qualification and ordering. The separately fitted evaluation LR uses achiral radius-3/2,048-bit inputs, L2 penalty, C=1, liblinear, balanced class weights, 2,000 maximum iterations and seed 2026100303. It converged in nine iterations. No hyperparameter search or probability-calibrator fit is performed.

For each fixed model, the calibration role selects a score threshold with empirical false-positive rate at most 0.05, maximizing true positives and breaking ties by fewer false positives and then a higher threshold. Distinct calibration scores and an above-range zero-positive sentinel are considered. The chosen RF and LR thresholds are 0.3068631933634494 and 0.42795301259064084. Both thresholds were saved before untouched-test prediction and new molecular generation. Threshold calibration does not establish that raw scores are calibrated probabilities.

On the untouched natural-distribution test set, RF/LR ROC-AUC is 0.632902/0.667869 and average precision is 0.123897/0.128222; prevalence is 217/6,169 (3.518%). RF/LR false-positive rates are 0.044523/0.063340. In particular, LR's test FPR exceeds its calibration FPR of 0.048883: the calibration constraint is not a guaranteed population-level FPR. Single-component test ROC-AUC is 0.607293/0.648691 and average precision 0.113329/0.107765. Performance is modest and is reported as a limitation, without refitting or loosening thresholds from test or generated-output results. The two models share the original assay definition and potential label bias despite disjoint training roles. Prior-training composition also remains unresolved.

### S18.2. Six new runs and nested permutation controls

Six independent REINVENT training runs use the same existing prior, 128 steps, batch size 16, DAP, sigma 128 and learning rate 0.0001, with the HIV RF as the sole task reward and no route feedback. Training seeds are 10305001–10305006; separate terminal sampling seeds are 10305101–10305106. Each terminal checkpoint supplies 64 raw draws with `unique_molecules=false` and randomized SMILES. All 384 positions are retained: 380 valid unique structures and four invalid outputs. No old kinase, TreeSHAP or development cohort is pooled with them.

Original, Morgan A, contribution B and three permutation realizations C0/C1/C2 use the same inventory, templates, source-relative QED condition, guiding threshold and bounded route policy. Each C realization is a complete separately budgeted policy; the comparison uses their arithmetic mean within each source/run. Their successes are neither unioned nor used to choose the best realization, and the realizations are not counted as eighteen independent generator runs. Permutation seeds add `realization_id*100000000` to `10304000 + task_index*1000000 + replicate*10000 + position*100 + root_index*10 + group_index`, with HIV task index 2. Parent planning seeds use `10302000 + task_index*1000000 + replicate*10000 + position`.

The primary endpoint is LR threshold qualification and audited route over the original 384 positions. It replaces guiding RF qualification in the endpoint, rather than adding a third RF requirement. LR never determines generation reward, trigger, precursor ranking, product qualification, final selection or retries. The primary contrast is B minus mean(C0,C1,C2). Descriptive uncertainty enumerates all 6^6=46,656 resamples of complete paired source runs.

### S18.3. Completed positive secondary results and null primary transfer

**Table S15. HIV outcomes with denominator 384 in every arm.** All 2,304 final rows have known endpoints.

| Endpoint | Original | A | B | C0 | C1 | C2 |
|---|---:|---:|---:|---:|---:|---:|
| Guiding RF qualification and audited route | 4 | 38 | 50 | 27 | 31 | 30 |
| Audited route | 168 | 187 | 191 | 181 | 182 | 183 |
| Evaluation LR qualification and audited route | 7 | 7 | 8 | 8 | 8 | 9 |
| Replaced | 0 | 34 | 46 | 23 | 27 | 26 |

The primary B−mean(C) difference is −1/3 of a position, or −0.086806 percentage points, with conditional interval [−0.781250, +0.868056]. One source run is positive, three negative and two tied. Thus this extension does not show a positive LR-evaluated transfer effect; its RF result must not replace that primary conclusion. The secondary RF difference is +20⅔ positions (+5.381944 points), positive in all six runs; the secondary route difference is +9 (+2.343750 points), positive in five and negative in one.

B has 48 distinct RF-joint structures among 50 successful positions and eight distinct LR-joint structures among eight. Forty-seven of the fifty RF-joint positions fail LR qualification. The two known exact labeled matches in all arms are retained LR-training negatives; removing their success contribution while keeping denominator 384 does not change the LR counts or primary contrast. Exact-match sensitivity is distinct from scaffold/nearest-training-similarity descriptions, which do not select or relabel outputs. The generated structures need not be scaffold-disjoint from model-development roles merely because the roles themselves are disjoint. Two reference strings that do not round-trip identically on reparsing are preserved and disclosed; no alias or role reassignment is introduced.

The result identifies predictor disagreement as an important limitation of the observed RF improvement. It is neither evidence of measured biological inactivity nor proof that future task-guided replacement cannot improve an independent endpoint. Any changed evaluator, threshold or candidate rule would be a separately planned intervention, not a reinterpretation of this completed primary test.

**Evidence:** [data audit](G:/逆合成结构化反馈/hiv_transfer_assets_20261003/asset_audit.json), [fixed model configuration](G:/逆合成结构化反馈/hiv_transfer_assets_20261003/models_v1/preparation_config.json), [role assignment](G:/逆合成结构化反馈/hiv_transfer_assets_20261003/models_v1/role_assignment.json), [measured models and untouched test](G:/逆合成结构化反馈/hiv_transfer_assets_20261003/models_v1/model_measurements.json), [prospective design preserved in candidate configuration](G:/逆合成结构化反馈/hiv_transfer384_role_separated_predictors_candidate_v1/candidate_configuration.json), [complete generation](G:/逆合成结构化反馈/hiv_transfer384_role_separated_predictors_candidate_v1/generation_complete.json), and [completed external analysis and overlap diagnostics](G:/逆合成结构化反馈/hiv_transfer384_role_separated_predictors_candidate_v1/external_analysis/analysis.json). Historical preparation-state labels describe when a configuration was written; terminal observations establish current completion.

## S19. Task reranking within retained native candidate pools

### S19.1. Separating task information from candidate-pool formation

R is an exploratory matched extension on the same New512 source positions. It uses the full saved pools from 431 completed native depth3/k10 enumerations, totaling 87,679 outputs. The other 73 positions retain their original sources, and eight native execution limitations remain scientific unknowns. Native enumeration is not rerun. The purpose is to test task information within the retained native opportunity set, rather than to reproduce the original SynTwins MPO/GA pipeline.

R computes the same positive source-environment weights as B, then ranks native products by retained positive source weight, saved native similarity and canonical SMILES. Candidate RF values do not enter this initial ranking; all plans are saved before RF-history lookup. The first eight enrolled products use the existing source-relative RF/QED qualification and the common 6+6/shared-45-second route rule. Original native, B and C outputs are retained references.

### S19.2. Outcomes, uncertainty and incurred work

**Table S16. Native-pool task reranking on the unchanged 512 positions.** Bounds preserve the eight terminal native limitations. The route lower bound can include an already known fallback success even when the final intervention is unknown.

| Endpoint | Native similarity ranking | Task reranking R | Saved B |
|---|---:|---:|---:|
| RF qualification and audited route | 153–161 | 193–201 | 311 |
| Audited-route bounds | 218–225 | 254–261 | 361 |
| External qualification and audited route | 58–66 | 90–98 | 176 |
| Logical candidate-planning model calls | 939 | 1,267 | 812 |
| Summed recorded candidate search seconds | 1,670.700 | 2,307.850 | 1,619.591 |

R has 253 directly known route successes; its conservative lower bound of 254 additionally accounts for known fallback information at an uncertain intervention. These are not conflicting totals. Paired external R−native bounds are +24 to +40 positions, accounting for the shared uncertain positions rather than subtracting only point estimates. B remains ahead of R by 78–86 external successes. Unknowns are not counted as chemical failures, and no fully determined run-bootstrap interval is supplied for these unresolved endpoints.

R adds 431 source decompositions and 1,728 RF predictions, with zero new native enumeration, forward replay, generator training or terminal sampling. It executes 224 new Graph queries, 951 new model calls and 1,734.086 new search seconds; 72 compatible historical queries carry 315 saved calls and 572.258 saved search seconds. The logical policy total is therefore greater than the saved native comparator's total. An additional 1,070 checkpoint/structure vectors were acquired for the fixed-output external evaluation; later analysis-only invocations reporting zero new predictions do not erase that earlier work.

All 140 distinct successful R query witnesses were subsequently loaded and audited without new model inference, and 72 physically reused graphs were exported with their explicit mappings. Declared-origin recovery retained exact requests and saved seeds; it did not repeat the native pools or substitute newly favorable searches. The positive result supports task information within this native pool at additional scoring/search work, rather than cost-free improvement or superiority to the full original native optimizer.

**Evidence:** [R design](G:/逆合成结构化反馈/syntwins_task_rerank512_v1/design.json), [candidate stage](G:/逆合成结构化反馈/syntwins_task_rerank512_v1/candidate_stage.json), [RF/route analysis](G:/逆合成结构化反馈/syntwins_task_rerank512_v1/analysis.json), [fixed-output external analysis](G:/逆合成结构化反馈/syntwins_task_rerank512_v1/external_fixed_outputs_v1/analysis.json), [declared-origin recovery](G:/逆合成结构化反馈/syntwins_task_rerank512_v1/declared_alias_recovery_20261004.json), and [all successful witness reaudits](G:/逆合成结构化反馈/syntwins_task_rerank512_v1/saved_route_reaudit_20261004.json).

## S20. Actual method execution, relocated analysis and private reproduction materials

The v0.2.0 method package exposes B/C execution with explicit RF and Graph interpreters and model/template/stock inputs. Its CLI defaults to fixed-task qualification, using the exact validated adapter from the prospective study; `--qualification source_preserving` explicitly reproduces historical strict qualification. The original chemistry and scoring kernels remain unchanged. Standard-library decision replay uses saved candidate/query observations and does not deserialize route graphs or execute models. Historical strict software fixtures remain execution evidence, not fresh demonstrations of fixed-mode efficacy.

An earlier selected A6000 fixture executed source preparation, candidate chemistry/RF, parent/candidate searches and saved-graph auditing. It recorded seventeen RF predictions, thirteen **candidate-formation** forward applications and four Graph queries totaling twenty-five calls in 87.671 controller seconds. Its legacy forward aggregate was candidate-only; it must not be described as complete chemistry expenditure. A separate two-source history execution added zero RF/forward/Graph calls but read and audited five saved queries, carrying thirty-seven historical calls and 73.250 historical search seconds. Its 93.581-second controller interval includes audit/retrieval overhead. These inputs differ, so the durations do not measure a fresh-versus-cached speedup.

Clean RF/chemistry and Graph environments were subsequently installed on the existing A6000 host. The RF environment supports the required chemistry, fingerprints and retained sklearn forest; the separate Graph environment supports the native planner without importing the incompatible RF runtime. Using retained assets and no history adapter, the same existing software fixture completed native source scoring, formation and planning in 86.316 controller seconds. It executed one source and sixteen candidate RF predictions; 121 source-replay plus thirteen candidate-formation forward invocations, giving **134** complete forward invocations; 1,071 reverse-template invocations; one source RF-path attribution; and four Graph queries using twenty-five calls and 56.090 resident search seconds. The thirteen-only legacy field is retained separately. The stored source cost counters and call-site mapping independently establish the complete 134 total; missing legacy source counters elsewhere remain unknown rather than silently filled from this fixture.

The clean fixture's source RF was 0.34. B's two bounded candidate searches did not yield an audited solution, while C's selected candidate did. These observations demonstrate branch coverage and actual model execution; this pre-existing selected fixture is excluded from all scientific-effect denominators. Its independent saved-output audit constructed no models and repeated no RF, reaction or route calculation. Clean installation on an existing asset host does not certify a completely independent new-machine scientific regeneration.

A caller-relative path defect was reproduced and repaired without changing scientific settings. Explicit input/asset/output paths are resolved against the caller before changing child-process directories; bare interpreter names retain PATH semantics. Ordinary software tests exercise real controller/subprocess boundaries with injected workers. Subsequent cost-reporting corrections expose source replay, candidate formation, reverse templates and attribution separately, while preserving the science kernels and candidate-only legacy field. Relocated saved-counter analysis reproduced 134 = 121 + 13 with site packages disabled, without scientific module imports. The private updated package preserves the previously executed code/receipts rather than advertising packaging validation as a new model experiment.

Private method, applicable retained assets and saved-output analysis materials are separate archives. The earlier analysis archive contains 4,608 kinase external rows, 2,304 HIV external rows, ten two-task vectors for 1,332 structures and 12,338 retained HIV test predictions. Its standard-library commands reproduce the corresponding kinase/HIV/native-reranking counts, paired contrasts and held-out predictor metrics. Both commands and two decision replays completed after extraction to different paths containing spaces and Chinese characters with Python 3.12.14 and site packages disabled, without count differences or inference. Earlier GraphGA, later costs and prospective GraphGA1024 are retained as separate additions; the older archive is not claimed to contain every later supplementary comparison.

Dependencies and asset permissions remain explicit. The tested RF deposit declares CC0, Graph2Edits checkpoint deposit MIT and PaRoutes stock CC BY 4.0. The REINVENT v1 public record carries Apache 2.0, but the retained prior's original acquisition and exact deposited correspondence remain unresolved. SynTwins/template and CombiMOTS weight redistribution permission, ultimate historical HIV mapping and detailed predictor/prior training composition remain limited. Upstream acquisition instructions do not authorize mirroring uncertain assets. No private ZIP is a public release, archival DOI or final availability statement.

**Evidence:** [earlier selected execution](G:/逆合成结构化反馈/method_roundtrip_branch_20261003_completed_20261004.tar.gz), [earlier relocated package receipt](G:/逆合成结构化反馈/submission_preparation_20261004/private_package_execution_report.json), [clean installation and complete-work receipt](G:/逆合成结构化反馈/submission_preparation_20261004/software_setup/clean_fresh_one_source_software_local_receipt_20261005.json), [independent fixture audit](G:/逆合成结构化反馈/submission_preparation_20261004/software_setup/clean_fresh_one_source_software_receipts_20261005/fresh_one_source_candidate_20261005/independent_saved_output_audit_20261005.json), [relocated complete-cost reproduction](G:/逆合成结构化反馈/submission_preparation_20261004/software_setup/candidate_package_validation_20261006/relocated_saved_CPU_costs_20261006.json), [private updated package verification](G:/逆合成结构化反馈/submission_preparation_20261004/software_setup/candidate_package_validation_20261006/private_candidate_package_validation_20261006.json), [analysis quickstart](G:/逆合成结构化反馈/submission_preparation_20261004/analysis_reproduction/README.md), [portable fresh quickstart](G:/逆合成结构化反馈/submission_preparation_20261004/software_setup/portable_fresh_quickstart.md), and [asset access/redistribution record](G:/逆合成结构化反馈/submission_preparation_20261004/asset_access_and_redistribution.md).

The new 16-run saved-analysis package is separate from the earlier archive. Actual extraction and execution from a different directory containing Chinese characters and spaces, with site packages disabled, reconstructed all 5,120 policy rows and 13,757 scalar fields. No substantive differences were found; sixteen final-bit floating differences were at most 1.28 × 10⁻¹⁵. Eight malformed-input tests passed without scientific imports. The canonical software suite passed 73 tests, with one expected platform skip; the retained fixed B/C decision replay matched without new scientific calls. These checks establish implementation and saved-data portability, not additional efficacy replication.

OpenAI Codex assisted with substantive software implementation and debugging, inspection of saved experimental records, arithmetic and source reconciliation, and editorial drafting of the English manuscript and Supplement. This assistance included generating prose from the recorded computational evidence. The resulting text remains a working draft for researcher rewriting and review; scientific interpretation, authorship and final factual verification remain the researchers' responsibility. Disclosure of this assistance does not itself establish compliance with a journal's authorship or manuscript-generation policy.

## S21. Earlier prospective transfer under source-relative qualification

### S21.1. Prospective population, operators and source sampling

This earlier block tested a second generator under the strict B/C method, with its design recorded before generation and route outcomes. All eight specified optimizer runs and all 512 terminal positions were completed without task omission, outcome-based resampling, threshold adjustment or evaluator replacement during that comparison. Its B−C contrast remains the original primary for that study. Subsequent fixed-qualification development used these same eight runs (S24); the independent qualification validation used sixteen separately generated runs (S23). Neither development reuse nor a new interpretation enlarges the earlier primary sample. The older excluded eight GraphGA molecules remain separate historical evidence.

The traced GraphGA uses GuacaMol graph crossover/mutation operators (Jensen, 2019; Brown et al., 2019) at commit `ae43219c89d5db134028336243f508606d81995e`. Each run has population 100, mating pool 200, one hundred offspring attempts per generation, 250 generations, mutation rate 0.01, no patience stop and two CPU jobs. Fitness is only the fixed task RF class-one probability, with invalid/disconnected inputs scored zero. There is no quality multiplier, constant synthesis feedback, route reward or external Chemprop reward.

Initialization independently permutes row positions in the retained GuacaMol v1 full corpus and accepts the first one hundred valid connected unique canonical structures, without RF, QED or route filtering. The corpus has 1,591,378 rows. Its official access record is GuacaMol All SMILES, Figshare 7322252 v2, file 13612745, DOI 10.6084/m9.figshare.7322252.v2. The processed deposit declares CC BY 4.0, while its underlying ChEMBL 24/24.1 release declares CC BY-SA 3.0; these data declarations are separate from GuacaMol/GraphGA software's MIT license. Filename and byte size agree with the deposit, but the original download receipt and exact deposited-binary identity were not verified. The full corpus includes its published training, validation and test partitions. New optimizer seeds do not establish exclusion from guiding RF, external Chemprop or prior training data.

All 25,000 offspring attempts per run are retained in lineage records. Each complete ordered terminal population of one hundred is saved before sixty-four uniform draws with replacement using a separate sampling seed. These are terminal draws, not mutation proposals, RF-top-k exports or post-filtered successful molecules. Duplicates remain separate positions.

**Table S17. GraphGA optimizer, initialization, sampling and source search seeds.** Search seed for a position is the listed base plus its position 0–63.

| Cohort | Optimizer seed | Initialization seed | Terminal sampling seed | Search seed base |
|---|---:|---:|---:|---:|
| GSK3B rep0 | 2026100401 | 2026100501 | 2026100601 | 10404000 |
| JNK3 rep0 | 2026100402 | 2026100502 | 2026100602 | 11404000 |
| GSK3B rep1 | 2026100403 | 2026100503 | 2026100603 | 10414000 |
| JNK3 rep1 | 2026100404 | 2026100504 | 2026100604 | 11414000 |
| GSK3B rep2 | 2026100405 | 2026100505 | 2026100605 | 10424000 |
| JNK3 rep2 | 2026100406 | 2026100506 | 2026100606 | 11424000 |
| GSK3B rep3 | 2026100407 | 2026100507 | 2026100607 | 10434000 |
| JNK3 rep3 | 2026100408 | 2026100508 | 2026100608 | 11434000 |

B and C retain the same guiding cutoffs, source-relative RF/QED qualification, stock, templates, sixteen forward/eight scoring-slot ceilings and parent12/45/candidate6+6/shared45 planning configuration. Their matched cohort uses the same GPU. Two queues distribute complete runs across the two authorized A6000 cards; no molecule-level pooling defines replication. The original method package's historical phrase “selected-source software round trip” is not interpreted as this experiment's scope: the controller records the new prospective GraphGA cohort explicitly.

All Original/B/C final structures were persisted before external prediction. The unchanged released ten two-task Chemprop models evaluate the fixed outputs; the primary endpoint is mean task probability at least 0.5 and audited route, without an added RF requirement. Its B−C contrast was recorded before outcomes. Conditional uncertainty enumerates all 65,536 ordered paired-run resamples, four complete optimizer runs within each fixed task with equal task weight. It is descriptive uncertainty for these tasks/assets, not a randomization test, a power guarantee or measured biological validation.

### S21.2. Completed joint yield, task heterogeneity and distinct outputs

All 512 draws are valid and RF-qualified. They contain 375 distinct task/structures (191 GSK3B, 184 JNK3), and none is shared between optimizer runs. The source RF min/median/max is 0.78/0.915/1.00. Sixty-four draws in GSK3B replicate 0 are exactly RF=1; all one hundred terminal population members of that run also score one. The remaining runs are not saturated. Overall source QED min/median/max is approximately 0.0430/0.1982/0.6811. These descriptors and predictor scores do not establish biological activity or failure.

**Table S18. GraphGA primary external/route outcomes.** RF/route and route-only counts happen to be identical because every route-positive selected output also qualifies under its guiding and external task-predictor criteria. This agreement is observed, not an endpoint redefinition. All 1,536 final rows are complete with zero unknowns.

| Task | Original raw (distinct) | B raw (distinct) | C raw (distinct) | B−C raw |
|---|---:|---:|---:|---:|
| GSK3B | 2 (1) | 98 (9) | 33 (13) | +65 |
| JNK3 | 67 (48) | 67 (48) | 67 (48) | 0 |
| Total | 69 (49) | 165 (57) | 100 (61) | +65 |

Primary raw rates are 13.477%, 32.227% and 19.531% for Original, B and C. B−C is +12.695 percentage points, with conditional interval [0.391, 27.148]; matched positions show seventy-nine wins and fourteen losses. Run differences are [0,+46,+18,+1] for GSK3B and [0,0,0,0] for JNK3, yielding three positive and five tied runs. The entire gain occurs in GSK3B. The saturated GSK3B run and all JNK3 runs do not improve.

B's raw success advantage does not extend to distinct discovery diversity: it has 57 distinct joint successes compared with C's 61 and 48 versus 52 nonempty successful scaffolds. B's ninety-six GSK3B replacements comprise only eight distinct products; C's thirty-one comprise twelve. The largest B product occupies fifty-nine positions from forty-seven distinct sources across two runs; the second occupies twenty-one from seventeen sources across two runs. These two products account for 80/98 (81.6%) of GSK3B B successes and 80/96 (83.3%) of its replacements. Maximum original-source multiplicity is five, so input repeats alone do not explain the observed convergence. The primary denominator remains raw original positions, and distinct counts are reported alongside it rather than substituted after outcomes.

### S21.3. Qualification and chemical-alert diagnostics

GSK3B B/C have 1,976/2,009 scored slots, of which 120/36 preserve the source RF score and QED condition. They give at least one qualified candidate on 106/33 sources and ninety-six/thirty-one solved replacements. Ten B and two C sources with eligible candidates exhaust their bounded candidate searches in replicate 2; they are known within-budget nonsolutions. In GSK3B replicate 0, every scored candidate falls below the saturated source score and neither arm replaces a source.

JNK3 B/C have 1,434/1,426 scored slots. Of these, 1,417/1,355 satisfy the fixed task cutoff, but none maintains its source RF score. Every B triggered source and all but one C triggered source have scored candidates; all scored slots in both tasks pass the source-relative QED condition. The best candidate-minus-source RF gap remains at most −0.23 for B and −0.17 for C over all JNK3 runs. Qualification, rather than demonstrated absence of feasible candidate routes, prevents replacement here. Rejected candidates were not searched and their route yield is unknown. The saved partial first-JNK3 diagnostic is preserved as history; the complete GSK3B results prevent generalizing that early zero-replacement observation to the whole experiment.

**Table S19. Descriptive screening filters on the unchanged GraphGA external-joint outputs.** These are saved-descriptor sensitivities retaining denominator 512, without candidate reselection or alteration of the primary endpoint.

| Counted external/route successes | Original | B | C |
|---|---:|---:|---:|
| No additional filter | 69 | 165 | 100 |
| No PAINS alert | 46 | 142 | 77 |
| No PAINS and no BRENK alert | 28 | 29 | 29 |
| Molecular weight ≤500 | 37 | 133 | 67 |

The B−C raw gain persists under no-PAINS and molecular-weight filters but disappears under the combined no-PAINS/BRENK condition. Among external-joint successes, B/C have 128/63 BRENK-flagged positions; across all final positions they have 280/248, versus 224 in Original. Mean final QED rises from 0.223458 to 0.285398 under B, while screening-alert concentration worsens. A positive QED or predicted task score is therefore not a substitute for broader screening suitability. No final output matches the available public task-labeled table, but that incomplete table cannot establish original training-set or scaffold exclusion. Neither a computational route nor an absence of these alerts is laboratory validation.

### S21.4. New work, route evidence, recovery and preservation

Generation performs 108,657 newly acquired task RF predictions; its foreground controller interval is 2,105.131 seconds. The completed method stage reuses the retained generation scores for all sources and adds 2,083 candidate RF predictions, 2,574 forward applications and 399 Graph queries totaling 4,364 model calls. The query set comprises 375 parents and twenty-four candidate queries. New sibling-cohort candidate/chemistry observations are not automatically imported across the unchanged cohort readers, so any repeated actual work remains charged. These totals are shared actual acquisition, not an isolated per-arm deployment benchmark.

All sixty-nine solved saved Syntheseus query graphs were subsequently reloaded and audited using the retained stock and audit algorithm, with no model constructors or additional inference. Every graph passed, and all 1,536 final-row query references agreed with their declared molecule and route status. The witnesses comprise forty-nine solved parent and twenty solved candidate queries. Different query records can concern the same final structure, so these numbers are distinct from the successful structure counts in Table S18.

External evaluation requires 391 unique structures including sources, or 3,910 ten-checkpoint structure vectors. Twenty compatible historical vectors are reused, and 3,890 missing vectors are newly computed on CPU. Saved inference batches total 14.025 seconds, excluding model loading, descriptors and historical cost; the successful external workflow invocation lasts 45.943 seconds. It uses Python 3.11.15, Chemprop 1.5.2, Torch 2.4.1+cu121, RDKit 2023.09.6 and NumPy 1.26.4, with GPU visibility empty and OMP/OpenBLAS/MKL thread settings of two. No old source generation, candidate RF or route search is repeated for external recovery.

The initial external invocation stopped before model construction because one historical structure absent from the fixed required set had a fold-vector difference of 1.1920928955078125×10⁻⁷. A read-only audit found zero conflicts among required structures. Recovery preserved full-file validation, existing tolerance and required-structure conflict rejection, while recording and excluding unused conflicts from returned evaluation vectors. It did not choose the more favorable value, loosen thresholds or rerun old inference. The two required structures represented in history retained their compatible vectors.

The full generation/method archive, including lineage, observations and native graphs, and the completed external archive are saved privately, in addition to compact metadata snapshots. A first detached launch that ended before producing scientific observations is retained as an infrastructure attempt; it is not counted as a failed molecular replicate. Foreground execution and terminal records establish completion. These backups and the reaudits support private preservation, not public availability.

**Evidence:** [prospective design](G:/逆合成结构化反馈/graphga_transfer512_v1/remote_snapshot_20261004/design.json), [completed generation](G:/逆合成结构化反馈/graphga_transfer512_v1/remote_snapshot_20261004/generation_complete.json), [completed method summary](G:/逆合成结构化反馈/graphga_transfer512_v1/remote_snapshot_20261004/method_stage_result.json), [full candidate/source diagnostics](G:/逆合成结构化反馈/graphga_transfer512_v1/source_diagnostics/complete_method_diagnostics_20261004.json), [primary analysis](G:/逆合成结构化反馈/graphga_transfer512_v1/remote_snapshot_20261004/external_fixed_outputs_v1/analysis.json), [unchanged final rows and descriptors](G:/逆合成结构化反馈/graphga_transfer512_v1/remote_snapshot_20261004/external_fixed_outputs_v1/enriched_rows.jsonl), [external runtime](G:/逆合成结构化反馈/graphga_transfer512_v1/remote_snapshot_20261004/external_fixed_outputs_v1/external_runtime.json), [all native-witness reaudits](G:/逆合成结构化反馈/graphga_transfer512_v1/witness_reaudit/witness_reaudit_20261004.json), [historical conflict audit](G:/逆合成结构化反馈/graphga_transfer512_v1/external_adapter/historical_conflict_audit_20261004.json), [corpus access/provenance](G:/逆合成结构化反馈/submission_preparation_20261004/graphga_corpus_access_addendum.md), [private generation/method backup](G:/逆合成结构化反馈/graphga_transfer512_v1/graphga_private_generation_method_backup_20261004.zip), and [completed external backup](G:/逆合成结构化反馈/graphga_transfer512_v1/graphga_completed_external_backup_20261004.zip).

## S22. Saved evidence for the first12 versus 6+6 allocation component

The completed New512 experiment already contains the relevant allocation comparison. B_joint12 and B_joint6plus6 use the same joint-failure trigger and candidate pool, carrying identical 2,398 RF score slots, 4,663 forward attempts and 439 source attributions. Both have a ceiling of twelve candidate model calls and one shared forty-five-second soft allowance. The 6+6 policy distributes attempts over at most two qualified structures rather than concentrating all twelve on the first; unused first-candidate calls do not transfer.

**Table S20. Completed New512 allocation results and recorded candidate-route work.**

| Quantity | B first12 | B 6+6 | 6+6 minus first12 |
|---|---:|---:|---:|
| RF qualification and audited route | 282 | 311 | +29 |
| External qualification and audited route | 155 | 176 | +21 |
| Audited route | 335 | 361 | +26 |
| Logical candidate-planning model calls | 1,013 | 812 | −201 (−19.84%) |
| Summed recorded candidate search seconds | 2,002.219 | 1,619.591 | −382.628 (−19.11%) |

Candidate RF/forward/attribution exposure is identical within every one of the eight source runs. All eight runs have positive RF and external allocation differences and fewer logical calls and summed search seconds under 6+6. The RF effect is +5.664 percentage points, with conditional interval [3.906, 7.617]; the external effect is +4.102 points [1.953, 7.031]. These remain recorded secondary component contrasts alongside the negative Mix primary result.

The per-source CSV independently reproduces all thirty-three audited selected cost totals across the seven policies, with maximum floating-point discrepancy 1.592×10⁻¹² seconds. A separate saved-record reanalysis also reproduces the two allocation arms and their effects. No old task or route was repeated to produce this evidence.

Recorded seconds include compatible measured-query work carried into each logical policy. They are not separate uncached end-to-end executions. Query setup, graph saving/audit, unequal target count, shared source acquisition and other incomplete isolated overhead prevent claiming complete deployment speed superiority or equal-total-compute optimality. The supported component conclusion is that, within the matched B pool and common candidate-route ceiling, allocating calls to a second candidate improves joint yield while consuming less recorded candidate-search work.

The earlier T1/T2 evidence is preserved with its original boundary. An old twenty-six-source first12/first6+second6 plan was not executed because compatible measurements existed. The Old256 diagnostic measured eighteen missing second targets using sixty-five new Graph calls and increased RF joint yield from 126 to 135 and external joint yield from 84 to 91. Its first-six prefix lacks a measured pure timer, and its two candidates each had a forty-five-second allowance, so it is not equal-time prospective evidence. The completed New512 comparison directly improves that allocation test without repeating the old sources. An older sixteen-position equal-incremental-action pilot favored continuation by one GSK3B joint success and produced zero edited joint successes; it remains negative development evidence rather than being omitted or repurposed as the current method test.

**Evidence:** [complete New512 costs](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/joint_validation512_v1/costs_20261002.json), [per-source logical work](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/joint_validation512_v1/costs_20261002_policy_sources.csv), [RF allocation contrast](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/joint_validation512_v1/analysis.json), [external allocation contrast](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/joint_validation512_v1/external_fixed_outputs_v1/analysis.json), [saved-cost audit](G:/逆合成结构化反馈/submission_preparation_20261004/saved_cost_evidence_audit_20261004.json), [independent allocation reanalysis](G:/逆合成结构化反馈/submission_preparation_20261004/allocation_cost_independent_reanalysis_20261004.json), [Old256 allocation diagnostic](G:/逆合成结构化反馈/bkunyun_migration_20261001/remote_snapshot/recovered_a6000/development/budget_replay256_v1/analysis.json), and [older equal-incremental pilot](G:/逆合成结构化反馈/fresh_equal_budget_analysis.json).

## S23. Prospective independent validation of fixed-task qualification

### S23.1. Question, enrollment and retained intervention

The prospective question was whether a fixed task cutoff improves externally predicted task/audited-route joint yield relative to source-relative RF qualification within the retained reaction-constrained replacement operation. Candidate development had already revealed a qualification bottleneck and substantial common gains in B and C on earlier sources. The new study specified sixteen independent task-only GraphGA optimizer runs, eight per fixed task, before their generation/outcomes. It retained all sixty-four sampled positions per run, for 1,024 raw positions, and all five output policies. No old optimizer run, favorable source subset, extra source refill, selected checkpoint or post-outcome size extension supplies primary inference.

The generator, corpus initialization and GraphGA operators are those described in S21, using the unchanged traced implementation and upstream commit. Each run uses population 100, mating pool 200, one hundred offspring attempts per generation, 250 generations, mutation rate 0.01, no patience stop and two CPU jobs. Optimization uses only the fixed task RF probability. Initialization accepts the first one hundred valid connected unique canonical structures from an independently seeded permutation of corpus rows, without task, QED or route filtering. Each full terminal population is saved before sixty-four uniform draws with replacement using a separate sampling seed. The 400,000 attempted offspring and 4,000 generation lineage files are retained. They are optimizer proposals, whereas the 1,024 evaluation entries are terminal sampled positions.

**Table S21. Seeds for the sixteen independent qualification-validation runs.** Parent seed for a source position is its listed base plus position 0–63. Exact-compatible reused queries retain requested and adopted seeds separately.

| Task/run | Optimizer seed | Initialization seed | Terminal sampling seed | Parent seed base |
|---|---:|---:|---:|---:|
| GSK3B rep0 | 2026101501 | 2026101601 | 2026101701 | 12500000 |
| JNK3 rep0 | 2026101502 | 2026101602 | 2026101702 | 12510000 |
| GSK3B rep1 | 2026101503 | 2026101603 | 2026101703 | 12520000 |
| JNK3 rep1 | 2026101504 | 2026101604 | 2026101704 | 12530000 |
| GSK3B rep2 | 2026101505 | 2026101605 | 2026101705 | 12540000 |
| JNK3 rep2 | 2026101506 | 2026101606 | 2026101706 | 12550000 |
| GSK3B rep3 | 2026101507 | 2026101607 | 2026101707 | 12560000 |
| JNK3 rep3 | 2026101508 | 2026101608 | 2026101708 | 12570000 |
| GSK3B rep4 | 2026101509 | 2026101609 | 2026101709 | 12580000 |
| JNK3 rep4 | 2026101510 | 2026101610 | 2026101710 | 12590000 |
| GSK3B rep5 | 2026101511 | 2026101611 | 2026101711 | 12600000 |
| JNK3 rep5 | 2026101512 | 2026101612 | 2026101712 | 12610000 |
| GSK3B rep6 | 2026101513 | 2026101613 | 2026101713 | 12620000 |
| JNK3 rep6 | 2026101514 | 2026101614 | 2026101714 | 12630000 |
| GSK3B rep7 | 2026101515 | 2026101615 | 2026101715 | 12640000 |
| JNK3 rep7 | 2026101516 | 2026101616 | 2026101716 | 12650000 |

All source positions are valid. They contain 744 distinct task/structures, with 280 repeated positions: 272 within-run duplicate positions and eight additional cross-run duplicates. Thirteen distinct source structures overlap the earlier GraphGA512 draw set, representing eighteen current positions. These duplicates remain enrolled, and exact compatible chemistry measurements can be reused; neither duplicate positions nor reused measurements become independent repetitions. The inference units are complete optimizer runs conditional on the two retained tasks and assets.

### S23.2. Five matched policies and temporal separation of evaluation

Original retains each source. B strict and C strict retain the source-relative RF condition; B fixed and C fixed use the fixed task cutoff. B uses aligned source-path contributions, while C permutes the same within-role environment weight multisets. The intervention trigger remains valid source failure of the guiding-RF cutoff and/or initial audited route. All four repair arms share the source QED minus 0.05 floor, selected exact-replay roots, template and stock restrictions, source/stock product exclusions and original two-missing-role formation operation. Each source/formation arm has at most sixteen forward proposals and eight unique QED-qualified score slots. Fixed and strict policies read the same completed, saved-before-RF candidate pool for that formation arm.

Candidate ordering is descending guiding RF, descending QED and canonical SMILES tie-breaking. Both qualifications take at most the top two eligible structures. The strict list must be the prefix of the corresponding fixed list; all 2,048 source/arm prefixes were independently checked. Parent planning uses twelve calls and forty-five soft seconds. Candidate planning uses six calls per candidate and one shared forty-five-second soft allowance, with first-solution stopping, nontransfer of unspent calls and target stock exclusion. The independent final route can differ from the candidate's proposing template root.

All 5,120 final policy-position objects were persisted before external model evaluation. The existing released ten-model two-task Chemprop ensemble then evaluated fixed structures, with mean intended-task probability at least 0.5 plus an audited route. **No guiding RF conjunct is added to this external endpoint.** The predictor architecture differs from the optimized RF, but both have a broad ExCAPE origin and incomplete original training composition; this is computational cross-architecture evaluation, not independent biological validation. The external values do not change formation, qualification, ranking, final selection or retries. All five-arm endpoints are known, with zero technical unknown rows.

Source optimization and formation use the retained RF runtime: Python 3.11.15, sklearn 1.2.2 and RDKit 2023.09.6. The external analysis records Python 3.11.15, Chemprop 1.5.2, Torch 2.4.1+cu121, RDKit 2023.09.6 and NumPy 1.26.4, with CPU-only inference. The clean software environments in S20 belong to the separate execution fixture and do not replace these recorded scientific runtimes. Compatible older vectors retain their named runtime evidence and exact structure/model identities.

### S23.3. Prespecified primary inference and complete outcomes

The sole primary is B fixed minus B strict on external/route joint yield. Each optimizer run contributes the paired count difference divided by sixty-four. Effects are averaged within each fixed task and the two task means receive equal weight. The main ninety-five-percent interval and two-sided test use the task-equal paired-run mean with Welch–Satterthwaite degrees of freedom. A prespecified 100,000-resample percentile bootstrap, seed 2026101801, resamples complete paired runs within task with equal task weights as a sensitivity. It does not resample positions, duplicate chemistry or checkpoints as independent experiments. Secondary contrasts are exploratory with unadjusted intervals; no favorable secondary replaces the primary. Technical unknowns would retain identification bounds, and uncertain fallback structures would not become known final-policy endpoints.

**Table S22. Complete independent qualification-validation outputs.** Every count uses the unchanged raw denominator 1,024. Distinct task/structure counts are descriptive unions rather than additional statistical units. RF/route and route-only counts happen to agree here; external/route counts are separate.

| Policy | External/route successes | Guiding-RF/route successes | Audited-route successes | Replaced | Distinct external-success task/structures | All distinct final task/structures |
|---|---:|---:|---:|---:|---:|---:|
| Original | 150 | 169 | 169 | 0 | 110 | 744 |
| B strict | 202 | 221 | 221 | 52 | 113 | 706 |
| C strict | 196 | 215 | 215 | 46 | 113 | 713 |
| B fixed (adopted) | 883 | 913 | 913 | 744 | 222 | 329 |
| C fixed (control) | 880 | 935 | 935 | 766 | 412 | 529 |

B fixed increases external joint successes from 202 to 883, a gain of 681 positions and a task-equal effect of **66.504 percentage points** (95% interval 50.497–82.511; two-sided P=4.577×10⁻⁷; Welch degrees of freedom 13.621). The bootstrap sensitivity interval is 52.051–79.297 points. GSK3B contributes 63.477 points and JNK3 69.531 points; fifteen runs improve and one ties, with no negative primary run. The unchanged Original external count is 150, compared with 883 for B fixed and 880 for C fixed. Route-only successes rise from 169 Original to 913 B fixed and 935 C fixed. The supported objective concerns routes to the final selected identities, not improved routes to the unchanged sources.

**Table S23. External/route counts in every prespecified run.** Each row retains sixty-four original positions, including duplicate chemistry. The last column is the primary paired count difference.

| Task/run | Original | B strict | C strict | B fixed | C fixed | B fixed minus B strict |
|---|---:|---:|---:|---:|---:|---:|
| GSK3B rep0 | 11 | 11 | 11 | 61 | 59 | 50 |
| JNK3 rep0 | 2 | 2 | 2 | 64 | 59 | 62 |
| GSK3B rep1 | 0 | 0 | 0 | 55 | 41 | 55 |
| JNK3 rep1 | 9 | 9 | 9 | 55 | 50 | 46 |
| GSK3B rep2 | 25 | 25 | 25 | 64 | 62 | 39 |
| JNK3 rep2 | 8 | 8 | 8 | 64 | 64 | 56 |
| GSK3B rep3 | 4 | 4 | 4 | 4 | 15 | 0 |
| JNK3 rep3 | 3 | 3 | 3 | 63 | 64 | 60 |
| GSK3B rep4 | 3 | 3 | 3 | 64 | 60 | 61 |
| JNK3 rep4 | 39 | 39 | 39 | 63 | 63 | 24 |
| GSK3B rep5 | 0 | 11 | 19 | 62 | 61 | 51 |
| JNK3 rep5 | 40 | 40 | 40 | 56 | 63 | 16 |
| GSK3B rep6 | 2 | 43 | 29 | 63 | 64 | 20 |
| JNK3 rep6 | 0 | 0 | 0 | 36 | 39 | 36 |
| GSK3B rep7 | 2 | 2 | 2 | 51 | 62 | 49 |
| JNK3 rep7 | 2 | 2 | 2 | 58 | 54 | 56 |

C fixed minus C strict also improves by 66.797 points (exploratory 95% interval 54.212–79.382). B minus C under strict qualification is only +0.586 points (−3.360 to +4.532), and under fixed qualification +0.293 points (−5.312 to +5.898). The qualification-by-formation interaction is −0.293 points (−7.303 to +6.717). These intervals include zero and do not independently establish B superiority; they are not equivalence tests. The earlier strict B−C gain in S21 remains its measured historical result, while its broad attribution-superiority interpretation is limited by this new comparison. The large current qualification gain is shared by aligned and permuted formation and cannot be attributed solely to the alignment mechanism.

### S23.4. Qualification coverage, score loss and output diversity

Current GSK3B source RF scores have mean 0.983906, median 0.98 and range 0.96–1.00; 131 positions score exactly one. Among 465 intervention positions, B has a fixed-qualified pool for all 465 but a strict-qualified pool for only 52; C has corresponding coverage 462 and 46. Current JNK3 source RF has mean 0.871797, median 0.90 and range 0.76–0.98. Among 390 intervention positions, B has 388 fixed-qualified pools and C 387, but neither formation arm has any strict-qualified JNK3 pool. These saved-pool diagnostics describe the operational source-relative ceiling. They neither claim that rejected candidates lack feasible routes nor establish a new biological mechanism.

The adopted fixed mode explicitly gives up source-relative RF nondecrease. Across all 1,024 positions, mean guiding RF changes are −0.182207 for B fixed and −0.234385 for C fixed; RF decreases occur at 692 and 720 positions. Mean external-probability changes are −0.000683 and −0.026173, while mean QED changes are +0.157189 and +0.154285. Among the selected 744 B replacements and 766 C replacements, the RF changes are −0.250780 and −0.313329, external changes −0.000940 and −0.034989, QED changes +0.216346 and +0.206251, and molecular-weight changes −82.094 and −79.920 Da. Replacement-subset values condition on an achieved decision and must not be presented as full-denominator effects. B's smaller observed mean external probability loss is a descriptive tradeoff, not demonstrated general B superiority.

Successful external task/structures increase from 110 Original to 222 B fixed and 412 C fixed; nonempty successful task/scaffolds are 83, 174 and 319. In contrast, the union of **all** final task/structures falls from 744 to 329 and 529. Thus the gain in qualified-route coverage coexists with output convergence, and C retains more distinct successful outputs than B. B's top two successful products occupy 105/883 successes (11.9%), compared with C's 75/880 (8.5%). Input duplicate retention and output diversity are reported separately; successful-position yield is not substituted for distinct chemical discovery.

**Table S24. Structural sensitivities of the fixed final outputs.** Counts retain denominator 1,024 and already selected outputs. Filters are descriptive and do not reselect candidates or define additional primary claims.

| Saved-output filter | Original | B strict | C strict | B fixed | C fixed |
|---|---:|---:|---:|---:|---:|
| No additional filter | 150 | 202 | 196 | 883 | 880 |
| No PAINS | 128 | 180 | 174 | 837 | 840 |
| No PAINS/BRENK | 45 | 45 | 45 | 506 | 469 |
| MW <=500 | 80 | 132 | 126 | 606 | 632 |
| Source QED tolerance | 150 | 202 | 196 | 883 | 880 |

The B fixed versus strict external gain remains under these saved-output filters, including 506 versus 45 with no PAINS/BRENK alerts. Its unfiltered 883 successes do not certify screening suitability or potency, and the different B/C alert and diversity distributions prevent a universal superiority claim. QED is an eligibility descriptor; it is not an activity, toxicity or synthesis measurement.

**Evidence and reproduction:** [prospective design](G:/逆合成结构化反馈/graphga_qualification_independent1024_v1/design.json), [completed generation and terminal sampling](G:/逆合成结构化反馈/graphga_qualification_independent1024_v1/remote_completed_fivearm_snapshot_20261006/generation_complete.json), [shared saved formation/scoring](G:/逆合成结构化反馈/graphga_qualification_independent1024_v1/remote_completed_fivearm_snapshot_20261006/candidate_stage.json), [five-arm results](G:/逆合成结构化反馈/graphga_qualification_independent1024_v1/remote_completed_fivearm_snapshot_20261006/method_result.json), [unchanged external analysis](G:/逆合成结构化反馈/graphga_qualification_independent1024_v1/remote_completed_fivearm_snapshot_20261006/external_fixed_outputs_v1/analysis.json), and [full audit completion](G:/逆合成结构化反馈/graphga_qualification_independent1024_v1/witness_reaudit/all_four_independent_audits_completed_20261006.json).

## S24. Development alternatives, full computation and independent audits

### S24.1. Development decisions without retrospective confirmatory relabeling

Fixed qualification was first assessed on existing GraphGA sources after the strict study revealed high source scores and restricted candidate coverage. The initial 128-position study reused two complete earlier runs; its B/C external joint counts were 122/126 versus 37/37 strict. Extending qualification development to all eight earlier runs gave 462/469 versus 165/100 strict, with unchanged Original 69/512. This common gain motivated a prospective independent qualification comparison; it did not establish aligned B superiority. The old fixed outputs and their optimizer runs are excluded from the new primary analysis, including when their exact scores, chemistry or route measurements supply compatible cache observations.

Two further development ideas were measured rather than hidden. A proposal-allocation rule prioritized products used fewer times earlier in the cohort, then RF/QED/SMILES. Both proposed candidates were reserved before route outcomes, including unused second proposals; it had no outcome-based count reversal or frequency parameter search. It reused the same two earlier 64-position runs and fixed candidate pools. The rule increased distinct successful outputs but reduced external joint yield. A retained-precursor operation used positions 0–15 from each of eight earlier runs. At two-missing-stock roots it combined one original nonstock precursor with a retrieved alternative on the other side, rather than treating the original precursor as stocked or solved. Root choices/rankings and budgets were retained. This also reduced the fixed-mode yield. Neither variant entered the independent recipe or was expanded after its nonpositive result.

**Table S25. Completed same-source development decisions.** Each row is a separate population or formation comparison. The allocation128 and retained128 selections differ and must not be pooled. No confirmatory P value is assigned to these adaptive development panels.

| Development comparison | Raw positions / earlier runs | B external/route count | C external/route count | B/C distinct successful task/structures | Interpretation |
|---|---:|---:|---:|---:|---|
| Earlier strict GraphGA | 512 / 8 | 165 | 100 | 57 / 61 | Original strict comparison retained |
| Fixed qualification on those runs | 512 / 8 | 462 | 469 | 123 / 245 | Shared qualification gain; motivates independent test |
| Fixed qualification controls for allocation | 128 / 2 | 122 | 126 | 32 / 56 | Same saved pools and positions |
| Proposal-count allocation variant | 128 / 2 | 111 | 118 | 52 / 97 | Diversity gain with lower joint yield; not adopted |
| Fixed qualification controls for precursor retention | 128 / 8 | 116 | 116 | 50 / 82 | First sixteen positions per earlier run |
| One-sided retained-nonstock-precursor variant | 128 / 8 | 74 | 65 | 39 / 48 | Joint-yield loss; not adopted |

The fixed-mode decision retains the original RF/QED ordering and two-stock formation operation. It does not combine B and C successes, choose a formation arm according to external results or adopt the negative Mix, fitted TreeSHAP switch, proposal-allocation or retained-precursor branches. Source-relative strict mode remains available for historical reproduction. The implemented recipe has empirical support for qualified-route yield at fixed cutoffs with explicit score/diversity tradeoffs; dependable cross-task activity transfer or attribution optimality has not been established.

### S24.2. Cumulative physical costs and measured timing scopes

Costs distinguish new native work, historical acquisition and logical references. The sixteen source runs execute 204,823 new guiding RF predictions in 4,016 completed batches. Their 756,762 logical RF requests include 551,939 cache consumptions; cache consumptions are not new predictions. Missing-only candidate scoring adds 3,054 unique task/structure predictions in 708 batches, bringing new guiding RF work to 207,877. The 12,961 scored slot references bind 3,349 unique task/structure pairs: 3,054 newly scored and 295 compatible older pairs. Source RF values are retained generation measurements, so zero new source RF in the CPU stage does not mean zero whole-pipeline scoring.

**Table S26. New physical work in the independent validation.** Shared work is counted once, rather than summing five policy outcome references or sixteen copied cohort aggregate reports.

| Stage/component | Newly executed physical work | Historical/cache scope and timer boundary |
|---|---|---|
| Source optimization | 16 runs; 400,000 offspring attempts; 204,823 RF predictions | 551,939 logical score-cache consumptions; controller 5,297.689 s |
| Root discovery and attribution | 660,807 reverse-template calls; 855 source path attributions | 617 new root discoveries ×1,071 reverse templates; 238 exact root-reuse memberships |
| Source-replay forward chemistry | 45,113 native calls | 45,173 logical requests; cumulative difference of 60 is reuse without new execution |
| Candidate forward chemistry | 4,142 native calls | 16,216 logical attempts: 9,044 current-study cache references and 3,030 older-origin references |
| Complete forward chemistry | 49,255 native calls | Source replay plus candidate formation; zero pending observations |
| Candidate RF | 3,054 new unique task/structure predictions | 12,961 score-slot references; 295 globally unique required pairs use older scores |
| Initial parent routes | 731 new queries; 8,273 model calls; 14,562.683 resident search s | 13 old queries carrying 106 calls and 168.251 search s |
| Candidate routes | 512 new queries; 1,433 model calls; 2,888.440 resident search s | 40 old queries carrying 104 calls and 197.821 search s |
| External fixed-output evaluation | 10,690 new checkpoint/structure vectors for 1,069 structures | 490 old vectors for 49 structures; full requirement 1,118 structures ×10 checkpoints |

New Graph work totals 1,243 queries, 9,706 model calls and 17,451.123 resident search seconds. The 53 exactly reused physical origins carry 210 older model calls and 366.072 older search seconds; those are not added to current new acquisition. The route audit inspected all 1,296 physical query records, including copied historical records. Candidate queries and five-policy logical decisions are different denominators. RF and Chemprop checkpoints do not define independent statistical replicates.

The CPU phase spans 70,673.065 calendar seconds from the original controller start to its final completion, including preparation, source replay, formation, RF scoring, persistence and a 874.754-second recovery/transition window. The latest invocation alone records 41,221.430 seconds; the raw `candidate_stage.elapsed_seconds` of 41,220.742 ends before final state/report saves and is not complete cumulative CPU time. Nested plan/root/forward/attribution timers include wrapper/cache/journal/product processing and overlap. Pure native reaction seconds and separate persistence seconds were not independently measured. The sixteen cohort exports copy full-study cost/elapsed fields, which must not be summed. The overall generation-to-external calendar span is 101,291.181 seconds, including interstage readiness, backup and audit gaps; it is not active compute or per-policy deployment latency. No isolated end-to-end speedup is inferred from the shared physical totals.

### S24.3. Recovery, independent auditing and scope of verification

The completed experiment retains native journals before calls, complete source plans and original graph/vector observations. When a saved plan already established completed source formation, continuation consumed that plan rather than rerunning the source. Compact JSON serialization reduced persistence overhead while preserving parsed values, finite-number enforcement, atomic write behavior and observation order. The actual old code, checkpoint and comparison receipts remain available; this infrastructure change is not an algorithmic efficacy intervention. After WSL loss during an initial relaunch, an independently observed hidden Windows process retained the foreground WSL controller until the scientific stages ended. All generation, parent, CPU, candidate-route and external execution stages have terminal successful receipts.

The first saved-CPU audit rejected the provenance label of exact root lists resumed from named current-study plans; those plans and pre-resume completion records supplied the missing origin distinction. The correction changes the auditor's validated origin handling only, with twenty-four model-free regression tests, and preserves the failed report. The native-property auditor separately exposed a borrowed RDKit reactant-query lifetime fault after its reaction owner was released. Retaining query ownership repairs that auditor; three native ownership regression tests passed before the complete recheck. The science kernels, candidate pools, scores, selected outputs and original failed records were preserved. Neither failure was converted into a molecular success, and no old RF/reaction/Graph/Chemprop measurement was repeated to obtain a favorable outcome.

Four full independent checks subsequently completed. The saved-CPU check reconciled all 1,024 sources, 12,961 score references, pending/new work bounds, permutations, strict/fixed prefixes and exact physical origins. The native descriptor check independently recalculated identity, standard stock keys, QED, Morgan environment grouping and complete stock-ranking arithmetic under saved weights: 14,939 QED calculations, 32,733,784 stock substructure tests and 10,710 full side/arm rankings. It repeated no RunReactants, forest prediction, Graph search or external prediction. The route check loaded all 1,296 native graphs and reaudited 514 solved graphs, 5,120 final rows and 2,048 qualification prefixes. The endpoint check validated 11,180 saved fold vectors, all 5,120 endpoints and twenty-one paired-run contrasts against independent arithmetic and statistical references. All four final reports bind to actual exit-zero execution receipts; generation observations and terminal sampling were independently checked beforehand. The final verification completed on 6 October 2026.

These audits verify saved computational evidence and decision semantics. They do not independently recompute forest explanation magnitudes, authenticate all checkpoint content through new inference, establish training-set exclusion, certify reaction prediction accuracy, validate toxicity or demonstrate biological activity/laboratory synthesis. Successful source/target identity, stock-supported saved routes and predictor endpoints remain conditional on the specified retained assets. Those boundaries persist after formal adoption of the fixed qualification recipe.

**Evidence:** [fixed512 development](G:/逆合成结构化反馈/graphga_fixedtask512_v1/remote_snapshot_20261004/external_fixed_outputs_v1/analysis.json), [allocation128 design](G:/逆合成结构化反馈/graphga_candidateallocation128_v1/remote_snapshot_20261004/design.json), [allocation128 results](G:/逆合成结构化反馈/graphga_candidateallocation128_v1/remote_snapshot_20261004/external_fixed_outputs_v1/analysis.json), [retained-precursor design](G:/逆合成结构化反馈/graphga_retainedprecursor128_v1/remote_snapshot_20261005/design.json), [retained-precursor results](G:/逆合成结构化反馈/graphga_retainedprecursor128_v1/remote_snapshot_20261005/external_fixed_outputs_v1/analysis.json), [complete physical-cost scopes](G:/逆合成结构化反馈/graphga_qualification_independent1024_v1/final_physical_cost_scope_summary_20261006.json), [complete CPU timing scope](G:/逆合成结构化反馈/graphga_qualification_independent1024_v1/logging_performance_repair_20261006/CPU_timing_scope_final_20261006.json), and [four final independent audit bindings](G:/逆合成结构化反馈/graphga_qualification_independent1024_v1/witness_reaudit/all_four_independent_audits_completed_20261006.json).

## References cited in this supplement

References below credit the component algorithms/resources. Project effects and implementation-specific settings derive from saved study records, not from these component papers. The existing checked [bibliography](G:/逆合成结构化反馈/bibliography_20261002.json) retains complete authors and source-verification limitations.

- Baell, J.B. and Holloway, G.A. (2010) New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening Libraries and for Their Exclusion in Bioassays. *J. Med. Chem.*, 53, 2719–2740. [DOI](https://doi.org/10.1021/jm901137j).
- Bemis, G.W. and Murcko, M.A. (1996) The Properties of Known Drugs. 1. Molecular Frameworks. *J. Med. Chem.*, 39, 2887–2893. [DOI](https://doi.org/10.1021/jm9602928).
- Bickerton, G.R. et al. (2012) Quantifying the chemical beauty of drugs. *Nat. Chem.*, 4, 90–98. [DOI](https://doi.org/10.1038/nchem.1243).
- Brenk, R. et al. (2008) Lessons Learnt from Assembling Screening Libraries for Drug Discovery for Neglected Diseases. *ChemMedChem*, 3, 435–444. [DOI](https://doi.org/10.1002/cmdc.200700139).
- Brown, N.; Fiscato, M.; Segler, M.H.S.; Vaucher, A.C. (2019) GuacaMol: Benchmarking Models for de Novo Molecular Design. *Journal of Chemical Information and Modeling*, **59**(3), 1096–1108. [doi:10.1021/acs.jcim.8b00839](https://doi.org/10.1021/acs.jcim.8b00839)
- Chen, B. et al. (2020) `Retro*: Learning Retrosynthetic Planning with Neural Guided A* Search.` *PMLR*, 119, 1608–1616. [Official proceedings](https://proceedings.mlr.press/v119/chen20k.html).
- Chen, S. et al. (2026) SynTwins: a retrosynthesis-guided framework for synthesizable molecular analog generation. *Chem. Sci.*, 17, 2255–2262. [DOI](https://doi.org/10.1039/d5sc05225d).
- Ertl, P.; Schuffenhauer, A. (2009) Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. *Journal of Cheminformatics*, **1**, 8. [doi:10.1186/1758-2946-1-8](https://doi.org/10.1186/1758-2946-1-8)
- Genheden, S. and Bjerrum, E. (2022) PaRoutes: towards a framework for benchmarking retrosynthesis route predictions. *Digital Discovery*, 1, 527–539. [DOI](https://doi.org/10.1039/d2dd00015f).
- Genheden, S., Thakkar, A., Chadimová, V., Reymond, J.-L., Engkvist, O. and Bjerrum, E. (2020) AiZynthFinder: a fast, robust and flexible open-source software for retrosynthetic planning. *Journal of Cheminformatics*, 12, 70. [Official article](https://link.springer.com/article/10.1186/s13321-020-00472-1); [DOI](https://doi.org/10.1186/s13321-020-00472-1).
- Huang, K. et al. (2021) Therapeutics Data Commons: Machine Learning Datasets and Tasks for Drug Discovery and Development. *NeurIPS Datasets and Benchmarks*, 1. [Official proceedings](https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/hash/4c56ff4ce4aaf9573aa5dff913df997a-Abstract-round1.html).
- Jensen, J.H. (2019) A graph-based genetic algorithm and generative model/Monte Carlo tree search for the exploration of chemical space. *Chemical Science*, **10**(12), 3567–3572. [doi:10.1039/C8SC05372C](https://doi.org/10.1039/C8SC05372C)
- Loeffler, H.H. et al. (2024) Reinvent 4: Modern AI–driven generative molecule design. *J. Cheminform.*, 16, 20. [DOI](https://doi.org/10.1186/s13321-024-00812-5).
- Lundberg, S.M., Erion, G.G. and Lee, S.-I. (2018) Consistent Individualized Feature Attribution for Tree Ensembles. *arXiv*, 1802.03888 (version 3, revised 2019). [Original preprint](https://arxiv.org/abs/1802.03888v3).
- Maziarz, K. et al. (2025) Re-evaluating retrosynthesis algorithms with Syntheseus. *Faraday Discuss.*, 256, 568–586. [DOI](https://doi.org/10.1039/d4fd00093e).
- Palczewska, A. et al. (2014) Interpreting Random Forest Classification Models Using a Feature Contribution Method. In *Integration of Reusable Systems*, AISC 263, 193–218. [DOI](https://doi.org/10.1007/978-3-319-04717-1_9).
- Rogers, D. and Hahn, M. (2010) Extended-Connectivity Fingerprints. *J. Chem. Inf. Model.*, 50, 742–754. [DOI](https://doi.org/10.1021/ci100050t).
- Southiratn, T. et al. (2025) CombiMOTS: Combinatorial Multi-Objective Tree Search for Dual-Target Molecule Generation. *PMLR*, 267, 56650–56691. [Official proceedings](https://proceedings.mlr.press/v267/southiratn25a.html).
- Sun, J. et al. (2017) ExCAPE-DB: an integrated large scale dataset facilitating Big Data analysis in chemogenomics. *J. Cheminform.*, 9, 17. [DOI](https://doi.org/10.1186/s13321-017-0203-5).
- Wettenstein, R., Mitchell, R. and Yu, P. (2026) Quadrature-TreeSHAP: Depth-Independent TreeSHAP and Shapley Interactions. *arXiv*, 2605.04497v1. [Original preprint](https://arxiv.org/abs/2605.04497v1); [DOI](https://doi.org/10.48550/arXiv.2605.04497).
- Wu, Z.; Ramsundar, B.; Feinberg, E.N.; Gomes, J.; Geniesse, C.; Pappu, A.S.; Leswing, K.; Pande, V. (2018) MoleculeNet: a benchmark for molecular machine learning. *Chemical Science*, **9**(2), 513–530. [doi:10.1039/C7SC02664A](https://doi.org/10.1039/C7SC02664A)
- Yang, K. et al. (2019) Analyzing Learned Molecular Representations for Property Prediction. *J. Chem. Inf. Model.*, 59, 3370–3388. [DOI](https://doi.org/10.1021/acs.jcim.9b00237).
- Zhong, W. et al. (2023) Retrosynthesis prediction using an end-to-end graph generative architecture for molecular graph editing. *Nat. Commun.*, 14, 3009. [DOI](https://doi.org/10.1038/s41467-023-38851-5).
