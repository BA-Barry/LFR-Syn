# Task-guided reaction-constrained replacement improves route-qualified molecular yield from fixed generator outputs


## Abstract

**Motivation:** Molecular optimization requires both predicted task performance and a supported synthesis route. We evaluated task-guided, reaction-constrained replacement on fixed generation slots.

**Results:** Four development and eight validation REINVENT runs supplied 256 and 512 GSK3B/JNK3 slots. Contribution guidance with joint-failure triggering and shared 6+6 route calls achieved 311/512 random-forest/route successes, versus 219 for Morgan and 249 for permuted contributions. Triggering and allocation added 28 and 29 successes. Fixed-output Chemprop/route counts were 176, 89 and 112. Prespecified contribution/Morgan alternation underperformed contribution alone by 15 RF and 20 Chemprop successes. Endpoints remain predictor-, stock- and planner-dependent.


## 1. Introduction

Generative molecules with desirable property scores can remain unresolved by synthesis planners; route assessment adds information beyond property and synthetic-complexity scores (Gao and Coley, 2020). We ask how a bounded post-generation workflow improves the number of outputs satisfying both a predicted task objective and an audited planning criterion. Retaining the raw sampling denominator keeps invalid and unsuccessful positions visible. Distinct-output counts also matter: repeatedly returning a previously known molecule can raise successful-slot yield without adding distinct structures.

Each original sampling position retains at most one final output. An eligible source undergoes additional planning or reaction-constrained replacement using a replayable disconnection, compatible inventory precursors, forward reconstruction, property qualification and saved-route verification. Failed interventions retain the original source. Complete policies are compared on identical positions, rather than only on successful candidates.

The mechanism studied is source-conditioned task information in precursor retrieval. A probability-form adaptation of random-forest decision-path contributions weights source fingerprint features (Palczewska et al., 2014). Contribution retrieval prioritizes compatible precursors retaining positively weighted features. A permutation control disrupts weight–environment correspondence while preserving weights and reaction-role membership; Morgan retrieval prioritizes structural similarity. Common downstream rules make these controlled workflow comparisons.

Related methods establish reaction-constrained and attribution-guided optimization precedents. SynTwins combines retrosynthesis, similar building-block retrieval and virtual synthesis (Chen et al., 2026). Syn-MolOpt combines substructure-mask explanations with functional reaction templates and planning-driven optimization (Yin et al., 2025), using graph-model fragment explanations (Wu et al., 2023). CombiMOTS explores reaction-based combinatorial space for dual-target generation (Southiratn et al., 2025).

Our contribution is the measured fixed-slot joint-yield effect of task-informed retrieval, intervention coverage and bounded candidate-route allocation, together with a matched information perturbation, distinct-output reporting and actual-work accounting. A disclosed SynTwins-compatible generator provides a close-method comparison. Published CombiMOTS models supply fixed-output cross-architecture sensitivity; this does not reproduce CombiMOTS search or establish biological activity.

## 2. Materials and methods

### 2.1. Cohorts, generation and experimental units

Development used two task-only REINVENT4 training runs per task (Loeffler et al., 2024), starting from the same prior with 128 optimization steps, batch size 16, DAP learning, sigma 128 and learning rate 0.0001. Every terminal checkpoint supplied 64 randomized-SMILES draws without a uniqueness requirement. Training seeds were 9301001/9301002 for GSK3B and 9301101/9301102 for JNK3; sampling seeds were 9301201/9301202 and 9301301/9301302. No route feedback trained these generators.

Further validation used four new training runs per task with unchanged generation settings: training/sampling seeds 10102001–10102004/10103001–10103004 for GSK3B and 10102101–10102104/10103101–10103104 for JNK3. Its seven policies and 512-slot denominator were recorded before new training. Development selected this extension but was not pooled into its estimates. Both designs are local recorded plans, not publicly registered protocols.

All scheduled positions remain in their denominators, including invalid structures and duplicates; extra sampling does not replace failures. A complete training run is the stochastic replicate; its 64 outputs share a checkpoint. Tasks share the prior and have fixed task-specific property assets. Development policies intervene at route-unresolved sources; validation additionally tests route-solved sources below the task cutoff.

### 2.2. Property qualification and reaction-root construction

Existing fixed PyTDC GSK3B/JNK3 random-forest assets supplied class-one probabilities (Breiman, 2001; Huang et al., 2021). The resources identify ExCAPE-DB as their data source (Sun et al., 2017). The implementation uses radius-two, 2,048-bit Morgan fingerprints, with joint thresholds 0.10 for GSK3B and 0.05 for JNK3. These are computational qualification thresholds, not measured potency. No property predictor is trained or independently benchmarked here; the fixed assets define operational endpoints. Accepted candidates require T(candidate) ≥ max[T(source), task cutoff] and QED(candidate) ≥ QED(source) − 0.05 (Bickerton et al., 2012).

The canonical inventory contains 13,432 structures derived from PaRoutes stock (Genheden and Bjerrum, 2022); planning membership uses InChIKey representation. SynTwins reverse templates generate precursor pairs for two-role forward reactions (Chen et al., 2026). A root is eligible only if forward replay reproduces the canonical source. Selection retains up to two distinct-template roots with one missing inventory precursor and two with two missing precursors, filling unused capacity from remaining ordered roots to at most four. This is a computational template relation, not experimental reaction validation.

Stocked roles retain original precursors. Missing roles retrieve inventory structures matching the reactant query: up to four alternatives for a one-missing-role root, or two per role and four combinations for a two-missing-role root. Enumeration visits candidate rank before root, with canonical product ordering and 16 forward attempts per arm/source. Duplicate products, the source and stocked products are excluded. The first eight canonical-unique QED-qualified products consume score allowances before task lookup; cached scores consume the same logical allowance. Qualified products are ordered by task score, QED and canonical SMILES.

### 2.3. Retrieval and the information control

Morgan ranks compatible precursors by radius-two, 2,048-bit Tanimoto similarity to the original precursor, then canonical SMILES (Rogers and Hahn, 2010). Contribution retrieval decomposes class-one probability along each source's actual decision paths: child-minus-parent increments are assigned to the splitting bit and averaged over trees. Root bias plus contributions is checked against the forest probability. This adapts Palczewska et al. (2014), whose classification-vote equivalence assumes pure leaves. This decomposition in the original seven-policy validation does not estimate SHAP values or chemical causal effects; a separate same-source TreeSHAP extension is described in Supplement S14.

Let H(x) contain unfolded radius-two source environments and φ_b(x) the contribution of folded bit b. Environment e receives w_e = max[0, φ_(e mod 2048)(x)]/n_(e mod 2048), where n_b counts distinct source environments mapping to b. For missing precursor p and replacement z, retrieval maximizes the sum of w_e for environments shared by x, p and z. Morgan similarity and canonical SMILES break ties. Apportionment does not eliminate fingerprint collision ambiguity.

Permutation shuffles weights within groups defined by presence in both missing roles, the first only or the second only. Seeds include task, training replicate, position, root and group. Development used 1,842 distinct seeds for 1,842 groups; one realization per group does not characterize permutation-to-permutation variation. All arms share reaction roots, qualification, final selection and fallback. Their successes are not retrospectively unioned.

### 2.4. Planning, allocation and native comparator

Syntheseus 0.7.2 (Maziarz et al., 2025) runs Graph2Edits (Zhong et al., 2023) through RetroStarSearch (Chen et al., 2020), with ConstantNodeEvaluator(0.0), ReactionModelLogProbCost, 20 predictions per expansion and maximum expansion depth 12. No learned Retro* value model is loaded. Initial searches allow 12 model calls and 45 soft seconds. The target's standard InChIKey is removed from stock, and search state resets between targets. Time is checked between expansions, permitting a single-call overrun. A positive requires a saved graph with the correct target, solved status and audited stock-supported leaves; it is conditional on these assets and settings.

Development candidate policies search the first qualified candidate. Continuation resumes the source's saved initial graph for up to 12 additional calls and 45 further seconds, rather than restarting a 24-call search. Validation contrasts first-candidate 12 calls with sequential six-call searches on the first two qualified candidates. Both share 45 soft seconds per policy/source; the second receives remaining time and cannot inherit unused calls. Six-call measurements are newly executed or exactly compatible saved queries, not truncated 12-call trajectories. Audited success permits replacement; otherwise the source is retained.

The SynTwins-compatible generator uses commit 53d6bbcfc540971460278ad85e31d312b749fa55, depth 3, k = 10, 1,024-bit fingerprints, fast enumeration, the canonical inventory and a fresh source tree. It receives common eight-score task/QED qualification and downstream planning. Enumeration has a 300-second source limit but no RDKit-template-call ceiling. Two disclosed compatibility edits fix a negative root slice and empty stock-role indexing. This adaptation does not reproduce the full inventory, native multiproperty-optimization experiments or oracle budget. Execution-limited decisions retain operational fallback and conservative scientific bounds.

### 2.5. Contrasts and statistical interpretation

Joint yield is final task qualification with an audited route divided by all raw slots. Route yield and distinct (task, canonical structure) joint outputs are also reported. The development mechanism contrast was contribution minus permutation. Validation's primary contrast is Mix minus contribution with joint triggering and 6+6 allocation; contribution-minus-permutation, joint-minus-route triggering, allocation and Morgan/native contrasts remain recorded secondary comparisons.

Mix alternates contribution and Morgan precursor-combination streams before task lookup, starting with contribution, under common 16-forward/eight-score limits. Duplicate reaction–reactant combinations are skipped while advancing the stream. It does not union two eight-score pools or use scores to reallocate allowance.

Paired wins/losses compare identical slots. Validation intervals enumerate all 65,536 empirical bootstrap combinations of four paired training-run effects within each fixed task. They describe conditional variation, not independent-molecule inference, randomization-test evidence, multiplicity-adjusted confirmation or biological generalization. Unknown native outcomes receive fallback-guaranteed lower and conservative upper bounds, which are not confidence intervals. Cohort size, policies and contrast labels were not changed in response to effects.

### 2.6. Fixed-output sensitivity and historical information

Published CombiMOTS fold_0–fold_9/model_0/model.pt classifiers average two-task Chemprop 1.5.2 predictions (Southiratn et al., 2025; Yang et al., 2019). These scores never guide selection. Development evaluated fixed thresholds 0.05–0.75; validation uses 0.5 descriptively and exploratory 0.3/0.7. RF threshold grids, PAINS/BRENK filters, descriptors and Murcko scaffolds are also applied without reselection (Baell and Holloway, 2010; Brenk et al., 2008; Bemis and Murcko, 1996). Exact RF-grid values were specified for postprocessing, not claimed as pretraining registration.

RF and Chemprop share broad ExCAPE origins. Chemprop arguments reference data/excape+aixfuse_curated.csv and 77,581–77,583 training entries per fold, whereas the public GSK3B_JNK3.csv has 95,159 rows and 95,157 unique canonical structures. Complete training composition and ordering are unavailable. Exact-match checks do not establish full training/scaffold exclusion or reconstructed out-of-fold predictions. Ensemble dispersion is model variation; thresholds are uncalibrated outputs. This is cross-architecture sensitivity, not biological validation.

A later fixed-output sensitivity compared original/B/C outputs from raw positions 0–15 of each validation run under one AiZynthFinder protocol (Genheden et al., 2020; Supplement S16).

A subsequent mechanism control, S, ranked reaction-compatible inventory precursors by ascending RDKit SA score (Ertl and Schuffenhauer, 2009), then molecular weight and canonical SMILES. Scores were calculated once for all 13,432 inventory structures. S retained the same saved roots, proposal order, 16-forward/eight-score limits, joint triggering, qualification and 6+6/shared-45-second planning. B−S was specified before S execution, after observing the original results; all 512 positions remained. It tests this particular precursor preference on existing sources.

A separately specified historical boundary supplies 117 distinct task–structures with prior task/QED and audited Graph12/45 evidence, excluding continuation-only solutions. Development compares repetition-allowing lookup, deterministic non-repeating lookup and maximum-cardinality assignment. Lookup lacks source-reaction reachability constraints and receives unequal historical information, so it diagnoses distinct discovery rather than fair cold-start efficiency. Validation reports fixed-output novelty relative to the same pool using task–structure and canonical identity.

### 2.7. Saved-record checks and actual-work accounting

Recorded scores, forward products and exactly compatible Graph queries are reused. Historical route reuse selects the first lexicographically compatible path independently of outcome and records adopted seeds. Development's seed-adoption exclusions supplement the unchanged denominator. Validation Graph caches are shared within each worker; concurrent workers do not constitute a globally deduplicated cache.

Postprocessing checks source correspondence, qualification, fallback, route-call/time limits and saved audits. Accounting separates generation, RF inference, CPU reaction operations, Graph inference, initialization and external postprocessing, distinguishing logical work from actual queries and cache hits. Concurrent component times and interrupted controller windows cannot be summed into continuous end-to-end wallclock. Shared preparation and unequal candidate coverage preclude equal-total-compute superiority. Commands, versions, seeds, reuse and component ledgers are provided in the Supplementary Methods.

## 3. Results

### 3.1. Separate 256-slot development evidence

All four development runs completed: 247 of 256 sources were valid, nine invalid. Contribution retrieval yielded 126 task/route successes and 125 distinct task–structures, compared with 93/92 for Morgan and 97/96 for permutation (Table 1). Contribution-minus-permutation had 49 paired wins and 20 losses (+11.33 percentage points); contribution-minus-Morgan had 53/20 (+12.89 points), positive in every training run. Four execution-limited native sources bound its joint count at 70–74. These results selected further validation and are not pooled with it.

**Table 1. Development outcomes on all 256 slots.** External joint uses fixed Chemprop mean ≥0.5. Native ranges reflect execution limitations, not sampling uncertainty.

|Policy|Audited route|RF/route joint|Distinct RF joint task–structures|Chemprop/route joint|
|---|---:|---:|---:|---:|
|Original|63|21|21|11|
|Saved-tree continuation +12|94|36|36|21|
|Morgan|135|93|92|55|
|Contribution|168|126|125|84|
|Permuted contribution|139|97|96|56|
|SynTwins-compatible|112 known|70–74|70 known|35–39|

Among 184 route-unresolved sources, Morgan/contribution/permutation supplied 79/126/98 qualified candidates and 72/105/76 audited replacements. Contribution's greater coverage is an operational pathway, not proof of greater intrinsic planning ease. Its 105 replacements included 69 already-qualified source rescues and 36 low-score conversions. Final joint count equalled the 126 original RF-qualified structures; development did not exceed unchanged-source qualification's theoretical upper bound. Permutation altered weights and selected ordering at 173 sources; all contribution–permutation wins/losses occurred within altered-weight sources.

### 3.2. Development sensitivities, historical lookup and costs

Saved-record checks reproduced 1,536 policy rows and 786 unique result files without errors or new model calls. Contribution retained 121 joint slots without PAINS and 88 without PAINS/BRENK, versus Morgan 89/62 and permutation 93/69. Its joint outputs occupied 120 nonempty Murcko scaffolds (maximum group three). Of 105 replacements, 104 changed scaffold; median similarity was 0.238, QED change +0.119 and molecular-weight change −93.7 Da. At GSK3B RF ≥0.5, contribution/Morgan/permutation joint counts were 29/17/11; at JNK3 ≥0.2, 11/3/7; at JNK3 ≥0.5, all were zero.

External Chemprop/route counts were 84/55/56 for contribution/Morgan/permutation, with contribution positive against both controls in all four runs. External scores nevertheless decreased at 20 of 105 replacements. The 525 newly evaluated structures had no exact public-table matches; 5,250 fold–structure predictions were newly executed and 176 earlier complete structure vectors reused. Neither result certifies exclusion from complete original model-training data.

Repetition-allowing historical lookup yielded 202 joint slots but only 27 distinct outputs and 175 duplicate occupancies; two structures occupied 85 and 84 successful slots. Non-repeating lookup and matching reached 109 distinct joint outputs. Contribution's 125 distinct joint structures, including 104 distinct replacements, all lay outside the specified pool. Thus lookup can exceed raw throughput while adding no pool-new structures. Excluding eight seed-adoption slots retained contribution/permutation/Morgan counts 120/93/87 over 248; excluding two cross-cohort assay-reuse slots gave 125/96/92 over 254.

Preparation used 8,008 training RF rows; parent planning required 2,558 new Graph calls, internal candidates 869, continuation 2,025 and native candidates 547. Logical Morgan/contribution/permutation candidate calls were 195/392/352. Shared chemistry executed 197,064 reverse-template calls, 11,188 source replays, 4,726 candidate forward attempts and 3,067 candidate RF predictions. Native enumeration recorded 7,012,890 reverse and 3,666,004 forward calls. Different allowances and shared caches prevent isolated-runtime or equal-cost conclusions.

### 3.3. Development budget diagnostic

A saved-trajectory diagnostic considered six calls on the first candidate followed by six on the second at 102 sources with at least two qualified candidates. Only missing measurements were added: one historical query reused and 18 new searches requiring 65 Graph calls, with nine audited solves and no execution/audit issues. Ten paired wins and one loss increased full-cohort RF joint yield 126→135 and distinct yield 125→133; external ≥0.5 joint yield increased 84→91 using 90 new fold–structure predictions. Run-level RF differences were +2, +5, −1 and +3. Logical calls decreased 328→280 in the eligible subset. This was neither fresh independent replication nor a common-wallclock comparison: each candidate had its own 45-second limit. It motivated the prospectively measured shared-time validation and does not replace its results.

### 3.4. Completed eight-run, seven-policy validation

All eight runs and seven policies completed. The 512 slots contained 491 valid canonical-unique sources and 21 invalids. Initial evaluation yielded 133 audited routes and 52 RF/route successes. Other valid strata comprised 196 task-qualified route-unresolved sources, 162 below-cutoff route-unresolved sources and 81 below-cutoff route-solved sources. Contribution/permutation/Morgan/Mix supplied eligible candidates for 286/233/185/276 sources. Native enumeration supplied 172 eligible decisions and eight execution-limited decisions, retained as unknown scientific endpoints. Deterministic record review covered all 3,584 policy rows without source/qualification/fallback mismatches or saved-route audit disagreements, using no new model calls. Final joint yields and paired policy effects are summarized in Figure 1.

**Table 2. Validation outcomes on all 512 slots.** Native ranges preserve fallback-guaranteed bounds; distinct native counts are known successes. External scores apply to unchanged final outputs. No row is an outcome-driven union.

|Policy|RF/route joint|RF joint (%)|Audited route|Distinct RF joint task–structures|Chemprop ≥0.5/route joint|
|---|---:|---:|---:|---:|---:|
|Original|52|10.16|133|52|23|
|Contribution, route trigger, first 12|254|49.61|335|253|147|
|Contribution, joint trigger, first 12|282|55.08|335|281|155|
|Contribution, joint trigger, shared 6+6|311|60.74|361|310|176|
|Mix, joint trigger, shared 6+6|296|57.81|347|295|156|
|Permutation, joint trigger, shared 6+6|249|48.63|303|247|112|
|Morgan, joint trigger, shared 6+6|219|42.77|276|219|89|
|SynTwins-compatible, joint trigger, shared 6+6|153–161|29.88–31.45|218–225|153 known|58–66|

The primary Mix-minus-contribution contrast produced 18 RF wins and 33 losses: −15/512, or −2.93 percentage points, with conditional run interval [−6.45, 0.59]. Three runs favored Mix and five contribution. At fixed-output Chemprop ≥0.5, Mix had 18 wins and 38 losses: −20/512, −3.91 points [−6.45, −0.98], positive in one run and negative in seven. Proposal alternation failed to improve the recorded primary endpoint and also underperformed in this property sensitivity. The primary result remains primary; positive secondary contrasts do not substitute for it.

**Table 3. Recorded paired component contrasts in validation.** Effects are percentage-point differences over 512 raw slots; brackets are descriptive 95% run-bootstrap intervals conditional on two fixed tasks. Secondary comparisons are not multiplicity-adjusted confirmation.

|Contrast|RF wins/losses|RF net slots; effect [interval]|Chemprop ≥0.5 wins/losses|External net slots; effect [interval]|
|---|---:|---|---:|---|
|Primary: Mix−contribution 6+6|18/33|−15; −2.93 [−6.45, 0.59]|18/38|−20; −3.91 [−6.45, −0.98]|
|Contribution−permutation 6+6|94/32|+62; +12.11 [8.98, 15.43]|87/23|+64; +12.50 [10.55, 14.84]|
|Joint−route triggering, first 12|28/0|+28; +5.47 [3.52, 7.23]|8/0|+8; +1.56 [0.78, 2.15]|
|Shared 6+6−first 12, contribution|33/4|+29; +5.66 [3.91, 7.62]|23/2|+21; +4.10 [1.95, 7.03]|
|Contribution−Morgan 6+6|119/27|+92; +17.97 [12.11, 23.83]|106/19|+87; +16.99 [12.89, 21.29]|

Contribution assignment and allocation favored contribution in all eight runs under both property definitions. Joint triggering had seven positive RF runs and one tie; external qualification had five positives and three ties. All 28 RF trigger gains came from already route-solved, below-cutoff sources, explaining unchanged route counts. Allocation increased route count by 26 and RF joint count by 29. Contribution's advantage over the native adaptation remained 150–158 RF slots and 110–118 external slots under conservative bounds; these comparisons retain unequal enumeration work and compatibility limitations.

### 3.5. What the final outputs establish

Contribution 6+6 yielded 164/256 GSK3B and 147/256 JNK3 RF joint outputs, versus Morgan 122/256 and 97/256; external joint counts were 123/256 and 53/256, versus 71/256 and 18/256. It recovered 143/196 task-qualified unresolved sources, 85/162 below-cutoff unresolved sources and 31/81 below-cutoff route-solved sources. Its 259 replacement slots represented 258 distinct replacement structures; 310 distinct final joint structures, including retained original successes, lay outside the unchanged 117-entry pool under task–structure and canonical identity. This is relative to that pool, not global chemical novelty.

Only 248 original sources meet the fixed RF thresholds. Even if every original structure acquired a route, their fixed-structure, fixed-predictor joint-count ceiling would be 248/512. The observed 311 joint outputs exceed it by 63. This deterministic cohort-specific ceiling distinguishes the attainable value of structural changes from planning more extensively on unchanged sources. It is not a newly executed continuation control or a ceiling under another predictor.

The 259 replacement routes contain one reaction step at 221 slots, two at 17, three at 15, four at four and five at two. Most successes are short inventory-supported routes; the study does not establish a benefit for difficult long synthesis routes. Of these replacements, 258 change the source Murcko scaffold. Median source similarity is 0.252, QED change +0.110 and molecular-weight change −77.7 Da. Retaining a prediction while changing scaffold does not demonstrate preserved biological function or pharmacophore.

### 3.6. Fixed-output sensitivity and actual validation work

At common RF cutoffs 0.25/0.50/0.75, contribution 6+6 yielded 127/48/9 joint slots, Mix 111/36/4, permutation 82/26/3 and Morgan 63/20/3; originals yielded 14/2/0. At external cutoffs 0.3/0.5/0.7, contribution yielded 197/176/150. External ≥0.5 distinct task–structures were 175 for contribution, 155 for Mix, 110 for permutation and 89 for Morgan. These are fixed-output sensitivities, not calibrated activity probabilities.

Contribution's RF joint set retained 306/311 without PAINS, 214 without PAINS/BRENK and 287 at molecular weight ≤500 Da. External ≥0.5 retained 173/176, 132 and 167 under those filters, versus Morgan 87/89, 68 and 85. These computational quality checks do not establish safety or activity. One contribution final output exactly matches the public labeled table, but is a non-joint fallback. Excluding public-matching outputs from success without reselection or denominator changes leaves contribution's 311 RF and 176 external successes unchanged. Complete training-set coverage remains unavailable.

Generation used 15,993 training RF predictions. Initial parent searches executed 4,987 Graph calls; candidate stages executed 3,918 new calls across 1,342 new queries, with 11 historical12 query reuses. Logical contribution 6+6/first12 candidate calls were 812/1,013, versus Morgan 405 and Mix 782: equal ceilings do not imply equal realized work or isolated runtimes. Internal chemistry executed 470,169 reverse-template calls and 10,779 new candidate forward attempts. Native enumeration separately recorded 14,973,625 reverse and 7,496,772 forward calls.

Candidate RF execution totalled 9,784 predictions: 6,646 internal and 3,138 native. Concurrent stages predicted 315 identical task–structure pairs twice, leaving 9,469 unique pairs with no disagreements. A separate audit found one internal prediction repeated from new training, also with identical score. All executed work remains charged; these omissions prevent a global zero-repetition claim. External postprocessing evaluated 1,225 unique structures including sources: 12,110 new fold–structure predictions and 140 historical predictions reused, with no historical recomputation.

The original controller encountered a broken output pipe after internal policies and native enumeration completed. Recovery ran only missing native candidate planning, final analysis and external postprocessing, retaining saved training, parent searches and internal results. Interrupted controller observations and accumulated component times remain separate from continuous end-to-end elapsed time; no overall speedup or equal-total-cost superiority is claimed.

### 3.7. Exploratory attribution extension on the same sources

After the original validation, exact path-dependent TreeSHAP was substituted for decision-path contributions under matched candidate and planning rules (Supplement S14). On the same 512 positions it produced 302 RF joint and 175 external-0.5 joint outputs, versus saved contribution counts of 311 and 176. Paired RF wins/losses were 27/36; external wins/losses were 27/28. These complete observations support workflow benefit under another attribution implementation but do not establish its superiority or equivalence. A source-only attribution-concentration routing replay, with thresholds fitted while leaving out each complete training run, retained 311 RF joint and yielded 179 external joint outputs (Supplement S15). Its feature was selected after observing the full cohort, so this replay is exploratory, not an independent validation or an executed combined policy. The original primary contrast remains unchanged.

### 3.8. Fixed outputs under a second planner

The 128-slot AiZ sensitivity retained six invalid sources and required 232 distinct target searches. Original/B/C yielded 43/80/65 audited routes, 21/68/49 RF joint outputs and 9/42/20 external-0.5 joint outputs. B−C had 23/4 RF wins/losses and 24/2 external wins/losses; RF effects were positive in all eight runs, and external effects positive in seven and tied in one. No route endpoint was unknown. These shared-source observations support fixed-output gains under a second planner. They do not test candidate selection within AiZ or equal computation across planners (Supplement S16).

### 3.9. Task information versus a precursor-complexity preference

Contribution B yielded 311 RF joint, 361 route and 176 external-0.5 joint outputs, versus S counts of 83, 156 and 35 on the same 512 slots. Paired B−S RF wins/losses were 234/6 (+44.53 percentage points; conditional interval 39.26–49.41); external wins/losses were 147/6 (+27.54 points; 24.61–30.47). Both joint contrasts were positive in all eight runs, with no unknown endpoint. Distinct joint task–structures numbered 310 versus 73 under RF and 175 versus 29 externally.

B supplied a qualified replacement candidate at 286 sources, versus 31 for S; all 31 S replacements acquired audited routes. Greater task-qualified candidate coverage explains much of this workflow difference. The follow-up required 1,446 new forward attempts, 992 RF predictions, 23 Graph calls and 180 external fold–structure predictions; compatible saved work was reused. It supports an increment over this specified SA preference, without excluding every complexity-related explanation or providing independent-source replication.

## 4. Discussion

Completed validation provides specific positive evidence for improving joint yield from fixed generator outputs. Task-informed precursor assignment improves joint yield relative to permutation and Morgan; joint-failure triggering converts low-score structures despite existing routes; sequential candidate allocation improves yield under a common route-call/time ceiling. Both the RF definition and fixed-output Chemprop sensitivity preserve these directions. The 311 RF joint outputs also exceed the cohort's unchanged-source qualification ceiling. These findings identify complementary ways to improve task-qualified, audited-route output yield, beyond demonstrating executable replacement software.

The negative primary Mix result constrains one combination scheme. Alternating score-blind retrieval streams reduced pooled yield relative to contribution alone, particularly under the external threshold. Candidate diversity can help when planning a second ranked candidate, yet alternating precursor streams within the same allowance need not supply better qualified candidates. This distinction supports the measured contribution/joint-trigger/allocation policy while leaving the tested mixing hypothesis unsupported. It does not justify redefining the primary contrast or claiming that all combinations improve yield.

Contribution–permutation tests information assignment within shared chemistry and qualification rules. The SA precursor control further supports task-informed retrieval over one simple complexity preference, with candidate qualification coverage as an operational pathway. These comparisons do not identify chemical causality, universally easier planning or an optimal ranking across targets. Primary positives require saved Graph2Edits witnesses. Fixed-output AiZ sensitivity adds second-planner support under different search units; neither establishes laboratory feasibility or biological potency.

Distinct-output reporting prevents repeated historical lookup from masquerading as discovery. Development lookup achieved higher raw yield by repetition, while contribution supplied structures outside the specified pool; validation supplied 310 distinct pool-new joint outputs. This is relative to a defined historical source. Public matches and incomplete checkpoint-training provenance preclude complete out-of-training claims, although removing public-matching successes does not explain the validation advantage.

Substantial scaffold changes and predominantly one-step replacement routes shape applicability. The workflow improves predicted-task and inventory-supported-route yield by changing structures substantially. It does not demonstrate conserved function for an active source or gains on difficult long routes. PAINS/BRENK, QED and weight sensitivities are computational descriptors, not toxicological or medicinal-chemistry validation.

Only two fixed kinase predictors and two development/four validation training runs per task were examined. External models share broad data origins; complete training provenance and biological measurements are absent. The native adaptation has disclosed compatibility changes, execution-limited unknowns and different enumeration work. The matched TreeSHAP extension addresses one alternative attribution definition; broader attribution comparisons and task transfer remain untested. Shared caches, RF repetition and unequal candidate coverage limit efficiency interpretation: bounded planning supports the measured policy contrasts, not equal-end-to-end-compute superiority.

## 5. Conclusion

Task-guided, reaction-constrained replacement increased the number of predicted-task and audited-route successes returned from fixed generator slots for the studied GSK3B/JNK3 assets. Recorded secondary contrasts support source-conditioned contribution retrieval, joint-failure triggering and sequential candidate-route allocation; the prespecified proposal-mixing comparison did not improve yield. The 311 RF joint outputs exceeded the 248-slot fixed-source qualification ceiling, while fixed-output Chemprop sensitivity retained the component advantages. These results establish a conditional computational benefit under the specified predictors, inventory and planner. They do not establish biological activity, conserved function, laboratory synthesis or equal-total-compute superiority.

## References

Baell, J.B.; Holloway, G.A. (2010) New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening Libraries and for Their Exclusion in Bioassays. *Journal of Medicinal Chemistry*, **53**(7), 2719-2740. [doi:10.1021/jm901137j](https://doi.org/10.1021/jm901137j)

Bemis, G.W.; Murcko, M.A. (1996) The Properties of Known Drugs. 1. Molecular Frameworks. *Journal of Medicinal Chemistry*, **39**(15), 2887-2893. [doi:10.1021/jm9602928](https://doi.org/10.1021/jm9602928)

Bickerton, G.R.; Paolini, G.V.; Besnard, J.; Muresan, S.; Hopkins, A.L. (2012) Quantifying the chemical beauty of drugs. *Nature Chemistry*, **4**(2), 90-98. [doi:10.1038/nchem.1243](https://doi.org/10.1038/nchem.1243)

Breiman, L. (2001) Random Forests. *Machine Learning*, **45**(1), 5-32. [doi:10.1023/a:1010933404324](https://doi.org/10.1023/a:1010933404324)

Brenk, R.; Schipani, A.; James, D.; Krasowski, A.; Gilbert, I.H.; Frearson, J.; Wyatt, P.G. (2008) Lessons Learnt from Assembling Screening Libraries for Drug Discovery for Neglected Diseases. *ChemMedChem*, **3**(3), 435-444. [doi:10.1002/cmdc.200700139](https://doi.org/10.1002/cmdc.200700139)

Chen, B.; Li, C.; Dai, H.; Song, L. (2020) Retro*: Learning Retrosynthetic Planning with Neural Guided A* Search. *Proceedings of Machine Learning Research*, **119**, 1608-1616. [Official proceedings](https://proceedings.mlr.press/v119/chen20k.html)

Chen, S.; Nam, G.; Aspuru-Guzik, A.; Jung, Y. (2026) SynTwins: a retrosynthesis-guided framework for synthesizable molecular analog generation. *Chemical Science*, **17**(4), 2255-2262. [doi:10.1039/d5sc05225d](https://doi.org/10.1039/d5sc05225d)

Ertl, P.; Schuffenhauer, A. (2009) Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. *Journal of Cheminformatics*, **1**, 8. [doi:10.1186/1758-2946-1-8](https://doi.org/10.1186/1758-2946-1-8)

Gao, W.; Coley, C.W. (2020) The Synthesizability of Molecules Proposed by Generative Models. *Journal of Chemical Information and Modeling*, **60**(12), 5714-5723. [doi:10.1021/acs.jcim.0c00174](https://doi.org/10.1021/acs.jcim.0c00174)

Genheden, S.; Bjerrum, E. (2022) PaRoutes: towards a framework for benchmarking retrosynthesis route predictions. *Digital Discovery*, **1**(4), 527-539. [doi:10.1039/d2dd00015f](https://doi.org/10.1039/d2dd00015f)

Genheden, S.; Thakkar, A.; Chadimová, V.; Reymond, J.-L.; Engkvist, O.; Bjerrum, E. (2020) AiZynthFinder: a fast, robust and flexible open-source software for retrosynthetic planning. *Journal of Cheminformatics*, **12**, 70. [doi:10.1186/s13321-020-00472-1](https://doi.org/10.1186/s13321-020-00472-1)

Huang, K.; Fu, T.; Gao, W.; Zhao, Y.; Roohani, Y.; Leskovec, J.; Coley, C.; Xiao, C.; Sun, J.; Zitnik, M. (2021) Therapeutics Data Commons: Machine Learning Datasets and Tasks for Drug Discovery and Development. *Proceedings of the Neural Information Processing Systems Track on Datasets and Benchmarks*, **1**. [Official proceedings](https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/hash/4c56ff4ce4aaf9573aa5dff913df997a-Abstract-round1.html)

Loeffler, H.H.; He, J.; Tibo, A.; Janet, J.P.; Voronov, A.; Mervin, L.H.; Engkvist, O. (2024) Reinvent 4: Modern AI–driven generative molecule design. *Journal of Cheminformatics*, **16**(1), 20. [doi:10.1186/s13321-024-00812-5](https://doi.org/10.1186/s13321-024-00812-5)

Maziarz, K.; Tripp, A.; Liu, G.; Stanley, M.; Xie, S.; Gaiński, P.; Seidl, P.; Segler, M.H.S. (2025) Re-evaluating retrosynthesis algorithms with Syntheseus. *Faraday Discussions*, **256**, 568-586. [doi:10.1039/d4fd00093e](https://doi.org/10.1039/d4fd00093e)

Palczewska, A.; Palczewski, J.; Marchese Robinson, R.; Neagu, D. (2014) Interpreting Random Forest Classification Models Using a Feature Contribution Method. In: Bouabana-Tebibel, T. and Rubin, S.H. (eds), *Integration of Reusable Systems*. Advances in Intelligent Systems and Computing, **263**, 193–218. Springer, Cham. [doi:10.1007/978-3-319-04717-1_9](https://doi.org/10.1007/978-3-319-04717-1_9)

Rogers, D.; Hahn, M. (2010) Extended-Connectivity Fingerprints. *Journal of Chemical Information and Modeling*, **50**(5), 742-754. [doi:10.1021/ci100050t](https://doi.org/10.1021/ci100050t)

Southiratn, T.; Koo, B.; Lu, Y.; Kim, S. (2025) CombiMOTS: Combinatorial Multi-Objective Tree Search for Dual-Target Molecule Generation. *Proceedings of Machine Learning Research*, **267**, 56650-56691. [Official proceedings](https://proceedings.mlr.press/v267/southiratn25a.html)

Sun, J.; Jeliazkova, N.; Chupakhin, V.; Golib-Dzib, J.-F.; Engkvist, O.; Carlsson, L.; Wegner, J.; Ceulemans, H.; Georgiev, I.; Jeliazkov, V.; Kochev, N.; Ashby, T.J.; Chen, H. (2017) ExCAPE-DB: an integrated large scale dataset facilitating Big Data analysis in chemogenomics. *Journal of Cheminformatics*, **9**(1), 17. [doi:10.1186/s13321-017-0203-5](https://doi.org/10.1186/s13321-017-0203-5)

Wu, Z.; Wang, J.; Du, H.; Jiang, D.; Kang, Y.; Li, D.; Pan, P.; Deng, Y.; Cao, D.; Hsieh, C.-Y.; Hou, T. (2023) Chemistry-intuitive explanation of graph neural networks for molecular property prediction with substructure masking. *Nature Communications*, **14**(1), 2585. [doi:10.1038/s41467-023-38192-3](https://doi.org/10.1038/s41467-023-38192-3)

Yang, K.; Swanson, K.; Jin, W.; Coley, C.; Eiden, P.; Gao, H.; Guzman-Perez, A.; Hopper, T.; Kelley, B.; Mathea, M.; Palmer, A.; Settels, V.; Jaakkola, T.; Jensen, K.; Barzilay, R. (2019) Analyzing Learned Molecular Representations for Property Prediction. *Journal of Chemical Information and Modeling*, **59**(8), 3370-3388. [doi:10.1021/acs.jcim.9b00237](https://doi.org/10.1021/acs.jcim.9b00237)

Yin, X.; Wang, X.; Wu, Z.; Li, Q.; Kang, Y.; Deng, Y.; Luo, P.; Liu, H.; Shi, G.; Wang, Z.; Yao, X.; Hsieh, C.-Y.; Hou, T. (2025) Syn-MolOpt: a synthesis planning-driven molecular optimization method using data-derived functional reaction templates. *Journal of Cheminformatics*, **17**(1), 27. [doi:10.1186/s13321-025-00975-9](https://doi.org/10.1186/s13321-025-00975-9)

Zhong, W.; Yang, Z.; Chen, C.Y.-C. (2023) Retrosynthesis prediction using an end-to-end graph generative architecture for molecular graph editing. *Nature Communications*, **14**(1), 3009. [doi:10.1038/s41467-023-38851-5](https://doi.org/10.1038/s41467-023-38851-5)

## Figure legend

**Figure 1. Fixed-slot joint yield and paired policy effects in eight complete validation runs.** (a) Outputs satisfying the RF task cutoff (0.10 for GSK3B; 0.05 for JNK3) and a saved, audited Graph2Edits route. (b) The same fixed outputs satisfying published Chemprop ensemble mean ≥0.5 and the unchanged route criterion. All 512 sampling slots, including 21 invalid sources, remain in the denominator. A ranks reaction-compatible precursors by Morgan similarity; B uses source-specific RF path contributions; C permutes these contributions; Mix alternates B and A proposal streams before task-score lookup. Joint triggering additionally intervenes on route-solved sources below the RF cutoff. Sequential 6+6 searches evaluate at most two candidates within the first-candidate policy's 12-call ceiling and shared 45-second soft allowance. Original sources provide fallback. SynTwins-compatible bars show lower bounds and pale extensions for eight execution-limited decisions: 153–161 RF/route and 58–66 Chemprop/route successes. (c) Paired joint-yield differences over all raw slots. Blue circles denote RF; orange squares denote Chemprop. Lines show conditional 95% empirical percentile intervals, enumerating all 65,536 bootstrap combinations of four paired training-run effects within each fixed task. The recorded primary contrast is Mix−B under joint triggering and 6+6 allocation; other contrasts are recorded secondary comparisons without multiplicity adjustment. Mix did not improve the primary endpoint. Intervals describe run variation under these fixed tasks. Predictions and route audits do not establish biological activity or laboratory synthesis.
