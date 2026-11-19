# Task qualification improves audited route yield in local molecular repair

## Abstract

**Motivation:** Molecular optimization can produce high predictor scores without a route to available building blocks. Requiring repairs to preserve every source score may exclude useful task-qualified alternatives.

**Results:** LFR-Syn combines reaction-constrained precursor replacement, bounded route planning and auditable output selection. In a prospective comparison across 16 independent GraphGA runs, fixed task qualification increased externally predicted task-qualified and audited-route outputs from 202 to 883 of 1,024 retained positions compared with source-score-preserving qualification. The paired task-equal gain was 66.5 percentage points (95% confidence interval 50.5–82.5). A contribution-permuted control yielded 880 successes; contribution-alignment superiority was not established in this comparison. Successful distinct structures increased, while overall output diversity contracted and guiding scores decreased. Historical component tests and negative transfer results delimit the evidence. These are computational outcomes, not demonstrated activity or synthesis.

**Availability and implementation:** A documented private reproduction package is prepared; public deposition remains pending.

## 1 Introduction

Molecular generators are commonly optimized against inexpensive property predictors. A high predicted score does not establish a route from available starting materials, and chemical complexity scores need not agree with retrosynthetic planning (Gao and Coley, 2020). Benchmarks such as GuacaMol and Therapeutics Data Commons facilitate reproducible property optimization, but the choice of endpoint determines what improvement means (Brown et al., 2019; Huang et al., 2021). For a design workflow, the useful output is a molecule meeting a stated task requirement for which the specified planner can provide an auditable route. Raising the score of an unchanged, route-unsolved molecule does not necessarily increase this output.

Retrosynthesis-guided analog generation and planning-driven molecular optimization provide alternatives to unconstrained generation (Chen et al., 2026; Yin et al., 2025). Local precursor replacement offers a practical interface: retain a reaction context, retrieve stock precursors, form candidate products and check their routes. Task information can prioritize this search, using fingerprint similarity or model contributions (Rogers and Hahn, 2010; Palczewska et al., 2014). Such guidance should be tested against matched controls because model explanation is not evidence of chemical causality, and alternative explanation methods have different assumptions (Wu et al., 2023).

A separate choice concerns candidate qualification. Preserving the guiding score of an already optimized source is stronger than meeting the task cutoff. It may exclude route-accessible candidates even when another predictor still regards them as task-qualified. Relaxing this requirement could increase usable outputs, but could also lower scores, concentrate many source positions into a few products or expose predictor disagreement. These consequences must be measured together.

We study LFR-Syn, an auditable workflow for reaction-constrained local replacement. Its prospective primary comparison isolates fixed task qualification from source-score-preserving qualification on the same candidate pools and route budgets. Earlier component studies test task contributions, triggering and allocation, while a separate task and planner probe the limits of transfer. The claim concerns predicted-task-qualified, audited-route output yield under specified assets. It does not concern the route solvability of unchanged sources, preserved biological function or laboratory synthesis.

## 2 Methods

### 2.1 Prospective design and historical evidence

The prospective study used 16 independent GraphGA optimizer runs, eight each for GSK3B and JNK3 (Jensen, 2019). Each run used population size 100, mating pool size 200, 100 offspring attempts per generation, 250 generations and mutation rate 0.01, without early stopping. The objective was the fixed task random-forest (RF) probability alone. Independent initialization seeds selected the first 100 valid, connected, canonically unique structures from a uniformly permuted GuacaMol corpus, without task, quality or route filtering. A separate sampling seed drew 64 positions uniformly with replacement from the complete ordered terminal population. These 1,024 positions, including repeated structures, were retained without outcome filtering.

The five paired policies were Original, contribution-aligned B and contribution-permuted C, each evaluated under source-score-preserving and fixed-task qualification (Fig. 1). Candidate pools were formed once for B and C and shared between the two qualification rules. Initialization, optimization, sampling and search seeds, all offspring attempts, final populations and source positions were saved. The complete optimizer run, rather than a molecule, checkpoint or repeated position, was the independent unit. The protocol and primary comparison were specified before generating these runs (Supplement S23).

Earlier experiments were retained with their original roles. They included eight REINVENT validation runs, eight earlier GraphGA runs, fixed-task development analyses of those earlier runs, a six-run HIV transfer study and fixed-output planner and ranking sensitivities. Earlier development runs were excluded from the new primary analysis. REINVENT used the retained Reinvent 4 implementation and RF objective; its settings and unfiltered source sampling are given in Supplement S2–S3 (Loeffler et al., 2024). No old generator cohort was rerun to obtain a newer result.

### 2.2 Reaction-constrained candidate formation and qualification

The RF predictors used 2,048-bit radius-2 Morgan fingerprints. Task cutoffs were 0.10 for GSK3B and 0.05 for JNK3; these are operational predictor thresholds, not measured activity thresholds (Breiman, 2001; Rogers and Hahn, 2010). The PaRoutes N1 inventory contained 13,432 canonical structures, represented by 13,415 standard InChIKeys (Genheden and Bjerrum, 2022). The source target identity was excluded from the search stock.

The workflow first assessed whether the source met both guiding-RF qualification and the parent-route criterion. Failure of either condition triggered candidate formation. Retained SynTwins reaction tables supplied at most four selected two-reactant roots, balancing one- and two-missing-stock-role cases (Chen et al., 2026). Source replay was checked before substitution. Each arm allowed at most 16 forward-formation attempts and eight QED-qualified logical scoring slots. The two-nonstock-precursor case retained the original 2 × 2 stock substitution rule. This limits the local search and does not enumerate all possible reactions or analogs.

B retrieved stock precursors using source-specific positive RF path contributions. For each tree, class-1 probability changes along the source path were attributed to fingerprint bits and averaged. Positive bit contributions were apportioned across the source's unfolded environments to account for collisions. Compatible precursors were ranked by retained positive environment weight, then ordinary Morgan similarity and canonical SMILES. This is a path-contribution method, not SHAP or a causal explanation. C permuted the same contribution values within the retained precursor-role groups, preserving chemistry, qualification and budgets. Each indexed control used one reproducible permutation, not a sampled permutation null distribution (Supplement S3–S6).

Both qualification rules required candidate QED to be at least source QED minus 0.05 (Bickerton et al., 2012). The historical source-preserving rule additionally required candidate RF to be at least the larger of source RF and task cutoff. The adopted fixed-task rule required candidate RF to meet the task cutoff without a source-relative score guarantee. Qualifying candidates were ranked by RF, then QED and canonical SMILES using the retained implementation. All score changes were recorded. This change admits additional candidates within the same scored pools; it does not change the reaction roots, stock, forward budget or ranker.

### 2.3 Bounded planning and external evaluation

Syntheseus 0.7.2 ran Graph2Edits with RetroStar search, a constant zero node evaluator and reaction-model log-probability cost (Chen et al., 2020; Maziarz et al., 2025; Zhong et al., 2023). The recipe requested 20 reaction predictions per expansion and maximum search depth 12. It did not use a learned Retro* value network. The parent received at most 12 reaction-model calls and 45 soft seconds. The top two qualifying replacement candidates received six calls each and shared a 45-second soft limit; unused calls were not transferred. Search stopped at the first audited solution, otherwise retaining the source. A soft deadline can be crossed by an in-flight model call.

A route counted only after reloading the saved native graph, extracting its route, checking stock leaves and excluding the target from stock. This is a planner- and inventory-specific computational criterion. Failed searches and technically unknown execution were distinguished. Compatible historical observations were reused with their original query identity, seed, budget and earliest saved origin, including unsuccessful observations. A shared policy prefix could reuse the same physical query but could not select a more favorable historical result. Logical policy work and newly executed work were recorded separately (Supplement S6, S23–S24).

External evaluation used all ten retained, released two-task Chemprop checkpoints from CombiMOTS only after final output selection (Yang et al., 2019; Southiratn et al., 2025). Their intended-task probabilities were averaged. The primary endpoint was mean probability at least 0.5 together with an audited route, without an additional guiding-RF conjunct. No checkpoint, cutoff or output was chosen using these external predictions. RF and Chemprop differ in architecture and evaluation role, but both kinase assets have ExCAPE provenance; incomplete original training composition prevents a claim of training-set independence (Sun et al., 2017).

### 2.4 Estimands and uncertainty

The sole prospective primary contrast was B fixed-task minus B source-preserving on the external joint endpoint. Each run supplied a paired difference in successes divided by its 64 retained positions. The estimand equally weighted the two fixed tasks. A Welch–Satterthwaite interval and two-sided test used the within-task variation of the eight paired run effects. A prespecified sensitivity bootstrap resampled complete runs within each task, with equal task weights, 100,000 replicates and seed 2026101801. Other policy contrasts and quality, diversity and mechanism summaries were exploratory; no multiplicity-adjusted confirmatory claims were made for them.

Uncertainty is conditional on the two tasks, predictors, stock and reaction assets. Repeated source positions, shared products and checkpoint predictions were not treated as independent replicates. Technical unknowns would retain lower and upper identification bounds; none remained in the completed prospective outputs. Global unions of successful task–canonical-SMILES pairs and nonempty Bemis–Murcko scaffolds were reported descriptively, together with the union of all final structures and alert sensitivities (Bemis and Murcko, 1996; Baell and Holloway, 2010; Brenk et al., 2008). Distinct-output counts do not provide an additional independent sample size.

### 2.5 Component, transfer and execution checks

Historical kinase comparisons retained source-score-preserving qualification. T1 compared joint-failure and route-failure triggering; T2 compared sequential 6 + 6 allocation with a first-candidate 12-call allocation on matched pools. The original REINVENT primary comparison, mixed B/A proposals versus B, retained its negative result. A fixed-output AiZynthFinder sensitivity changed the planner for 128 positions, without claiming an AiZ-guided workflow or equal total work (Genheden et al., 2020).

The HIV study used six independent source-generation runs and role-separated RF guidance and logistic-regression (LR) evaluation, trained on separate full-fingerprint/scaffold-group partitions of the shared NCI HIV labels (Wu et al., 2018). Thresholds were selected on calibration data before test prediction. The primary was B minus the mean of three nested C controls on LR qualification plus route, not on RF. These controls were not additional independent generator runs (Supplement S18).

Additional saved-pool analyses tested SA-based retrieval, exact TreeSHAP guidance and task reranking of retained native SynTwins candidate pools (Ertl and Schuffenhauer, 2009). Native compatibility changes, bounded execution unknowns and historical reuse remain documented. These analyses do not implement native multi-property optimization or establish state-of-the-art superiority. Four independent saved-record audits checked CPU records, native identity/QED/stock/ranking, route semantics and external endpoint/statistical reconstruction. Software checks and portable saved-vector reconstruction evaluated execution fidelity rather than adding efficacy replicates.

## 3 Results

### 3.1 Fixed task qualification increased joint output yield

All 16 specified runs, 1,024 retained positions and 5,120 policy outputs completed. The sources comprised 744 distinct task–structure pairs; duplicates were retained. Original outputs provided 169 audited routes and 150 external joint successes. B source-preserving produced 202 external joint successes, whereas B fixed-task produced 883 (Table 1). The primary gain was 681 positions, or 66.5 percentage points (95% confidence interval 50.5–82.5; two-sided P = 4.58 × 10⁻⁷). Fifteen paired run effects were positive and one was tied. Mean gains were 63.5 points for GSK3B and 69.5 for JNK3. The bootstrap sensitivity interval was 52.1–79.3 points (Fig. 2).

The comparison held candidate chemistry and search allowances constant. It therefore supports changing the candidate eligibility requirement for these optimized sources, rather than increasing reaction enumeration or route budgets. Relative to Original, B fixed-task supplied 733 more external joint successes; this was a secondary comparison, not the primary contrast.

**Table 1. Prospective five-policy outputs from 1,024 retained positions.** Joint success requires the ten-checkpoint intended-task mean to be at least 0.5 and an audited route. Distinct counts are global task–canonical-SMILES unions, not independent replicates.

| Policy | Audited routes | Joint successes | Distinct successful structures | All final structures |
| --- | ---: | ---: | ---: | ---: |
| Original | 169 | 150 | 110 | 744 |
| B source-preserving | 221 | 202 | 113 | 706 |
| C source-preserving | 215 | 196 | 113 | 713 |
| B fixed-task | 913 | 883 | 222 | 329 |
| C fixed-task | 935 | 880 | 412 | 529 |

### 3.2 Contribution alignment and qualification are different questions

C also benefited from fixed task qualification, increasing external joint successes from 196 to 880. Its exploratory gain was 66.8 points (95% interval 54.2–79.4). Under fixed qualification, B minus C was only 0.3 points (−5.3 to 5.9); under source preservation it was 0.6 points (−3.4 to 4.5). Neither result establishes B superiority or B/C equivalence (Table 2). The similar qualification gains indicate a shared eligibility bottleneck, without attributing that gain to correct assignment of contributions.

The distinction matters for the earlier GraphGA evidence. In eight earlier runs, B exceeded C by 65 external joint positions under source preservation, with an original run-bootstrap interval of 0.4–27.1 points. That result was preserved as the earlier primary, but the contribution advantage did not robustly recur in the fresh 16-run comparison. Fixed qualification on those earlier development inputs had yielded B/C counts of 462/469. These older inputs motivated the new study and were not pooled with it.

**Table 2. Paired run-level external joint contrasts.** Effects and 95% intervals are in percentage points, with 16 complete paired runs and equal task weights. Only the first row is the prespecified primary. Other rows are exploratory, without multiplicity adjustment.

| Contrast | Effect | 95% interval | Positive / tied / negative runs |
| --- | ---: | --- | --- |
| B fixed-task minus B source-preserving | 66.50 | 50.50 to 82.51 | 15 / 1 / 0 |
| C fixed-task minus C source-preserving | 66.80 | 54.21 to 79.38 | 16 / 0 / 0 |
| B minus C, fixed-task | 0.29 | −5.31 to 5.90 | 8 / 2 / 6 |
| B minus C, source-preserving | 0.59 | −3.36 to 4.53 | 1 / 14 / 1 |

### 3.3 Yield increased with score and diversity tradeoffs

B fixed-task replaced 744 positions, compared with 52 under source preservation. Its mean guiding-RF change over all 1,024 positions was −0.1822; the external mean change was −0.0007 and mean QED change was +0.1572. Corresponding C changes were −0.2344, −0.0262 and +0.1543. These full-position summaries include unchanged fallbacks; replacement-only changes were reported separately. They describe the observed tradeoff and are not additional confirmatory comparisons.

Successful distinct task–structure pairs increased from 113 to 222 for B and from 113 to 412 for C. Successful nonempty scaffolds increased from 86 to 174 for B and to 319 for C. However, the union of all final structures contracted from 706 to 329 for B and from 713 to 529 for C. A method can increase the diversity of its successful subset while reducing overall output diversity. Thus 883 successful positions do not represent 883 distinct discoveries, and B cannot be described as generally more diverse than C.

Alert-restricted sensitivities kept the full 1,024-position denominator. Excluding final structures flagged by either PAINS or BRENK left 45 Original and source-preserving B successes, versus 506 B fixed-task and 469 C fixed-task successes. These descriptive checks show that the new yield gain was not confined to flagged outputs; they establish neither compound safety nor assay activity. Thirty B fixed-task and 55 C fixed-task audited-route outputs remained below the external cutoff. A route alone did not guarantee external task qualification (Supplement S23–S24).

### 3.4 Historical components provide support with explicit limits

In the historical REINVENT validation, B and C provided 176 and 112 external joint successes across 512 raw positions. The recorded secondary B minus C effect was 12.5 points (empirical run-bootstrap interval 10.5–14.8). T1 added eight external joint successes without changing the 335-route count, supporting repair of below-cutoff sources even when a route already existed. T2 added 21 external joint successes and 26 routes; recorded logical candidate Graph calls decreased from 1,013 to 812 and resident search time from 2,002.2 to 1,619.6 seconds. This is component work, not an end-to-end speedup claim.

The original mixed-proposal primary was negative: Mix supplied 156 successes versus B's 176, a −3.9-point external effect with interval −6.4 to −1.0. Mixing proposal streams was therefore not retained as evidence of improvement. The six-run HIV transfer primary was also unsupported: B supplied eight LR joint successes versus control counts of eight, eight and nine. The paired effect was −0.087 points (−0.781 to 0.868), despite positive RF and route secondary effects. Guiding-predictor gains did not establish external-task transfer (Table 3).

**Table 3. Selected historical results and their original roles.** All counts refer to each panel's retained raw denominator. Kinase panels used source-preserving qualification. HIV used its separately calibrated RF/LR thresholds. Historical results were not pooled into the new primary.

| Panel and comparison | Positions / independent runs | External joint counts | Evidence role |
| --- | --- | --- | --- |
| REINVENT B versus C | 512 / 8 | 176 versus 112 | Recorded secondary |
| T1 joint versus route trigger | 512 / 8 | 155 versus 147 | Recorded secondary |
| T2 6 + 6 versus first-candidate 12 | 512 / 8 | 176 versus 155 | Recorded secondary |
| REINVENT Mix versus B | 512 / 8 | 156 versus 176 | Original primary, negative |
| Earlier GraphGA B versus C | 512 / 8 | 165 versus 100 | Earlier primary, not robustly replicated |
| HIV B versus three C controls | 384 / 6 | 8 versus 8, 8, 9 | LR primary, no positive support |

Fixed-output AiZynthFinder gave 42/20 external joint successes for B/C on 128 positions, versus nine Original successes. SA preference produced 35 successes versus B's 176 on the historical 512-position pool; TreeSHAP guidance produced 175. These are sensitivities, not proof that all complexity explanations are excluded or that guidance methods are equivalent. Reranking retained native SynTwins pools yielded conservative external gains of 24–40 positions over native ranking, with eight technical unknowns retained. Native-compatible execution does not constitute native multi-property optimization or a state-of-the-art benchmark (Supplement S14, S16–S22).

### 3.5 Execution fidelity and measured work

The prospective study executed 49,255 new native forward calls: 45,113 source-replay calls and 4,142 candidate-formation calls. Candidate scoring executed 3,054 new RF predictions; optimizer work was recorded separately. Route planning executed 1,243 new physical queries and 9,706 model calls, while 53 compatible historical queries were reused. External evaluation executed 10,690 new checkpoint–structure predictions and reused 490, covering all 1,118 distinct evaluated structures. These totals distinguish physical execution from repeated policy references.

All four saved-record audits completed without unresolved scientific discrepancies. Native checking recomputed 14,939 QED values and independently checked stock identities and weighted rankings. The endpoint audit reconstructed all policy outcomes, ten-fold prediction means and 21 run-level endpoint contrasts. A portable, saved-data-only reconstruction also reproduced 13,757 scalar fields, with no mismatches beyond final-bit floating differences. It required no model inference or GPU. The canonical executable includes explicit fixed-task and historical source-preserving modes, tested qualification behavior and a matching retained fixed-output replay.

CPU formation calendar time was 19.63 hours; generation-to-external completion was 28.14 hours. These include persistence, recovery and interstage handling and are not pure chemical-kernel runtime or an algorithm speed comparison. Audit repairs concerned resume-origin representation and native reaction-object lifetime in the auditor, without changing study scores, roots, candidates or searches. Saved failed attempts and repair checks were retained (Supplement S24).

## 4 Discussion

The prospective result supports a specific design choice: qualify local replacements against the task requirement rather than insisting that they preserve the score of a source optimized by the same predictor. On the studied sources, reaction pools and budgets, this increased output positions meeting a separate predictor criterion and an audited-route criterion. The large run-level gain occurred in both fixed tasks and was supported by an increased number of successful distinct structures. It therefore cannot be reduced to more repeated occupancy alone. At the same time, overall diversity decreased and guiding scores fell, so the objective remains a tradeoff rather than unrestricted improvement.

The evidence separates eligibility from information assignment. Historical kinase comparisons supported contribution guidance, joint triggering and sequential planning, while the new B/C comparison did not establish a contribution advantage. C's broader successful-structure union is relevant when diverse output is the priority. B's small observed raw-yield advantage cannot justify declaring it best. The primary contribution of this study is the validated qualification change within an auditable local-repair workflow; contribution alignment remains a conditional component rather than a universally necessary explanation.

Negative results constrain transfer without negating the primary comparison. Mixed proposals reduced yield under a common proposal allowance, and HIV failed its LR joint primary despite RF/route improvements. Retained-precursor and proposal-count-first development variants also failed to improve their matched controls and were not adopted (Supplement S24). These outcomes argue against simply expanding algorithm complexity or interpreting every positive guided score as transferable utility.

Several limitations remain. Both kinase predictors originate from ExCAPE-related assets, with incomplete training composition; architectural separation does not establish biochemical independence. HIV used a different data source and separate predictor roles, but its modest test discrimination and possible distribution shift limit interpretation. Fixed thresholds specify this experiment's requirements rather than validated activity probabilities. Reaction templates, a finite stock and short bounded searches restrict the route claim. Native route auditing validates stored computational witnesses, not reaction yields, conditions or laboratory execution. QED and alert filters are screening descriptors, not safety or function guarantees.

Local replacement changes the molecule. It can increase task-qualified route output while leaving the original target unsolved and changing scaffold or function. Repeated terminal draws, shared products and reusable queries also make raw yield, distinct successful output and physical compute different quantities. We report all three rather than inferring more discoveries or equal-compute superiority from a higher position count. Wider task, stock and long-route applicability requires separate studies; the completed study does not estimate those effects.

## 5 Conclusion

Fixed task qualification increased predicted-task-qualified, audited-route output yield in a prospective 16-run local-repair comparison under unchanged reaction and route allowances. The gain was accompanied by more distinct successful structures, lower guiding scores and reduced overall diversity. Historical component support, unreplicated contribution superiority and negative task transfer delimit the result. These findings establish a computational workflow improvement under specified assets, not demonstrated biological activity or synthesis.

## Data and code availability

Documented source code, test data, seed tables, policy outputs, saved route records and prediction vectors are retained in private reproduction packages. The portable primary-analysis package reconstructs the new study from saved data without model inference. Third-party assets retain their original provenance and access conditions. Stable public software availability and archival deposition remain pending; no public release or repository accession is claimed here.

## AI assistance disclosure

OpenAI Codex assisted with implementation, experiment orchestration, saved-record checking, analysis documentation and editorial drafting. Scientific outcomes reported here come from retained executed programs and saved observations, rather than language-model predictions. This English text is an AI-assisted working manuscript for researcher revision; author-written submission text and final author verification remain pending. Detailed assistance and verification records are retained with the supplementary preparation materials.

## References

Baell, J.B.; Holloway, G.A. (2010) New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening Libraries and for Their Exclusion in Bioassays. *Journal of Medicinal Chemistry*, **53**(7), 2719-2740. [doi:10.1021/jm901137j](https://doi.org/10.1021/jm901137j)

Bemis, G.W.; Murcko, M.A. (1996) The Properties of Known Drugs. 1. Molecular Frameworks. *Journal of Medicinal Chemistry*, **39**(15), 2887-2893. [doi:10.1021/jm9602928](https://doi.org/10.1021/jm9602928)

Bickerton, G.R.; Paolini, G.V.; Besnard, J.; Muresan, S.; Hopkins, A.L. (2012) Quantifying the chemical beauty of drugs. *Nature Chemistry*, **4**(2), 90-98. [doi:10.1038/nchem.1243](https://doi.org/10.1038/nchem.1243)

Breiman, L. (2001) Random Forests. *Machine Learning*, **45**(1), 5-32. [doi:10.1023/a:1010933404324](https://doi.org/10.1023/a:1010933404324)

Brenk, R.; Schipani, A.; James, D.; Krasowski, A.; Gilbert, I.H.; Frearson, J.; Wyatt, P.G. (2008) Lessons Learnt from Assembling Screening Libraries for Drug Discovery for Neglected Diseases. *ChemMedChem*, **3**(3), 435-444. [doi:10.1002/cmdc.200700139](https://doi.org/10.1002/cmdc.200700139)

Brown, N.; Fiscato, M.; Segler, M.H.S.; Vaucher, A.C. (2019) GuacaMol: Benchmarking Models for de Novo Molecular Design. *Journal of Chemical Information and Modeling*, **59**(3), 1096–1108. [doi:10.1021/acs.jcim.8b00839](https://doi.org/10.1021/acs.jcim.8b00839)

Chen, B.; Li, C.; Dai, H.; Song, L. (2020) Retro\*: Learning Retrosynthetic Planning with Neural Guided A\* Search. *Proceedings of Machine Learning Research*, **119**, 1608-1616. [Official proceedings](https://proceedings.mlr.press/v119/chen20k.html)

Chen, S.; Nam, G.; Aspuru-Guzik, A.; Jung, Y. (2026) SynTwins: a retrosynthesis-guided framework for synthesizable molecular analog generation. *Chemical Science*, **17**(4), 2255-2262. [doi:10.1039/d5sc05225d](https://doi.org/10.1039/d5sc05225d)

Ertl, P.; Schuffenhauer, A. (2009) Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. *Journal of Cheminformatics*, **1**, 8. [doi:10.1186/1758-2946-1-8](https://doi.org/10.1186/1758-2946-1-8)

Gao, W.; Coley, C.W. (2020) The Synthesizability of Molecules Proposed by Generative Models. *Journal of Chemical Information and Modeling*, **60**(12), 5714-5723. [doi:10.1021/acs.jcim.0c00174](https://doi.org/10.1021/acs.jcim.0c00174)

Genheden, S.; Bjerrum, E. (2022) PaRoutes: towards a framework for benchmarking retrosynthesis route predictions. *Digital Discovery*, **1**(4), 527-539. [doi:10.1039/d2dd00015f](https://doi.org/10.1039/d2dd00015f)

Genheden, S.; Thakkar, A.; Chadimová, V.; Reymond, J.-L.; Engkvist, O.; Bjerrum, E. (2020) AiZynthFinder: a fast, robust and flexible open-source software for retrosynthetic planning. *Journal of Cheminformatics*, **12**, 70. [doi:10.1186/s13321-020-00472-1](https://doi.org/10.1186/s13321-020-00472-1)

Huang, K.; Fu, T.; Gao, W.; Zhao, Y.; Roohani, Y.; Leskovec, J.; Coley, C.; Xiao, C.; Sun, J.; Zitnik, M. (2021) Therapeutics Data Commons: Machine Learning Datasets and Tasks for Drug Discovery and Development. *Proceedings of the Neural Information Processing Systems Track on Datasets and Benchmarks*, **1**. [Official proceedings](https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/hash/4c56ff4ce4aaf9573aa5dff913df997a-Abstract-round1.html)

Jensen, J.H. (2019) A graph-based genetic algorithm and generative model/Monte Carlo tree search for the exploration of chemical space. *Chemical Science*, **10**(12), 3567–3572. [doi:10.1039/C8SC05372C](https://doi.org/10.1039/C8SC05372C)

Loeffler, H.H.; He, J.; Tibo, A.; Janet, J.P.; Voronov, A.; Mervin, L.H.; Engkvist, O. (2024) Reinvent 4: Modern AI–driven generative molecule design. *Journal of Cheminformatics*, **16**(1), 20. [doi:10.1186/s13321-024-00812-5](https://doi.org/10.1186/s13321-024-00812-5)

Maziarz, K.; Tripp, A.; Liu, G.; Stanley, M.; Xie, S.; Gaiński, P.; Seidl, P.; Segler, M.H.S. (2025) Re-evaluating retrosynthesis algorithms with Syntheseus. *Faraday Discussions*, **256**, 568-586. [doi:10.1039/d4fd00093e](https://doi.org/10.1039/d4fd00093e)

Palczewska, A.; Palczewski, J.; Marchese Robinson, R.; Neagu, D. (2014) Interpreting Random Forest Classification Models Using a Feature Contribution Method. In: Bouabana-Tebibel, T. and Rubin, S.H. (eds), *Integration of Reusable Systems*. Advances in Intelligent Systems and Computing, **263**, 193–218. Springer, Cham. [doi:10.1007/978-3-319-04717-1_9](https://doi.org/10.1007/978-3-319-04717-1_9)

Rogers, D.; Hahn, M. (2010) Extended-Connectivity Fingerprints. *Journal of Chemical Information and Modeling*, **50**(5), 742-754. [doi:10.1021/ci100050t](https://doi.org/10.1021/ci100050t)

Southiratn, T.; Koo, B.; Lu, Y.; Kim, S. (2025) CombiMOTS: Combinatorial Multi-Objective Tree Search for Dual-Target Molecule Generation. *Proceedings of Machine Learning Research*, **267**, 56650-56691. [Official proceedings](https://proceedings.mlr.press/v267/southiratn25a.html)

Sun, J.; Jeliazkova, N.; Chupakhin, V.; Golib-Dzib, J.-F.; Engkvist, O.; Carlsson, L.; Wegner, J.; Ceulemans, H.; Georgiev, I.; Jeliazkov, V.; Kochev, N.; Ashby, T.J.; Chen, H. (2017) ExCAPE-DB: an integrated large scale dataset facilitating Big Data analysis in chemogenomics. *Journal of Cheminformatics*, **9**(1), 17. [doi:10.1186/s13321-017-0203-5](https://doi.org/10.1186/s13321-017-0203-5)

Wu, Z.; Ramsundar, B.; Feinberg, E.N.; Gomes, J.; Geniesse, C.; Pappu, A.S.; Leswing, K.; Pande, V. (2018) MoleculeNet: a benchmark for molecular machine learning. *Chemical Science*, **9**(2), 513–530. [doi:10.1039/C7SC02664A](https://doi.org/10.1039/C7SC02664A)

Wu, Z.; Wang, J.; Du, H.; Jiang, D.; Kang, Y.; Li, D.; Pan, P.; Deng, Y.; Cao, D.; Hsieh, C.-Y.; Hou, T. (2023) Chemistry-intuitive explanation of graph neural networks for molecular property prediction with substructure masking. *Nature Communications*, **14**(1), 2585. [doi:10.1038/s41467-023-38192-3](https://doi.org/10.1038/s41467-023-38192-3)

Yang, K.; Swanson, K.; Jin, W.; Coley, C.; Eiden, P.; Gao, H.; Guzman-Perez, A.; Hopper, T.; Kelley, B.; Mathea, M.; Palmer, A.; Settels, V.; Jaakkola, T.; Jensen, K.; Barzilay, R. (2019) Analyzing Learned Molecular Representations for Property Prediction. *Journal of Chemical Information and Modeling*, **59**(8), 3370-3388. [doi:10.1021/acs.jcim.9b00237](https://doi.org/10.1021/acs.jcim.9b00237)

Yin, X.; Wang, X.; Wu, Z.; Li, Q.; Kang, Y.; Deng, Y.; Luo, P.; Liu, H.; Shi, G.; Wang, Z.; Yao, X.; Hsieh, C.-Y.; Hou, T. (2025) Syn-MolOpt: a synthesis planning-driven molecular optimization method using data-derived functional reaction templates. *Journal of Cheminformatics*, **17**(1), 27. [doi:10.1186/s13321-025-00975-9](https://doi.org/10.1186/s13321-025-00975-9)

Zhong, W.; Yang, Z.; Chen, C.Y.-C. (2023) Retrosynthesis prediction using an end-to-end graph generative architecture for molecular graph editing. *Nature Communications*, **14**(1), 3009. [doi:10.1038/s41467-023-38851-5](https://doi.org/10.1038/s41467-023-38851-5)

## Figure legends

**Figure 1. Auditable local replacement and paired qualification comparison.** The same retained sources generate B and contribution-permuted C pools under common reaction, stock, QED and formation allowances. Each pool is read under source-preserving or fixed-task RF qualification. Parent planning uses at most 12 calls; at most two candidates receive six calls each and share 45 soft seconds. Saved native routes are audited before counting. External ten-checkpoint predictions are obtained after final output selection and do not guide selection. Qualification changes eligibility, not the molecular generator or reaction budget.

**Figure 2. Prospective run-level qualification gains and exploratory control comparison.** Panels a and b show external joint success proportions for source-preserving and fixed-task B, joined within each of eight GSK3B and eight JNK3 optimizer runs. Each run retains 64 terminal draws, including duplicates. Panel c shows the primary equal-task mean B qualification effect and exploratory C qualification and fixed-task B minus C effects, with 95% Welch–Satterthwaite intervals. The primary effect is +66.5 percentage points (50.5–82.5); the B/C interval crosses zero. Runs, rather than positions or checkpoint predictions, supply replication.
