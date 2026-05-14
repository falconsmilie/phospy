# Scientific Coverage

PhosPy implements selected PhosR-style phosphoproteomics workflows. It does not
claim full package equivalence with PhosR.

## Current Supported Public Lanes

PhosPy currently supports these public analysis lanes:

1. build an `AnalysisReadyPhosphoDataset`
2. run differential phosphorylation analysis
3. run kinase scoring and prediction
4. optionally run signalome analysis from the kinase result

Differential analysis requires analysis-ready numeric inputs plus valid
`ExperimentalDesign` and `Contrast` metadata. It does not infer design from
sample names and does not replace upstream preprocessing requirements.

### Differential Parity Envelope (Current Release)

`DifferentialAnalysisWorkflow` parity claims are currently scoped to:

- two-condition unpaired designs with biological-replicate rows
- simple condition-vs-condition contrasts with one `+1` and one `-1` term
- empirical-Bayes modes: `method="standard"` and `method="robust"` with
  optional `trend=True`
- Benjamini-Hochberg multiple-testing adjustment (`adj.P.Val`)

Explicitly unsupported in this release:

- batch-aware differential modelling (`batch`)
- block/paired/repeated-measure differential modelling (`block`)

Contract difference vs limma/PhosR surface:

- analysis-ready inputs must be complete at boundary; missing values are
  rejected before differential execution instead of being handled inside
  differential model fitting.

Bundled runtime references in the current release are rat-only. Human and mouse
analysis can be run by passing an explicit `ReferenceBundle` in Python.

## Scientific Confidence Labels

Use these labels when discussing coverage confidence:

| Label | Meaning |
| --- | --- |
| `PARITY_GATED_ACTIVE_SCIENCE` | protected by active fixture-backed parity tests |
| `PHOSPY_VALIDATED_SCIENCE` | validated by PhosPy contract, unit, and integration tests |
| `SUPPORTED_CONTRACT_CHANGED` | intentionally supported with a changed public contract |
| `OPEN_GAP` | not yet covered or not yet claimed |

## Intended PhosR Parity Scope Labels

Use these labels to describe parity intent, separate from current confidence:

| Scope label | Meaning |
| --- | --- |
| `required parity` | PhosPy is expected to match PhosR-style behavior for the scoped lane |
| `deliberate scope difference` | PhosPy intentionally keeps a different surface or contract |
| `useful future extension` | scientifically useful but not currently in the supported core lane |
| `not planned` | intentionally outside PhosPy scope |

## PhosR Parity Scope Matrix

Parity is feature-specific, not global. "Supported" means supported in the
documented PhosPy contract, not full PhosR equivalence.

| Area | Current confidence | Intended PhosR parity scope | Current PhosPy support | Notes / limits | Test expectation |
| --- | --- | --- | --- | --- | --- |
| Input formats | `SUPPORTED_CONTRACT_CHANGED` | `deliberate scope difference` | Tabular/DataFrame-first boundary with typed request models, not Bioconductor object classes | PhosPy keeps explicit Python-native contracts (`AnalysisReadyPhosphoDataset`, typed workflow requests/results) | Contract tests for request/result models and boundary validation (`tests/unit/test_public_contract_*.py`) |
| Phosphosite representation | `PHOSPY_VALIDATED_SCIENCE` | `required parity` | Strict site-ID and site-metadata validation at dataset boundary | Site identity must be explicit and aligned; peptide/site ambiguity policies remain bounded by contract | Unit/integration contract tests plus parity tests for site-sensitive lanes (`tests/parity/test_preprocessing_science_parity.py`) |
| Site/flanking sequence | `PARITY_GATED_ACTIVE_SCIENCE` | `required parity` | `site_sequence` is required at analysis-ready boundary and used in scoring/prediction | Sequence quality is an upstream dependency; parity claims are lane-scoped | Fixture-backed parity checks for sequence-sensitive scoring/prediction outputs (`tests/parity/test_l6_prediction_parity.py`, `tests/parity/test_prediction_science_parity.py`) |
| Localisation confidence | `PHOSPY_VALIDATED_SCIENCE` | `useful future extension` | Optional localisation confidence fields and validation policies are supported | No first-class localisation-filter workflow parity claim in current lane | Unit tests for localisation validation boundaries (`tests/unit/test_localisation_policy_preprocessing.py`, `tests/unit/test_validator_boundaries.py`) |
| Replicate/condition modelling | `PHOSPY_VALIDATED_SCIENCE` | `required parity` | Typed `ExperimentalDesign` and `Contrast` contract for differential workflows | No implicit condition inference from sample names; technical replicates require explicit policy (`reject`, `mean`, `median`) | Contract and boundary tests for design/contrast and technical-replicate validation (`tests/unit/test_experimental_design_contract.py`, `tests/unit/test_differential_analysis.py`, `tests/integration/test_differential_with_technical_replicates.py`) |
| Missing-value handling | `PARITY_GATED_ACTIVE_SCIENCE` | `required parity` | Missing-data preprocessing policies plus downstream score preconditioning are supported | Policy choice changes retained rows and downstream behavior | Unit policy tests plus preprocessing parity fixture checks (`tests/unit/test_missing_data_split_modules.py`, `tests/parity/test_preprocessing_science_parity.py`) |
| Normalisation | `PARITY_GATED_ACTIVE_SCIENCE` | `required parity` | Supported preprocessing methods include `none`, `median_center`, `quantile` with stage-order provenance | Method-by-method scope; not a blanket normalization equivalence claim | Unit preprocessing tests plus preprocessing parity tests (`tests/unit/test_dataset_preprocessing_subsystem.py`, `tests/parity/test_preprocessing_science_parity.py`) |
| Imputation | `PHOSPY_VALIDATED_SCIENCE` | `required parity` | Explicit preprocessing imputation policies are supported (`row_median`, `minprob`, `knn`) | Imputation is policy-dependent and scientifically consequential; parity claims stay lane-scoped | Unit tests for imputation behavior and invariants (`tests/unit/test_dataset_preprocessing_subsystem.py`, `tests/unit/test_scientific_invariants.py`) |
| Batch correction | `OPEN_GAP` | `useful future extension` | RUV readiness diagnostics exist; executable SPS/RUV batch correction is not in the supported public lane | Must not be described as supported until executable workflow + parity evidence exists | Keep as open gap; add fixture-backed parity tests only if/when executable lane is implemented |
| Differential phosphorylation | `PARITY_GATED_ACTIVE_SCIENCE` | `required parity` | Public `DifferentialAnalysisWorkflow` route with design/contrast validation, per-site statistics, and multiple-testing adjustment | Parity envelope is explicitly scoped to two-condition unpaired simple contrasts, with batch/block/paired modes rejected in this release; analysis-ready boundary rejects missing values before differential execution | Unit and integration tests for design validation, unsupported-mode rejection, contrast orientation/sign, adjusted p-values, result alignment, and fixture-backed limma parity comparison (`tests/unit/test_differential_analysis.py`, `tests/unit/test_differential_unsupported_designs.py`, `tests/unit/test_differential_result_contract.py`, `tests/parity/test_differential_analysis_parity.py`, `tests/parity/test_differential_limma_parity.py`) |
| Kinase/substrate analysis | `PARITY_GATED_ACTIVE_SCIENCE` | `required parity` | Public kinase scoring/prediction lane with active fixture-backed parity checks | Parity is scoped to documented fixture lanes and published comparison rules | Active parity tests in `tests/parity/test_kinase_workflow_parity.py`, `tests/parity/test_l6_prediction_parity.py`, `tests/parity/test_public_predmat_parity.py` |
| Motif/sequence-aware analysis | `PARITY_GATED_ACTIVE_SCIENCE` | `required parity` | Rank-weighted motif/profile evidence fusion and sequence-aware scoring are supported | Relative-support ranking, not calibrated causal inference | Unit scientific-policy tests and fixture-backed parity checks (`tests/unit/test_prediction_science.py`, `tests/parity/test_prediction_science_parity.py`) |
| Enrichment analysis | `SUPPORTED_CONTRACT_CHANGED` | `deliberate scope difference` | KSEA-style substrate-set enrichment is supported in kinase activity lane | Broader pathway/gene-set enrichment utilities are not currently a first-class core lane | Unit tests for KSEA-style activity behavior (`tests/unit/test_activity_science.py`) |
| Clustering/time-series | `PARITY_GATED_ACTIVE_SCIENCE` | `useful future extension` | Signalome clustering/network construction is supported for documented lanes; no first-class time-series workflow contract | Time-series analysis remains out of current public workflow contract | Signalome parity and backend parity tests (`tests/parity/test_signalome_workflow_parity.py`, `tests/parity/test_signalome_clustering_backend_parity.py`) |
| Visualisation | `OPEN_GAP` | `deliberate scope difference` | No first-class visualization API in core PhosPy workflows | Visualization is intentionally not treated as core scientific parity in this release | Keep documented as scope difference unless a separate visualization scope is approved |
| Reproducibility/reporting | `PHOSPY_VALIDATED_SCIENCE` | `required parity` | Provenance payloads, fixture provenance docs, and reloadable output bundles are supported | Reproducibility claims depend on documented provenance and fixture governance | Integration and provenance regression tests (`tests/integration/test_kinase_workflow_integration.py`, `tests/integration/test_signalome_workflow_integration.py`, `tests/unit/test_provenance_regressions.py`) |
| Workflow composition/extensibility | `SUPPORTED_CONTRACT_CHANGED` | `deliberate scope difference` | Governed public composition lanes are explicit (`DatasetBuilder`, `DifferentialAnalysisWorkflow`, `KinaseWorkflow`, `SignalomeWorkflow`) | PhosPy does not mirror full PhosR package layering or object model | Public-contract and integration workflow tests (`tests/unit/test_public_contract_workflows.py`, `tests/integration/test_dataset_builder_integration.py`, `tests/integration/test_kinase_workflow_integration.py`, `tests/integration/test_signalome_workflow_integration.py`) |
| Full PhosR package equivalence | `OPEN_GAP` | `not planned` | Not claimed | PhosPy parity is lane-scoped and evidence-scoped; full-package equivalence is intentionally outside scope | Guardrail documentation plus parity suite scoping in `docs/parity.md` and `tests/parity/` |

Required-parity rows must maintain explicit test expectations and should not be
promoted to parity-equivalent language without fixture provenance and active
comparison tests.

## Interpretation Limits

- Weighted activity output (`simplified_weighted_substrate_activity_v1`) is a
  heuristic summary over predicted substrates above threshold/top-N support.
- KSEA-style activity output (`ksea_zscore_activity_v1`) applies unweighted
  substrate-set enrichment z-scores after evidence thresholding and reports
  p-values (and q-values when enabled).
- KSEA-style activity is not equivalent to full PhosR kinase activity inference.
- Rank-weighted fusion scores combine profile-correlation and motif-frequency
  evidence using rank-derived weights.
- Signalome module/network scores are derived summaries, not probabilities,
  calibrated confidence values, or causal proof.
- Missing kinase correlations stay missing. `0.0` means a finite near-zero
  correlation was estimated.
- Differential phosphorylation results depend on valid design matrices, contrast
  definitions, replicate structure, and upstream preprocessing quality.
- Differential analysis does not resolve peptide/site ambiguity, localisation
  confidence, imputation, normalisation, or batch correction unless those steps
  were already performed or explicitly configured in the route.
- Adjusted p-values control false discovery rate according to the implemented
  correction method; they do not validate biological causality.

## Scientific Policy Records

Workflow provenance includes machine-readable `scientific_policies` records.
Each record carries:

- stable policy ID
- name and version
- plain-language description
- active parameters
- scientific assumptions
- output scale/meaning

Ownership of scientific policy modules is domain-scoped:

- shared models: `phospy.provenance.scientific_policy_models`
- prediction: `phospy.science.prediction.scientific_policies`
- activities: `phospy.science.activities.scientific_policies`
- preprocessing: `phospy.science.datasets.preprocessing.scientific_policies`
- signalome workflow: `phospy.workflows.signalome.scientific_policies`
- signalome clustering: `phospy.science.signalomes.clustering.scientific_policies`
- differential aggregation: `phospy.science.differential.aggregation.scientific_policies`

Differential outputs now expose structured policy provenance through
`DifferentialAnalysisResult.policy_provenance`, including:

- design formula and design-matrix summary
- explicit contrast definitions
- replicate/group requirements and technical-replicate lineage
- empirical-Bayes moderation settings
- p-value and adjusted p-value methods
- missing-value handling policy
- intentionally rejected unsupported design features (batch/block/paired)

### `profile_correlation_shifted_unit_v1`

- What it does:
  transforms profile correlations from `[-1, 1]` to `[0, 1]` using `(r + 1) / 2`.
- Assumptions:
  positive correlation increases support.
- Parameters:
  transform formula, clipping to `[0, 1]`, preserve undefined values as missing.
- Output meaning:
  relative support score; larger means stronger positive agreement.
- Output does not mean:
  calibrated probability or direct evidence of inhibition/activation.
  Negative correlations are treated as lower support, not explicit inhibitory
  evidence.

### `kinase_profile_scoring_v1`

- What it does:
  records kinase profile-construction and scoring behavior, including
  self-inclusion vs leave-one-out semantics.
- Assumptions:
  profile rows can include the same substrate site later scored unless a
  leave-one-out policy is explicitly enabled.
- Parameters:
  profile missing-value strategy, self-inclusion behavior, leave-one-out flag,
  and scoring substrate floors.
- Output meaning:
  explicit provenance of the profile-scoring policy context used for
  downstream support scores.

### `motif_profile_rank_fusion_v1`

- What it does:
  fuses motif-frequency and profile-correlation evidence using rank-derived
  logarithmic weights.
- Assumptions:
  motif-library size and quantified-substrate count proxy evidence strength.
- Parameters:
  motif/profile weight formulas and fallback/diagnostic flags.
- Output meaning:
  relative downstream support for kinase-site ranking.
- Output does not mean:
  statistical enrichment p-value or calibrated confidence.

### `simplified_weighted_substrate_activity_v1`

- What it does:
  computes prediction-weighted activity and thresholded substrate-mean activity.
- Assumptions:
  predicted substrate support can summarize relative kinase activity in-run.
- Parameters:
  threshold, `min_substrates`, `top_n_substrates`, and explicit scoring rules.
- Output meaning:
  relative sample-by-kinase activity summaries.
- Output does not mean:
  full KSEA-style enrichment statistics.

### `ksea_zscore_activity_v1`

- What it does:
  computes KSEA-style z-score substrate-set enrichment activity.
- Assumptions:
  kinase substrate membership is unweighted after evidence thresholding.
- Parameters:
  evidence threshold, minimum substrates, z-score formula, p-value method, and
  optional q-value adjustment.
- Output meaning:
  statistically interpretable substrate-set enrichment activity z-scores with
  accompanying p-values.
- Output does not mean:
  PhosR-equivalent kinase activity inference.

### `candidate_substrate_selection_v1`

- What it does:
  records candidate substrate filtering for kinase prediction.
- Assumptions:
  top-k ranking, score-threshold filtering, and inclusion floor jointly define
  usable candidate support.
- Parameters:
  `top_k`, threshold rule, threshold value, inclusion floor, and site
  restriction behavior.
- Output meaning:
  explicit provenance of the candidate-selection rule that gates kinase ranking
  and prediction outputs.

### `signalome_module_candidate_score_v1`

- What it does:
  ranks candidate module counts using within-cluster correlation summaries.
- Assumptions:
  stronger within-cluster profile coherence indicates better candidate module
  structure.
- Parameters:
  requested/resolved candidate-scoring policies, mode, guards, and skip/evaluated
  diagnostics.
- Output meaning:
  candidate module-count support score used for ranking/selection.
- Output does not mean:
  biological certainty or causal regulation evidence.

### `signalome_missing_value_clustering_v1`

- What it does:
  records missing-value handling for clustering distance/tree inputs.
- Assumptions:
  non-finite values are normalized to missing; missing values are imputed for
  clustering internals.
- Parameters:
  missing-value policy name, applicability scope, and whether imputed values
  appear in output tables.
- Output meaning:
  explicit provenance for clustering-matrix preparation rules.

### `signalome_score_preconditioning_v1`

- What it does:
  records row-retention policy for downstream score preconditioning before
  signalome execution.
- Assumptions:
  all-missing rows are unsupported and can be dropped or treated as boundary
  errors depending on policy.
- Parameters:
  preconditioning policy, row-retention rule, and input/dropped/retained row
  counts.
- Output meaning:
  explicit provenance for score-row retention behavior that can change site
  coverage and assignments.

### `protein_module_from_site_membership_v1`

- What it does:
  derives protein module IDs from site-cluster membership incidence patterns.
- Assumptions:
  shared site-cluster membership reflects shared protein-level module context.
- Parameters:
  membership-vector representation and module ID assignment rule.
- Output meaning:
  integer protein module IDs for grouping.
- Output does not mean:
  direct mechanistic proof of shared regulation.

### `preprocessing_stage_order_v1`

- What it does:
  records explicit preprocessing stage order used to construct analysis-ready
  dataset inputs.
- Assumptions:
  stage order is scientifically meaningful and can change transformed values,
  row retention, and derived comparison outputs.
- Parameters:
  configured stage order, default order, and supported stage order metadata.
- Output meaning:
  explicit provenance for preprocessing execution order.

### `peptide_to_site_aggregation_v1`

- What it does:
  records how peptide-level differential statistics are aggregated to site-level
  summaries.
- Assumptions:
  aggregation strategy and variance rules change site-level uncertainty and
  significance behavior.
- Parameters:
  aggregation strategy, minimum peptides per site, missing-variance policy, and
  weighting mode.
- Output meaning:
  explicit provenance for site-level differential summary construction.

## Where Details Live

- [Scientific Coverage](scientific-coverage.md) is the maintained scope and
  coverage matrix.
- [Parity](parity.md) tracks PhosR comparison evidence, fixture locations, and
  parity test references.
- [Performance Contracts](performance.md) covers scale limits.
- [ADR Index](adr/index.md) stores maintainer decision records.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi,
R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.

Xiao, D., Yang, P., & Kim, H. J. (2026). *PhosR* (R package manual).
Bioconductor. https://bioconductor.org/packages/release/bioc/manuals/PhosR/man/PhosR.pdf

Xiao, D., Yang, P., & Kim, H. J. (2026). *An introduction to PhosR package*
(Bioconductor vignette). https://bioconductor.org/packages/release/bioc/vignettes/PhosR/inst/doc/PhosR.html
