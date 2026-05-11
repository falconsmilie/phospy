# Scientific Coverage

PhosPy implements selected PhosR-style phosphoproteomics workflows. It does not
claim full package equivalence with PhosR.

## Current Supported Lane

The supported public lane is:

1. build an `AnalysisReadyPhosphoDataset`
2. run kinase scoring and prediction
3. optionally run signalome analysis from the kinase result

Bundled runtime references in the current release are rat-only. Human and mouse analysis can
be run by passing an explicit `ReferenceBundle` in Python.

## Scientific Confidence Labels

Use these labels when discussing coverage:

| Label | Meaning |
| --- | --- |
| `PARITY_GATED_ACTIVE_SCIENCE` | protected by active fixture-backed parity tests |
| `PHOSPY_VALIDATED_SCIENCE` | validated by PhosPy contract, unit, and integration tests |
| `SUPPORTED_CONTRACT_CHANGED` | intentionally supported with a changed public contract |
| `OPEN_GAP` | not yet covered or not yet claimed |

## Active Coverage

| Area | Current status |
| --- | --- |
| Dataset boundary | strict PhosPy dataset contract |
| Kinase scoring and prediction | active fixture-backed parity and workflow tests |
| Activity output | supports weighted heuristic activity and KSEA-style z-score substrate-set enrichment activity |
| Signalome workflow | supported from kinase result with explicit `protein_id` |
| Output publishing | supported simple publishers and reloadable bundle services |
| Human/mouse bundled references | open gap for bundled runtime data in this release |

## PhosR compatibility and scope

This matrix is maintained to prevent "PhosR-inspired" from being read as global
PhosR parity.

- Parity is feature-specific, not global.
- "Supported" means supported in the documented PhosPy contract, not full PhosR
  equivalence.
- PhosR-equivalent claims require source-labelled parity fixtures and parity
  tests for that feature lane.

| Area | PhosR capability | PhosPy current support | Current status | Gap type | Gap class | Priority | Planned action | Evidence / tests | References |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Input object/data model | Bioconductor object-oriented workflow (`SummarizedExperiment` / `PhosphoExperiment`) and function-driven processing/downstream interfaces. | `AnalysisReadyPhosphoDataset` plus typed workflow requests/results and explicit validation boundary. | deliberate scope difference | Deliberate scope difference | architectural gap | maintain | Keep explicit contract-bound Python-native model; do not mirror R object stack. | `docs/workflow_contracts.md`; `docs/adr/adr_0001_public_api_contract.md`; `tests/unit/test_public_contract_*.py` | PhosR vignette; PhosR package page |
| Phosphosite representation | Site-centric workflows with explicit phosphosite IDs and sequence context used downstream. | Standard site ID contract with strict site metadata and boundary validation. | partially supported | Required parity | scientific gap | high | Maintain current site-ID contract and add/keep fixture-backed comparisons for site identity-sensitive lanes. | `docs/adr/adr_0018_phosphosite_identity_and_localisation_policy.md`; `tests/parity/test_preprocessing_science_parity.py` | PhosR vignette; ADR-0018 references |
| Site/flanking sequence | Sequence-aware scoring lanes (`kinaseSubstrateScore`) and motif-aware downstream analysis. | `site_sequence` is required at analysis-ready boundary and consumed by kinase scoring/prediction. | supported | Required parity | implementation gap | high | Keep parity checks on sequence-sensitive prediction/scoring outputs. | `docs/workflow_contracts.md`; `tests/parity/test_l6_prediction_parity.py`; `tests/parity/test_prediction_science_parity.py` | PhosR vignette; PhosR manual |
| Localisation probability | Localisation confidence can be carried with phosphosite records and used in downstream filtering/QC decisions. | Optional `localisation_probability` with validator policies (`allow_unknown`, `require_present`, `require_threshold`). | partially supported | Useful future extension | scientific gap | medium | Keep validator-level controls; add fixture-backed parity only if localisation-based filtering becomes a first-class lane. | `docs/adr/adr_0018_phosphosite_identity_and_localisation_policy.md`; `tests/unit/test_validator_boundaries.py` | ADR-0018 references |
| Replicate/condition modelling | Experimental design drives downstream comparisons and modelling choices. | Typed `ExperimentalDesign` and `Contrast` contract; no condition inference from sample names. | partially supported | Required parity | architectural gap | high | Keep typed contract; add parity fixtures per supported design class only. | `docs/api/differential-workflow.md`; `docs/adr/adr_0019_experimental_design_and_contrast_contract.md`; `tests/parity/test_differential_analysis_parity.py` | PhosR vignette; limma references in ADR-0019 |
| Filtering | Processing-stage filtering and downstream selection steps (for example dynamic-site selection in tutorials/workflows). | Supported preprocessing filters and workflow-specific selection policies; no claim of full PhosR filtering surface. | partially supported | Required parity | implementation gap | medium | Add fixture-backed parity only for explicitly documented filtering rules. | `docs/workflow_contracts.md`; `tests/parity/test_preprocessing_science_parity.py`; `tests/unit/test_dataset_preprocessing_subsystem.py` | PhosR vignette |
| Missing value handling | Multiple imputation/handling options used as part of preprocessing and batch-correction preparation. | Explicit missing-data policies (`forbid`, row median, MinProb, KNN) plus signalome score preconditioning diagnostics. | supported | Required parity | implementation gap | high | Keep policy-level tests; expand parity fixtures only for lanes claimed equivalent. | `docs/workflow_contracts.md`; `docs/validation.md`; `tests/parity/test_preprocessing_science_parity.py` | PhosR package page; PhosR vignette |
| Normalisation | Processing includes normalisation methods and integrated preprocessing workflows. | Dataset preprocessing supports `none`, `median_center`, `quantile`; explicit transform/normalisation order provenance. | partially supported | Required parity | implementation gap | high | Maintain current methods; add parity fixtures method-by-method before equivalence claims. | `docs/workflow_contracts.md`; `docs/scientific-coverage.md` (`preprocessing_stage_order_v1`); `tests/parity/test_preprocessing_science_parity.py` | PhosR package page; PhosR vignette |
| Batch correction / RUV / SPS | Dedicated SPS-guided `RUVphospho` batch-correction workflow. | RUV readiness reporting exists; no executable RUV/SPS batch correction in current supported lane. | not supported | Useful future extension | scientific gap | high | Treat as future extension; do not describe as supported until executable workflow + fixture-backed parity exists. | `docs/validation.md` (RUV readiness is report-only); `docs/workflow_contracts.md` | PhosR vignette (batch correction); PhosR manual (`RUVphospho`) |
| Differential phosphorylation | Differential analysis lane with empirical Bayes and contrast-based outputs. | First-class `DifferentialAnalysisWorkflow` with parity tests for selected limma-style fixtures; batch/block modelling remains non-executable. | partially supported | Required parity | implementation gap | high | Keep lane release-ready only per documented contract + fixture coverage; expand feature support deliberately. | `docs/api/differential-workflow.md`; `tests/parity/test_differential_analysis_parity.py`; `tests/fixtures/rewrite_parity/differential_r_reference/PROVENANCE.md` | limma references in ADR-0019 |
| Kinase/substrate analysis | Kinase-substrate scoring and prediction from motif and profile evidence. | Supported scoring/prediction lane with active parity gates on L6/public fixtures. | supported | Required parity | implementation gap | high | Continue parity-gated maintenance for scoring/prediction seams. | `docs/parity.md`; `tests/parity/test_kinase_workflow_parity.py`; `tests/parity/test_l6_prediction_parity.py`; `tests/parity/test_public_predmat_parity.py` | PhosR vignette; PhosR manual |
| Signalome construction | Signalome construction from kinase-substrate scoring/prediction outputs and network/module summarisation. | Supported `SignalomeWorkflow` lane with fixture-backed output checks for defined L6/public reference lanes only. | requires parity fixture | Required parity | implementation gap | high | Do not claim full PhosR signalome parity outside source-labelled fixture lanes; expand only with explicit fixture provenance and parity tests. | `docs/parity.md`; `tests/parity/test_signalome_workflow_parity.py`; `tests/fixtures/public_workflow_reference/PROVENANCE.md` | PhosR vignette; PhosR manual |
| Motif/sequence-aware analysis | Motif-aware kinase-substrate relationship scoring integrated with dynamic profiles. | Rank-weighted motif/profile fusion is supported; parity is lane-specific, not full-method global equivalence. | partially supported | Required parity | scientific gap | medium | Keep current motif/profile lane claims bounded to fixture-backed seams. | `docs/scientific-coverage.md` (`motif_profile_rank_fusion_v1`); `tests/parity/test_prediction_science_parity.py`; `tests/fixtures/rewrite_parity/fragile_support_reference/PROVENANCE.md` | PhosR vignette; PhosR manual |
| Enrichment analysis | Gene- and phosphosite-centric enrichment utilities are part of broader PhosR downstream analysis. | Kinase activity includes KSEA-style substrate-set enrichment only; broader pathway enrichment utilities are not in core PhosPy lane. | partially supported | Deliberate scope difference | scientific gap | medium | Keep KSEA-style support explicit; classify broader enrichment as future extension, not parity-by-default. | `docs/workflow_contracts.md` (KSEA limitations); `src/phospy/activities/methods/ksea_zscore.py`; `tests/unit/test_activity_science.py` | PhosR vignette (pathway enrichment); Yang et al., 2019 |
| Clustering / time-series | Signalome clustering plus additional temporal/pathway-oriented analyses in broader package workflows. | Signalome clustering/network construction is supported for documented lanes; no first-class time-series workflow contract. | partially supported | Useful future extension | architectural gap | medium | Keep clustering lane parity-scoped; treat time-series as separate future workflow design. | `docs/workflow_contracts.md`; `tests/parity/test_signalome_workflow_parity.py`; `tests/parity/test_signalome_clustering_backend_parity.py` | PhosR vignette |
| Visualisation | Built-in plotting/visualisation helpers for pathway and signalome interpretation. | No first-class visualization API/plotting contract in core PhosPy workflows. | not planned | Unnecessary feature creep | UX/documentation gap | low | Keep out of core scientific contract unless a separate visualization scope is approved. | `docs/workflow_contracts.md`; `README.md` (supported lane focus) | PhosR vignette |
| Reproducibility / provenance | Reproducible analysis context and interpretable outputs are emphasized in documented workflows. | Strong provenance contracts, deterministic fixture lanes, exact/tolerance hashes, and bundle manifests. | supported | Required parity | implementation gap | high | Keep provenance goldens and parity fixture provenance docs release-gated. | `docs/parity.md`; `tests/fixtures/public_workflow_reference/PROVENANCE.md`; `tests/integration/test_kinase_workflow_integration.py`; `tests/integration/test_signalome_workflow_integration.py` | PhosR vignette; Yang et al., 2019 |
| Workflow composition / extensibility | End-to-end processing and downstream workflow composition in package ecosystem. | Public composition lane exists (`DatasetBuilder -> KinaseWorkflow -> SignalomeWorkflow`) with separate differential lane and typed requests. | supported | Deliberate scope difference | architectural gap | medium | Continue explicit contracts and avoid inferring support from internal helpers alone. | `README.md`; `docs/workflow_contracts.md`; `docs/api/guide.md`; `tests/integration/test_*workflow*_integration.py` | PhosR package page |

### Gap Class Legend

- `scientific gap`: missing science method/assumption coverage
- `architectural gap`: contract/model shape differs from PhosR interface style
- `implementation gap`: intended scope exists but is not fully delivered or parity-gated
- `UX/documentation gap`: presentation or discoverability gap without core scientific engine impact

### PhosR Equivalence Guardrail

Rows marked `requires parity fixture` (or planned to require one before parity
language) must not be described as PhosR-equivalent until:

1. fixture source is labelled,
2. comparison rule/tolerance is documented, and
3. parity tests are active in `tests/parity/`.

## Interpretation Limits

- Weighted activity output (`simplified_weighted_substrate_activity_v1`) is a
  heuristic summary over predicted substrates above threshold/top-N support.
- KSEA-style activity output (`ksea_zscore_v1`) applies unweighted substrate-set
  enrichment z-scores after evidence thresholding and reports p-values (and
  q-values when enabled).
- KSEA-style activity is not equivalent to full PhosR kinase activity inference.
- Rank-weighted fusion scores combine profile-correlation and motif-frequency
  evidence using rank-derived weights.
- Signalome module/network scores are derived summaries, not probabilities,
  calibrated confidence values, or causal proof.
- Missing kinase correlations stay missing. `0.0` means a finite near-zero
  correlation was estimated.

## Scientific Policy Records

Workflow provenance includes machine-readable `scientific_policies` records.
Each record carries:

- stable policy ID
- name and version
- plain-language description
- active parameters
- scientific assumptions
- output scale/meaning

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

## Where Details Live

- [Parity](parity.md) tracks PhosR comparison evidence and fixture locations.
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
