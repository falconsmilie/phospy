# Scientific Coverage

PhosPy implements selected PhosR-style phosphoproteomics workflows.
"PhosR-style" and "PhosR-inspired" in this repository always mean
feature-scoped, evidence-scoped comparison lanes. They do not mean full PhosR
package equivalence.

## Current Supported Public Lanes

PhosPy currently supports these public analysis lanes:

1. build an `AnalysisReadyPhosphoDataset`
2. run differential phosphorylation analysis
3. run kinase scoring and prediction
4. optionally run signalome analysis from the kinase result

The prediction science layer includes a pure Kinase Library-style motif scorer.
`KinaseWorkflow` exposes it through explicit scoring modes only; the default
kinase lane remains the PhosR-style rank-weighted scoring mode. Kinase Library
workflow scoring requires a caller-supplied `KinaseLibraryResource` with
compatible organism, residue-class lanes, score matrices, sequence-window
definition, and provenance.

Differential analysis requires analysis-ready numeric inputs plus valid
`ExperimentalDesign` and `Contrast` metadata. It does not infer design from
sample names and does not replace upstream preprocessing requirements.

At dataset-construction boundary, PhosPy uses a protein-scoped analysis-ready
row key (`site_key`) and keeps `display_id` (for example `GENE;SITE;`) as a
human-readable label. `site_key` is required to be unique, while `display_id`
may repeat once `site_key` is the row identity. Direct analysis-ready datasets
must use `site_key` indexes and include auditable protein context metadata
(`organism`, `protein_namespace`, `protein_identifier`, and `site`); they must
not silently fall back to display-site identity. Builder input may accept legacy
display-indexed shape only when enough protein context exists to derive
`site_key`. Workflows operate on `site_key`, and site-level outputs that
materialize row identity include both `site_key` and `display_id`. `display_id`
is a human-readable label and may repeat. Rows that resolve to the same
`site_key` are a scientific ambiguity; the default duplicate-site policy fails,
and non-error policies are deliberate preprocessing choices that change which
evidence enters the analysis-ready dataset. Differential result tables and
direct public `DifferentialAnalysisResult` construction require encoded
protein-scoped `site_key` indexes plus explicit `site_key`, `display_id`,
`organism`, `protein_namespace`, `protein_identifier`, `gene_symbol`, and `site`
columns. Workflow-created differential results preserve that required protein
context and optional protein metadata such as `protein_id` when present.
Validation fails rather than repairing display-indexed, display-keyed,
arbitrary-keyed, or stat-only public result tables. Kinase references may use
display IDs only through the explicit
reference-to-dataset mapping layer; references remain reference/display
identifiers. Analysis-ready datasets may carry `site_metadata.protein_id` as
optional metadata, and it may be absent or incomplete at the base dataset
boundary. Signalome requires complete `site_metadata.protein_id` values as
algorithm-specific protein grouping metadata; that field is not the dataset row
identity and is not encoded in `site_key`. See
[ADR-0024: Protein-Scoped Phosphosite Row Identity](adr/adr_0024_protein_scoped_phosphosite_row_identity.md).

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
analysis can be run by passing an explicit `ReferenceBundle` in Python. No
packaged human or mouse reference lane is committed in this release because no
approved redistributable source bundle with complete license and provenance
metadata is included.

## Roadmap Visibility and Guardrails

[ADR-0025: Competitive Phosphoproteomics Workflow Coverage Roadmap](adr/adr_0025_competitive_phosphoproteomics_workflow_coverage.md)
records intended future direction for references, kinase inference, importers,
richer differential designs, enrichment, visualisation, and possible CLI
support.

Roadmap entries are not current feature claims. A roadmap item becomes supported
only when executable implementation, public contracts, documentation, and tests
exist, and this page is updated to the correct scope category.

| Roadmap area | Current status | Direction, not current support |
| --- | --- | --- |
| References | Bundled runtime references are rat-only. Human and mouse workflows require an explicit caller-supplied `ReferenceBundle`. | Broader reference handling should use explicit provenance, compatibility checks, and external bundle validation. New bundled data requires redistribution permission, provenance, docs, and tests before `_BUNDLED_DEFAULTS` is updated. |
| Kinase inference | Kinase scoring/prediction and two explicit activity methods are executable. Scores are relative support or substrate-set summaries, not calibrated causal inference. | Additional kinase inference or activity methods should be added one method at a time with stable scientific policy records and method-specific validation. |
| Importers | PhosPy supports analysis-ready tables and generic table I/O contracts used by Python workflows. It does not currently provide broad semantic importers for vendor, search-engine, or upstream statistical outputs. | Semantic importers should produce typed tables or requests that still pass builder and workflow validation; they must not bypass site identity or provenance contracts. |
| Richer differential designs | Current parity-protected differential lane is two-condition unpaired simple contrasts. Batch-aware, block, paired, and repeated-measure modelling are not executable in this release. | Richer designs require explicit design/result contracts, provenance, validation, and parity or method-specific evidence before any support claim. |
| Enrichment | KSEA-style substrate-set activity exists only inside the kinase activity lane. Broader pathway or gene-set enrichment is not a current core workflow lane. | Enrichment should be a separately contracted workflow or method, with clear null model, input universe, and output-scale documentation. |
| Visualisation | Core PhosPy has no first-class visualisation workflow/API. | Visualisation should consume validated result objects and must not become a hidden analysis engine or source of scientific truth. |
| CLI workflow support | Scientific workflow execution through a CLI is not currently supported; the Python API is the supported interface. | Any future CLI must be a thin wrapper over Python API requests/workflows and satisfy ADR-0022 reintroduction criteria before support is claimed. |

## Scope Categories

Every scientific scope claim in public docs must map to one category:

| Category | Meaning |
| --- | --- |
| `parity-gated` | Executable lane with active fixture-backed parity checks in release gates |
| `validated PhosPy implementation` | Executable lane validated by PhosPy contract/unit/integration tests; not a PhosR-equivalence claim by itself |
| `experimental` | Executable but intentionally provisional/approximate behavior with explicit caveats |
| `open gap` | Not currently executable in the supported public workflow lane |
| `deliberate scope difference` | Intentionally different from PhosR surface or intentionally narrowed contract |
| `not planned` | Intentionally outside supported scope |

## Scientific Scope Matrix (Single Source Of Truth)

This matrix is the maintained user-facing scope source for PhosPy. Parity is
feature-specific and evidence-scoped. Full PhosR package equivalence is not
claimed.

| Area | Scope category | Current executable support | Evidence and release checks | Limits and non-claims |
| --- | --- | --- | --- | --- |
| Differential analysis | `parity-gated` | `DifferentialAnalysisWorkflow` for two-condition unpaired simple contrasts with empirical-Bayes `standard`/`robust` and optional `trend` | `tests/parity/test_differential_analysis_parity.py`, `tests/parity/test_differential_limma_parity.py`, plus unit/integration design and result-contract tests | Batch-aware, block/paired, and repeated-measure designs are rejected in this release. Missing values are rejected at analysis-ready boundary before model fitting. |
| Kinase scoring | `parity-gated` | `KinaseWorkflow` default `scoring_mode="phosr_rank_weighted"` profile/motif scoring and rank-weighted fusion | `tests/parity/test_kinase_workflow_parity.py`, `tests/parity/test_prediction_science_parity.py`, `tests/parity/test_l6_prediction_parity.py` | Relative support scoring only; not calibrated causal inference. Kinase Library scoring is not the default parity lane. |
| Kinase Library motif scoring | `validated PhosPy implementation` | Pure science-layer `KinaseLibraryMotifScorer` / `score_kinase_library_motifs`, plus opt-in `KinaseWorkflow` modes `kinase_library_motif` and `combined_profile_motif` for supplied Kinase Library-style resources | `tests/unit/test_kinase_library_motif_scoring.py`, `tests/science/test_kinase_library_motif_scoring_science.py`, `tests/integration/test_kinase_library_workflow_scoring.py` | Requires explicit compatible resource. Workflow motif scores are normalized to unit interval per kinase matrix for within-run ranking support; raw science-layer motif scores preserve provider scale. Scores are not probabilities. Ser/Thr and Tyr matrix lanes are not interchangeable. |
| Kinase prediction | `parity-gated` | Deterministic and adaptive kinase prediction in `KinaseWorkflow` | `tests/parity/test_public_predmat_parity.py`, `tests/parity/test_l6_prediction_parity.py`, `tests/parity/test_adaptive_prediction_parity.py`, `tests/parity/test_adaptive_replay_parity.py` | Prediction scores are ranking support, not probabilities. |
| Kinase activity scoring | `validated PhosPy implementation` | Supported activity methods: `simplified_weighted_substrate_activity_v1` and `ksea_zscore_activity_v1` | Unit activity tests (`tests/unit/test_activity_science.py`) and parity activity gate (`tests/parity/test_activity_stage_parity.py`) | KSEA-style activity is not a claim of full PhosR kinase activity equivalence. |
| Signalome analysis | `parity-gated` | `SignalomeWorkflow` module assignment, network outputs, and protein-site context | `tests/parity/test_signalome_workflow_parity.py`, `tests/parity/test_signalome_clustering_backend_parity.py` | Derived summaries, not causal proof. Requires explicit signalome protein grouping metadata in `site_metadata.protein_id`. |
| Signalome sampled candidate scoring policy | `experimental` | `SignalomeConfig.sampled_candidate_scoring()` approximates candidate module-count scoring | Parity/contract coverage through signalome parity tests and workflow contract checks | Approximation applies to candidate scoring only; tree generation remains exact-policy governed. |
| Sequence context | `parity-gated` | `site_sequence` required at analysis-ready boundary and used in kinase scoring/prediction | `tests/parity/test_l6_prediction_parity.py`, `tests/parity/test_prediction_science_parity.py` | Sequence quality remains an upstream dependency. |
| Localisation handling | `validated PhosPy implementation` | Localisation confidence validation and fail-fast threshold policies are supported | `tests/unit/test_localisation_policy_preprocessing.py`, `tests/unit/test_validator_boundaries.py` | No full localisation-filter workflow parity claim in this release. |
| Missing values | `parity-gated` | Missing-data policy execution in preprocessing and downstream score preconditioning | `tests/parity/test_preprocessing_science_parity.py`, unit missing-data tests | Policy choice changes retained rows and downstream behavior. |
| Imputation | `validated PhosPy implementation` | Supported policies include `row_median`, `minprob`, `knn` | Unit preprocessing/scientific invariant tests | Policy-dependent behavior; not blanket PhosR-equivalent imputation. |
| Normalisation | `parity-gated` | Supported methods: `none`, `median_center`, `quantile` with stage-order provenance | `tests/parity/test_preprocessing_science_parity.py`, unit preprocessing tests | Method-specific claims only; no blanket normalisation equivalence claim. |
| Batch correction / RUV | `open gap` | No executable SPS/RUV correction lane in current public workflow | N/A for execution; readiness diagnostics documented in workflow contracts | `ruv_readiness` is diagnostic/report-only and must not be interpreted as correction support. |
| Enrichment | `deliberate scope difference` | KSEA-style substrate-set enrichment exists within kinase activity lane | Unit activity tests and scientific policy provenance | Broader pathway/gene-set enrichment lane is not part of current core workflow contract. |
| Visualisation | `deliberate scope difference` | No first-class visualization workflow/API in core PhosPy | N/A | Visualization is intentionally out of current scientific parity scope. |
| Supported bundled organisms and references | `deliberate scope difference` | Bundled runtime references are rat-only for `ReferencePreset.AUTO` in this release | Runtime behavior, reference compatibility tests, manifest approval checks, and workflow docs | Human/mouse are valid organisms but require explicit caller-supplied `ReferenceBundle` unless a future release commits approved redistributable packaged data. |
| Full PhosR package equivalence claim | `not planned` | Not claimed | Guardrail documentation in this matrix and `docs/parity.md` | Any implication of global PhosR parity is out of scope. |

## Release-Gated Scientific Checks

Release-bearing scientific checks are documented and executable with these exact
commands/workflows:

- Local release gate command: `make test-release-gate`
- `make test-release-gate` executes:
  - `pytest tests/unit tests/integration -m "not parity and not performance and not release_gate"`
  - `pytest tests/unit/test_provenance_regressions.py tests/integration/test_kinase_workflow_integration.py::test_kinase_public_predmat_provenance_matches_golden_contract tests/integration/test_signalome_workflow_integration.py::test_signalome_l6_provenance_matches_golden_contract -m "release_gate and (reproducibility or golden)"`
  - `pytest tests/parity -m "parity and not parity_diagnostic" -s`
  - `pytest tests/performance -m "performance or release_gate" -q`
- Publish pipeline release gate workflow:
  - `.github/workflows/publish.yml` job `release-gate` runs `make test-release-gate`
- CI parity workflows:
  - `.github/workflows/ci.yml` job `activity-parity-gate` runs `pytest tests/parity/test_activity_stage_parity.py -m "parity and activity_parity" -s`
  - `.github/workflows/ci.yml` job `parity-tests` runs `pytest tests/parity -m parity -s`

## Interpretation Limits

- Weighted activity output (`simplified_weighted_substrate_activity_v1`) is a
  heuristic summary over predicted substrates above threshold/top-N support.
- KSEA-style activity output (`ksea_zscore_activity_v1`) applies unweighted
  substrate-set enrichment z-scores after evidence thresholding and reports
  p-values (and q-values when enabled).
- KSEA-style activity is not equivalent to full PhosR kinase activity inference.
- Rank-weighted fusion scores combine profile-correlation and motif-frequency
  evidence using rank-derived weights.
- Kinase Library-style science-layer motif scores are raw position-specific
  matrix sums on the caller-supplied score scale. Optional percentiles/ranks are
  empirical summaries against caller-supplied reference distributions only.
- Kinase Library workflow motif scores are normalized per kinase matrix to a
  unit interval for within-run ranking support. They are not calibrated
  probabilities and do not imply activity without an explicit activity method.
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

### `kinase_library_motif_scoring_v1`

- What it does:
  scores phosphosite sequence windows against Kinase Library-style
  position-specific motif matrices.
- Assumptions:
  Ser/Thr and Tyr residue-class lanes are distinct and must not be substituted
  for each other.
- Parameters:
  score scale, residue classes, upstream/downstream window size, sequence-window
  semantics, and whether reference distributions were supplied.
- Output meaning:
  raw provider-scale motif score sums in the science layer; workflow integration
  records the provider scale and exposes normalized unit-interval motif support
  scores with site and kinase diagnostics.
- Output does not mean:
  calibrated probability, causal kinase-substrate proof, or activity inference
  by itself.

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
- [ADR-0025](adr/adr_0025_competitive_phosphoproteomics_workflow_coverage.md)
  records future coverage direction and guardrails; it is not a current
  support claim by itself.
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
