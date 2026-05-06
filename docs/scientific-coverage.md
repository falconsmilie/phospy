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
