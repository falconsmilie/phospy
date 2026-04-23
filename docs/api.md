# API Guide

This guide describes the current supported public contract only.

> Audience: advanced users and maintainers who need exact public contract detail.
> If you are new to PhosPy, start with
> [What is PhosPy?](getting-started/what-is-phospy.md) and
> [Quickstart](getting-started/quickstart-first-workflow.md) first.

## Supported Lanes

PhosPy has one public dataset boundary and two public workflow stories:

- Dataset construction:
  `DatasetBuildRequest -> AnalysisReadyDatasetBuilder.run(request) -> AnalysisReadyPhosphoDataset`
- Kinase workflow:
  `KinaseWorkflow.run(KinaseWorkflowRequest(...)) -> KinaseWorkflowResult`
- Signalome workflow:
  `SignalomeWorkflow.run(SignalomeWorkflowRequest(...)) -> SignalomeWorkflowResult`

All public executors use `run(request)`.

## Package Boundary

- `src/phospy/`: supported package

## Installation Contract

Release-facing install:

```bash
pip install phospy
```

Contributor install from a local clone:

```bash
pip install -e ".[dev]"
```

Standard installation includes all dependencies required by supported prediction
lanes, including `mode="adaptive_ensemble"` (scikit-learn is part of base
dependencies). No extra install step is required for adaptive mode.

## Recommended First-Run Lane

The supported first-run path is:

1. build `AnalysisReadyPhosphoDataset` with `dataset.organism=Organism.RAT`
2. run `KinaseWorkflowRequest(..., references=ReferencePreset.AUTO)`
3. optionally run `SignalomeWorkflowRequest(kinase_result=...)`

Reference expectations for this lane:

- `ReferencePreset.AUTO` requires `dataset.organism`
- bundled runtime references are rat-only in this release
- human/mouse execution requires an explicit caller-supplied `ReferenceBundle`

## Import Contract

`phospy.api` is the canonical namespace where public API types are defined and
organised in source.

This split is intentional:

- `phospy.api` is the authoritative full public contract namespace.
- top-level `phospy` is a curated convenience surface for only:
  `AnalysisReadyDatasetBuilder`, `AnalysisReadyPhosphoDataset`,
  `KinaseWorkflow`, `SignalomeWorkflow`.
- requests, configs, results, support types, and exceptions are imported from
  `phospy.api`.

Simple product-entrypoint usage:

```python
from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
```

Full contract usage:

```python
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    KinaseWorkflowRequest,
    PhosPyValidationError,
)
```

## Public Types

Use `phospy.api` for the supported contract surface.

- Main product entrypoints (convenience imports available from top-level
  `phospy` and from `phospy.api`):
  `AnalysisReadyDatasetBuilder`, `AnalysisReadyPhosphoDataset`,
  `KinaseWorkflow`, `SignalomeWorkflow`
- Request/config/result/reference/enum types:
  import from `phospy.api`
- Error taxonomy:
  import from `phospy.api`

## Builder Contract

There is one public builder story: `AnalysisReadyDatasetBuilder.run(DatasetBuildRequest(...))`.

`DatasetBuildRequest` supports both public input routes:

- pandas `DataFrame` values
- file paths (`str`, `pathlib.Path`, `os.PathLike`) for supported table formats

Required request fields:

- `phospho`
- `site_metadata`

Optional request fields:

- `sample_metadata`
- `total`
- `organism`
- `preprocessing_config` (`DatasetPreprocessingConfig`)

After loading, both routes share the same normalization and validation path.

Supported site-metadata alias mapping is explicit and narrow:

- `gene_symbol`: `gene_symbol`, `gene_name`
- `site`: `site`
- `site_sequence`: `site_sequence`, `centralized_sequence`
- `protein_id`: `protein_id`

Unsupported legacy aliases (`gene`, `residue`, `phosphosite`, `site_position`,
`sequence`, `protein`) are rejected instead of guessed.

If `gene_symbol` and/or `site` are absent, one derivation convention is supported:
`site_metadata.index` values formatted exactly as `"<gene_symbol>;<site>;"`.
This derivation populates only `gene_symbol` and `site`; it does not infer
`protein_id`.

Current transformation establishment policy in this public builder lane is
intentionally narrow:

- builder execution establishes only the supported pass-through `linear` state
- quantitative matrix values are preserved as provided after builder normalization,
  except when an explicit preprocessing policy modifies them
- no additional log or heuristic transformation path is exposed in the public builder

`DatasetPreprocessingConfig` is grouped and builder-owned:

- `missing_data` (`DatasetMissingDataConfig`)
- `total_protein_correction` (`DatasetTotalProteinCorrectionConfig`)
- `site_matrix` (`DatasetSiteMatrixConfig`)
- `comparisons` (`DatasetComparisonBuildingConfig`)

Current supported policies:

- `missing_data.policy="forbid"` (default): no missing-value preprocessing.
- `missing_data.policy="impute_row_median"`: drop rows below
  `missing_data.min_observed_values`, then row-median imputation.
- `total_protein_correction.policy="none"` (default): no total/protein correction.
- `total_protein_correction.policy="ratio_to_total"`: subtract matched
  `total` abundance from `phospho` abundance per sample column in builder
  preprocessing. This policy requires `total` input, exact
  `total.columns == phospho.columns`, and complete
  `site_metadata.gene_symbol` vs `total.index` matching.
- `site_matrix.policy="as_input"` (default): preserve interpreted phospho/site rows.
- `site_matrix.policy="build_from_metadata"`: construct site-matrix-ready rows
  from `site_metadata.gene_symbol`/`site_metadata.site` after upstream missing-data and
  total/protein-correction stages.
  This path is sequence-dependent at preprocessing time: usable sequence support
  is established row-wise from supplied `site_metadata.site_sequence` values
  and/or bundled derivation when available. Rows lacking usable sequence support
  are excluded from this path, so retained rows can be narrower than the
  original metadata table.
  Bundled derivation resolves per-row site identity from `gene_symbol` + `site`
  (with row-index fallback when needed), so mixed-support inputs can retain
  derivable/supplied rows together.
  Mixed-support inputs therefore keep resolvable rows instead of collapsing the
  entire derived-sequence path.
  This requirement is specific to the selected preprocessing policy and does
  not change the final dataset boundary where `site_sequence` remains optional.
  Additional supported policy controls under `DatasetSiteMatrixConfig`:
  - `missing_data_policy="drop_any_missing"` (default): keep only complete rows.
    This is the supported public complete-case lane for
    `AnalysisReadyPhosphoDataset` construction.
  - retained-missingness site-matrix modes
    (`retain_missing`, `require_min_observed_values`) are internal-only
    compatibility behavior and are rejected at public request validation.
  - `minimum_observed_values` is internal-only compatibility state and must
    remain unset in the supported public builder lane.
  - `duplicate_site_strategy="max_mean_signal"` (default): keep strongest row.
  - `duplicate_site_strategy="first"`: keep first duplicate row.
  - `duplicate_site_strategy="aggregate_mean"`: aggregate duplicate rows by mean.
  - `duplicate_site_strategy="aggregate_median"`: aggregate duplicate rows by median.
  - `duplicate_site_strategy="error"`: fail on duplicate constructed site IDs.
  Site IDs are constructed deterministically as canonical
  `GENE_SYMBOL;SITE_TOKEN;` identifiers.
  When `site_matrix.policy="as_input"`, these execution-only fields must remain
  at defaults and are rejected if overridden.
- `comparisons.policy="none"` (default): do not construct comparison columns.
- `comparisons.policy="sample_metadata_pairs"`: build dataset-level pairwise
  comparison columns in `dataset.comparisons` from grouped sample metadata.
  Required and supported inputs:
  - `sample_metadata` must be provided
  - `sample_metadata.index` must align to `phospho.columns`
  - `sample_metadata[comparisons.sample_group_column]` must contain one
    non-empty group label per sample
  - `comparisons.pairs` can pass explicit `(left_group, right_group)` pairs;
    when omitted, all unique pairwise combinations are inferred from observed
    groups.

## Final Dataset Boundary

`AnalysisReadyPhosphoDataset` is strict and workflow-facing.

- It owns validated tables, not input files.
- It requires DataFrame values for `phospho` and `site_metadata` at construction time.
- `site_metadata` must contain `gene_symbol`, `site` with non-empty strings.
- Site identity coherence is strict: each `phospho.index` row ID must be
  interpretable as `"<gene_symbol>;<site>;"`, and parsed values must exactly
  match the corresponding `site_metadata.gene_symbol` and `site_metadata.site`
  row values.
- `site_sequence` is optional at this boundary; when present it must be non-empty.
- Final-boundary optionality does not remove policy-specific preprocessing
  requirements: `site_matrix.policy="build_from_metadata"` still requires
  usable sequence-bearing rows for that construction path.
- Site identifiers must already be canonical and non-colliding.
- `sample_metadata` (if present) must align to `phospho.columns`.
- `total` (if present) must be numeric and column-aligned to `phospho`.
- `comparisons` (if present) must be numeric, non-missing, and aligned to
  `phospho.index`.
- `transformation_state` is mandatory and must be both coherent and
  established through a supported PhosPy path.

Supported establishment paths are:

- `AnalysisReadyDatasetBuilder.run(...)` (default supported lane)
- supported transformer execution through the dataset builder executor/resolver
- supported bundle reconstruction paths

For `AnalysisReadyDatasetBuilder.run(...)`, the established state is currently
the pass-through `linear` path.
When `total_protein_correction.policy="ratio_to_total"` is used, corrected
quantitative values are represented directly in `dataset.phospho`.

Directly declared transformation-state objects are not accepted as authoritative
at the dataset boundary unless they were established through one of the supported
PhosPy paths above.
Direct caller-side minting is explicitly rejected: `TransformationState.established_raw(...)`
and `establish_transformation_state(...)` are accepted only when an approved
internal establishment authority is provided by the supported
builder/transformer or bundle-reconstruction lanes.

Builder flexibility does not weaken this final dataset strictness.
Workflows consume only `AnalysisReadyPhosphoDataset`.

## Workflow Contract

`KinaseWorkflowRequest` fields:

- `dataset: AnalysisReadyPhosphoDataset`
- `references: ReferencePreset | ReferenceBundle`
- `scoring_config: KinaseScoringConfig`
- `prediction_config: KinasePredictionConfig`
- `activity_config: KinaseActivityConfig | None`

`SignalomeWorkflowRequest` fields:

- `kinase_result: KinaseWorkflowResult`
- `config: SignalomeConfig`
- supported prerequisite: `kinase_result.dataset.site_metadata.protein_id` must be
  present and non-empty for all interpreted sites

Signalome protein-identity contract (supported lane):

- signalome requires explicit `dataset.site_metadata.protein_id` values
- gene-symbol site IDs (for example `"<gene_symbol>;<site>;"`) are site identity,
  not protein identity, and are not used as fallback protein mapping
- this is an intentional scientific boundary: downstream signalome grouping and
  module assignment are protein-identity-aware
- builder flexibility at ingestion does not weaken this downstream workflow
  contract

`KinaseScoringConfig` fields:

- `min_substrates` (validated floor: `>= 2`)
- `include_diagnostic_scoring_tables`
- `profile_missing_value_strategy` (`"strict"` or `"median_skipna"`)

`KinasePredictionConfig` fields:

- `top_k` (validated floor: `>= 1`)
- `ensemble_size` (validated floor: `>= 1`)
- `mode` (`"deterministic_ranking"` or `"adaptive_ensemble"`)
- `adaptive_policy` (`"stable"` or `"r_parity"`)
- `n_iterations` (validated floor: `>= 1`; used by adaptive lane)
- `random_state` (`None` or integer `>= 0`)

## Kinase Contract Difference Inventory (2026-04-22)

This table is the supported contract truth source for kinase
scoring/prediction behavior where current and historical baselines differ.

| Behavior seam | Supported contract | Historical baseline | Classification |
| --- | --- | --- | --- |
| Profile-only fallback in profile+motif combine | Workflow scoring always calls profile+motif combine with profile-only fallback enabled (`allow_profile_only_fallback=True`) | Legacy public config default was `allow_profile_only_fallback=False` | Intentional and supported |
| Missing motif values when profile is present | Combined scores restore profile values for `(motif is missing) AND (profile is present)` cells instead of leaving missing combined values | Legacy combine path did not apply this rescue mask | Intentional and supported |
| Candidate filter defaults in workflow prediction lanes | Deterministic/adaptive workflow lanes use fixed candidate filters `score_threshold=0.0` and `inclusion=1`; only `top_k` is caller-facing | Legacy defaults were `score_threshold=0.8`, `inclusion=20`, `top=50` | Intentional and supported |
| Prediction/scoring public knobs | Historical knobs (`min_motif_size`, `allow_profile_only_fallback`, `score_threshold`, `inclusion`, `svm_mode`, `profile_policy`) are not public request fields in the supported lane | Those knobs were exposed in historical `PredictionRunConfig` | Intentional and supported |
| Public prediction mode surface | Public `mode` is explicit (`deterministic_ranking` default, `adaptive_ensemble` optional), with adaptive behavior further selected by `adaptive_policy` | Historical public lane was `svm_mode`-centric and did not expose deterministic ranking as the default contract lane | Intentional and supported |

Inventory status in this audit:

- Temporary parity differences to remove: none in runtime behavior; parity/reporting wording drift is tracked in `docs/parity.md`.
- Unresolved design decisions: none currently tracked for kinase scoring/prediction runtime behavior.

`"adaptive_ensemble"` is part of the normal supported contract and is expected
to work after a standard package install.

Kinase scoring authority in the supported lane:

- Scoring is a function of:
  - `request.dataset` (analysis-ready phospho matrix)
  - resolved `ReferenceBundle` content (`kinase_substrate_map`, `site_sequences`)
  - `request.scoring_config`
- Scoring is not changed by:
  - `prediction_config.mode` (`deterministic_ranking` vs `adaptive_ensemble`)
  - whether references were supplied via `ReferencePreset` or explicit
    `ReferenceBundle` when both resolve to equivalent reference content
- Adaptive prediction is downstream of scoring and consumes the same authoritative
  scoring outputs as deterministic prediction.
- Adaptive-policy settings (`adaptive_policy`) are prediction-stage controls and
  do not opt into a separate scoring contract.

`SignalomeConfig` fields:

- `substrate_support_cutoff`
- `network_correlation_threshold`
- `network_policy` (`"positive_only"`, `"absolute_threshold"`, or `"signed"`)
- `assignment_policy` (`"cutoff_binary"` or `"weighted_top"`)
- `score_preconditioning_policy` (`"allow_and_report"` default, or `"error_on_drop"`)
- `module_count` (`None` for automatic module-count selection)
- `module_selection_primary_correlation_threshold`
- `module_selection_fallback_correlation_threshold`
- `module_selection_max_clusters`

## Result Contract (Nested Stage Outputs)

`KinaseWorkflowResult`:

- `result.dataset`
- `result.references`
- `result.scoring_result.profile_scores`
- `result.scoring_result.combined_scores`
- `result.scoring_result.motif_scores` (optional diagnostic field)
- `result.scoring_result.weights` (optional diagnostic field)
- `result.prediction_result.pred_mat`
- `result.prediction_result.substrate_list` (optional)
- `result.activity_result` (`None` when activity is disabled)

`SignalomeWorkflowResult`:

- `result.dataset`
- `result.kinase_result` (full upstream nested lineage)
- `result.module_assignments.table`
- `result.signalome_modules.table`
- `result.kinase_network.edges`
- `result.kinase_network.nodes` (optional)
- `result.module_selection_diagnostics`
- `result.score_preconditioning_diagnostics`
- `result.expanded_signalome` (official supported output; optional by type for
  compatibility, populated in the supported executor lane)

`module_selection_diagnostics` fields:

- `strategy` (`"correlation_thresholds"` or `"explicit_module_count"`)
- `selected_module_count`
- `requested_module_count`
- `threshold_used`
- `max_clusters_evaluated`
- `candidate_scores` (per-candidate `min_median_correlation` and `mean_median_correlation`)
- `reason`
- `zero_variance_profile_count`
- `near_constant_profile_count`
- `excluded_from_correlation_count`

`score_preconditioning_diagnostics` fields:

- `policy` (`"allow_and_report"` or `"error_on_drop"` from
  `SignalomeConfig.score_preconditioning_policy`)
- `input_row_count` (aligned downstream-score rows before preconditioning)
- `dropped_all_missing_row_count` (rows removed because all kinase scores are missing)
- `retained_row_count` (rows retained for signalome score-driven stages)

`module_assignments.table` includes site-level assignment metadata:

- structural fields: `protein_id`, `module_id`, `top_kinase`, `top_score`
- tie/ambiguity fields: `top_kinase_candidates`, `top_kinase_weights`,
  `top_kinase_tie_count`, `top_kinase_is_ambiguous`,
  `top_kinase_selection_policy`
- module-level attribution fields: `module_top_kinase`,
  `module_top_kinase_candidates`, `module_top_kinase_tie_count`,
  `module_top_kinase_is_ambiguous`, `module_top_kinase_selection_policy`

In the supported signalome lane, `protein_id` is always taken from validated
`dataset.site_metadata.protein_id` (required and non-empty for interpreted
sites). Gene-symbol-prefixed site IDs are not used as protein-identity
fallbacks. This strictness is intentional in the supported scientific lane.

`assignment_policy` behavior:

- `"cutoff_binary"`: support is binary from `prediction_result.pred_mat >
  substrate_support_cutoff`
- `"weighted_top"`: support is fractional from per-site `top_kinase_weights`
  tie metadata propagated through module/expanded outputs

`expanded_signalome` flattened schema:

- `kinase`: focal kinase
- `row_kind`: `"site"` or `"summary"`
- `assignment_policy`: emitted assignment policy
- `linked_kinases`: JSON array string of focal + linked kinases
- `regulated_module_ids`: JSON array string of module IDs where focal-kinase
  share is strictly greater than `1.0` percent
- `site_id`: selected phosphosite ID (`""` on summary rows)
- `site_order`: zero-based position in the original module-assignment order
- `protein_id`, `module_id`, `top_kinase`, `top_score`: selected site metadata
- `support_kinases`: JSON array string of linked kinases supporting the site
- `support_weight`: site support weight under the selected assignment policy

`expanded_signalome` row-selection conditions:

- The executor emits one focal-kinase block per kinase in `signalome_modules.columns`.
- A `row_kind="site"` row is emitted when:
  - the site belongs to a regulated module for that focal kinase
    (`signalome_modules[module_id, kinase] > 1.0`), and
  - at least one linked kinase provides positive site support under the active
    `assignment_policy`.
- If no site rows qualify for a focal kinase, one `row_kind="summary"` row is
  emitted for that kinase to preserve linked-kinase and regulated-module metadata.

No top-level convenience mirrors flatten nested stage outputs.

## Supported Science vs Deferred Science

Supported public lane (stable and recommended):

- Kinase scoring stage always outputs `profile_scores` and `combined_scores`.
- Diagnostic scoring tables (`motif_scores`, `weights`) are opt-in via
  `KinaseScoringConfig(include_diagnostic_scoring_tables=True)`.
- Motif scoring in the supported kinase lane uses `references.site_sequences`
  from the resolved `ReferenceBundle`.
- Prediction stage uses a downstream score matrix that resolves to
  `combined_scores` when present and falls back to `profile_scores` only when
  combined scores are unavailable.
- Activity stage is supported and optional inside `KinaseWorkflow`.
- Signalome stage consumes the same downstream score-matrix lane as prediction
  and outputs module assignments, module matrix, kinase network, and
  `expanded_signalome`.
- Signalome assignment policy is explicit:
  `assignment_policy="cutoff_binary"` keeps cutoff/binary support semantics;
  `assignment_policy="weighted_top"` propagates fractional
  `top_kinase_weights` support into module shares.
- Signalome network policy is explicit:
  `network_policy="positive_only"` keeps only positive correlations above
  threshold; `network_policy="absolute_threshold"` keeps absolute correlations
  above threshold and emits unsigned edge correlations; `network_policy="signed"`
  keeps absolute correlations above threshold and emits signed edge
  correlations.
- Downstream score missingness is part of the supported scientific contract:
  all-missing score rows are explicitly policy-governed preconditioning
  candidates, partially missing rows are retained, and infinite values remain
  invalid.
- Preconditioning drop policy is caller-owned in the supported lane:
  `score_preconditioning_policy="allow_and_report"` drops all-missing rows and
  reports counts via `result.score_preconditioning_diagnostics`;
  `score_preconditioning_policy="error_on_drop"` fails interpretation when any
  all-missing row would be dropped.
- Active parity regression for this lane is execution against committed
  rewrite-owned fixtures; no archived runtime-module execution is part of the
  active parity gate.

Deferred/experimental/not yet ported into the public lane:
- Additional legacy science lanes listed as roadmap follow-ons.

## Reference Resolution

- `ReferencePreset.AUTO` requires `dataset.organism`.
- Preset/dataset organism compatibility is enforced.
- Bundled runtime references are currently rat-only.
- Enum breadth is broader than bundled runtime support by design:
  `ReferencePreset.HUMAN` and `ReferencePreset.MOUSE` remain public enum lanes,
  but bundled resolution for those presets is intentionally unsupported in this release.
- Non-rat execution uses explicit caller-provided `ReferenceBundle`.
- Public enum availability is syntactic contract surface, not by itself a
  scientific-support claim for bundled runtime execution.
- Any future broadening of bundled organism support must ship both:
  - bundled reference data for the new organism lane, and
  - parity-backed evidence for that lane in the active regression contract
    before docs can claim it as supported bundled runtime behavior.

## User-Handleable Exceptions

Import public exceptions from `phospy.api`.
`phospy.errors` remains available as the underlying taxonomy module.

## Quick Usage Pattern

```python
import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.api import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
    SignalomeWorkflowRequest,
)

phospho = pd.DataFrame(
    {
        "sample_a": [1.00, 0.70],
        "sample_b": [1.10, 0.80],
        "sample_c": [0.95, 0.75],
    },
    index=["TSC2;S939;", "GSK3B;S9;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["TSC2", "GSK3B"],
        "site": ["S939", "S9"],
        "site_sequence": [
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
            "_______MSGRPRTTSFAESCKPVQQPSAFG",
        ],
        "protein_id": ["TSC2", "GSK3B"],
    },
    index=phospho.index.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
    )
)

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(dataset=dataset, references=ReferencePreset.AUTO)
)

pred_mat = kinase_result.prediction_result.pred_mat
signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(kinase_result=kinase_result)
)
```

`SignalomeWorkflow` requires non-empty `dataset.site_metadata.protein_id` values
for interpreted sites.

If you choose `site_matrix.policy="build_from_metadata"`, inspect row-retention
counts after builder execution (`dataset.phospho.shape[0]` versus input row
count). Sequence support is resolved row-by-row from supplied and/or derivable
`site_sequence` values; unsupported rows do not participate in that
construction path and are excluded there. See
[`examples/dataset_builder_demo.py`](../examples/dataset_builder_demo.py) for a
concrete retained-vs-excluded example.

For CLI and bundle persistence details, see:

- [`cli.md`](cli.md)
- [`output_bundles.md`](output_bundles.md)

## Where Next

- Validation guarantees: [Validation Guide](validation.md)
- Workflow-oriented navigation: [Workflow guides](workflow-guides/index.md)
- Scientific confidence and parity policy: [Scientific parity and governance](science/index.md)
