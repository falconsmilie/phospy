# PhosPy

PhosPy exposes a focused public product in `src/phospy/`.

## Installation

Standard installation includes support for both public prediction modes,
including `mode="adaptive_ensemble"`. No extra dependency step is required for
the supported adaptive path.

```bash
pip install .
```

## Public Import Contract

`phospy.api` is the canonical namespace where public API types are defined and
organised in source.

Both namespaces are public, but they have different roles:

- `phospy.api` is the authoritative full public contract namespace.
- top-level `phospy` is a curated convenience surface for the four main product
  entrypoints only:
  `AnalysisReadyDatasetBuilder`, `AnalysisReadyPhosphoDataset`,
  `KinaseWorkflow`, `SignalomeWorkflow`.
- Import requests, configs, results, enums/references, and errors from
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

## Public Product Shape

- `AnalysisReadyPhosphoDataset` is the workflow dataset boundary.
- There is one public dataset builder story:
  `AnalysisReadyDatasetBuilder().run(DatasetBuildRequest(...))`.
- The builder supports both public input routes:
  pandas `DataFrame` values and table file paths (`.csv`, `.tsv`/`.txt`, `.parquet`).
- Builder preprocessing policy is configured through
  `DatasetBuildRequest.preprocessing_config` (`DatasetPreprocessingConfig`).
- Public workflows are:
  `KinaseWorkflow().run(KinaseWorkflowRequest(...))` and
  `SignalomeWorkflow().run(SignalomeWorkflowRequest(...))`.
- Public workflows are one request DTO in, one result DTO out.
- Result models stay nested by stage:
  `result.scoring_result.profile_scores`,
  `result.prediction_result.pred_mat`,
  `result.activity_result.weighted_activity` (when activity is enabled),
  `signalome_result.kinase_result.prediction_result.pred_mat`,
  `signalome_result.module_selection_diagnostics`,
  `signalome_result.expanded_signalome`.

## Boundary Contract

- Builder boundary is flexible about source type (in-memory frames or file paths).
- Final dataset boundary is strict:
  `AnalysisReadyPhosphoDataset` validates DataFrame structure/content, canonical site IDs,
  required metadata (`gene_symbol`, `site`), optional `site_sequence` quality when present,
  and transformation-state coherence.
- Public builder transformation-state establishment is intentionally narrow:
  it establishes the supported pass-through `linear` state only.
- Builder preprocessing can still apply explicit supported policies before
  state establishment (`impute_row_median`, `ratio_to_total`,
  `build_from_metadata`, and grouped comparison construction).
- `site_matrix.policy="build_from_metadata"` is sequence-dependent builder
  preprocessing: it requires usable `site_sequence` rows to construct
  site-level sequence context, and rows without usable sequence are excluded
  from that construction path.
- In this lane, sequence support is established row-by-row from supplied
  `site_metadata.site_sequence` values and/or bundled derivation when available.
  Bundled derivation resolves each row from `gene_symbol` + `site` identity
  (falling back to row index labels when needed).
  Mixed-support inputs retain resolvable rows and exclude only unresolved rows.
- Choosing `build_from_metadata` can therefore reduce row retention relative to
  the original metadata table.
- This is a policy-specific preprocessing requirement, not a contradiction of
  the final dataset boundary where `site_sequence` may still be optional.
- Missing-data preprocessing policy is explicit and grouped:
  `preprocessing_config.missing_data.policy="forbid"` (default) or
  `"impute_row_median"` (requires
  `preprocessing_config.missing_data.min_observed_values`).
- Workflows consume only `AnalysisReadyPhosphoDataset` (not raw input files/frames).

## Supported Science vs Deferred

Scientific confidence is tiered in this project. `implemented`, `supported`,
`parity-gated`, and `closed` are not interchangeable claims. See
[`docs/parity.md`](docs/parity.md) for tier definitions and
[`docs/architecture/legacy_science_gap_audit.md`](docs/architecture/legacy_science_gap_audit.md)
for area-level evidence.

Supported in the current public lane:

- Kinase scoring with authoritative downstream outputs:
  `profile_scores`, `combined_scores`.
- Optional diagnostic scoring tables:
  `motif_scores`, `weights` via
  `KinaseScoringConfig(include_diagnostic_scoring_tables=True)`.
- Profile-driven prediction ranking and matrix assembly (`prediction_result.pred_mat`).
- Kinase prediction config supports `mode` (`deterministic_ranking` or
  `adaptive_ensemble`), `adaptive_policy` (`stable` or `r_parity`), and
  `n_iterations` (adaptive resampling iterations).
- Kinase scoring config supports
  `profile_missing_value_strategy` (`strict` or `median_skipna`).
- Optional kinase activity stage inside `KinaseWorkflow` (`activity_config=None` or
  `enabled=False` disables it).
- Signalome workflow outputs:
  `module_assignments`, `signalome_modules`, `kinase_network`,
  `module_selection_diagnostics`, and `expanded_signalome`.
- Signalome config supports module-selection controls:
  `module_count`, `module_selection_primary_correlation_threshold`,
  `module_selection_fallback_correlation_threshold`,
  `module_selection_max_clusters`.
- Builder preprocessing supports explicit:
  `total_protein_correction.policy="ratio_to_total"`,
  `site_matrix.policy="build_from_metadata"`, and
  `comparisons.policy="sample_metadata_pairs"` lanes.

Contract-changed supported lanes include:

- adaptive public naming with `adaptive_policy` rather than legacy `svm_mode`
- signalome input contracted to `SignalomeWorkflowRequest(kinase_result=...)`
- signalome protein identity from explicit `site_metadata.protein_id` (no legacy
  site-id-prefix fallback)
- motif sequence authority from `references.site_sequences` in the resolved
  reference bundle

Deferred or out of the supported default lane:
- legacy or experimental science lanes not yet ported into the public path
- legacy-science surfaces not listed in the audited inventory documents above

## Current Limits

- Supported site-metadata aliases in builder normalization are intentionally narrow:
  - `gene_symbol`: `gene_symbol`, `gene_name`
  - `site`: `site`
  - `site_sequence`: `site_sequence`, `centralized_sequence`
  - `protein_id`: `protein_id`
- Unsupported legacy aliases (`gene`, `residue`, `phosphosite`, `site_position`,
  `sequence`, `protein`) are rejected with actionable errors.
- `ReferencePreset.AUTO` requires `dataset.organism`.
- Bundled runtime references are currently rat-only.
- `ReferencePreset.HUMAN` and `ReferencePreset.MOUSE` remain valid public enum values,
  but require caller-supplied `ReferenceBundle` in this release.
- Kinase scoring enforces `scoring_config.min_substrates >= 2`.

## Example

```python
from phospy import KinaseWorkflow
from phospy.api import KinaseWorkflowRequest

result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
    )
)

profile_scores = result.scoring_result.profile_scores
pred_mat = result.prediction_result.pred_mat
if result.activity_result is not None:
    weighted_activity = result.activity_result.weighted_activity
```

Signalome outputs now include module-selection diagnostics and expanded signalome rows:

```python
from phospy import SignalomeWorkflow
from phospy.api import SignalomeWorkflowRequest

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(kinase_result=result)
)

diagnostics = signalome_result.module_selection_diagnostics
expanded = signalome_result.expanded_signalome  # populated in supported executor lane
```

## Package Boundary

```text
src/phospy/                  # supported package
legacy_archive/phospy_legacy # migration reference only (not installed package content)
```

## Examples

- [`examples/dataset_builder_demo.py`](examples/dataset_builder_demo.py)
- [`examples/kinase_workflow_demo.py`](examples/kinase_workflow_demo.py)
- [`examples/signalome_workflow_demo.py`](examples/signalome_workflow_demo.py)

`examples/dataset_builder_demo.py` includes a `site_matrix.policy="build_from_metadata"`
row-retention walkthrough, including unusable-sequence handling and
retained-row diagnostics.

## Docs

- [`docs/api.md`](docs/api.md)
- [`docs/cli.md`](docs/cli.md)
- [`docs/output_bundles.md`](docs/output_bundles.md)
- [`docs/validation.md`](docs/validation.md)
- [`docs/roadmap.md`](docs/roadmap.md)
