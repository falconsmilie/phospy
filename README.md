# PhosPy

PhosPy exposes a focused rewrite public product in `src/phospy/`.

## Public Product Shape

- `AnalysisReadyPhosphoDataset` is the workflow dataset boundary.
- There is one public dataset builder story:
  `AnalysisReadyDatasetBuilder().run(DatasetBuildRequest(...))`.
- The builder supports both public input routes:
  pandas `DataFrame` values and table file paths (`.csv`, `.tsv`/`.txt`, `.parquet`).
- Public workflows are:
  `KinaseWorkflow().run(KinaseWorkflowRequest(...))` and
  `SignalomeWorkflow().run(SignalomeWorkflowRequest(...))`.
- Public workflows are one request DTO in, one result DTO out.
- Result models stay nested by stage:
  `result.scoring_result.profile_scores`,
  `result.prediction_result.pred_mat`,
  `result.activity_result.weighted_activity` (when activity is enabled),
  `signalome_result.kinase_result.prediction_result.pred_mat`.

## Boundary Contract

- Builder boundary is flexible about source type (in-memory frames or file paths).
- Final dataset boundary is strict:
  `AnalysisReadyPhosphoDataset` validates DataFrame structure/content, canonical site IDs,
  required metadata (`gene_symbol`, `site`), optional `site_sequence` quality when present,
  and transformation-state coherence.
- Workflows consume only `AnalysisReadyPhosphoDataset` (not raw input files/frames).

## Supported Science vs Deferred

Supported in the current rewrite lane:

- Kinase scoring with nested outputs:
  `profile_scores`, `motif_scores`, `combined_scores`, `weights`.
- Supported kinase motif scoring consumes `references.site_sequences` from the resolved
  reference bundle.
- Profile-driven prediction ranking and matrix assembly (`prediction_result.pred_mat`).
- Optional kinase activity stage inside `KinaseWorkflow` (`activity_config=None` or
  `enabled=False` disables it).
- Signalome workflow outputs:
  `module_assignments`, `signalome_modules`, `kinase_network`.

Deferred or not in the supported default lane:

- `SignalomeWorkflowResult.expanded_signalome` population (`expanded_signalome` is optional and currently `None` in the default route).
- Legacy or experimental science lanes not yet ported into the public rewrite path.

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
from phospy import KinaseWorkflow, KinaseWorkflowRequest

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

## Package Boundary

```text
src/phospy/                  # supported rewrite package
legacy_archive/phospy_legacy # migration reference only (not installed package content)
```

## Examples

- [`examples/dataset_builder_demo.py`](examples/dataset_builder_demo.py)
- [`examples/simple_workflow_demo.py`](examples/simple_workflow_demo.py)
- [`examples/signalome_workflow_demo.py`](examples/signalome_workflow_demo.py)

## Docs

- [`docs/api.md`](docs/api.md)
- [`docs/cli.md`](docs/cli.md)
- [`docs/output_bundles.md`](docs/output_bundles.md)
- [`docs/validation.md`](docs/validation.md)
- [`docs/roadmap.md`](docs/roadmap.md)
