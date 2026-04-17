# Validation Guide

This guide covers validation for the currently supported rewrite contract only.

For public types and signatures, see [`api.md`](api.md).

## Boundary

- Supported inputs are in-memory pandas `DataFrame` objects through `DatasetBuildRequest`.
- Supported workflow route is `SimpleKinaseWorkflow.run(SimpleKinaseWorkflowRequest(...))`.
- `SignalomeWorkflow.run(SignalomeWorkflowRequest(...))` is a placeholder shell.
- File-ingestion and legacy convenience routes are outside the current rewrite contract.

## Dataset Builder Validation

`DatasetBuildRequest` and `AnalysisReadyPhosphoDataset` enforce:

- `phospho` and `site_metadata` must be pandas `DataFrame` values
- `phospho` must be numeric, non-empty, with unique index and unique columns
- `site_metadata` must be non-empty and index-aligned to `phospho`
- `site_metadata` must include:
  `gene_symbol`, `site`, `site_sequence`
- `site_sequence` values must be non-empty strings
- `sample_metadata` (if present) must be index-aligned to `phospho.columns`
- `total` (if present) must be numeric and share columns with `phospho`
- `organism` (if present) must be an `Organism` enum

## Reference Validation

`ReferencePreset` and `ReferenceBundle` resolution enforces:

- `ReferencePreset.AUTO` requires `dataset.organism`
- explicit preset organism must match `dataset.organism` when both are set
- bundled references are currently packaged for rat only
- `ReferenceBundle.kinase_substrate_map` must be non-empty with:
  `kinase`, `substrate_site`
- `ReferenceBundle.site_sequences` must be non-empty with:
  `site_sequence` and unique index
- each `substrate_site` must exist in `site_sequences.index`

## Workflow Validation

`SimpleKinaseWorkflowRequest` enforces:

- `dataset` is `AnalysisReadyPhosphoDataset`
- `references` is `ReferencePreset` or `ReferenceBundle`
- `scoring_config` is `KinaseScoringConfig`
- `prediction_config` is `KinasePredictionConfig`
- `activity_config` is `KinaseActivityConfig` or `None`

Rewrite-era boundary diagnostics (raised as `WorkflowBoundaryError`) also enforce:

- interpreted reference coverage must overlap dataset phosphosites
- at least one kinase must meet `scoring_config.min_substrates` after overlap
- prediction ranking must produce at least one ensemble kinase
- activity (when enabled) must have at least one positive predicted site assignment

Boundary error messages include:

- the failing seam (for example `simple_kinase.interpreter.eligible_kinases`)
- concrete counts (`dataset_sites`, `overlap_sites`, `eligible_kinases`, etc.)
- active config values (`scoring_config_min_substrates`,
  `prediction_config_ensemble_size`, `activity_config_threshold`)
- a `next_action` hint for likely recovery

Stage result access is nested and stable:

- `result.scoring_result.profile_scores`
- `result.prediction_result.pred_mat`
- `result.activity_result.activity_scores` (when enabled)

## Signalome Placeholder Validation

`SignalomeWorkflowRequest` currently validates:

- `kinase_result` is `SimpleKinaseWorkflowResult`
- `config` is `SignalomeConfig`

Current signalome execution returns placeholder (empty) output tables while the
workflow implementation is pending.

## Quick Troubleshooting

| Problem | Usually means | Good next step |
| --- | --- | --- |
| Builder rejects input format | A request field is not a pandas `DataFrame` | Pass in-memory `DataFrame` objects for `phospho` and `site_metadata` |
| `ReferencePreset.AUTO` fails | `dataset.organism` is missing | Set `organism` on `DatasetBuildRequest` |
| Reference mismatch error | Dataset and selected preset organisms conflict | Align `dataset.organism` with `ReferencePreset` |
| Workflow request type error | Request field types are not the public models | Build the request from top-level `phospy` models |
| `simple_kinase.interpreter.reference_coverage` | None of the reference substrate sites overlap `dataset.phospho.index` | Use references for the same identifier scheme/organism and verify site IDs |
| `simple_kinase.interpreter.eligible_kinases` | Overlap exists, but no kinase reaches `scoring_config.min_substrates` | Lower `min_substrates` or use references with deeper site overlap |
| `simple_kinase.executor.prediction_ensemble` | Scoring completed, but no kinase had a finite prediction ranking | Use at least two informative samples and relax strict scoring thresholds |
| `simple_kinase.executor.activity_support` | Activity was enabled, but predictions had no positive site assignments | Increase `top_k` and/or lower `activity_config.threshold` for sparse data |
| Signalome result is empty | Signalome route is still a placeholder | Use simple kinase outputs as the supported scientific route for now |

Runnable rewrite examples:

- [`../examples/dataset_builder_demo.py`](../examples/dataset_builder_demo.py)
- [`../examples/simple_workflow_demo.py`](../examples/simple_workflow_demo.py)
- [`../examples/signalome_placeholder_demo.py`](../examples/signalome_placeholder_demo.py)
