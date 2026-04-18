# Validation Guide

This guide covers validation for the currently supported rewrite contract only.

For public types and signatures, see [`api.md`](api.md).

## Boundary

- Supported dataset-build inputs are pandas `DataFrame` objects or file paths through
  `DatasetBuildRequest`.
- Supported workflow route is `KinaseWorkflow.run(KinaseWorkflowRequest(...))`.
- Supported downstream route is
  `SignalomeWorkflow.run(SignalomeWorkflowRequest(...))`.
- File-ingestion is part of the dataset builder route only.

## Constructor Validation Policy

Public request/result/data containers validate type contracts at construction time
and raise PhosPy exception types (not raw built-in exceptions) for intended
caller-facing validation failures.

- Request DTOs (`DatasetBuildRequest`, `KinaseWorkflowRequest`,
  `SignalomeWorkflowRequest`) validate core type shape immediately.
- Result/data containers validate nested public types and pandas container fields
  at construction time.
- Deeper semantic policies (numeric ranges, overlap constraints, scientific floors)
  still run in validator/interpreter boundaries and raise domain-specific PhosPy
  exceptions.

## Dataset Builder Validation

`DatasetBuildRequest` and `AnalysisReadyPhosphoDataset` enforce:

- `phospho` and `site_metadata` must be pandas `DataFrame` values or file paths
- `phospho` must be numeric, non-empty, with unique index and unique columns
- `site_metadata` must be non-empty and index-aligned to `phospho`
- `site_metadata` must include:
  `gene_symbol`, `site`, `site_sequence`
- `site_metadata.protein_id` is optional but preserved when supplied
- `site_sequence` values must be non-empty strings
- `sample_metadata` (if present) must be index-aligned to `phospho.columns`
- `total` (if present) must be numeric and share columns with `phospho`
- `organism` (if present) must be an `Organism` enum
- `transformation_state` is required on direct dataset construction
- boundary site identifiers must already be canonical (non-empty stripped strings)
- boundary constructors do not canonicalize or repair incoming tables
- dataset builder collaborators are responsible for canonicalization and shaping

## Reference Validation

`ReferencePreset` and `ReferenceBundle` resolution enforces:

- `ReferencePreset.AUTO` requires `dataset.organism`
- explicit preset organism must match `dataset.organism` when both are set
- bundled references are currently packaged for rat only
- `ReferencePreset.HUMAN` and `ReferencePreset.MOUSE` remain valid public enum
  values but intentionally fail bundled resolution in the current release
- non-rat workflows must provide a caller-supplied `ReferenceBundle`
- `ReferenceBundle.kinase_substrate_map` must be non-empty with:
  `kinase`, `substrate_site`
- `kinase` and `substrate_site` values must already be canonical strings
- duplicate `(kinase, substrate_site)` pairs are rejected
- `ReferenceBundle.site_sequences` must be non-empty with:
  `site_sequence` and unique canonical index
- each `substrate_site` must exist in `site_sequences.index`
- boundary constructors do not trim/canonicalize/deduplicate reference inputs
- bundled-reference provider/loading paths perform bundled-data shaping

## Workflow Validation

`KinaseWorkflowRequest` enforces:

- `dataset` is `AnalysisReadyPhosphoDataset`
- `references` is `ReferencePreset` or `ReferenceBundle`
- `scoring_config` is `KinaseScoringConfig`
- `scoring_config.min_substrates` is an int and must be `>= 2`
- `prediction_config` is `KinasePredictionConfig`
- `activity_config` is `KinaseActivityConfig` or `None`
- `activity_config.min_substrates` is an int and must be `>= 1`
- `activity_config.top_n_substrates` is an int and must be `>= 1`

Rewrite-era boundary diagnostics (raised as `WorkflowBoundaryError`) also enforce:

- interpreted reference coverage must overlap dataset phosphosites
- at least one kinase must meet `scoring_config.min_substrates` after overlap
- prediction ranking must produce at least one ensemble kinase
- activity (when enabled) requires prediction/phospho overlap and at least one
  valid kinase candidate after activity-stage filters

Boundary error messages include:

- the failing seam (for example `kinase.interpreter.eligible_kinases`)
- concrete counts (`dataset_sites`, `overlap_sites`, `eligible_kinases`, etc.)
- active config values (`scoring_config_min_substrates`,
  `prediction_config_ensemble_size`, `prediction_config_top_k`)
- a `next_action` hint for likely recovery

Stage result access is nested and stable:

- `result.scoring_result.profile_scores`
- `result.prediction_result.pred_mat`
- `result.activity_result.weighted_activity` (when enabled)
- `result.activity_result.ksea_scores` (when enabled)
- `result.activity_result.ksea_counts` (when enabled)
- `result.activity_result.target_counts` (when enabled)
- `result.activity_result.target_table` (when enabled)

`SignalomeWorkflowRequest` enforces:

- `kinase_result` is `KinaseWorkflowResult`
- `config` is `SignalomeConfig`
- `config.substrate_support_cutoff` is numeric in `[0.0, 1.0]`
- `config.network_correlation_threshold` is numeric in `[0.0, 1.0]`
- `kinase_result` prediction/scoring matrices are non-empty numeric DataFrames
- threshold semantics are explicit:
  - substrate support cutoff controls kinase-to-substrate inclusion from prediction
    scores (biological support selection)
  - network correlation threshold controls kinase-kinase edge inclusion from score
    correlations (graph sparsity control)

Signalome rewrite boundary diagnostics (raised as `WorkflowBoundaryError`) enforce:

- interpreted site alignment has usable overlap across dataset/prediction/score
- interpreted prediction/scoring matrices have overlapping kinase sets
- interpreted sites resolve to usable protein identifiers from one explicit path:
  `dataset.site_metadata.protein_id` (preferred) or non-empty protein prefixes in
  `dataset.phospho.index` / interpreted site IDs
- at least one kinase has support above `substrate_support_cutoff`
- module construction produces non-degenerate usable outputs
- network construction has required kinases, usable score variance, and at least one
  pair above `network_correlation_threshold`

Signalome boundary error messages include:

- the failing seam (for example `signalome.executor.network`)
- concrete counts (`dataset_sites`, `shared_kinases`, `supported_kinases`,
  `score_variance_kinases`, etc.)
- active config values (`substrate_support_cutoff`,
  `network_correlation_threshold`)
- a `next_action` hint for likely recovery

## Quick Troubleshooting

| Problem | Usually means | Good next step |
| --- | --- | --- |
| Builder rejects input format | A request field is neither a `DataFrame` nor a supported file path | Pass `DataFrame` values or supported file paths (`.csv`, `.tsv`, `.parquet`) |
| `ReferencePreset.AUTO` fails | `dataset.organism` is missing | Set `organism` on `DatasetBuildRequest` |
| Reference mismatch error | Dataset and selected preset organisms conflict | Align `dataset.organism` with `ReferencePreset` |
| Workflow request type error | Request field types are not the public models | Build the request from top-level `phospy` models |
| `kinase.interpreter.reference_coverage` | None of the reference substrate sites overlap `dataset.phospho.index` | Use references for the same identifier scheme/organism and verify site IDs |
| `kinase.interpreter.eligible_kinases` | Overlap exists, but no kinase reaches `scoring_config.min_substrates` | Lower `min_substrates` (not below `2`) or use references with deeper site overlap |
| `kinase.executor.prediction_ensemble` | Scoring completed, but no kinase had a finite prediction ranking | Provide at least two non-constant sample columns in `dataset.phospho` and/or lower `scoring_config.min_substrates` (not below `2`) |
| `kinase.activity.input_overlap` | Prediction and activity phospho matrices do not share sufficient phosphosite rows | Ensure `prediction_result.pred_mat` and `dataset.phospho` come from the same run and share site IDs |
| `kinase.activity.valid_candidates` | Activity stage filtered all kinases out | Lower `activity_config.min_substrates`, raise `activity_config.top_n_substrates`, or lower `activity_config.threshold` |
| `signalome.interpreter.site_alignment` | Dataset sites and interpreted scoring/prediction site IDs do not overlap | Ensure score/prediction outputs were generated from this dataset and share site IDs |
| `signalome.interpreter.kinase_overlap` | Score and prediction kinase columns have no shared kinase set | Regenerate kinase outputs so both matrices come from the same lane |
| `signalome.interpreter.protein_mapping` | Interpreted sites do not resolve to usable proteins | Populate `dataset.site_metadata.protein_id` or provide site IDs with non-empty protein prefixes |
| `signalome.executor.kinase_support` | No kinase has prediction support above `substrate_support_cutoff` | Lower `substrate_support_cutoff` or provide stronger prediction support |
| `signalome.executor.module_construction` | Module table collapsed to empty/trivial output | Increase kinase diversity and ensure multiple supported kinases |
| `signalome.executor.network` | Required kinases are missing from scores, score variance is unusable, or no pair passes the network correlation filter | Align score/prediction kinases, provide variable scoring signal, or lower `network_correlation_threshold` |

Runnable rewrite examples:

- [`../examples/dataset_builder_demo.py`](../examples/dataset_builder_demo.py)
- [`../examples/simple_workflow_demo.py`](../examples/simple_workflow_demo.py)
- [`../examples/signalome_workflow_demo.py`](../examples/signalome_workflow_demo.py)
