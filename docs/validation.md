# Validation Guide

This guide documents validation in the supported rewrite contract.

For API shapes, see [`api.md`](api.md).

## Boundary Model

Validation is split across three boundaries:

- Builder request boundary:
  validates supported source types (`DataFrame` or file path) and builder conventions.
- Final dataset boundary:
  validates strict `AnalysisReadyPhosphoDataset` invariants.
- Workflow boundaries:
  validate request DTOs, config ranges, reference compatibility, and runtime/science seams.

## Builder Flexibility vs Dataset Strictness

`DatasetBuildRequest` accepts:

- `pandas.DataFrame`
- file path (`str`, `Path`, `PathLike`) to `.csv`, `.tsv`/`.txt`, `.parquet`

This flexibility exists only at the builder input boundary.

After source loading, both routes pass through the same normalization path:

- site ID canonicalization
- site-metadata alias resolution
- optional derivation of `gene_symbol` and `site` from strict index format
- fail-fast rejection of ambiguous/unsupported legacy aliases

`AnalysisReadyPhosphoDataset` itself is strict and DataFrame-only.
Workflows consume only this dataset type.

## AnalysisReady Dataset Validation

`AnalysisReadyPhosphoDataset` constructor composes:

- `AnalysisReadyDatasetValidator`
- `TransformationStateValidator`

Enforced dataset invariants:

- `phospho`: non-empty numeric DataFrame, unique index/columns,
  canonical site IDs.
- `site_metadata`: non-empty DataFrame, exact index alignment with `phospho.index`,
  required columns `gene_symbol`, `site`, `site_sequence` with non-empty strings.
- `sample_metadata` (if present): index aligns to `phospho.columns`.
- `total` (if present): non-empty numeric DataFrame, unique index,
  columns align to `phospho.columns`.
- `organism` (if present): `Organism` enum.
- `transformation_state`: typed `TransformationState`, coherent with presence/absence
  of `total`, and coherent transformation kind between matrices when both exist.

Boundary constructors validate; they do not silently repair invalid data.

## Builder Convention Rules

Supported site-metadata aliases:

- `gene_symbol`: `gene_symbol`, `gene_name`
- `site`: `site`
- `site_sequence`: `site_sequence`, `centralized_sequence`
- `protein_id`: `protein_id`

Unsupported legacy aliases are rejected:

- `gene`
- `residue`
- `phosphosite`
- `site_position`
- `sequence`
- `protein`

If `gene_symbol` and/or `site` are missing, builder derivation is allowed only from
index values exactly matching `"<gene_symbol>;<site>;"`.

If both `site_metadata.index` and `site_metadata.site_id` are present, they must match
after canonicalization.

## Reference Validation

Reference resolution/compatibility enforces:

- `ReferencePreset.AUTO` requires `dataset.organism`.
- Explicit preset must match dataset organism when both are set.
- Explicit `ReferenceBundle.organism` must match dataset organism when both are set.
- Bundled runtime references are currently rat-only.
- `ReferencePreset.HUMAN` and `ReferencePreset.MOUSE` are valid enum values but are
  not bundled runtime lanes in this release.

`ReferenceBundle` structure/content validation enforces non-empty, canonical,
internally consistent reference tables.

## Workflow Validation

`KinaseWorkflowValidator` enforces:

- request type is `KinaseWorkflowRequest`
- `dataset` is `AnalysisReadyPhosphoDataset`
- `references` is `ReferencePreset | ReferenceBundle`
- config types and numeric floors/ranges:
  `scoring_config.min_substrates >= 2`,
  `prediction_config.top_k >= 1`,
  `prediction_config.ensemble_size >= 1`,
  activity config bounds when activity is enabled/configured

`SignalomeWorkflowValidator` enforces:

- request type is `SignalomeWorkflowRequest`
- `kinase_result` is `KinaseWorkflowResult`
- signalome config bounds in `[0.0, 1.0]`
- upstream score/prediction matrices are usable numeric matrices for signalome execution

Interpreters/executors enforce seam-level scientific/runtime boundary checks and raise
`WorkflowBoundaryError` with seam names, concrete counts, and `next_action` hints.

## Nested Result Access (Validation-Relevant Contract)

Stable access paths are nested by stage:

- `result.scoring_result.profile_scores`
- `result.scoring_result.motif_scores`
- `result.scoring_result.combined_scores`
- `result.scoring_result.weights`
- `result.prediction_result.pred_mat`
- `result.activity_result.weighted_activity` (when activity is enabled)
- `signalome_result.kinase_result.prediction_result.pred_mat`

Optional outputs must be checked before dereference:

- `result.activity_result`
- `result.prediction_result.substrate_list`
- `signalome_result.kinase_network.nodes`
- `signalome_result.expanded_signalome`

## Validation Ownership Summary

| Invariant | Owner |
| --- | --- |
| Builder input source type checks | `DatasetBuildRequestValidator` + `DatasetInputSourceValidator` |
| Builder convention normalization/derivation | `DatasetBuildRequestInterpreter` collaborators |
| Analysis-ready dataset structure/content | `AnalysisReadyDatasetValidator` |
| Transformation-state coherence | `TransformationStateValidator` |
| Reference compatibility (dataset vs preset/bundle) | `ReferenceCompatibilityValidator` / `ReferenceResolver` |
| Reference bundle structure/content | `ReferenceBundleValidator` |
| Kinase workflow request/config validity | `KinaseWorkflowValidator` + `WorkflowConfigValidator` |
| Signalome workflow request/config validity | `SignalomeWorkflowValidator` + `WorkflowConfigValidator` |
| Kinase/signalome runtime seam diagnostics | workflow interpreters/executors (`WorkflowBoundaryError`) |

## Quick Troubleshooting

| Problem | Usually means | Good next step |
| --- | --- | --- |
| Builder rejects input format | Field is neither DataFrame nor supported file path | Pass DataFrame or path to `.csv`/`.tsv`/`.txt`/`.parquet` |
| Dataset constructor fails on site metadata | Required strict boundary columns/values are missing | Provide `gene_symbol`, `site`, `site_sequence` with non-blank strings |
| `ReferencePreset.AUTO` fails | Dataset organism is missing | Set `organism` in `DatasetBuildRequest` |
| Bundled human/mouse preset fails | Bundled references are rat-only in this release | Provide explicit non-rat `ReferenceBundle` |
| Kinase boundary seam fails | Overlap/support constraints were not met | Read seam details and adjust dataset/references/config |
| Signalome protein-mapping seam fails | Interpreted sites did not resolve to protein identity | Provide `site_metadata.protein_id` or resolvable site-ID protein prefixes |
