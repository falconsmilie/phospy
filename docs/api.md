# API Guide

PhosPy currently exposes one supported end-to-end rewrite route:

`DatasetBuildRequest -> AnalysisReadyDatasetBuilder.run(request) -> AnalysisReadyPhosphoDataset -> KinaseWorkflow.run(request)`

## Support Boundary

- `src/phospy/`: supported public rewrite package
- `legacy_archive/phospy_legacy/`: migration reference only; not an installed package target
- Dataset builder: supported
- Kinase workflow: supported
- Signalome workflow: first real vertical slice implemented

## Public Types

Import from top-level `phospy`:

- Dataset and references:
`AnalysisReadyPhosphoDataset`, `Organism`, `ReferencePreset`, `ReferenceBundle`
- Builder:
`DatasetBuildRequest`, `AnalysisReadyDatasetBuilder`
- Workflows and requests:
`KinaseWorkflow`, `KinaseWorkflowRequest`,
`SignalomeWorkflow`, `SignalomeWorkflowRequest`
- Config models:
`KinaseScoringConfig`, `KinasePredictionConfig`, `KinaseActivityConfig`,
`SignalomeConfig`
- Result models:
`KinaseWorkflowResult`, `SignalomeWorkflowResult`,
`KinaseScoringResult`, `KinasePredictionResult`, `KinaseActivityResult`

All public executors use `run(request)`.

## User-Handleable Exceptions

Public constructors and workflow/builder boundaries raise PhosPy exception types
as the intended failure story. User-handleable exceptions are exported from
top-level `phospy` and include:

- Base:
`PhosPyError`
- Input/build:
`PhosPyInputError`, `UnsupportedInputFormatError`, `PhosPyBuildError`,
`DatasetBuildError`
- Validation:
`PhosPyValidationError`, `DatasetValidationError`, `ReferenceValidationError`,
`TransformationValidationError`, `WorkflowValidationError`
- Reference:
`PhosPyReferenceError`, `ReferenceResolutionError`,
`ReferenceCompatibilityError`, `UnsupportedOrganismError`
- Transformation:
`PhosPyTransformationError`, `InvalidTransformationStateError`,
`TransformationStateEstablishmentError`, `TransformerExecutionError`
- Workflow:
`PhosPyWorkflowError`, `WorkflowBoundaryError`, `WorkflowStageError`

## Builder Contract

`DatasetBuildRequest` accepts in-memory pandas `DataFrame` values or file paths.

- Required: `phospho`, `site_metadata`
- Optional: `sample_metadata`, `total`, `organism`
- `site_metadata.protein_id` is an optional identity column and is preserved when
  provided (it is not treated as `gene_symbol`)
- Transformation state is established inside PhosPy through the supported
  transformer path at builder execution time.
- If callers construct `AnalysisReadyPhosphoDataset` directly (instead of using
  the builder), `transformation_state` must be provided explicitly.
- Public dataset/reference models are strict boundaries: they validate but do not
  trim/canonicalize/deduplicate inputs.
- Dirty or colliding site identifiers at the public boundary fail fast.
- Input shaping and cleanup responsibilities live below the boundary in builder
  and bundled-reference provider/loading paths.
- `ReferencePreset.AUTO` requires `dataset.organism` at workflow execution time
- Bundled references are currently available for rat only
- `ReferencePreset.HUMAN` and `ReferencePreset.MOUSE` are part of the public enum,
  but bundled runtime resolution for those presets is intentionally unsupported in
  this release and fails with `UnsupportedOrganismError`
- For non-rat workflows, supply an explicit `ReferenceBundle` instead of relying
  on bundled presets

## Frame Ownership

PhosPy applies a package-wide ownership rule for pandas frames:

- public boundary objects own copies of caller-provided frames
- internal stage DTOs may alias already-owned frames
- internal assembly paths may transfer ownership without re-copying

See [`frame_ownership.md`](frame_ownership.md) for the full policy.

## Result Contract

`KinaseWorkflowResult` keeps stage outputs nested:

- `result.scoring_result`
- `result.prediction_result`
- `result.activity_result`

Supported scoring behavior is currently profile-only:

- `result.scoring_result.profile_scores` is the score matrix used for prediction
  and downstream signalome execution.
- `result.scoring_result.motif_scores`, `combined_scores`, and `weights` are
  optional compatibility fields and are `None` on this supported route.

No mirrored top-level convenience aliases are part of the rewrite contract.
Reusable save/load output bundles are provided as external services in
`phospy.io` (`save_kinase_workflow_bundle`,
`load_kinase_workflow_bundle`) so result models remain plain containers.

`SignalomeWorkflowResult` currently returns real tables (`module_assignments`,
`signalome_modules`, `kinase_network`) from the first vertical slice, while
`expanded_signalome` remains optional and currently `None`.
`SignalomeConfig` now uses two explicit thresholds:

- `substrate_support_cutoff`: prediction-score cutoff for kinase substrate support
- `network_correlation_threshold`: absolute correlation cutoff for kinase network
  edges

Signalome protein grouping resolves from explicit protein identity:
`dataset.site_metadata.protein_id` when present, otherwise interpreted site-ID
protein prefixes. Missing protein identity fails with `WorkflowBoundaryError`.

Kinase activity output exposes the full downstream activity stage:

- `result.activity_result.weighted_activity`
- `result.activity_result.ksea_scores`
- `result.activity_result.ksea_counts`
- `result.activity_result.target_counts`
- `result.activity_result.target_table`

Activity configuration controls:

- `activity_config.threshold`
- `activity_config.min_substrates`
- `activity_config.top_n_substrates`

Kinase scoring support floor:

- `scoring_config.min_substrates` defaults to `2`
- values below `2` are rejected by validation
- kinases with fewer than `min_substrates` quantified overlap sites are excluded
  from profile/scoring stages

See [`output_bundles.md`](output_bundles.md) for manifest format and table
inventory.

## Example

```python
import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    KinaseScoringConfig,
    Organism,
    ReferenceBundle,
    KinaseWorkflow,
    KinaseWorkflowRequest,
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 0.7], "sample_b": [1.2, 0.8]},
            index=["MAPK14;Y182;", "GSK3B;S9;"],
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "GSK3B"],
                "site": ["Y182", "S9"],
                "site_sequence": [
                    "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                    "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
                ],
            },
            index=["MAPK14;Y182;", "GSK3B;S9;"],
        ),
        organism=Organism.RAT,
    )
)

references = ReferenceBundle(
    organism=Organism.RAT,
    kinase_substrate_map=pd.DataFrame(
        {
            "kinase": ["MAP2K6", "MAP2K6"],
            "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
        }
    ),
    site_sequences=pd.DataFrame(
        {"site_sequence": dataset.site_metadata.loc[:, "site_sequence"]},
        index=pd.Index(dataset.site_metadata.index, name="site_id"),
    ),
)

result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=2),
    )
)

pred_mat = result.prediction_result.pred_mat
```

## CLI

The rewrite CLI now supports the same narrow lane from files:

1. `dataset-build` (files -> analysis-ready dataset outputs)
2. `kinase` (files -> dataset build + kinase workflow outputs)

Examples:

```bash
phospy dataset-build \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --outdir ./out

phospy kinase \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --reference auto \
  --outdir ./out
```

Input files currently support `.csv`, `.tsv`/`.txt`, and `.parquet`.
