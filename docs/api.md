# API Guide

PhosPy is in rewrite cutover mode. The first supported end-to-end route is:

`DatasetBuildRequest -> AnalysisReadyDatasetBuilder.run() -> AnalysisReadyPhosphoDataset -> SimpleKinaseWorkflow.run()`

## Migration Boundary

- `src/phospy/` is the authoritative package location for all new work.
- `src/phospy_legacy/` is reference-only and not a supported public API target.
- Legacy structures must not be extended.
- Temporary functional incompleteness in `phospy` is expected during this stage.

## Current Public Surface

Import from the top-level package:

- dataset model: `AnalysisReadyPhosphoDataset`
- reference models: `Organism`, `ReferencePreset`, `ReferenceBundle`
- builder request + entry point: `DatasetBuildRequest`, `AnalysisReadyDatasetBuilder`
- workflow requests: `SimpleKinaseWorkflowRequest`, `SignalomeWorkflowRequest`
- workflow configs: `KinaseScoringConfig`, `KinasePredictionConfig`, `KinaseActivityConfig`, `SignalomeConfig`
- stage results: `KinaseScoringResult`, `KinasePredictionResult`, `KinaseActivityResult`
- workflow results: `SimpleKinaseWorkflowResult`, `SignalomeWorkflowResult`
- workflow entry points: `SimpleKinaseWorkflow`, `SignalomeWorkflow`

All public execution entry points use `run(request)`.

## Supported Builder Input Mode

`DatasetBuildRequest` currently supports in-memory pandas `DataFrame` inputs only.

- required: `phospho`, `site_metadata`
- optional: `sample_metadata`, `total`
- unsupported in this phase: file-path and string-based ingestion
- bundled presets currently available in the rewrite path: rat (`ReferencePreset.RAT`)

## Example

```python
import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    Organism,
    ReferencePreset,
    SimpleKinaseWorkflow,
    SimpleKinaseWorkflowRequest,
)

builder = AnalysisReadyDatasetBuilder()
dataset = builder.run(
    DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        organism=Organism.RAT,
    )
)

kinase_result = SimpleKinaseWorkflow().run(
    SimpleKinaseWorkflowRequest(dataset=dataset, references=ReferencePreset.AUTO)
)

profile_scores = kinase_result.scoring_result.profile_scores
pred_mat = kinase_result.prediction_result.pred_mat
```

`SimpleKinaseWorkflowResult` keeps stage outputs nested (`scoring_result`, `prediction_result`, `activity_result`) rather than exposing mirrored top-level convenience aliases.

`SignalomeWorkflow` remains a placeholder in the rewrite tree and is not yet a
scientifically complete execution path.
