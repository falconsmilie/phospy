# API Guide

PhosPy is in rewrite cutover mode.

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

## Example

```python
import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    ReferencePreset,
    SimpleKinaseWorkflow,
    SimpleKinaseWorkflowRequest,
)

builder = AnalysisReadyDatasetBuilder()
dataset = builder.run(
    DatasetBuildRequest(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["GENEA;S1;"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["GENEA"],
                "site": ["S1"],
                "site_sequence": ["AAAAAAA"],
            },
            index=["GENEA;S1;"],
        ),
    )
)

kinase_result = SimpleKinaseWorkflow().run(
    SimpleKinaseWorkflowRequest(dataset=dataset, references=ReferencePreset.AUTO)
)
```

`SimpleKinaseWorkflowResult` keeps stage outputs nested (`scoring_result`, `prediction_result`, `activity_result`) rather than exposing mirrored top-level convenience aliases.
