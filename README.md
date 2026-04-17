# PhosPy

PhosPy is in a rewrite cutover phase, with the first real vertical slice now
implemented for the simple kinase workflow.

## Package Boundary

```text
src/phospy/         # new architecture and supported public package
src/phospy_legacy/  # old implementation, reference-only during migration
```

- New implementation work must land under `src/phospy/`.
- `phospy_legacy` is internal migration reference material and is not a supported public API target.
- The simple kinase path is real for DataFrame-backed, analysis-ready inputs.
- Signalome in `src/phospy/` is still a placeholder shell.

## Current Public Shell API

Top-level `phospy` exports include:

- dataset and references: `AnalysisReadyPhosphoDataset`, `Organism`, `ReferencePreset`, `ReferenceBundle`
- builder lane: `DatasetBuildRequest`, `AnalysisReadyDatasetBuilder`
- workflows: `SimpleKinaseWorkflow`, `SignalomeWorkflow`
- workflow requests/results and stage config/result shells

All public execution entry points use `run(request)`.

## First Supported Route

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

result = SimpleKinaseWorkflow().run(
    SimpleKinaseWorkflowRequest(dataset=dataset, references=ReferencePreset.AUTO)
)

scoring = result.scoring_result.profile_scores
pred_mat = result.prediction_result.pred_mat
```

Notes:
- `DatasetBuildRequest` currently supports in-memory pandas `DataFrame` inputs only.
- `ReferencePreset.AUTO` requires `dataset.organism`.
- Built-in bundled references are currently supported for the rat lane.

## Docs

- [`docs/api.md`](docs/api.md)
- [`docs/architecture/rewrite_cutover_boundary.md`](docs/architecture/rewrite_cutover_boundary.md)
- [`docs/architecture/phospy_architecture_reset_notes.md`](docs/architecture/phospy_architecture_reset_notes.md)
