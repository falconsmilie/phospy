# PhosPy

PhosPy is in a rewrite cutover phase.

## Package Boundary

```text
src/phospy/         # new architecture and supported public package
src/phospy_legacy/  # old implementation, reference-only during migration
```

- New implementation work must land under `src/phospy/`.
- `phospy_legacy` is internal migration reference material and is not a supported public API target.
- Temporary functional incompleteness in the new package is expected at this stage.

## Current Public Shell API

Top-level `phospy` exports include:

- dataset and references: `AnalysisReadyPhosphoDataset`, `Organism`, `ReferencePreset`, `ReferenceBundle`
- builder lane: `DatasetBuildRequest`, `AnalysisReadyDatasetBuilder`
- workflows: `SimpleKinaseWorkflow`, `SignalomeWorkflow`
- workflow requests/results and stage config/result shells

All public execution entry points use `run(request)`.

## Quick Smoke Example

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

result = SimpleKinaseWorkflow().run(
    SimpleKinaseWorkflowRequest(dataset=dataset, references=ReferencePreset.AUTO)
)
```

## Docs

- [`docs/api.md`](docs/api.md)
- [`docs/architecture/rewrite_cutover_boundary.md`](docs/architecture/rewrite_cutover_boundary.md)
- [`docs/architecture/phospy_architecture_reset_notes.md`](docs/architecture/phospy_architecture_reset_notes.md)
