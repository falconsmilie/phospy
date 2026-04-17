# PhosPy

PhosPy is in rewrite cutover mode. Public support is intentionally narrow.

## Package Boundary

```text
src/phospy/         # supported rewrite package
src/phospy_legacy/  # migration reference only (not public rewrite contract)
```

## Supported Today

- Dataset builder route: supported
- Simple kinase workflow route: supported
- Signalome workflow route: first real vertical slice (module assignments, modules, network)

## Rewrite Contract

All public execution entry points use `run(request)`.

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

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=pd.DataFrame(
            {"sample_a": [1.0], "sample_b": [1.2]},
            index=["MAPK14;Y182;"],
        ),
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

profile_scores = result.scoring_result.profile_scores
pred_mat = result.prediction_result.pred_mat
activity_scores = result.activity_result.activity_scores
```

## Current Limits

- `DatasetBuildRequest` accepts pandas `DataFrame` inputs only.
- `ReferencePreset.AUTO` requires `dataset.organism`.
- Bundled references are currently available only for the rat lane.
- `SignalomeWorkflowResult.expanded_signalome` remains optional and is currently `None`.

## Examples

- [`examples/dataset_builder_demo.py`](examples/dataset_builder_demo.py)
- [`examples/simple_workflow_demo.py`](examples/simple_workflow_demo.py)
- [`examples/signalome_placeholder_demo.py`](examples/signalome_placeholder_demo.py)

## Docs

- [`docs/api.md`](docs/api.md)
- [`docs/validation.md`](docs/validation.md)
- [`docs/roadmap.md`](docs/roadmap.md)
