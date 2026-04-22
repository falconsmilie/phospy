# Quickstart: First Workflow

This guide gets you to a first successful PhosPy run with the supported happy path.

## 1. Install

```bash
pip install .
```

## 2. Build an Analysis-Ready Dataset

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBuildRequest,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    Organism,
)

phospho = pd.DataFrame(
    {
        "sample_a": [1.00, 0.70],
        "sample_b": [1.15, 0.80],
        "sample_c": [0.95, 0.75],
    },
    index=["MAPK14;Y182;", "GSK3B;S9;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "GSK3B"],
        "site": ["Y182", "S9"],
        "site_sequence": [
            "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
            "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
        ],
        "protein_id": ["MAPK14", "GSK3B"],
    },
    index=phospho.index.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        preprocessing_config=DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(policy="forbid"),
        ),
    )
)
```

## 3. Run Kinase Workflow

```python
from phospy import KinaseWorkflow
from phospy.api import (
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)

references = ReferenceBundle(
    organism=Organism.RAT,
    kinase_substrate_map=pd.DataFrame(
        {
            "kinase": ["MAP2K6", "MAP2K6", "AKT1", "AKT1"],
            "substrate_site": [
                "MAPK14;Y182;",
                "GSK3B;S9;",
                "MAPK14;Y182;",
                "GSK3B;S9;",
            ],
        }
    ),
    site_sequences=pd.DataFrame(
        {"site_sequence": dataset.site_metadata.loc[:, "site_sequence"]},
        index=pd.Index(dataset.site_metadata.index, name="site_id"),
    ),
)

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=2),
    )
)

pred_mat = kinase_result.prediction_result.pred_mat
```

## 4. (Optional) Run Signalome Workflow

```python
from phospy import SignalomeWorkflow
from phospy.api import SignalomeConfig, SignalomeWorkflowRequest

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=SignalomeConfig(
            substrate_support_cutoff=0.5,
            network_correlation_threshold=0.5,
        ),
    )
)
```

## 5. Use Official Examples

If you want runnable scripts instead of snippets:

- `python examples/dataset_builder_demo.py`
- `python examples/kinase_workflow_demo.py`
- `python examples/signalome_workflow_demo.py`

## What To Learn Next

- Concepts: [Core concepts](../concepts/core-concepts.md)
- Practical usage: [Tutorials and user guides](../user-guides/index.md)
- Full contract details: [API Guide](../api.md), [Validation Guide](../validation.md)
