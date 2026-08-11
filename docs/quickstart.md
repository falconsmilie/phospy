# Run Your First Analysis

This example builds a small log2 dataset and tests one treatment-versus-control
contrast. It is intentionally compact so you can check the installation and the
public API before using your own data.

## Build the Dataset

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.advanced import DatasetIntensityTransformConfig
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)

phospho = pd.DataFrame(
    {
        "control_1": [1000.0, 900.0],
        "control_2": [1050.0, 880.0],
        "treated_1": [1800.0, 930.0],
        "treated_2": [1750.0, 920.0],
    },
    index=["MAPK14;Y182;", "GSK3A;S21;"],
)

site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "GSK3A"],
        "site": ["Y182", "S21"],
        "site_sequence": [
            "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
            "PSGGGPGGSGRARTSSFAEPGGGGGGGGGGP",
        ],
        "protein_identifier": ["MAPK14", "GSK3A"],
        "localisation_confidence": [0.95, 0.94],
    },
    index=phospho.index,
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        preprocessing_config=DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            localisation=DatasetLocalisationConfig(
                confidence_column="localisation_confidence",
                min_confidence=0.75,
            ),
        ),
    )
)
```

## Describe the Experiment and Run the Contrast

```python
design = ExperimentalDesign(
    samples=(
        SampleDesignRecord("control_1", "control", "control_r1"),
        SampleDesignRecord("control_2", "control", "control_r2"),
        SampleDesignRecord("treated_1", "treated", "treated_r1"),
        SampleDesignRecord("treated_2", "treated", "treated_r2"),
    )
)

result = DifferentialAnalysisWorkflow().run(
    DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=(
            Contrast(
                name="treated_vs_control",
                numerator_condition="treated",
                denominator_condition="control",
            ),
        ),
    )
)

print(
    result.table_for("treated_vs_control").loc[
        :, ["display_id", "logFC", "P.Value", "adj.P.Val"]
    ]
)
```

Positive `logFC` values indicate higher fitted phosphorylation in the treated
condition. Adjusted *p* values describe statistical evidence within this
contrast; they are not effect sizes.

## Continue

- [Differential Analysis](api/differential-analysis.md)
- [Kinase Analysis](api/kinase.md)
- [Enrichment](api/enrichment.md)
- [Signalome Analysis](api/signalome.md)
- [Scientific Interpretation and Limitations](scientific-interpretation.md)
