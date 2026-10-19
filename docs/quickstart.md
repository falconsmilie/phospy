# Run Your First Analysis

This walkthrough prepares a small phosphoproteomics dataset and tests one
treatment-versus-control comparison. Run the Python blocks in order in one
notebook or script: each block continues from the values created before it.

## Start With the Input Tables

A phosphosite is an amino-acid position on a protein that may carry a phosphate
group. The intensity table below contains measured phosphosite abundance. Each
row is one phosphosite and each column is one sample. These example values are
positive intensities on a linear scale.

```python
import pandas as pd

phospho = pd.DataFrame(
    {
        "control_1": [1000.0, 900.0],
        "control_2": [1050.0, 880.0],
        "treated_1": [1800.0, 930.0],
        "treated_2": [1750.0, 920.0],
    },
    index=["MAPK14;Y182;", "GSK3A;S21;"],
)
```

The intensity values need biological context. `site_metadata` provides one
aligned row of descriptive information for each phosphosite:

- `gene_symbol` gives the readable gene name and `protein_id` supplies the
  protein identity used for this example.
- `site` identifies the modified residue and its position.
- `site_sequence` gives the surrounding amino-acid sequence.
- `localisation_confidence` records how confidently the modification was
  assigned to that position.

```python
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "GSK3A"],
        "site": ["Y182", "S21"],
        "site_sequence": [
            "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
            "PSGGGPGGSGRARTSSFAEPGGGGGGGGGGP",
        ],
        "protein_id": ["MAPK14", "GSK3A"],
        "localisation_confidence": [0.95, 0.94],
    },
    index=phospho.index,
)
```

The two tables use the same row labels so PhosPy can align their measurements
and metadata. The [dataset preparation guide](api/dataset-build-workflow.md)
explains how to prepare larger input tables and use them with the builder.

## Choose the Preprocessing

PhosPy's differential workflow expects an established log2 intensity scale.
Because this example begins with positive linear values, configure a log2
transformation. This turns multiplicative differences into additive
differences, making the direction and size of `logFC` easier to interpret.
`policy="log2"` asks PhosPy to calculate `log2(value + pseudocount)`. The
default pseudocount is `1.0`.

```python
from phospy.advanced import DatasetIntensityTransformConfig

intensity_transform = DatasetIntensityTransformConfig(policy="log2")
```

Localisation confidence indicates how confidently a phosphate group was
assigned to a particular residue. This configuration tells PhosPy to read the
`localisation_confidence` metadata column and retain sites with a confidence of
at least `0.75`. The threshold is illustrative, and both example sites pass it.
For a real analysis, choose a threshold appropriate to the upstream
localisation method and study protocol.

```python
from phospy.api import DatasetLocalisationConfig

localisation = DatasetLocalisationConfig(
    mode="require_threshold",
    confidence_column="localisation_confidence",
    min_confidence=0.75,
)
```

`DatasetPreprocessingConfig` groups the preprocessing choices that should be
applied together when the dataset is built.

```python
from phospy.api import DatasetPreprocessingConfig

preprocessing_config = DatasetPreprocessingConfig(
    intensity_transform=intensity_transform,
    localisation=localisation,
)
```

## Describe and Build the Dataset

A `DatasetBuildRequest` describes the dataset you want PhosPy to create. It
combines the input tables, their organism, and the preprocessing configuration.
Creating the request records these choices. It does not build the dataset yet.

```python
from phospy.api import DatasetBuildRequest, Organism

dataset_request = DatasetBuildRequest(
    phospho=phospho,
    site_metadata=site_metadata,
    organism=Organism.RAT,
    preprocessing_config=preprocessing_config,
)
```

The builder provides the supported path from the request to an analysis-ready
dataset. Its `run()` method validates the request, aligns the tables, applies
the configured preprocessing, and returns the result.

```python
from phospy import AnalysisReadyDatasetBuilder

builder = AnalysisReadyDatasetBuilder()
dataset = builder.run(dataset_request)
```

## Inspect the Dataset

Check the matrix shape and established scale before starting an analysis. Then
look at a small metadata summary. PhosPy preserves the readable input labels in
`display_id`.

```python
dataset_shape = dataset.phospho.shape
intensity_scale = dataset.intensity_scale_state.label
site_summary = dataset.site_metadata.loc[
    :, ["display_id", "gene_symbol", "site", "localisation_confidence"]
].reset_index(drop=True)

print("Phosphosites by samples:", dataset_shape)
print("Intensity scale:", intensity_scale)
print(site_summary)
```

The expected shape is `(2, 4)`: two retained phosphosites measured in four
samples. The scale is `log2`, confirming that the requested transformation was
applied. The summary shows the readable site label, gene name, modified
residue, and localisation confidence that was checked during the build.

!!! note "This quickstart uses non-imputed data"
    Imputation means replacing missing measurements with estimated values. This
    example has no missing values and performs no imputation. If you select any
    dataset imputation policy, PhosPy's differential workflow rejects an
    upstream-imputed dataset by default unless you deliberately select the
    documented alternative policy. Read
    [how differential analysis handles imputed values](api/differential-analysis.md#missing-imputed-and-authoritative-matrix-policy)
    before analysing an imputed dataset.

## Describe the Experiment

Differential analysis also needs to know which samples belong to each
condition. A biological replicate is an independently measured biological
sample belonging to an experimental condition. A `SampleDesignRecord` connects
one matrix column to its condition and biological replicate. The collection
below describes two control replicates and two treated replicates.

```python
from phospy.api import SampleDesignRecord

sample_records = (
    SampleDesignRecord(
        sample_id="control_1",
        condition="control",
        biological_replicate_id="control_r1",
    ),
    SampleDesignRecord(
        sample_id="control_2",
        condition="control",
        biological_replicate_id="control_r2",
    ),
    SampleDesignRecord(
        sample_id="treated_1",
        condition="treated",
        biological_replicate_id="treated_r1",
    ),
    SampleDesignRecord(
        sample_id="treated_2",
        condition="treated",
        biological_replicate_id="treated_r2",
    ),
)
```

`ExperimentalDesign` holds that explicit sample-to-condition description.
PhosPy does not infer the design from sample names.

```python
from phospy.api import ExperimentalDesign

design = ExperimentalDesign(samples=sample_records)
```

A `Contrast` names the comparison to test. Here, `treated` is the numerator
and `control` is the denominator. The resulting `logFC` therefore represents
treated minus control on the log2 scale: a positive value means higher fitted
phosphosite intensity in the treated condition. An intensity difference alone
does not establish altered phosphorylation occupancy or a biological
mechanism.

```python
from phospy.api import Contrast

treated_vs_control = Contrast(
    name="treated_vs_control",
    numerator_condition="treated",
    denominator_condition="control",
)
```

The request accepts a tuple because one analysis can contain several
contrasts. This walkthrough uses only the single named contrast.

```python
contrasts = (treated_vs_control,)
```

Now combine the analysis-ready dataset, experimental design, and contrast in a
`DifferentialAnalysisRequest`.

```python
from phospy.api import DifferentialAnalysisRequest

differential_request = DifferentialAnalysisRequest(
    dataset=dataset,
    design=design,
    contrasts=contrasts,
)
```

## Run and Inspect the Contrast

`DifferentialAnalysisWorkflow.run()` validates the complete request and runs
the differential analysis.

```python
from phospy.api import DifferentialAnalysisWorkflow

differential_workflow = DifferentialAnalysisWorkflow()
differential_result = differential_workflow.run(differential_request)
```

The result model stores one table for each named contrast. Select this
contrast's table and display its readable site labels, estimated effects, and
statistical evidence columns.

```python
result_table = differential_result.table_for("treated_vs_control")
result_columns = ["display_id", "logFC", "P.Value", "adj.P.Val"]
result_view = result_table.loc[:, result_columns].reset_index(drop=True)

print("Displayed result shape:", result_view.shape)
print(result_view)
```

The displayed result has shape `(2, 4)`: two tested phosphosites and the four
selected columns. For these values, the MAPK14 site has a `logFC` of about
`0.792`, while the GSK3A site has a `logFC` of about `0.056`. `logFC` is the
fitted effect in the direction defined by the contrast. `P.Value` and
`adj.P.Val` summarise statistical evidence. They are not effect sizes. This
two-site dataset is kept small for learning, so its statistics demonstrate the
workflow rather than supporting a biological conclusion.

## Complete Runnable Example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.advanced import DatasetIntensityTransformConfig
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
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
        "protein_id": ["MAPK14", "GSK3A"],
        "localisation_confidence": [0.95, 0.94],
    },
    index=phospho.index,
)

intensity_transform = DatasetIntensityTransformConfig(policy="log2")

localisation = DatasetLocalisationConfig(
    mode="require_threshold",
    confidence_column="localisation_confidence",
    min_confidence=0.75,
)

preprocessing_config = DatasetPreprocessingConfig(
    intensity_transform=intensity_transform,
    localisation=localisation,
)

dataset_request = DatasetBuildRequest(
    phospho=phospho,
    site_metadata=site_metadata,
    organism=Organism.RAT,
    preprocessing_config=preprocessing_config,
)

builder = AnalysisReadyDatasetBuilder()
dataset = builder.run(dataset_request)

dataset_shape = dataset.phospho.shape
intensity_scale = dataset.intensity_scale_state.label
site_summary = dataset.site_metadata.loc[
    :, ["display_id", "gene_symbol", "site", "localisation_confidence"]
].reset_index(drop=True)

print("Phosphosites by samples:", dataset_shape)
print("Intensity scale:", intensity_scale)
print(site_summary)

sample_records = (
    SampleDesignRecord(
        sample_id="control_1",
        condition="control",
        biological_replicate_id="control_r1",
    ),
    SampleDesignRecord(
        sample_id="control_2",
        condition="control",
        biological_replicate_id="control_r2",
    ),
    SampleDesignRecord(
        sample_id="treated_1",
        condition="treated",
        biological_replicate_id="treated_r1",
    ),
    SampleDesignRecord(
        sample_id="treated_2",
        condition="treated",
        biological_replicate_id="treated_r2",
    ),
)

design = ExperimentalDesign(samples=sample_records)

treated_vs_control = Contrast(
    name="treated_vs_control",
    numerator_condition="treated",
    denominator_condition="control",
)

contrasts = (treated_vs_control,)

differential_request = DifferentialAnalysisRequest(
    dataset=dataset,
    design=design,
    contrasts=contrasts,
)

differential_workflow = DifferentialAnalysisWorkflow()
differential_result = differential_workflow.run(differential_request)

result_table = differential_result.table_for("treated_vs_control")
result_columns = ["display_id", "logFC", "P.Value", "adj.P.Val"]
result_view = result_table.loc[:, result_columns].reset_index(drop=True)

print("Displayed result shape:", result_view.shape)
print(result_view)
```

## Continue

- [Prepare a Dataset](api/dataset-build-workflow.md)
- [Differential Analysis](api/differential-analysis.md)
- [Kinase Analysis](api/kinase.md)
- [Enrichment](api/enrichment.md)
- [Signalome Analysis](api/signalome.md)
- [Scientific Interpretation and Limitations](scientific-interpretation.md)
