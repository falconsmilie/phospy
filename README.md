# PhosPy

[![PyPI version](https://img.shields.io/pypi/v/phospy.svg)](https://pypi.org/project/phospy/)
[![Python versions](https://img.shields.io/pypi/pyversions/phospy.svg)](https://pypi.org/project/phospy/)
[![Tests](https://github.com/falconsmilie/phospy/actions/workflows/ci.yml/badge.svg)](https://github.com/falconsmilie/phospy/actions/workflows/ci.yml)
[![License](https://img.shields.io/pypi/l/phospy.svg)](https://github.com/falconsmilie/phospy/blob/main/LICENSE)

PhosPy helps you turn phosphosite intensity tables into analysis-ready datasets
and run focused phosphoproteomics workflows from Python.

Use PhosPy to prepare and validate phosphosite data, test differential
phosphorylation, run offline enrichment, explore kinase-substrate support, and
build signalome summaries. The supported user interface is the Python API.

## Install PhosPy

PhosPy supports Python 3.11 and 3.12.

```bash
pip install phospy
```

Add Parquet support when you need it:

```bash
pip install "phospy[parquet]"
```

## Try a Kinase Analysis

The public workflow follows a simple pattern: prepare a dataset, create a
request, and run the workflow. The complete example below uses a small rat
dataset and bundled rat references.

<details>
<summary><strong>View the Complete Kinase Example</strong></summary>

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.advanced import KinaseReliabilityProfile, KinaseScoringConfig
from phospy.api import (
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    IntensityScaleKind,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
)

phospho = pd.DataFrame(
    {
        "control_1": [8200.0, 9100.0, 6000.0],
        "control_2": [8000.0, 9000.0, 5900.0],
        "treated_1": [16200.0, 9150.0, 13000.0],
    },
    index=["MAPK14;Y182;", "GSK3A;S21;", "TSC2;S939;"],
)

site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "GSK3A", "TSC2"],
        "site": ["Y182", "S21", "S939"],
        "site_sequence": [
            "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
            "PSGGGPGGSGRARTSSFAEPGGGGGGGGGGP",
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
        ],
        "protein_identifier": ["MAPK14", "GSK3A", "TSC2"],
        "protein_group_id": ["MAPK14", "GSK3A", "TSC2"],
        "localisation_confidence": [0.95, 0.94, 0.96],
    },
    index=phospho.index,
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        input_intensity_scale=IntensityScaleKind.LINEAR,
        preprocessing_config=DatasetPreprocessingConfig(
            localisation=DatasetLocalisationConfig(
                confidence_column="localisation_confidence",
                min_confidence=0.75,
            )
        ),
    )
)

result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(
            reliability_profile=KinaseReliabilityProfile.CUSTOM,
        ),
        activity_config=None,
    )
)

print(result.prediction_result.pred_mat.head())
```

The values are illustrative. Use study-specific references, biological
replicates, and thresholds for real analyses.

</details>

The main dataset and workflow entry points are available from `phospy`:

```python
from phospy import AnalysisReadyDatasetBuilder
from phospy import DifferentialAnalysisWorkflow, KinaseWorkflow, SignalomeWorkflow
```

Enrichment is available from `phospy.api`:

```python
from phospy.api import EnrichmentConfig, EnrichmentWorkflow, EnrichmentWorkflowRequest
```

## Read the Guides

- [Prepare a dataset](docs/api/dataset-build-workflow.md)
- [Run differential analysis](docs/api/differential-analysis.md)
- [Run enrichment](docs/api/enrichment.md)
- [Run kinase analysis](docs/api/kinase.md)
- [Run signalome analysis](docs/api/signalome.md)

The full documentation is available at [PhosPy Docs](https://phospy.com/docs/).
See the [contributing guide](docs/contributing.md), [license](LICENSE), and
[`CITATION.cff`](CITATION.cff) for project and citation details.

## Release Checks

Maintainers should run the lightweight release gate before publishing:

```bash
make release-check
```

This provides normal CI/build confidence, not formal
exact-source/exact-artifact attestation.

## Scientific Scope

PhosPy is careful about what each result can support. Read the
[interpretation guide](docs/scientific-interpretation.md) before drawing
biological conclusions from an analysis.

<details>
<summary><strong>Review the Current Scientific Boundaries</strong></summary>

PhosPy is inspired by PhosR, but differential analysis is limited to tested
design and contrast envelopes; it is not full limma or PhosR parity.

Bundled runtime references are rat-only. `ReferencePreset.AUTO` is therefore
for rat datasets; human, mouse, and custom analyses must pass an explicit
`ReferenceBundle`. Kinase scores show relative support within a run. They are
not calibrated probabilities, and activity outputs are not direct proof of
kinase activation or causal pathway activity.

Enrichment provides offline over-representation analysis (ORA) with
caller-supplied identifier sets and an explicit background. ORA does not imply
GSEA or PTM-SEA support.

`linear_residualize_batch`, a limited fixed-effect residualisation, rejects
confounded batch/condition designs. It is not ComBat, not RUV, not limma
`removeBatchEffect` parity, and not mixed-effects modelling. Native
SPS/RUV-style correction through `SpsRuvBatchCorrectionConfig` is a separate,
explicit preprocessing method. It is not PhosR-equivalent SPS/RUV-III parity
and not PhosR-equivalent batch correction. Replicate metadata is checked and
recorded but is not used for numerical unwanted-factor estimation.
`ruv_readiness` values are report-only readiness signals and do not apply
correction; use native correction only after those prerequisites are
implemented in the request.

</details>
