# PhosPy

PhosPy is a Python package for selected phosphoproteomics workflows inspired by
PhosR. It is aimed at scientists who want a clear Python lane from phosphosite
intensity tables to differential phosphorylation analysis, kinase scoring,
kinase prediction, and optional signalome analysis.

"PhosR-inspired" in PhosPy docs means scoped, feature-level comparison lanes. It
does not imply full PhosR package parity or full PhosR API compatibility.

PhosPy does **not** provide HTTP endpoints or a web service. The supported user
interface is the Python API.

## Recommended Reading
You can view the full documentation here: [Phospy Docs](https://phospy.com/docs)

## Install

PhosPy requires Python 3.10 or newer.

```bash
pip install phospy
```

For `.parquet` input or output support:

```bash
pip install "phospy[parquet]"
```

For local development from a clone:

```bash
pip install -e ".[dev]"
pyright
pytest -m "not parity"
```

For reproducible scientific/regression runs aligned to CI:

```bash
pip install -c constraints/ci.txt -e ".[dev,test]"
pytest tests/parity -m parity -s
```

For full release-gate validation (unit/integration, reproducibility goldens,
parity, and performance):

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet]"
make test-release-gate
```

## Quick Start

1. Build an analysis-ready phosphoproteomics dataset.
2. Run a kinase workflow.
3. Explore full API workflow documentation:
   - [Dataset building](docs/api/dataset-build-workflow.md)
   - [Differential workflow](docs/api/differential-workflow.md)
   - [Kinase workflow](docs/api/kinase-workflow.md)
   - [Signalome workflow](docs/api/signalome-workflow.md)

Bundled runtime references in the current release are rat-only. For human or
mouse work, create and pass an explicit `ReferenceBundle` in Python instead of
using `ReferencePreset.AUTO`.

Scientific scope categories and parity/open-gap status are maintained in
[`docs/scientific-coverage.md`](docs/scientific-coverage.md). Parity fixture
evidence lives in [`docs/parity.md`](docs/parity.md).

## Kinase Workflow Example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.api import (
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
)

# Tiny synthetic example for workflow mechanics only (not biological discovery).
phospho = pd.DataFrame(
    {
        "control_rep1": [8200.0, 9100.0, 6000.0],
        "control_rep2": [8000.0, 9000.0, 5900.0],
        "treatment_rep1": [16200.0, 9150.0, 13000.0],
        "treatment_rep2": [15800.0, 9050.0, 12800.0],
    },
    index=["MAPK14;Y182;", "GSK3B;S9;", "TSC2;S939;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "GSK3B", "TSC2"],
        "site": ["Y182", "S9", "S939"],
        "site_sequence": [
            "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
        ],
        "protein_id": ["MAPK14", "GSK3B", "TSC2"],
        "localisation_confidence": [0.95, 0.94, 0.96],
    },
    index=phospho.index.copy(),
)
sample_metadata = pd.DataFrame(
    {
        "condition": ["control", "control", "treatment", "treatment"],
    },
    index=phospho.columns.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        organism=Organism.RAT,
        preprocessing_config=DatasetPreprocessingConfig(
            # Site-level workflows should fail fast when localisation is missing
            # or below threshold, because ambiguous site assignment can
            # mis-state kinase/substrate interpretation.
            localisation=DatasetLocalisationConfig(
                mode="require_threshold",
                confidence_column="localisation_confidence",
                min_confidence=0.75,
            )
        ),
    )
)

# Dataset construction validates required site metadata, including site_sequence.
print(dataset.site_metadata.loc[:, ["gene_symbol", "site", "site_sequence"]])

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        activity_config=None,
    )
)

print(kinase_result.prediction_result.pred_mat.round(3).iloc[:3, :5])
if kinase_result.prediction_result.substrate_list is not None:
    print(kinase_result.prediction_result.substrate_list.head(5))
```

## Import Contract

Use top-level `phospy` for the five main entrypoints:

```python
from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy import DifferentialAnalysisWorkflow, KinaseWorkflow, SignalomeWorkflow
```

Use `phospy.api` for requests, configs, results, enums, references, and public
exceptions.

## Documentation

1. [Quickstart](https://phospy.com/docs/quickstart/)
2. [API Guide](https://phospy.com/docs/api/)
3. [Workflow Contracts](https://phospy.com/docs/workflow_contracts/)
4. [Validation Guide](https://phospy.com/docs/validation/)
5. [Scientific Coverage Matrix](https://phospy.com/docs/scientific-coverage/)

## Citation

If you use PhosPy in scientific work, cite this software release using
[`CITATION.cff`](CITATION.cff) and also cite the upstream PhosR project and
publications described in [`NOTICE.md`](NOTICE.md).
