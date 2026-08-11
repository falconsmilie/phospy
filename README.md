# PhosPy

[![PyPI version](https://img.shields.io/pypi/v/phospy.svg)](https://pypi.org/project/phospy/)
[![Python versions](https://img.shields.io/pypi/pyversions/phospy.svg)](https://pypi.org/project/phospy/)
[![Tests](https://github.com/falconsmilie/phospy/actions/workflows/ci.yml/badge.svg)](https://github.com/falconsmilie/phospy/actions/workflows/ci.yml)
[![License](https://img.shields.io/pypi/l/phospy.svg)](https://github.com/falconsmilie/phospy/blob/main/LICENSE)

PhosPy is a Python package for selected phosphoproteomics workflows.

It helps you build strict phosphosite datasets, run differential phosphorylation
analysis, run offline over-representation analysis (ORA), score kinase support,
predict kinase-substrate support, and create optional signalome summaries. The
supported interface is the Python API.

PhosPy is PhosR-inspired, but current differential analysis is limited to tested
design and contrast envelopes; it is not full limma or PhosR parity.

## Installation

PhosPy requires Python 3.11 or 3.12.

```bash
pip install phospy
```

For Parquet input/output:

```bash
pip install "phospy[parquet]"
```

For local development:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet,docs]"
```

For maintained release checks, run `make release-check`. It provides normal
CI/build confidence, not formal exact-source/exact-artifact attestation.

## Compact example

This example builds a small rat dataset and runs the public kinase workflow.
Real analyses should use study-specific biological replicates, references, and
reporting thresholds.

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.api import (
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    IntensityScaleKind,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
)
from phospy.advanced import KinaseReliabilityProfile, KinaseScoringConfig

phospho = pd.DataFrame(
    {
        "control_rep1": [8200.0, 9100.0, 6000.0],
        "control_rep2": [8000.0, 9000.0, 5900.0],
        "treatment_rep1": [16200.0, 9150.0, 13000.0],
        "treatment_rep2": [15800.0, 9050.0, 12800.0],
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
        "display_id": ["MAPK14;Y182;", "GSK3A;S21;", "TSC2;S939;"],
        "organism": ["rat", "rat", "rat"],
        "protein_namespace": ["protein_id", "protein_id", "protein_id"],
        "protein_identifier": ["MAPK14", "GSK3A", "TSC2"],
        "protein_group_id": ["MAPK14", "GSK3A", "TSC2"],
        "localisation_confidence": [0.95, 0.94, 0.96],
    },
    index=phospho.index.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        input_intensity_scale=IntensityScaleKind.LINEAR,
        preprocessing_config=DatasetPreprocessingConfig(
            localisation=DatasetLocalisationConfig(
                mode="require_threshold",
                confidence_column="localisation_confidence",
                min_confidence=0.75,
            )
        ),
    )
)

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(
            reliability_profile=KinaseReliabilityProfile.CUSTOM,
        ),
        activity_config=None,
    )
)

print(kinase_result.prediction_result.pred_mat.head())
```

Use top-level imports for the main dataset and workflow entry points:

```python
from phospy import AnalysisReadyDatasetBuilder
from phospy import DifferentialAnalysisWorkflow, KinaseWorkflow, SignalomeWorkflow
```

`EnrichmentWorkflow` is public through `phospy.api`:

```python
from phospy.api import EnrichmentConfig, EnrichmentWorkflow, EnrichmentWorkflowRequest
```

## Supported workflows

- [Dataset preparation](docs/api/dataset-build-workflow.md)
- [Differential analysis](docs/api/differential-analysis.md)
- [Enrichment](docs/api/enrichment.md)
- [Kinase analysis](docs/api/kinase.md)
- [Signalome analysis](docs/api/signalome.md)

The full documentation starts at [PhosPy Docs](https://phospy.com/docs/).
See the [contributing guide](docs/contributing.md), [license](LICENSE), and
[`CITATION.cff`](CITATION.cff) for project and citation details.

## Scientific-use notes

Bundled runtime references are rat-only in this release. Use
`ReferencePreset.AUTO` only for rat datasets; for human or mouse work, pass an
explicit `ReferenceBundle`.

Kinase scores are relative support within a run, not calibrated probabilities.
Kinase activity outputs are not direct proof of kinase activation or causal
pathway activity. Enrichment ORA uses caller-supplied selected identifiers, set
collections, and a background universe; ORA does not imply GSEA or PTM-SEA
support.

Batch-correction preprocessing is explicit. `linear_residualize_batch` is
limited fixed-effect residualisation; it rejects confounded batch/condition
metadata and is not ComBat, not RUV, not limma `removeBatchEffect` parity, not
mixed-effects modelling, not native SPS/RUV-style correction, and not
PhosR-equivalent batch correction. Native SPS/RUV-style correction through
`SpsRuvBatchCorrectionConfig` has implemented prerequisites, but it is not
PhosR-equivalent SPS/RUV-III parity; those prerequisites are implemented as a
native PhosPy preprocessing lane, not current PhosR parity. Replicate metadata is validated and
recorded for provenance and diagnostics only; it is not used for numerical
unwanted-factor estimation. `ruv_readiness` diagnostics are report-only
RUV-readiness metadata; they are readiness signals and do not apply correction.
