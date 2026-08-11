# Running your first analysis

This quickstart runs one small kinase analysis. It is intentionally small so you
can confirm the installation, dataset preparation, references, and public API
imports before working with a real experiment.

## Install

```bash
pip install phospy
```

## Build a small dataset and run kinase analysis

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.advanced import (
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
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
        "sample_a": [1.00, 0.70, 0.85],
        "sample_b": [1.10, 0.80, 0.88],
        "sample_c": [0.95, 0.75, 0.92],
    },
    index=["TSC2;S939;", "GSK3A;S21;", "MAPK14;Y182;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["TSC2", "GSK3A", "MAPK14"],
        "site": ["S939", "S21", "Y182"],
        "site_sequence": [
            ("A" * 15) + "S" + ("A" * 15),
            ("A" * 15) + "S" + ("A" * 15),
            ("A" * 15) + "Y" + ("A" * 15),
        ],
        "display_id": ["TSC2;S939;", "GSK3A;S21;", "MAPK14;Y182;"],
        "organism": ["rat", "rat", "rat"],
        "protein_namespace": ["protein_id", "protein_id", "protein_id"],
        "protein_identifier": ["TSC2", "GSK3A", "MAPK14"],
        "protein_group_id": ["TSC2", "GSK3A", "MAPK14"],
        "localisation_confidence": [0.96, 0.93, 0.95],
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
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        activity_config=None,
    )
)

print(dataset.site_metadata.loc[:, ["site_key", "display_id", "site_sequence"]])
print(kinase_result.prediction_result.pred_mat.head())
```

`ReferencePreset.AUTO` uses bundled rat-only references in this release. Use an
explicit `ReferenceBundle` for human, mouse, or custom reference contexts.

## What to read next

- [Preparing a dataset](api/dataset-build-workflow.md)
- [Differential analysis](api/differential-analysis.md)
- [Kinase analysis](api/kinase.md)
- [Signalome analysis](api/signalome.md)
- [Scientific interpretation and limitations](scientific-interpretation.md)
