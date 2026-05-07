# API Guide

Welcome to the PhosPy API guide. This page gives you the public import contract,
the supported workflow shape, and links to the workflow-specific API pages in
`docs/api/`.

PhosPy does not expose HTTP endpoints. The supported programmatic interface is
the Python API.

## Workflow Pages

The workflow documentation is split into dedicated pages:

| Workflow | Page | Description |
| --- | --- | --- |
| Dataset | [Dataset Workflow](dataset-build-workflow.md) | Start here when you have phosphosite intensity data and want a strict `AnalysisReadyPhosphoDataset` for kinase and signalome analysis.|
| Kinase | [Kinase Workflow](kinase-workflow.md) | `KinaseWorkflow` resolves references, scores kinase-substrate evidence, predicts candidate kinase regulation, and can optionally compute kinase activity tables. |
| Signalome | [Signalome Workflow](signalome-workflow.md) | `SignalomeWorkflow` interprets kinase score profiles into module assignments, signalome module summaries, kinase networks, and protein-site context tables |

The usual order is:

```python
dataset = AnalysisReadyDatasetBuilder().run(dataset_request)
kinase_result = KinaseWorkflow().run(kinase_request)
signalome_result = SignalomeWorkflow().run(signalome_request)
```

## Import Contract

Use top level `phospy` for the main entrypoints:

```python
from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
    SignalomeWorkflow,
)
```

Use `phospy.api` for requests, configs, results, enums, references, and public
exceptions:

```python
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflowRequest,
    UnsupportedInputFormatError,
    WorkflowValidationError,
)
```

All public executors use `run(request)`.

## Public Workflow Shape

1. `DatasetBuildRequest` -> `AnalysisReadyDatasetBuilder.run(...)` -> `AnalysisReadyPhosphoDataset`
2. `KinaseWorkflowRequest` -> `KinaseWorkflow.run(...)` -> `KinaseWorkflowResult`
3. `SignalomeWorkflowRequest` -> `SignalomeWorkflow.run(...)` -> `SignalomeWorkflowResult`

The beginner lane is rat first because bundled runtime references in the current
release are rat only. Human and mouse workflows need an explicit
`ReferenceBundle`.

The dataset that leaves the builder must be missing-value-free. This strict
boundary keeps kinase scoring, prediction, and signalome interpretation easier
to audit.

For concise scientist facing assumptions and interpretation notes, see
[Workflow Contracts](../workflow_contracts.md).

## Result Models

### Analysis-Ready Phospho Dataset

Important fields on `AnalysisReadyPhosphoDataset` include:

- `phospho`
- `site_metadata`
- `sample_metadata`
- `total`
- `comparisons`
- `organism`
- `intensity_scale_state`
- `processing_state`
- `preprocessing_report`
- `provenance`

Read `intensity_scale_state.label` together with
`intensity_scale_state.quantity`. For example, `log2` describes numeric scale,
while `phospho_total_log_ratio` describes what the values mean scientifically.

Use `dataset.to_dataframe()` for a safe phospho snapshot:

```python
phospho_snapshot = dataset.to_dataframe()
```

### Kinase Workflow Result

Important fields on `KinaseWorkflowResult` include:

- `dataset`
- `references`
- `scoring_result`
- `prediction_result`
- `activity_result`
- `provenance`

Common tables include `profile_scores`, `rank_weighted_fusion_scores`,
`pred_mat`, and activity tables when activity is enabled. Use
`activity_result.activity_scores` as the primary activity-score matrix, and
`weighted_activity` as a compatibility alias.

Use export helpers on scoring, prediction, and activity result objects for safe
snapshot copies:

```python
profile_scores = kinase_result.scoring_result.to_dataframe()
prediction_matrix = kinase_result.prediction_result.to_dataframe()
```

`kinase_result.provenance.scientific_policies` lists the active scientific
scoring policies with stable IDs, assumptions, parameters, and output scale
notes for auditability.

### Signalome Workflow Result

Important fields on `SignalomeWorkflowResult` include:

- `dataset`
- `kinase_result`
- `module_assignments`
- `signalome_modules`
- `kinase_network`
- `module_selection_diagnostics`
- `score_preconditioning_diagnostics`
- `expanded_signalome`
- `site_membership`
- `protein_site_context`
- `provenance`

Undefined kinase correlations are preserved as missing values. A correlation of
`0.0` means a correlation was estimated and is near zero.

Use public export helpers for safe sidecar snapshots:

```python
expanded_signalome = signalome_result.to_dataframe()
site_membership = signalome_result.site_membership_dataframe()
protein_context = signalome_result.protein_site_context_dataframe()
```

## References

`Organism` values are:

```python
Organism.HUMAN
Organism.MOUSE
Organism.RAT
```

Their string values are:

```python
"human"
"mouse"
"rat"
```

`ReferencePreset` values are:

```python
ReferencePreset.AUTO
ReferencePreset.HUMAN
ReferencePreset.MOUSE
ReferencePreset.RAT
```

Their string values are:

```python
"auto"
"human"
"mouse"
"rat"
```

Enum presence does not mean bundled runtime data exists for every organism in
this release. Use `ReferenceBundle` for custom references. It requires:

- `organism`
- `kinase_substrate_map` with `kinase` and `substrate_site`
- `site_sequences` indexed by site ID with `site_sequence`

Example:

```python
from phospy.api import Organism, ReferenceBundle

references = ReferenceBundle(
    organism=Organism.HUMAN,
    kinase_substrate_map=kinase_substrate_map,
    site_sequences=site_sequences,
)
```

## Public Exceptions

All user facing exception types are available from `phospy.api`. Common ones are:

- `PhosPyInputError`
- `UnsupportedInputFormatError`
- `PhosPyValidationError`
- `ReferenceResolutionError`
- `ReferenceCompatibilityError`
- `WorkflowValidationError`
- `WorkflowBoundaryError`
- `SignalomeScaleError`

Example:

```python
from phospy.api import PhosPyValidationError, WorkflowValidationError

try:
    kinase_result = KinaseWorkflow().run(kinase_request)
except WorkflowValidationError as error:
    print(f"Please check the workflow configuration: {error}")
except PhosPyValidationError as error:
    print(f"Please check the input tables: {error}")
```

## Small Working Example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.api import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
)

phospho = pd.DataFrame(
    {
        "sample_a": [1.00, 0.70],
        "sample_b": [1.10, 0.80],
        "sample_c": [0.95, 0.75],
    },
    index=["TSC2;S939;", "GSK3B;S9;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["TSC2", "GSK3B"],
        "site": ["S939", "S9"],
        "site_sequence": [
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
            "_______MSGRPRTTSFAESCKPVQQPSAFG",
        ],
        "protein_id": ["TSC2", "GSK3B"],
    },
    index=phospho.index.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
    )
)
kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        activity_config=None,
    )
)
print(kinase_result.prediction_result.pred_mat)
```
