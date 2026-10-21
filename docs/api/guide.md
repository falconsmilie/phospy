# Public Python API

PhosPy has a small public Python API built around a consistent pattern:

1. describe the work with a request;
2. pass the request to a builder or workflow;
3. inspect the typed result.

Start with a workflow guide if you are learning PhosPy. Use this page when you
need to confirm where a public class is imported from or how the public API is
organised.

## Start With a Workflow Guide

Choose the guide that matches your task:

- [Build an Analysis-Ready Dataset](dataset-build-workflow.md)
- [Differential Analysis](differential-analysis.md)
- [Enrichment](enrichment.md)
- [Kinase Analysis](kinase.md)
- [Signalome Analysis](signalome.md)
- [SPS Discovery and RUV-Style Correction](sps-ruv.md)

The guides introduce requests, configuration, workflows, and results only when
they are needed. You do not need to learn the complete API before running an
analysis.

## Understand the Three Public Surfaces

PhosPy is currently beta software for Python 3.11 and 3.12. APIs may evolve
during beta. ADR-0054 defines three exact aggregate surfaces.

### `phospy`

Use the package root for the dataset builder and the Kinase and Signalome
workflow conveniences. It deliberately contains only a small set of common
entry points.

### `phospy.api`

`phospy.api` is the stable user-facing route during beta for ordinary requests,
configurations, results, references, enums, errors, and workflow classes.

### `phospy.advanced`

`phospy.advanced` is a supported advanced route for documented specialist
policies, diagnostics, reference helpers, and publishers. It is expected to
evolve more readily than `phospy.api`.

Internal modules are unsupported import targets for external callers.

Python may allow imports from implementation modules, but that does not make
them public. There are no supported forwarding submodules or runtime
deprecation aliases.

[ADR-0054](../adr/adr_0054_public_api_reset_and_compatibility_removal.md)
contains the formal public-surface decision.

## Match a Task to Its Public Objects

| Task | Request | Builder or Workflow | Result |
| --- | --- | --- | --- |
| Prepare data | `DatasetBuildRequest` | `AnalysisReadyDatasetBuilder` | `AnalysisReadyPhosphoDataset` |
| Compare conditions | `DifferentialAnalysisRequest` | `DifferentialAnalysisWorkflow` | `DifferentialAnalysisResult` |
| Test set overlap | `EnrichmentWorkflowRequest` | `EnrichmentWorkflow` | `EnrichmentWorkflowResult` |
| Score kinase support | `KinaseWorkflowRequest` | `KinaseWorkflow` | `KinaseWorkflowResult` |
| Build Signalome summaries | `SignalomeWorkflowRequest` | `SignalomeWorkflow` | `SignalomeWorkflowResult` |

A request records what you want to run. The builder or workflow validates the
request against the supplied data when `run(...)` is called.

## Import Common Entry Points

The package root contains only the three curated conveniences:

```python
from phospy import (
    AnalysisReadyDatasetBuilder,
    KinaseWorkflow,
    SignalomeWorkflow,
)
```

`phospy.__all__` contains exactly these three names. Differential analysis,
enrichment, and all request and result contracts use `phospy.api`.

## Import Stable API Objects

Stable user-facing names live in `phospy.api`:

```python
from phospy.api import (
    DatasetBuildRequest,
    ExperimentalDesign,
    Contrast,
    SampleDesignRecord,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SignalomeWorkflowRequest,
    UnsupportedInputFormatError,
    WorkflowValidationError,
)
```

Enrichment is also available through `phospy.api`:

```python
from phospy.api import (
    EnrichmentConfig,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
)
```

The [Enrichment guide](enrichment.md) documents additional provenance objects
for identifier sets derived from quantitative results.

## Import Advanced Supported Objects

Use `phospy.advanced` when a workflow guide or reference page asks for a
specialist object:

```python
from phospy.advanced import (
    SignalomeConfig,
    publish_dataset,
    publish_kinase_workflow,
    publish_signalome_workflow,
)
```

Examples of advanced supported features include:

- SPS discovery and RUV-style correction;
- specialist differential and Signalome configuration;
- local reference-resource helpers;
- diagnostic models; and
- publisher helpers for writing result tables.

Advanced support does not make private validators, workflow executors,
underscored helpers, or nearby implementation modules public APIs.

## Build Datasets Through the Builder

Use `AnalysisReadyDatasetBuilder` to turn input tables into a validated,
analysis-ready dataset:

```python
from phospy import AnalysisReadyDatasetBuilder

builder = AnalysisReadyDatasetBuilder()
dataset = builder.run(dataset_request)
```

The direct `AnalysisReadyPhosphoDataset(...)` constructor raises immediately.
The advanced `AnalysisReadyPhosphoDataset.from_trusted_tables(...)` route is
for callers who already own fully prepared tables and all required typed
evidence. It is not a shortcut around validation.

Most users should follow [Build an Analysis-Ready Dataset](dataset-build-workflow.md).

## Follow the Request-and-Run Pattern

This shortened example assumes that `phospho` and `site_metadata` have already
been created. The [Quickstart](../quickstart.md) defines the complete input
tables and explains each step.

First, establish the log2 transformation and localisation requirement:

```python
from phospy.advanced import DatasetIntensityTransformConfig
from phospy.api import DatasetLocalisationConfig, DatasetPreprocessingConfig

intensity_transform = DatasetIntensityTransformConfig(policy="log2")

localisation = DatasetLocalisationConfig(
    confidence_column="localisation_confidence",
    min_confidence=0.75,
)

preprocessing = DatasetPreprocessingConfig(
    intensity_transform=intensity_transform,
    localisation=localisation,
)
```

The value `0.75` is illustrative. Choose a localisation threshold that matches
the upstream evidence and study policy.

Next, describe and build the dataset:

```python
from phospy import AnalysisReadyDatasetBuilder
from phospy.api import DatasetBuildRequest, Organism

dataset_request = DatasetBuildRequest(
    phospho=phospho,
    site_metadata=site_metadata,
    organism=Organism.RAT,
    preprocessing_config=preprocessing,
)

dataset = AnalysisReadyDatasetBuilder().run(dataset_request)
```

Then describe the experiment and its contrast:

```python
from phospy.api import Contrast, ExperimentalDesign, SampleDesignRecord

design = ExperimentalDesign(
    samples=(
        SampleDesignRecord("control_1", "control", "control_r1"),
        SampleDesignRecord("control_2", "control", "control_r2"),
        SampleDesignRecord("treated_1", "treated", "treated_r1"),
        SampleDesignRecord("treated_2", "treated", "treated_r2"),
    )
)

contrast = Contrast(
    "treated_vs_control",
    "treated",
    "control",
)
```

Finally, construct and run the differential request:

```python
from phospy.api import DifferentialAnalysisRequest, DifferentialAnalysisWorkflow

differential_request = DifferentialAnalysisRequest(
    dataset=dataset,
    design=design,
    contrasts=(contrast,),
)

result = DifferentialAnalysisWorkflow().run(differential_request)
```

## Work With Results

Result models are typed containers. Documented table accessors return
independent pandas snapshots, so changing a returned DataFrame does not change
the stored result.

Use the workflow-specific accessors:

```python
differential_table = result.table_for("treated_vs_control")
kinase_scores = kinase_result.scoring_result.authoritative_scores
predictions = kinase_result.prediction_result.pred_mat
activity = kinase_result.activity_result.activity_matrix
enrichment_table = enrichment_result.table
assignments = signalome_result.module_assignments.table
```

Results can also contain diagnostics, caveats, provenance, and attrition
information. Review these alongside the main tables before interpreting a
missing or filtered row as biological absence.

Use [Output Bundles](../output_bundles.md) when you need to write supported
workflow tables to disk.

## Handle Validation Errors

Configuration objects may reject invalid local values when they are created.
Scientific and cross-table checks usually run at the builder or workflow
boundary.

Common public exception families are available from `phospy.api`:

```python
from phospy import KinaseWorkflow
from phospy.api import PhosPyValidationError, WorkflowValidationError

try:
    result = KinaseWorkflow().run(kinase_request)
except WorkflowValidationError as error:
    print(f"Check the workflow request: {error}")
except PhosPyValidationError as error:
    print(f"Check the input data: {error}")
```

See the [Validation Guide](../validation.md) for the main validation boundaries
and common remedies.

## Understand the API Scope

PhosPy's API is Python-native. It does not provide HTTP endpoints, REST routes,
or a hosted service.

The stable and advanced surfaces define supported imports. They do not promise
that every object will remain unchanged throughout beta. Review the
[release notes](../release-notes.md) when upgrading.
