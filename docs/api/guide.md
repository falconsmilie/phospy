# Public Python API

PhosPy uses a small public Python API. Start with the workflow guides for
end-to-end examples; use this page when you need to confirm where a public class
should be imported from.

## Workflow Map

| Task | Request | Workflow | Result |
| --- | --- | --- | --- |
| Prepare data | `DatasetBuildRequest` | `AnalysisReadyDatasetBuilder` | `AnalysisReadyPhosphoDataset` |
| Test condition contrasts | `DifferentialAnalysisRequest` | `DifferentialAnalysisWorkflow` | `DifferentialAnalysisResult` |
| Run offline enrichment | `EnrichmentWorkflowRequest` | `EnrichmentWorkflow` | `EnrichmentWorkflowResult` |
| Score kinase support | `KinaseWorkflowRequest` | `KinaseWorkflow` | `KinaseWorkflowResult` |
| Build signalome summaries | `SignalomeWorkflowRequest` | `SignalomeWorkflow` | `SignalomeWorkflowResult` |

The complete request and response contract for each workflow lives on its own
page:

- [Dataset Preparation](dataset-build-workflow.md)
- [Differential Analysis](differential-analysis.md)
- [Enrichment](enrichment.md)
- [Kinase Analysis](kinase.md)
- [Signalome Analysis](signalome.md)

## Import From `phospy`

Use the package root for the dataset builder and the main workflow entry points:

```python
from phospy import (
    AnalysisReadyDatasetBuilder,
    DifferentialAnalysisWorkflow,
    KinaseWorkflow,
    SignalomeWorkflow,
)
```

## Import From `phospy.api`

Stable public API names live in `phospy.api`.

Use `phospy.api` for stable requests, common configuration objects, references,
results, enums, and user-facing exceptions:

```python
from phospy.api import (
    DatasetBuildRequest,
    ExperimentalDesign,
    Contrast,
    SampleDesignRecord,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SignalomeWorkflowRequest,
    UnsupportedInputFormatError,
    WorkflowValidationError,
)
```

Enrichment is also public through `phospy.api`:

```python
from phospy.api import (
    EnrichmentConfig,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
)
```

## Import From `phospy.advanced`

Advanced supported API names live in `phospy.advanced`.

Use `phospy.advanced` only when a workflow guide asks for a specialized policy,
configuration object, or table publisher:

```python
from phospy.advanced import (
    SignalomeConfig,
    publish_dataset,
    publish_kinase_workflow,
    publish_signalome_workflow,
)
```

The stable and advanced surfaces are intentional. Do not build user code around
private validators, internal workflow executors, underscored helpers, or nearby
implementation modules simply because Python can import them.

Internal / experimental API names are not supported import targets.

## Build Datasets Through the Builder

The supported construction path is:

```python
dataset = AnalysisReadyDatasetBuilder().run(dataset_request)
```

The direct `AnalysisReadyPhosphoDataset(...)` constructor raises immediately.
The advanced/trusted route,
`AnalysisReadyPhosphoDataset.from_trusted_tables(...)`, is for callers who
already own fully prepared, `site_key`-indexed tables and the required typed
evidence. It is not a shortcut around validation. Most users should use
`AnalysisReadyDatasetBuilder`.

## A Small Public-API Pattern

The example below shows the request-and-run pattern for differential analysis.
The [quickstart](../quickstart.md) provides complete input tables.

```python
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
            Contrast("treated_vs_control", "treated", "control"),
        ),
    )
)
```

Constructing a request records your intent. The builder or workflow validates
that request when `run(...)` is called.

## Work With Results

Result models are typed containers. Their documented table helpers return
independent pandas snapshots, so changing a returned DataFrame does not change
the result object.

Use the workflow-specific helpers whenever possible:

```python
differential_table = result.table_for("treated_vs_control")
kinase_scores = kinase_result.scoring_result.authoritative_scores
predictions = kinase_result.prediction_result.pred_mat
enrichment_table = enrichment_result.table
assignments = signalome_result.module_assignments.table
```

Result objects also carry diagnostics, provenance, and caveats. Review these
before interpreting filtered or missing rows as biological absence.

## Handle Validation Errors

Common user-facing exceptions are available from `phospy.api`:

```python
from phospy.api import PhosPyValidationError, WorkflowValidationError

try:
    result = KinaseWorkflow().run(kinase_request)
except WorkflowValidationError as error:
    print(f"Check the workflow request: {error}")
except PhosPyValidationError as error:
    print(f"Check the input data: {error}")
```

Configuration objects may reject invalid local values when they are created.
Scientific and cross-table checks usually run at the builder or workflow
boundary.

## API Scope

PhosPy's API is Python-native; it does not provide HTTP endpoints, REST routes,
or a hosted service. Detailed implementation boundaries and compatibility
policy are documented in [ADR-0031](../adr/adr_0031_public_api_stability_tiers.md).
