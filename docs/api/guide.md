# Public Python API

PhosPy uses a small public Python API. Start with the workflow guides for
end-to-end examples; use this page when you need to confirm where a public class
should be imported from.

## Stability Tiers During Beta

PhosPy is currently beta software for Python 3.11 and 3.12. APIs may evolve
during beta, but compatibility expectations follow the governed stability tier
for the import route you use:

`phospy.api` is the stable user-facing route for ordinary workflow code.
`phospy.advanced` is a supported advanced route for documented specialist
configuration, diagnostics, references, and publishing helpers.
Implementation modules are unsupported import targets for external callers.
Non-default specialist policy objects, including
`DifferentialAnalysisConfig` and paired-design policy constants such as
`PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION`, are supported from
`phospy.advanced`; they are intentionally outside the stable aggregate
`phospy.api` export set during beta.
The beta compatibility expectations are:

| Stability tier | Supported import route | Compatibility expectation |
| --- | --- | --- |
| Stable public API | Use `phospy.api` for ordinary workflow code. Selected root-package convenience imports, such as `from phospy import KinaseWorkflow`, alias stable facade names. | Contains the primary beta-user contracts. Compatibility and migration treatment follow the project's stable public API policy. "Stable" means policy-governed during beta, not frozen forever. |
| Advanced supported API | Use `phospy.advanced` when a guide asks for specialized configuration, diagnostic models, reference helpers, or publishing helpers. | These APIs are supported and documented, but may evolve more readily. Release notes or migration guidance accompany material changes where project policy requires it. |
| Internal / experimental API | Do not use direct imports from implementation modules as public contracts. | Python allowing an import does not make it public. Internal modules may change without compatibility treatment, and private validators must not be used as external APIs. |

[ADR-0031](../adr/adr_0031_public_api_stability_tiers.md) governs the policy
behind these tiers. This guide points to supported facades and workflow pages
instead of listing every exported symbol.

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

Use the package root only for curated convenience aliases to stable facade
names, such as the dataset builder and selected workflow entry points:

```python
from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
    KinaseWorkflow,
    SignalomeWorkflow,
)
```

These imports are conveniences for stable public API objects. Do not infer that
every symbol visible from the package root is stable. For request, result,
configuration, reference, enum, and exception contracts, prefer the explicit
`phospy.api` imports below.

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
    EnrichmentDerivedQuantitativeSetProvenance,
    EnrichmentDerivedSetMissingValueRule,
    EnrichmentDerivedSetSourceResultKind,
    EnrichmentDerivedSetThresholdDirection,
    EnrichmentDerivedSetValueMeaning,
    EnrichmentDerivedSetValueScale,
    EnrichmentConfig,
    EnrichmentIdentifierSetProvenance,
    EnrichmentIdentifierSetSourceType,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
    InputIntensityScaleEvidence,
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

SPS discovery is also supported from this facade. Request construction only
stores caller intent. The workflow boundary validates the reference matrices
as one multi-dataset input while keeping the validator itself private, then
derives a consensus control set from those explicit references:

```python
from phospy.advanced import (
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsDiscoveryWorkflow,
    SpsDiscoveryValidationError,
    SpsReferenceDataset,
)

request = SpsDiscoveryRequest(
    reference_datasets=(reference_a, reference_b),
    config=SpsDiscoveryConfig(top_n=100),
)

try:
    result = SpsDiscoveryWorkflow().run(request)
except SpsDiscoveryValidationError as exc:
    validation = exc.validation_result
```

Each `SpsReferenceDataset` supplies its own sample-to-condition mapping, so
reference datasets need not share samples or condition labels. Cross-dataset
phosphosite identity is the governed `site_key` index. The matrix should contain
condition-relative phosphorylation measurements, such as log-ratios. Discovery
averages biological replicates within each condition and scores a site in each
reference by its largest absolute condition mean. It ranks smaller change
magnitudes as more stable and combines reference ranks with a Fisher-style
consensus score. `result.selected_site_keys` contains the requested leading
controls when that many remain, and `result.control_site_set` feeds the existing
SPS/RUV control pathway directly.

Missing replicate cells are ignored only when at least one finite value remains
for that site and condition. A site with a completely missing condition does not
contribute a rank for that reference; infinite measurements reject the request.
The site must still contribute to at least
`minimum_datasets_per_site` independent references; discovery never relaxes
that threshold. If fewer than `top_n` eligible controls remain, all eligible
controls are returned and the requested and actual counts are recorded. If
fewer than `minimum_shared_sites` are rankable, execution fails with structured
`insufficient_rankable_overlap` validation. Equal change magnitudes receive an
equal within-reference rank and final consensus ties use ascending `site_key`.

Reference matrices are caller-supplied evidence; the target experiment is not
silently reused for control discovery. SPS discovery is PhosR-inspired and is
designed for later external parity validation, but it does not claim PhosR or
RUV-III parity.

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
        contrasts=(Contrast("treated_vs_control", "treated", "control"),),
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
