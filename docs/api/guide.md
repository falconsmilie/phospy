# API Guide

Welcome to the PhosPy API guide. This page shows the supported Python import
contract, the workflow map, and the main boundaries to keep in mind when writing
analysis code.

## API Endpoint Status

PhosPy does not expose HTTP endpoints, REST routes, or a web service. The supported user interface is the Python API. In this guide, "API" means importable Python classes, request objects, configuration objects, workflow entrypoints, and result objects.

This guide describes executable interfaces, not broad PhosR-equivalence claims.
Scope categories and parity/open-gap status live in
[Scientific Coverage](../scientific-coverage.md), with fixture comparison
details in [Parity](../parity.md).

## Workflow Pages

Each public workflow has its own API page:

| Area | Page | Public API |
| --- | --- | --- |
| Dataset builder | [Dataset Build API](dataset-build-workflow.md) | `AnalysisReadyDatasetBuilder.run(DatasetBuildRequest)` returns `AnalysisReadyPhosphoDataset`. |
| Differential analysis | [Differential Analysis Workflow](differential-analysis.md) | `DifferentialAnalysisWorkflow.run(DifferentialAnalysisRequest)` returns `DifferentialAnalysisResult`. |
| Enrichment | [Enrichment Workflow](enrichment.md) | `EnrichmentWorkflow.run(EnrichmentWorkflowRequest)` returns `EnrichmentWorkflowResult`. |
| Kinase | [Kinase Workflow](kinase.md) | `KinaseWorkflow.run(KinaseWorkflowRequest)` returns `KinaseWorkflowResult`. |
| Signalome | [Signalome Workflow](signalome.md) | `SignalomeWorkflow.run(SignalomeWorkflowRequest)` returns `SignalomeWorkflowResult`. |
| Importers | [Phosphosite Importers](../importers.md) | Importers produce builder input candidates; they do not bypass dataset validation. |

The usual workflow shape is:

```python
dataset = AnalysisReadyDatasetBuilder().run(dataset_request)
differential_result = DifferentialAnalysisWorkflow().run(differential_request)
kinase_result = KinaseWorkflow().run(kinase_request)
signalome_result = SignalomeWorkflow().run(signalome_request)
```

Enrichment is independent of the dataset/kinase/signalome chain:

```python
enrichment_result = EnrichmentWorkflow().run(enrichment_request)
```

## Import Contract

Use top-level `phospy` for the main convenience entrypoints:

```python
from phospy import (
    AnalysisReadyDatasetBuilder,
    DifferentialAnalysisWorkflow,
    KinaseWorkflow,
    SignalomeWorkflow,
)
```

`AnalysisReadyPhosphoDataset` remains exported for advanced/trusted callers who
already own fully prepared `site_key`-indexed tables and processing-state
provenance. Ordinary dataset construction should use
`AnalysisReadyDatasetBuilder().run(DatasetBuildRequest(...))`.

The `phospy.api.datasets` submodule is stable-only and exports
`AnalysisReadyPhosphoDataset`. It does not export processing-state internals,
missing-data or normalisation state records, batch-correction report classes, or
other nested diagnostic model classes. Inspect diagnostics through returned
objects such as `dataset.preprocessing_report`,
`dataset.preprocessing_report.batch_correction`, and workflow result properties
when those reports are present; importing the underlying diagnostic classes is
reserved for PhosPy internals through their owning implementation modules.

Use `phospy.api` for stable request, workflow, primary result, reference, enum,
and common exception contracts. The aggregate facade is intentionally smaller
than the implementation modules:

```python
from phospy.api import (
    DatasetBuildRequest,
    ExperimentalDesign,
    Contrast,
    SampleDesignRecord,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
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

`EnrichmentWorkflow` is a supported public workflow from `phospy.api`:

```python
from phospy.api import (
    EnrichmentConfig,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
)
```

### API Stability Tiers

Stable public API is the default user-facing surface. It includes the dataset
builder, `AnalysisReadyPhosphoDataset`, core workflow request objects, workflow
classes with `run(...)`, primary workflow result objects, reference-bundle
entrypoints, example-level enums, and common exception types.

Advanced supported API is public but should be imported deliberately. It
includes selected specialized configuration objects, control-site policy
helpers, local Kinase Library-style resource loaders, and explicit result-table
inspection helpers such as `filter_differential_results` and
`rank_differential_results`.

Internal / experimental API is not exported through `phospy.api`. This includes
validators, private result assemblers, internal scoring helpers, private
provenance serialization functions, low-level workflow interpreters and
executors, reference manifest validation internals, processing-state internals,
nested diagnostic records, and compatibility constants. Submodule wildcard
surfaces under `phospy.api` follow the same inventory: names classified as
internal / experimental are not exported from `__all__`.

See [ADR-0031](../adr/adr_0031_public_api_stability_tiers.md) for the current
inventory and promotion policy.

### Public and Semi-Public Routes

The supported public API remains:

- `phospy` for the small top-level workflow convenience surface.
- `phospy.api` for public request, config, result, enum, reference, workflow,
  and exception names listed in `phospy.api.__all__`. This list contains stable
  and explicitly advanced supported names only.

Selected `phospy.science.*` routes are semi-public compatibility routes for
advanced extension, parity, and backend-contract use. They are not promoted to
`phospy.api`, and neighbouring private helpers are not public.

| Status | Supported route | Supported names |
| --- | --- | --- |
| Public | `phospy` | Top-level convenience entrypoints listed above. |
| Public | `phospy.api` | Public names listed in `phospy.api.__all__`. |
| Public | `phospy.api.datasets` | `AnalysisReadyPhosphoDataset` only. |
| Semi-public | `phospy.science.datasets.preprocessing.stage_registry` | `PreprocessingStageMetadata` and exported registry helpers. |
| Semi-public | `phospy.science.signalomes.clustering.protocol` | `ClusterTreeEngine`, `SignalomeClusteringEngine`. |
| Semi-public | `phospy.science.signalomes.clustering.exact_python` | Exact-Python clustering compatibility facade names exported in `__all__`. |
| Semi-public | `phospy.science.prediction.scoring` | `fuse_profile_and_motif_scores_by_rank_weight` for parity and advanced scoring checks. |

Unsupported import routes include underscored helpers and root/API imports for
semi-public science helpers. For example,
`from phospy.api import PreprocessingStageMetadata` is not supported.

See [ADR-0028](../adr/adr_0028_semi_public_science_import_policy.md) for the
semi-public route policy.

## Request Validation Boundary

Public request dataclasses are lightweight command payloads. Constructing a
request records user intent, but it does not mean the request is scientifically
valid.

Validation happens when the relevant builder or workflow is run:

- `AnalysisReadyDatasetBuilder.run(request)` validates inputs, preprocessing
  compatibility, site-resolution state, and the strict analysis-ready dataset
  boundary.
- `DifferentialAnalysisWorkflow.run(request)` validates the dataset, explicit
  design, contrasts, replicate requirements, and differential config.
- `EnrichmentWorkflow.run(request)` validates explicit identifier semantics,
  background universe, and supplied set collections before ORA execution.
- `KinaseWorkflow.run(request)` validates the dataset, references, workflow
  configs, localisation requirements, site sequences, and reference projection.
- `SignalomeWorkflow.run(request)` validates the upstream kinase result, matrix
  alignment, site identity, protein grouping metadata, and signalome config.

Config objects may reject invalid local policy values at construction time
because those invariants belong to the config itself.

## Dataset Boundary

`AnalysisReadyPhosphoDataset` is the strict analysis-ready dataset boundary.
Downstream workflows expect it to be complete, auditable, and keyed by
`site_key`.

Important identity rules:

- `dataset.phospho.index` is `site_key`.
- `dataset.site_metadata.index` is `site_key`.
- `dataset.site_metadata["site_key"]` matches the index.
- `dataset.site_metadata["display_id"]` is a human-readable label.
- Duplicate `display_id` values remain valid when the corresponding `site_key`
  values differ.
- Duplicate rows that resolve to the same `site_key` are a scientific ambiguity
  and fail by default unless an explicit non-error duplicate-site preprocessing
  policy is chosen.

For preprocessing options, including localisation, missing data, total-protein
correction, protein-aware preparation, batch residualisation, native
SPS/RUV-style correction, and RUV-readiness reporting, see
[Dataset Build API](dataset-build-workflow.md).

`ruv_readiness` is report-only. It does not select SPS controls, run
correction, or imply PhosR-equivalent batch correction. Native SPS/RUV-style
correction is available only through explicit `SpsRuvBatchCorrectionConfig`.

## Result Snapshots

Result models are typed containers. Public helpers such as `to_dataframe()`,
`*_dataframe()`, `table`, `result_table`, and `to_payload()` return defensive
in-memory snapshots for inspection or handoff.

They are not exporters, plotting helpers, report generators, or places to run
additional scientific post-processing.

## References

`ReferencePreset.AUTO` is intended for the bundled rat beginner lane in this
release. Human and mouse workflows should pass an explicit `ReferenceBundle`.

Use `ReferenceBundleBuilder` when building references from local source files so
provenance and validation are recorded consistently. The builder reads local
files only; it does not scrape web resources or invent missing sequence windows.

## Public Exceptions

Common user-facing exception types are available from `phospy.api`:

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

This tiny example builds a rat analysis-ready dataset and runs the kinase
workflow with activity disabled. The numbers are synthetic and demonstrate API
wiring only.

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
            "ATMSGRPRTTSFAESSSPVQQPSAFGQAAAL",
        ],
        "display_id": ["TSC2;S939;", "GSK3B;S9;"],
        "organism": ["rat", "rat"],
        "protein_namespace": ["protein_id", "protein_id"],
        "protein_identifier": ["TSC2", "GSK3B"],
        "protein_id": ["TSC2", "GSK3B"],
        "localisation_confidence": [0.95, 0.92],
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
        activity_config=None,
    )
)

print(kinase_result.prediction_result.pred_mat)
```
