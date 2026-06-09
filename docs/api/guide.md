# API Guide

Welcome to the PhosPy API guide. This page gives you the public import contract,
the supported workflow shape, and links to the workflow-specific API pages in
`docs/api/`.

PhosPy does not expose HTTP endpoints. The supported programmatic interface is
the Python API.

This API guide describes executable interfaces, not global PhosR-equivalence
claims. Scope categories and parity/open-gap status are maintained in
[`docs/scientific-coverage.md`](../scientific-coverage.md).

## Workflow Pages

The workflow documentation is split into dedicated pages:

| Workflow | Page | Description |
| --- | --- | --- |
| Dataset | [Dataset Workflow](dataset-build-workflow.md) | Start here when you have phosphosite intensity data and want a strict `AnalysisReadyPhosphoDataset` for kinase and signalome analysis.|
| Differential | [Differential Workflow](differential-workflow.md) | `DifferentialAnalysisWorkflow` runs moderated differential analysis over an `AnalysisReadyPhosphoDataset` using explicit design and contrast definitions. |
| Kinase | [Kinase Workflow](kinase-workflow.md) | `KinaseWorkflow` resolves references, scores kinase-substrate evidence, predicts candidate kinase regulation, and can optionally compute kinase activity tables. |
| Signalome | [Signalome Workflow](signalome-workflow.md) | `SignalomeWorkflow` interprets kinase score profiles into module assignments, signalome module summaries, kinase networks, and protein-site context tables |

The usual order is:

```python
dataset = AnalysisReadyDatasetBuilder().run(dataset_request)
differential_result = DifferentialAnalysisWorkflow().run(differential_request)
kinase_result = KinaseWorkflow().run(kinase_request)
signalome_result = SignalomeWorkflow().run(signalome_request)
```

## Import Contract

Use top level `phospy` for the main entrypoints:

```python
from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
    KinaseWorkflow,
    SignalomeWorkflow,
)
```

Use `phospy.api` for requests, configs, results, enums, references, and public
exceptions:

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

All public executors use `run(request)`.

Most short snippets below show one concept at a time. For a copy/paste run,
use the full examples in [Quickstart](../quickstart.md) or each workflow page.

## Scientific Policy Module Ownership

Scientific policy records are owned by domain modules, not a root dumping-ground
module.

- Shared policy record models:
  `phospy.provenance.scientific_policy_models`
- Prediction scientific policies:
  `phospy.science.prediction.scientific_policies`
- Activity scientific policies:
  `phospy.science.activities.scientific_policies`
- Preprocessing scientific policies:
  `phospy.science.datasets.preprocessing.scientific_policies`
- Signalome workflow scientific policies:
  `phospy.workflows.signalome.scientific_policies`
- Signalome clustering scientific policies:
  `phospy.science.signalomes.clustering.scientific_policies`
- Differential aggregation scientific policies:
  `phospy.science.differential.aggregation.scientific_policies`

`phospy.scientific_policies` is intentionally not part of the import contract.

## Public Workflow Shape

1. `DatasetBuildRequest` -> `AnalysisReadyDatasetBuilder.run(...)` -> `AnalysisReadyPhosphoDataset`
2. `DifferentialAnalysisRequest` -> `DifferentialAnalysisWorkflow.run(...)` -> `DifferentialAnalysisResult`
3. `KinaseWorkflowRequest` -> `KinaseWorkflow.run(...)` -> `KinaseWorkflowResult`
4. `SignalomeWorkflowRequest` -> `SignalomeWorkflow.run(...)` -> `SignalomeWorkflowResult`

The beginner lane is rat first because bundled runtime references in the current
release are rat only. Human and mouse workflows need an explicit
`ReferenceBundle`.

The dataset that leaves the builder must be missing-value-free. This strict
boundary keeps kinase scoring, prediction, and signalome interpretation easier
to audit. At this boundary, `site_key` is the unique analysis-ready row
identity and `display_id` is the human-readable `GENE;SITE;` label. The public
dataset indexes are:

- `AnalysisReadyPhosphoDataset.phospho.index`: `site_key`
- `AnalysisReadyPhosphoDataset.site_metadata.index`: `site_key`
- `AnalysisReadyPhosphoDataset.site_metadata["site_key"]`: same values as the
  index
- `AnalysisReadyPhosphoDataset.site_metadata["display_id"]`: display label

`display_id` may repeat when distinct `site_key` values preserve distinct
protein context. Direct `AnalysisReadyPhosphoDataset` construction requires
encoded `site_key` indexes plus auditable protein context metadata
(`organism`, `protein_namespace`, `protein_identifier`, `gene_symbol`, `site`,
and `site_sequence`). It does not fall back to display-site identity. Builder
ingestion may accept legacy display-indexed input only when protein context is
sufficient to derive `site_key`. See
[ADR-0024](../adr/adr_0024_protein_scoped_phosphosite_row_identity.md).

Workflows operate on `site_key`. User-facing site-level outputs that materialize
row identity include both `site_key` and `display_id`. Differential result
tables are stricter public scientific outputs: direct
`DifferentialAnalysisResult` construction requires encoded `site_key` indexes
and non-empty `site_key`, `display_id`, `gene_symbol`, and `site` columns.
Workflow-created differential results preserve available protein context such as
`organism`, `protein_namespace`, `protein_identifier`, and `protein_id`.
Display-indexed or stat-only differential result tables are not valid public
inputs.

The lower-level differential statistical executor may produce an internal
stat-only computation payload for workflow assembly. The public API result is
only `DifferentialAnalysisResult`, after the workflow has attached dataset
identity metadata.

## Result Construction Contracts

Public-looking result classes do not all have the same construction contract:

| Result object | Direct construction contract | Identity guarantee |
| --- | --- | --- |
| `DifferentialAnalysisResult` | Strict user-constructible public result. Use direct construction only with complete public contrast tables. | Requires encoded `site_key` index, matching `site_key` column, non-empty `display_id`, `gene_symbol`, and `site`, coherent display/site metadata, and contrast tables aligned to residual-statistic indexes. |
| `KinaseScoringResult`, `KinasePredictionResult`, `KinaseActivityResult` | Directly constructible stage result tables with schema validation. | Their own public table schemas are validated. Cross-object workflow coherence is guaranteed only when produced by `KinaseWorkflow.run(...)`. |
| `KinaseWorkflowResult` | Workflow-owned container with intentionally minimal direct construction. | Direct construction does not revalidate nested object types, reference compatibility, dataset alignment, scoring, prediction, activity, eligibility, or provenance coherence. Use `KinaseWorkflow.run(...)` for scientifically coherent results. |
| `SignalomeWorkflowResult` | Workflow-owned result. Direct construction is supported for reconstruction/tests and validates owned public sidecar table contracts. | Site-level public sidecars that claim analysis-ready phosphosite rows must use encoded `site_key`, non-empty `display_id`, and align to `result.dataset`. Full module/network/scoring coherence is guaranteed only when produced by `SignalomeWorkflow.run(...)`. |

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
- `site_sequences` indexed by display site ID with `site_sequence`

Kinase references may use display IDs at the reference boundary. During workflow
interpretation, those display IDs are matched against dataset `display_id`
metadata and projected to internal `site_key` rows through an explicit mapping
layer. References remain display-ID keyed at the reference boundary and are not
converted into analysis-ready row identity.

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
    IntensityScaleKind,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
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
            # Fail fast on missing/low-confidence localisation so site-level
            # kinase interpretation does not rely on ambiguous phosphosite mapping.
            localisation=DatasetLocalisationConfig(
                mode="require_threshold",
                confidence_column="localisation_confidence",
                min_confidence=0.75,
            )
        ),
    )
)
print(
    dataset.site_metadata.loc[
        :,
        [
            "site_key",
            "display_id",
            "gene_symbol",
            "site",
            "protein_namespace",
            "protein_identifier",
            "protein_id",
        ],
    ]
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
