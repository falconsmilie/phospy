# API Guide

Welcome to the PhosPy API guide. This page gives you the public import contract,
the supported workflow shape, and links to the workflow-specific API pages in
`docs/api/`.

PhosPy does not expose HTTP endpoints. The supported programmatic interface is
the Python API.

This API guide describes executable interfaces, not global PhosR-equivalence
claims. Scope categories and parity/open-gap status are maintained in
[`docs/scientific-coverage.md`](../scientific-coverage.md).
Future coverage direction is governed by
[ADR-0025](../adr/adr_0025_competitive_phosphoproteomics_workflow_coverage.md)
and does not expand the current API support contract by itself.

## Workflow Pages

The workflow documentation is split into dedicated pages:

| Workflow | Page | Description |
| --- | --- | --- |
| Dataset | [Dataset Workflow](dataset-build-workflow.md) | Start here when you have phosphosite intensity data and want a strict `AnalysisReadyPhosphoDataset` for kinase and signalome analysis.|
| Importers | [Phosphosite Importers](../importers.md) | Translate upstream search-engine outputs into dataset-builder input candidates without bypassing builder validation. |
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

## Request Validation Boundary

Public request dataclasses are lightweight command payloads. Constructing
`DatasetBuildRequest`, `DifferentialAnalysisRequest`, `KinaseWorkflowRequest`,
or `SignalomeWorkflowRequest` does not mean the request is scientifically valid.

Scientific validation happens when the relevant builder or workflow is run:

- `AnalysisReadyDatasetBuilder.run(request)` validates dataset-build request
  fields, input sources, preprocessing compatibility, and site-resolution state
  before building a dataset.
- `DifferentialAnalysisWorkflow.run(request)` validates the dataset, explicit
  design, contrasts, replicate requirements, and differential config before
  statistical execution.
- `KinaseWorkflow.run(request)` validates the dataset, references,
  workflow configs, localisation requirements, and reference-projection policy
  before kinase interpretation and scoring.
- `SignalomeWorkflow.run(request)` validates the upstream kinase result,
  score/prediction matrices, site identity, protein grouping metadata, and
  signalome config before signalome interpretation.

Config objects may be stricter than request objects. For example, config
dataclasses can reject invalid local policy values at construction time because
those invariants belong to the config itself. Request dataclasses should not be
treated as mini-workflow validators.

## Enrichment Contract Boundary

`EnrichmentWorkflowRequest`, `EnrichmentConfig`, and
`EnrichmentWorkflowResult` provide a typed foundation for future native
enrichment support. They are contracts only; PhosPy does not yet expose an
`EnrichmentWorkflow` executor and does not calculate enrichment statistics from
these objects.

Enrichment requests must provide exactly one identifier source:
`input_table` or `selected_identifiers`. They also require an
`identifier_column`, an explicit `identifier_kind`, a homogeneous
`EnrichmentSetCollection`, `GeneSetCollection`, or `PtmSetCollection`, an
explicit non-empty `background_universe`, and an `EnrichmentConfig`. Supported
identifier kinds are `gene_symbol`, `protein_id`, `site_key`, `display_id`,
and `phosphosite`. Gene-set collections are separate from PTM-set collections
so gene-level and phosphosite-level semantics are not collapsed into a generic
string flag.

The initial config supports `method="over_representation"` and
`multiple_testing_correction` values `"none"` or `"benjamini_hochberg"`.
Background universes are never inferred by these contracts, and no online
resources are loaded.

```python
from phospy.api import (
    EnrichmentConfig,
    EnrichmentSet,
    EnrichmentSetCollection,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
    PtmSetCollection,
)
```

Collections can be created directly from set records. Duplicate identifiers
inside one set are de-duplicated in first-seen order, and a collection must use
one explicit `identifier_kind` throughout.

```python
collection = EnrichmentSetCollection(
    sets=(
        EnrichmentSet(
            set_id="MAPK_PATHWAY",
            name="MAPK pathway",
            identifiers=("AKT1", "MAPK1", "AKT1"),
            identifier_kind="gene_symbol",
            source_name="local curated sets",
            source_version="2026.06",
            description="Example offline gene set",
        ),
    )
)

assert collection.members_by_set_id["MAPK_PATHWAY"] == ("AKT1", "MAPK1")
```

Local enrichment-set readers live under `phospy.io.readers`:

```python
from phospy.io.readers import (
    read_enrichment_sets_gmt,
    read_enrichment_sets_table,
)

gmt_collection = read_enrichment_sets_gmt(
    "gene_sets.gmt",
    identifier_kind="gene_symbol",
    source_name="local GMT",
)

table_collection = read_enrichment_sets_table("sets.csv")
```

GMT-like files are interpreted as
`set_id<TAB>description<TAB>identifier...` and require the caller to pass
`identifier_kind` because GMT does not carry identifier semantics. CSV/TSV
tables require `set_id`, `name`, and `identifier` columns, may include
`identifier_kind`, `source_name`, `source_version`, and `description`, and also
require an explicit `identifier_kind` argument when the file has no
`identifier_kind` column. These readers only parse local files. They do not
fetch GO, KEGG, Reactome, PTMsigDB, Enrichr, gseapy, clusterProfiler, or any
online database.

## Differential Design Declarations

`ExperimentalDesign` can explicitly declare fixed-effect covariates
without inferring anything from `dataset.sample_metadata`:

```python
from phospy.api import (
    BatchCovariate,
    CategoricalCovariate,
    ContinuousCovariate,
    ExperimentalDesign,
    SampleDesignRecord,
)

design = ExperimentalDesign(
    samples=(
        SampleDesignRecord(
            sample_id="control_rep1",
            condition="control",
            batch="run_1",
            covariates={"sex": "F", "dose": 0.0},
        ),
        SampleDesignRecord(
            sample_id="control_rep2",
            condition="control",
            batch="run_2",
            covariates={"sex": "M", "dose": 1.0},
        ),
        SampleDesignRecord(
            sample_id="control_rep3",
            condition="control",
            batch="run_1",
            covariates={"sex": "M", "dose": 2.0},
        ),
        SampleDesignRecord(
            sample_id="treated_rep1",
            condition="treated",
            batch="run_2",
            covariates={"sex": "F", "dose": 0.5},
        ),
        SampleDesignRecord(
            sample_id="treated_rep2",
            condition="treated",
            batch="run_1",
            covariates={"sex": "M", "dose": 1.5},
        ),
        SampleDesignRecord(
            sample_id="treated_rep3",
            condition="treated",
            batch="run_2",
            covariates={"sex": "F", "dose": 2.5},
        ),
    ),
    fixed_effects=(
        BatchCovariate(),
        CategoricalCovariate("sex"),
        ContinuousCovariate("dose"),
    ),
)
```

Each declaration records the covariate `name`, `kind`, whether it is
`required`, and whether it is intended to `include_in_model`. Modelled fixed
effects are validated before differential interpretation: missing covariates,
invalid levels, non-finite continuous values, rank-deficient designs, and
non-estimable contrasts are rejected. Result provenance records the resolved
design formula or description, condition columns, fixed-effect covariate
columns and kinds, contrast vectors, and validation status. Fixed batch terms
are model covariates; they are not batch correction.

Supported fixed-effect covariates are ordinary fixed terms in the differential
linear model. This is not ComBat, RUV, `removeBatchEffect`,
`duplicateCorrelation`, or mixed-effects modelling.

Paired or blocked designs use a single public sample metadata name:
`SampleDesignRecord.block_id`. `block_id` must be supplied explicitly by the
caller; PhosPy does not infer it from sample names, column order, or
`dataset.sample_metadata`. The differential config exposes
`paired_design_policy`, which defaults to `"reject"`. Setting
`paired_design_policy="fixed_block"` validates an explicit fixed-effect block
design and adds block terms to the differential design matrix. Each sample must
have `block_id`; each block must contain at least two samples and both sides of
every requested contrast; the resolved design must be full rank and contrasts
must be estimable. The block terms are fixed effects. This is not limma
`duplicateCorrelation`, not mixed-effects modelling, and not random subject
modelling; no mixed effects are fitted. Incomplete or partially covered blocks
are rejected before execution; PhosPy does not drop them to continue. Simple
unpaired workflows are unchanged unless the caller explicitly opts into
`paired_design_policy="fixed_block"`.

The lower-level design-matrix builder can represent included fixed effects for
design-domain inspection and validation. Categorical and batch covariates are
dummy-encoded with deterministic level order. Continuous covariates are emitted
as one raw numeric column named by the covariate; sample order is preserved,
values must be numeric and finite, string values such as `"2.5"` or
`"unknown"` are rejected instead of parsed or treated as missing, and no
centering or scaling is applied.

## Importer Boundary

`PhosphositeImportRequest` and `PhosphositeImportResult` support upstream table
translation before dataset building. Importers produce:

- `phospho_matrix_candidate`
- `site_metadata_candidate`
- optional `peptide_evidence`
- explicit `sample_column_mapping`
- `localisation_confidence_column`
- `warnings` and `diagnostics`

Importer output still feeds `AnalysisReadyDatasetBuilder`:

```python
from phospy import AnalysisReadyDatasetBuilder
from phospy.api import Organism, PhosphositeImportRequest
from phospy.io.readers import MappedPhosphositeTableImporter

import_result = MappedPhosphositeTableImporter().run(import_request)

dataset = AnalysisReadyDatasetBuilder().run(
    import_result.to_dataset_build_request(
        organism=Organism.RAT,
        input_intensity_scale="linear",
    )
)
```

Importers do not infer sample groups from column names and do not infer
differential design. Peptide-evidence handoff requires an explicit
`multi_site_policy`; ambiguous localisation and multi-site rows are retained
with diagnostics instead of being silently dropped.

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
and non-empty `site_key`, `display_id`, `organism`, `protein_namespace`,
`protein_identifier`, `gene_symbol`, and `site` columns. Workflow-created
differential results preserve that required protein context and optional
workflow-relevant protein metadata such as `protein_id` when present.
Display-indexed or stat-only differential result tables are not valid public
inputs.

The lower-level differential statistical executor may produce an internal
stat-only computation payload for workflow assembly. The public API result is
only `DifferentialAnalysisResult`, after the workflow has attached dataset
identity metadata.

## Total Protein And Protein-Aware Preparation

Total-protein correction and protein-aware preparation are separate
preprocessing contracts.

`DatasetTotalProteinCorrectionConfig(policy="subtract_log_total")` subtracts
matched log-scale total-protein abundance from log-scale phosphosite abundance:
`log2_phospho - log2_total`. This changes the phosphosite matrix values and the
dataset quantitative meaning. It requires total-protein input data and
log2-scale phospho/total values. It is not joint PTM/protein modelling.

`DatasetProteinAwarePreparationConfig(policy="prepare_model_inputs")` prepares
aligned phosphosite/protein input contracts and diagnostics. It does not change
the phosphosite matrix, does not subtract total protein, does not normalise
intensities, and does not run differential analysis. The default policy is
`"disabled"`.

Prepared protein-aware inputs are represented by
`ProteinAwarePreparationResult` and `ProteinAwarePreparationReport`. These are
preparation/audit objects, not model results. They expose matched
phosphosite/protein pairs, a sample-aligned protein covariate matrix,
per-site eligibility rows, missing protein-abundance diagnostics, ambiguous
mapping diagnostics, sample-alignment diagnostics, and policy/provenance fields.
The covariate matrix is a future modelling input contract; current
`DifferentialAnalysisWorkflow` execution does not consume it.

The public builder preparation stage uses explicit protein identifiers from
`site_metadata` (`protein_accession`, `protein_id`, or `protein_group_id`) and
matches them to `total.index`. Gene-symbol matching is not the public default.
`protein_mapping_policy="require_unambiguous"` makes missing or ambiguous
mappings ineligible for preparation. `protein_mapping_policy="allow_missing_with_report"`
allows missing site-protein identifiers or missing total-protein rows to remain
as phospho-only fallback rows in the report. Ambiguous site-to-protein mappings
or ambiguous total-protein row mappings are still excluded and reported.

Protein-aware preparation requires `total` input data. Missing total-protein
rows are per-site diagnostics; missing total input data is a build error when
`policy="prepare_model_inputs"` is selected. Phospho and total sample columns
must match in the same order at the builder boundary. Reordered, missing, or
extra total-protein sample columns are reported and make sites ineligible; the
builder does not reorder matrices for this preparation stage. Phospho and total
transformation states must also be compatible, meaning the same scale kind and
transformed flag. Incompatible transformation state is reported and excludes
sites from preparation.

Full joint PTM/protein modelling is not a dataset-preprocessing policy. It is
not enabled by total-protein subtraction and is not executed by protein-aware
preparation config. PhosPy does not claim MSstatsPTM-style inference or
MSstatsPTM equivalence for this preparation stage. Protein-aware preparation
does not run joint PTM/protein differential modelling.

```python
from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DatasetProteinAwarePreparationConfig,
    DatasetTotalProteinCorrectionConfig,
    Organism,
)

preprocessing = DatasetPreprocessingConfig()
assert preprocessing.total_protein_correction.policy == "none"
assert preprocessing.protein_aware_preparation.policy == "disabled"

subtraction = DatasetPreprocessingConfig(
    total_protein_correction=DatasetTotalProteinCorrectionConfig(
        policy="subtract_log_total"
    )
)

preparation_intent = DatasetPreprocessingConfig(
    protein_aware_preparation=DatasetProteinAwarePreparationConfig(
        policy="prepare_model_inputs",
        protein_mapping_policy="require_unambiguous",
    )
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        total=total,
        organism=Organism.RAT,
        input_intensity_scale="log2",
        preprocessing_config=preparation_intent,
    )
)

preparation = dataset.protein_aware_preparation
assert preparation is not None

matched_pairs = preparation.matched_pairs_dataframe()
protein_covariates = preparation.protein_covariate_matrix_dataframe()
report = dataset.preprocessing_report.protein_aware_preparation
site_eligibility = report.site_eligibility_dataframe()
missing_total_rows = report.missing_protein_abundance_diagnostics
ambiguous_mappings = report.ambiguous_mapping_diagnostics
sample_alignment = report.sample_alignment_diagnostics.to_payload()
```

## Batch Correction Preprocessing

Dataset preprocessing exposes an explicit batch-correction contract. The
default remains disabled:

```python
from phospy.api import DatasetBatchCorrectionConfig, DatasetPreprocessingConfig

preprocessing = DatasetPreprocessingConfig(
    batch_correction=DatasetBatchCorrectionConfig()
)
assert preprocessing.batch_correction.method == "none"
```

Users who want fixed-effect residualisation can opt in by name:

```python
preprocessing = DatasetPreprocessingConfig(
    batch_correction=DatasetBatchCorrectionConfig(
        method="linear_residualize_batch",
        batch_column="batch",
        condition_column="condition",
        preserve_condition_effects=True,
    )
)
```

`linear_residualize_batch` means fixed-effect residualisation of batch terms
while preserving condition effects by design. It is not ComBat, not RUV, and not
limma `removeBatchEffect` parity, and not mixed-effects modelling. During
dataset build, preprocessing resolves the configured batch and condition columns
from `sample_metadata`, validates design adequacy, and applies correction before
total-protein correction, site-matrix construction, normalisation, and
comparison building. If a log2 intensity transform is configured, correction
runs after that transform. Perfectly confounded batch/condition designs are
rejected because removing batch would also remove protected condition signal.

Dataset preprocessing reports can include a typed `batch_correction` sidecar:
`BatchCorrectionReport` with `BatchCorrectionPolicy` and
`BatchCorrectionDiagnostics`. The report records method, status (`"disabled"`,
`"applied"`, or `"rejected"`), batch and condition columns, observed batch and
condition levels when available, matrix shapes before and after the
batch-correction boundary, the design-preservation policy, confounding-check
status, warnings, and limitations. The builder reports default `method="none"`
as `"disabled"`. Requested correction either returns `"applied"` with the
corrected phosphosite matrix in the analysis-ready dataset, or fails clearly
when metadata or design adequacy is invalid.

Minimal dataset-build example with batch and condition metadata:

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBatchCorrectionConfig,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    IntensityScaleKind,
    Organism,
)

phospho = pd.DataFrame(
    {
        "sample_1": [10.0, 2.0],
        "sample_2": [15.0, 7.0],
        "sample_3": [14.0, 1.0],
        "sample_4": [19.0, 6.0],
    },
    index=["MAPK14;Y182;", "AKT1;T308;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "AKT1"],
        "site": ["Y182", "T308"],
        "site_sequence": [
            "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
        ],
        "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
        "organism": ["rat", "rat"],
        "protein_namespace": ["protein_id", "protein_id"],
        "protein_identifier": ["MAPK14", "AKT1"],
        "protein_id": ["MAPK14", "AKT1"],
        "localisation_confidence": [0.95, 0.92],
    },
    index=phospho.index.copy(),
)
sample_metadata = pd.DataFrame(
    {
        "batch": ["run_1", "run_2", "run_1", "run_2"],
        "condition": ["control", "control", "treated", "treated"],
    },
    index=phospho.columns.copy(),
)
preprocessing = DatasetPreprocessingConfig(
    batch_correction=DatasetBatchCorrectionConfig(
        method="linear_residualize_batch",
        batch_column="batch",
        condition_column="condition",
        preserve_condition_effects=True,
    )
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        organism=Organism.RAT,
        input_intensity_scale=IntensityScaleKind.LOG2,
        preprocessing_config=preprocessing,
    )
)

report = dataset.preprocessing_report.batch_correction
assert report is not None
print(report.status)
print(report.method)
print(report.confounding_check_status)
print(report.batch_levels)
print(report.condition_levels)
print(report.limitations)
```

The `condition` column is used here only to protect condition effects during
batch residualisation. It does not replace the explicit
`ExperimentalDesign` required by differential analysis.

## Result Construction Contracts

Public-looking result classes do not all have the same construction contract:

| Result object | Direct construction contract | Identity guarantee |
| --- | --- | --- |
| `DifferentialAnalysisResult` | Strict user-constructible public result. Use direct construction only with complete public contrast tables. | Requires encoded `site_key` index, matching `site_key` column, non-empty `display_id`, `organism`, `protein_namespace`, `protein_identifier`, `gene_symbol`, and `site`, coherent protein-scoped display/site metadata, and contrast tables aligned to residual-statistic indexes. |
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
`pred_mat`, and activity tables when activity is enabled. Opt-in Kinase
Library scoring modes additionally expose `kinase_library_motif_scores` and,
for combined scoring, `combined_profile_motif_scores`. Use
`activity_result.activity_matrix` as the primary activity-score matrix.
`activity_scores`, `weighted_activity`, and `to_dataframe()` are compatibility
accessors for that same primary matrix.

Stable kinase activity result fields are:

- `activity_matrix`: primary kinase-by-condition activity scores for the selected method.
- `substrate_count_matrix`: kinase-by-condition substrate counts used by the selected method when defined.
- `method_diagnostics`: typed method diagnostics, for example weighted-substrate or KSEA diagnostics.
- `policy_provenance`: scientific policy records attached to the activity method.

Optional activity fields are present only when the selected method produces
them:

- `p_value_matrix`
- `q_value_matrix`
- `confidence_interval_low`
- `confidence_interval_high`

Legacy activity sidecars remain available for existing users:
`thresholded_substrate_mean_activity`, `thresholded_substrate_counts`,
`activity_substrate_counts`, `target_counts`, `target_table`, and
`statistics_table`.

Use export helpers on scoring, prediction, and activity result objects for safe
snapshot copies:

```python
profile_scores = kinase_result.scoring_result.to_dataframe()
prediction_matrix = kinase_result.prediction_result.to_dataframe()
```

`kinase_result.provenance.scientific_policies` lists the active scientific
scoring policies with stable IDs, assumptions, parameters, and output scale
notes for auditability.

The default kinase scoring mode remains `"phosr_rank_weighted"`. To score with
Kinase Library-style sequence motifs, pass `KinaseScoringConfig` with
`scoring_mode="kinase_library_motif"` and provide a compatible
`kinase_library_resource` on `KinaseWorkflowRequest`. Kinase Library workflow
scores are relative motif support scores normalized to a unit interval per
kinase matrix; they are not calibrated probabilities or direct activity
evidence.

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
this release. Rat has a bundled beginner lane in this release. Human and mouse
workflows should use an explicit `ReferenceBundle`; the supported way to build
one from local source files is `ReferenceBundleBuilder`.

Use `ReferenceBundle` directly only when you already have validated DataFrames.
It requires:

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

For human or mouse local source files, prefer the builder so file-to-reference
normalisation and provenance are consistent:

```python
from phospy.api import (
    Organism,
    ReferenceBundleBuildRequest,
    ReferenceBundleBuilder,
)

references = ReferenceBundleBuilder().run(
    ReferenceBundleBuildRequest(
        organism=Organism.MOUSE,
        kinase_substrate_path="mouse_kinase_substrates.csv",
        site_sequence_path="mouse_site_sequences.csv",
        source_name="local curated kinase reference",
        source_version="2026-06-11",
        retrieved_at="2026-06-11",
        license="record the source license here",
        redistribution_status="record redistribution status here",
        identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
    )
)
```

The builder reads only local files. It does not scrape web resources, does not
invent missing sequence windows, and fails if any kinase-substrate site lacks
sequence context.

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
