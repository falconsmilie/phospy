# Signalome Analysis

`SignalomeWorkflow` turns a `KinaseWorkflowResult` into site modules, module-level
kinase summaries, and kinase score-profile associations.

Use it after kinase analysis when you want an exploratory view of how
kinase-supported phosphosites group within the current dataset. The result keeps
site, protein, and diagnostic context alongside the summary tables.

!!! note "Interpret summaries carefully"
    Signalome outputs are descriptive, score-derived summaries. They are not
    pathway probabilities, causal evidence, proof of physical interaction, or
    experimental validation of signaling relationships.

This walkthrough starts from raw teaching tables so that `kinase_result` is
never hidden notebook state. It builds a five-site rat dataset, runs the kinase
workflow, and then summarises it with a configuration derived from the
production preset with one deliberate reference-context compatibility-policy
override. Run the Python blocks in order in one fresh process.

## Prepare the Phosphosite Measurements

Signalome groups phosphosites by their upstream kinase score profiles. Five
sites are used because the production network policy requires at least five
paired finite observations for a kinase pair.

```python
import pandas as pd

phospho = pd.DataFrame(
    {
        "sample_a": [1.00, 0.70, 0.85, 0.92, 0.66],
        "sample_b": [1.10, 0.80, 0.88, 0.96, 0.69],
        "sample_c": [0.95, 0.75, 0.92, 0.90, 0.72],
    },
    index=[
        "TSC2;S939;",
        "GSK3A;S21;",
        "MAPK14;Y182;",
        "AKT1;T308;",
        "SRC;Y416;",
    ],
)
```

`protein_group_id` is the explicit protein grouping identity used by Signalome.
It is separate from readable gene and site labels and must be non-empty for
every interpreted site.

```python
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["TSC2", "GSK3A", "MAPK14", "AKT1", "SRC"],
        "site": ["S939", "S21", "Y182", "T308", "Y416"],
        "site_sequence": [
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
            "PSGGGPGGSGRARTSSFAEPGGGGGGGGGGP",
            "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
            ("A" * 15) + "T" + ("A" * 15),
            ("A" * 15) + "Y" + ("A" * 15),
        ],
        "protein_identifier": ["TSC2", "GSK3A", "MAPK14", "AKT1", "SRC"],
        "protein_namespace": ["protein_id"] * phospho.shape[0],
        "protein_group_id": ["TSC2", "GSK3A", "MAPK14", "AKT1", "SRC"],
        "localisation_confidence": [0.95] * phospho.shape[0],
    },
    index=phospho.index,
)
```

## Require Localisation Evidence

Production Signalome requires site-level localisation evidence with probability
at least `0.75`. The dataset builder checks that evidence before either
workflow runs.

```python
from phospy.api import DatasetLocalisationConfig

localisation = DatasetLocalisationConfig(
    mode="require_threshold",
    confidence_column="localisation_confidence",
    min_confidence=0.75,
)
```

```python
from phospy.api import DatasetPreprocessingConfig

preprocessing_config = DatasetPreprocessingConfig(
    localisation=localisation,
)
```

## Build the Analysis-Ready Dataset

The input values are already suitable linear abundances for the upstream
profile-scoring method, so the request declares the linear scale without a
transformation.

```python
from phospy.api import DatasetBuildRequest, Organism

dataset_request = DatasetBuildRequest(
    phospho=phospho,
    site_metadata=site_metadata,
    organism=Organism.RAT,
    input_intensity_scale="linear",
    preprocessing_config=preprocessing_config,
)
```

```python
from phospy import AnalysisReadyDatasetBuilder

dataset_builder = AnalysisReadyDatasetBuilder()
dataset = dataset_builder.run(dataset_request)
```

## Configure the Upstream Kinase Analysis

This small dataset uses a custom reliability profile and deliberately permits
unknown dataset/reference context with a recorded caveat. The deterministic
prediction policy makes the upstream result repeatable.

```python
from phospy.advanced import (
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)

kinase_scoring_config = KinaseScoringConfig(
    reliability_profile=KinaseReliabilityProfile.CUSTOM,
    reference_context_compatibility_policy=(
        ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
    ),
)
```

```python
from phospy.advanced import KinasePredictionConfig

kinase_prediction_config = KinasePredictionConfig.deterministic()
```

Construct the kinase request with the bundled rat references.

```python
from phospy.api import KinaseWorkflowRequest, ReferencePreset

kinase_request = KinaseWorkflowRequest(
    dataset=dataset,
    references=ReferencePreset.AUTO,
    scoring_config=kinase_scoring_config,
    prediction_config=kinase_prediction_config,
    activity_config=None,
)
```

```python
from phospy import KinaseWorkflow

kinase_workflow = KinaseWorkflow()
kinase_result = kinase_workflow.run(kinase_request)
```

## Configure Signalome

`SignalomeConfig.production()` supplies the production localisation and network
guardrails. The teaching dataset has unknown reference-context metadata, so
copy that configuration and deliberately change only the compatibility policy.
This allows the unknown context and records a caveat; it does not establish
reference compatibility. Ordinary production use should require a known match.

```python
from dataclasses import replace

from phospy.advanced import SignalomeConfig

production_config = SignalomeConfig.production()
signalome_config = replace(
    production_config,
    validation=replace(
        production_config.validation,
        reference_context_compatibility_policy=(
            ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
        ),
    ),
)
```

## Construct and Run the Request

The request consumes the complete `KinaseWorkflowResult`, preserving its
dataset, score matrix, prediction matrix, references, diagnostics, and
provenance.

```python
from phospy.api import SignalomeWorkflowRequest

signalome_request = SignalomeWorkflowRequest(
    kinase_result=kinase_result,
    config=signalome_config,
)
```

```python
from phospy import SignalomeWorkflow

signalome_workflow = SignalomeWorkflow()
signalome_result = signalome_workflow.run(signalome_request)
```

## Inspect and Interpret the Result

Inspect site assignments, module-level kinase support, and retained kinase
score-profile associations through their public tables.

```python
assignment_view = signalome_result.module_assignments.table.loc[
    :, ["display_id", "module_id", "top_kinase"]
]
module_table = signalome_result.signalome_modules.table
network_edges = signalome_result.kinase_network.edges

print("Assignment shape:", assignment_view.shape)
print(assignment_view)
print("Module table shape:", module_table.shape)
print("Network edge shape:", network_edges.shape)
```

The assignment view has shape `(5, 3)`, the module table has shape `(1, 2)`,
and the retained edge table has shape `(1, 4)`. `module_id` labels
score-derived groups within this run. `top_kinase` identifies the strongest
support under the upstream scoring policy. Neither field establishes a pathway,
physical interaction, or causal mechanism.

Module-table cells are kinase support shares from 0 to 100. Network correlations
are associations across paired finite site scores. A missing edge means the
pair was ineligible or did not pass the selected filter; it is not evidence of
biological absence.

## At a Glance

| Question | Answer |
| --- | --- |
| **Input** | A completed `KinaseWorkflowResult` with aligned score and prediction matrices. |
| **Main request** | `SignalomeWorkflowRequest` |
| **Run** | `SignalomeWorkflow().run(request)` |
| **Main result** | `SignalomeWorkflowResult` |
| **Site output** | `result.module_assignments.table` |
| **Module output** | `result.signalome_modules.table` |
| **Network output** | `result.kinase_network.edges` |

## Input Requirements

The upstream `KinaseWorkflowResult` must provide:

- an `AnalysisReadyPhosphoDataset` with unique `site_key` rows;
- `prediction_result.pred_mat`;
- `scoring_result.authoritative_scores`;
- overlapping site rows and kinase columns across those matrices;
- usable `site_sequence` values from the dataset; and
- non-empty `site_metadata.protein_group_id` values for interpreted sites.

`protein_group_id` is the grouping identity used by Signalome. Core protein
identity still comes from `organism`, `protein_namespace`, and
`protein_identifier`. Signalome does not infer `protein_group_id` from
`gene_symbol` or `display_id`.

Production mode requires site-level localisation evidence with a minimum
probability of 0.75. Configure this requirement when building the upstream
dataset. Production Signalome fails when localisation is missing or below the
threshold.

```python
from phospy.api import DatasetLocalisationConfig, DatasetPreprocessingConfig

preprocessing = DatasetPreprocessingConfig(
    localisation=DatasetLocalisationConfig(
        mode="require_threshold",
        confidence_column="localisation_confidence",
        min_confidence=0.75,
    )
)
```

Reference-context compatibility is conservative by default. An unknown context
fails unless `"allow_unknown_with_caveat"` is selected deliberately in both the
upstream kinase settings and the Signalome validation settings.

## Reuse Another Kinase Result

The primary journey shows how to pass a completed `kinase_result` into a
`SignalomeWorkflowRequest`. Use `SignalomeConfig.production()` for routine
analysis when reference context is known to match. If the upstream result has
unknown context, select `"allow_unknown_with_caveat"` deliberately as shown in
the primary journey. Compatibility and permissive presets are intended for
explicit exploratory work and should not be mistaken for production evidence.

## Request

Create a `SignalomeWorkflowRequest`.

| Parameter | Type | Required or Default | Description | Main Constraint |
| --- | --- | --- | --- | --- |
| `kinase_result` | `KinaseWorkflowResult` | Required | Upstream kinase result to summarize. | Must contain aligned prediction and authoritative score matrices, plus `protein_group_id` for interpreted sites. |
| `config` | `SignalomeConfig` | `SignalomeConfig()` | Scientific, clustering, validation, output, and performance policies. | The default mode is `"production"`. |

### `SignalomeConfig`

| Parameter | Default | Description |
| --- | --- | --- |
| `scientific` | `SignalomeScientificConfig()` | Substrate support and assignment policy. |
| `clustering` | `SignalomeClusteringConfig()` | Module selection and clustering policy. |
| `validation` | `SignalomeValidationConfig()` | Score, localisation, reference-context, and quantitative guardrails. |
| `output` | `SignalomeOutputConfig()` | Network threshold and paired-observation settings. |
| `performance` | `SignalomePerformanceConfig()` | Scale guardrails for exact and candidate-scoring paths. |
| `mode` | `"production"` | Use `"exploratory_compatibility"` only for a named compatibility run. |

Useful constructors include:

- `SignalomeConfig.production()` for the recommended production policy;
- `SignalomeConfig.strict()` for strict score preconditioning;
- `SignalomeConfig.permissive_missing_scores()` to report and drop all-missing
  score rows;
- `SignalomeConfig.sampled_candidate_scoring()` for sampled module-count
  candidate scoring; and
- `SignalomeConfig.compatibility()` for explicit exploratory compatibility.

<details markdown="1">
<summary><strong>Complete Configuration Parameters</strong></summary>

### `SignalomeScientificConfig`

| Parameter | Default | Description | Main Constraint |
| --- | --- | --- | --- |
| `substrate_support_cutoff` | `0.5` | Prediction support cutoff used to assign kinase-supported substrates. | Must be between 0 and 1. |
| `assignment_policy` | `"cutoff_binary"` | Controls how site support contributes to modules. | `"weighted_top"` uses top-kinase weights when available. |

### `SignalomeClusteringConfig`

| Parameter | Default | Description | Main Constraint |
| --- | --- | --- | --- |
| `module_count` | `None` | Requested number of modules. | `None` enables automatic selection; an integer must be at least 1. |
| `module_selection_primary_correlation_threshold` | `0.5` | Primary threshold for automatic module-count selection. | Must be between 0 and 1. |
| `module_selection_fallback_correlation_threshold` | `0.1` | Fallback threshold for automatic module-count selection. | Must be between 0 and 1. |
| `module_selection_max_clusters` | `10` | Maximum candidate cluster count. | Must be at least 1. |
| `candidate_scoring_policy` | `"full"` | Candidate module-count scoring strategy. | `"sampled"` approximates candidate scoring. |
| `clustering_engine` | `"scipy_hierarchical"` | Tree-building engine. | `"exact_python"` is also supported. |

### `SignalomeValidationConfig`

| Parameter | Default | Description | Main Constraint |
| --- | --- | --- | --- |
| `score_preconditioning_policy` | `"error_on_drop"` | Handles score rows that are entirely missing. | `"allow_and_report"` drops those rows and records the decision. |
| `localisation_requirement` | `LocalisationRequirement.production_site_level()` | Requires site-level localisation evidence. | Production requires present evidence with probability at least 0.75. |
| `allow_mixed_total_protein_quantitative_meaning` | `False` | Allows mixed total-corrected and uncorrected meaning. | Keep `False` unless the mixture is intentional and documented. |
| `reference_context_compatibility_policy` | `"require_known_match"` | Handles unknown reference context. | `"allow_unknown_with_caveat"` records a caveat. |

`LocalisationRequirement` uses `require_present=True` and
`minimum_probability=0.75` in the production preset.

### `SignalomeOutputConfig`

| Parameter | Default | Description | Main Constraint |
| --- | --- | --- | --- |
| `network_correlation_threshold` | `0.5` | Threshold for retaining kinase score-profile associations. | Must be between 0 and 1. |
| `network_policy` | `"signed"` | Controls sign handling. | Also supports `"positive_only"` and `"absolute_threshold"`. |
| `network_min_paired_finite_observations` | `5` | Minimum paired finite sites for a kinase pair. | Public floor is 3; production requires an effective value of at least 5. |

### `SignalomePerformanceConfig`

| Parameter | Default | Description |
| --- | --- | --- |
| `max_exact_tree_sites` | `2000` | Maximum site count for exact tree construction. |
| `max_full_candidate_scoring_sites` | `2000` | Maximum site count for full candidate scoring. |

</details>

## Execution Behaviour and Errors

The primary journey above shows the public request-and-run boundary.

For the same upstream result, configuration, and package version, the production
path is deterministic. A seed used by an upstream adaptive kinase run still
affects the input to Signalome.

`WorkflowValidationError` is raised before result construction when required
matrices are missing, site or kinase axes do not overlap, `protein_group_id` is
missing, localisation requirements fail, strict score preconditioning would
remove rows, or network settings fall below the supported floor.

## Response

`SignalomeWorkflow.run(...)` returns a `SignalomeWorkflowResult`.

### Main Outputs

| Attribute | Format | Meaning |
| --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | The analysis-ready dataset carried from the upstream kinase result. |
| `kinase_result` | `KinaseWorkflowResult` | The complete upstream kinase result used by Signalome. |
| `module_assignments` | `SignalomeAssignments` | Site-level assignments; use `.table` or `.to_pandas()`. |
| `signalome_modules` | `SignalomeModules` | Module-by-kinase support shares; use `.table` or `.to_pandas()`. |
| `kinase_network` | `KinaseNetwork` | Edge, node, and candidate-correlation tables. |
| `module_selection_diagnostics` | Typed diagnostics or `None` | Why a module count was selected. |
| `clustering_preparation_diagnostics` | Typed diagnostics or `None` | How the score matrix was prepared for clustering. |
| `score_preconditioning_diagnostics` | Typed diagnostics or `None` | Rows retained or removed before clustering. |
| `alignment_diagnostics` | Typed diagnostics or `None` | Dataset, prediction, and score-matrix alignment. |
| `expanded_signalome` | `pandas.DataFrame` or `None` | Optional flattened reporting view. |
| `site_membership` | `pandas.DataFrame` or `None` | Optional site-level inclusion sidecar. |
| `protein_site_context` | `pandas.DataFrame` or `None` | Optional protein-level context sidecar. |
| `provenance` | `RunProvenance` or `None` | Resolved settings and input fingerprints. |
| `caveats` | `tuple[ResultCaveat, ...]` | Structured scientific limitations and warnings. |

### Module Assignments

`module_assignments.table` is indexed by `site_key`. The most useful columns are:

| Column | Meaning |
| --- | --- |
| `site_key`, `display_id` | Stable row identity and readable site label. |
| `gene_symbol`, `site`, `protein_group_id` | Site and protein context. |
| `module_id` | Module assigned within this run. |
| `top_kinase`, `top_score` | Highest-supported kinase label and score for the site. |
| `top_kinase_candidates`, `top_kinase_weights` | Candidate labels and weights when support is tied or weighted. |
| `top_kinase_tie_count`, `top_kinase_is_ambiguous` | Tie and ambiguity diagnostics. |
| `module_top_kinase` | Highest-supported kinase summarised across the module. |

`protein_accession`, `isoform_id`, selection-policy fields, and module-level tie
diagnostics are also retained when available.

### Module Summary

`signalome_modules.table` uses `module_id` rows and kinase columns. Each cell is
the kinase support share within that module, expressed from 0 to 100. A row
totals approximately 100 when support is assigned, or 0 when no support is
available.

### Kinase Network

`kinase_network.edges` contains retained pairwise associations:

| Column | Meaning |
| --- | --- |
| `source_kinase`, `target_kinase` | The two kinase labels. These are table labels, not inferred direction. |
| `correlation` | Pairwise score-profile association after the selected network policy. |
| `valid_observations` | Paired finite site count used for the estimate. |

`kinase_network.nodes`, when present, reports each kinase's retained edge degree
and substrate count. `candidate_correlations`, when present, adds
`correlation_status` and `correlation_reason` so that a near-zero finite
association can be distinguished from an undefined or ineligible pair.

<details markdown="1">
<summary><strong>Sidecar and Reporting Schemas</strong></summary>

### `site_membership`

| Column | Meaning |
| --- | --- |
| `site_key`, `display_id`, `site_id`, `site` | Site identity and readable labels. |
| `protein_group_id`, `protein_accession`, `isoform_id` | Protein context. |
| `site_cluster`, `protein_module_id` | Site- and protein-level assignments. |
| `included_in_module_table`, `excluded_reason` | Inclusion status and reason. |
| `top_kinase`, `top_kinase_score`, `top_kinase_weight` | Site-level kinase support. |
| `n_supported_kinases` | Number of supported kinases for the site. |

### `protein_site_context`

| Column | Meaning |
| --- | --- |
| `protein_group_id`, `n_sites` | Protein-group identity and site count. |
| `site_ids`, `site_keys`, `display_ids`, `site_clusters` | Sites and cluster labels in the group. |
| `n_distinct_site_clusters`, `protein_module_id` | Module summary. |
| `multi_site_protein`, `ambiguous_module_context` | Multi-site and ambiguity diagnostics. |
| `top_kinases_by_site`, `module_ids_by_site` | Per-site kinase and module mappings. |
| `site_key_to_display_id` | Mapping from stable identity to readable labels. |

### `expanded_signalome`

This optional table is a flattened reporting view. It contains at least
`site_key` and `display_id` for site rows. Treat the typed assignment, module,
and network objects as the primary result contracts.

</details>

## Interpret the Result

`module_id` labels describe score-derived groups for this dataset and
configuration. They are not stable pathway identifiers and should not be
compared across runs without a deliberate matching strategy.

`top_kinase` and `module_top_kinase` identify the strongest support under the
selected scoring policy. They do not establish a causal mechanism.

Network correlations summarize kinase score-profile associations across paired
finite sites. A missing correlation means the pair was not estimable or did not
pass filtering; a value near 0 means a finite, near-zero association was
estimated. Edges do not establish causality, physical interaction, or direction.

Automatic module-count diagnostics are descriptive. Stability frequencies and
similarities are not *p* values, confidence probabilities, or biological
validation.

## Common Issues

| Issue | What to Check |
| --- | --- |
| `protein_group_id` is missing. | Add non-empty values before kinase and Signalome analysis. Do not substitute gene symbols or display labels. |
| Localisation validation fails. | Build the upstream dataset with localisation evidence; production Signalome fails when localisation is missing. |
| No site rows overlap. | Confirm that `pred_mat` and `authoritative_scores` share `site_key` rows. |
| No kinase columns overlap. | Confirm that prediction and score matrices use matching kinase labels. |
| All-missing rows fail. | Keep strict mode, or choose `"allow_and_report"` only when the drop is acceptable and reviewed. |
| The network has no edges. | Review the threshold, sign policy, paired finite observations, constant profiles, and candidate-correlation diagnostics. |
| Production has too few sites. | Production edges require at least five paired finite observations. Add evidence rather than weakening the setting without justification. |

## Complete Runnable Example

```python
import pandas as pd

phospho = pd.DataFrame(
    {
        "sample_a": [1.00, 0.70, 0.85, 0.92, 0.66],
        "sample_b": [1.10, 0.80, 0.88, 0.96, 0.69],
        "sample_c": [0.95, 0.75, 0.92, 0.90, 0.72],
    },
    index=[
        "TSC2;S939;",
        "GSK3A;S21;",
        "MAPK14;Y182;",
        "AKT1;T308;",
        "SRC;Y416;",
    ],
)

site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["TSC2", "GSK3A", "MAPK14", "AKT1", "SRC"],
        "site": ["S939", "S21", "Y182", "T308", "Y416"],
        "site_sequence": [
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
            "PSGGGPGGSGRARTSSFAEPGGGGGGGGGGP",
            "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
            ("A" * 15) + "T" + ("A" * 15),
            ("A" * 15) + "Y" + ("A" * 15),
        ],
        "protein_identifier": ["TSC2", "GSK3A", "MAPK14", "AKT1", "SRC"],
        "protein_namespace": ["protein_id"] * phospho.shape[0],
        "protein_group_id": ["TSC2", "GSK3A", "MAPK14", "AKT1", "SRC"],
        "localisation_confidence": [0.95] * phospho.shape[0],
    },
    index=phospho.index,
)

from phospy.api import DatasetLocalisationConfig

localisation = DatasetLocalisationConfig(
    mode="require_threshold",
    confidence_column="localisation_confidence",
    min_confidence=0.75,
)

from phospy.api import DatasetPreprocessingConfig

preprocessing_config = DatasetPreprocessingConfig(
    localisation=localisation,
)

from phospy.api import DatasetBuildRequest, Organism

dataset_request = DatasetBuildRequest(
    phospho=phospho,
    site_metadata=site_metadata,
    organism=Organism.RAT,
    input_intensity_scale="linear",
    preprocessing_config=preprocessing_config,
)

from phospy import AnalysisReadyDatasetBuilder

dataset_builder = AnalysisReadyDatasetBuilder()
dataset = dataset_builder.run(dataset_request)

from phospy.advanced import (
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)

kinase_scoring_config = KinaseScoringConfig(
    reliability_profile=KinaseReliabilityProfile.CUSTOM,
    reference_context_compatibility_policy=(
        ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
    ),
)

from phospy.advanced import KinasePredictionConfig

kinase_prediction_config = KinasePredictionConfig.deterministic()

from phospy.api import KinaseWorkflowRequest, ReferencePreset

kinase_request = KinaseWorkflowRequest(
    dataset=dataset,
    references=ReferencePreset.AUTO,
    scoring_config=kinase_scoring_config,
    prediction_config=kinase_prediction_config,
    activity_config=None,
)

from phospy import KinaseWorkflow

kinase_workflow = KinaseWorkflow()
kinase_result = kinase_workflow.run(kinase_request)

from dataclasses import replace

from phospy.advanced import SignalomeConfig

production_config = SignalomeConfig.production()
signalome_config = replace(
    production_config,
    validation=replace(
        production_config.validation,
        reference_context_compatibility_policy=(
            ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
        ),
    ),
)

from phospy.api import SignalomeWorkflowRequest

signalome_request = SignalomeWorkflowRequest(
    kinase_result=kinase_result,
    config=signalome_config,
)

from phospy import SignalomeWorkflow

signalome_workflow = SignalomeWorkflow()
signalome_result = signalome_workflow.run(signalome_request)

assignment_view = signalome_result.module_assignments.table.loc[
    :, ["display_id", "module_id", "top_kinase"]
]
module_table = signalome_result.signalome_modules.table
network_edges = signalome_result.kinase_network.edges

print("Assignment shape:", assignment_view.shape)
print(assignment_view)
print("Module table shape:", module_table.shape)
print("Network edge shape:", network_edges.shape)
```

## Continue

- [Kinase analysis](kinase.md)
- [Prepare a dataset](dataset-build-workflow.md)
- [Scientific interpretation and limitations](../scientific-interpretation.md)
