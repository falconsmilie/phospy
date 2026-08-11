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

## Before You Start

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

## Example

Start with `kinase_result` returned by the [kinase workflow](kinase.md), then run
Signalome with the production policy:

```python
from phospy import SignalomeWorkflow
from phospy.advanced import SignalomeConfig
from phospy.api import SignalomeWorkflowRequest

request = SignalomeWorkflowRequest(
    kinase_result=kinase_result,
    config=SignalomeConfig.production(),
)

result = SignalomeWorkflow().run(request)

assignments = result.module_assignments.table
modules = result.signalome_modules.table
edges = result.kinase_network.edges

print(assignments.loc[:, ["display_id", "module_id", "top_kinase"]])
print(modules)
print(edges.head())
```

Use `SignalomeConfig.production()` for routine analysis. Compatibility and
permissive presets are intended for explicitly exploratory work and should not
be mistaken for production evidence.

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

## Run the Workflow

```python
result = SignalomeWorkflow().run(request)
```

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
| `module_top_kinase` | Highest-supported kinase summarized across the module. |

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

## Related Guides

- [Kinase analysis](kinase.md)
- [Prepare a dataset](dataset-build-workflow.md)
- [Scientific interpretation and limitations](../scientific-interpretation.md)
