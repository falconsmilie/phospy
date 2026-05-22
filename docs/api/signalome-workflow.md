# Signalome Workflow

This page explains the signalome workflow API. Use it after the kinase workflow
has produced a `KinaseWorkflowResult`.

## Purpose

`SignalomeWorkflow` interprets kinase score profiles into module assignments,
signalome module summaries, kinase networks, and protein-site context tables.
It requires explicit protein identity for interpreted sites.

```python
signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(kinase_result=kinase_result)
)
```

## Localisation Prerequisite

Signalome runs on site-level kinase outputs and should inherit explicit
localisation-confidence policy from dataset building:

- metadata column: `site_metadata["localisation_confidence"]`
- minimum threshold: `0.75`
- failure behaviour: dataset build fails if localisation metadata is missing,
  invalid, missing per-row, or below threshold
- why it matters: ambiguous site localisation can distort module assignments and
  kinase-network interpretation

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

The signalome request examples below focus on signalome API wiring and assume
this upstream dataset-localisation policy is already in place.

## Imports

```python
from phospy import SignalomeWorkflow
from phospy.api import (
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeOutputConfig,
    SignalomePerformanceConfig,
    SignalomeScientificConfig,
    SignalomeValidationConfig,
    SignalomeWorkflowRequest,
)
```

## Request Parameters

| Parameter | Type | Default | Required | How to Use It |
| --- | --- | --- | --- | --- |
| `kinase_result` | `KinaseWorkflowResult` | None | Yes | Result returned by `KinaseWorkflow.run(...)`. Its dataset must include explicit, non-empty `site_metadata.protein_id` for all interpreted sites. |
| `config` | `SignalomeConfig` | `SignalomeConfig()` | No | Grouped signalome settings for scientific interpretation, clustering, validation, output, and scale guardrails. |

Minimal call:

```python
result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(kinase_result=kinase_result)
)
```

Explicit config call:

```python
result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=SignalomeConfig(),
    )
)
```

## Protein Identity Requirement

Signalome requires explicit, non-empty `dataset.site_metadata.protein_id` for all
interpreted sites. The gene-symbol prefix in a site ID such as `TSC2;S939;` is
not treated as a protein-identity fallback.

```python
protein_ids = kinase_result.dataset.site_metadata["protein_id"]
if not protein_ids.astype("string").str.strip().ne("").all():
    raise ValueError("Signalome requires protein_id for every interpreted site.")
```

## Signalome Configuration

`SignalomeConfig` groups options by user intent.

| Parameter | Type | Default | How to Use It |
| --- | --- | --- | --- |
| `scientific` | `SignalomeScientificConfig` | `SignalomeScientificConfig()` | Controls substrate support and assignment policy. |
| `clustering` | `SignalomeClusteringConfig` | `SignalomeClusteringConfig()` | Controls module count, automatic module selection, candidate scoring, and clustering engine. |
| `validation` | `SignalomeValidationConfig` | `SignalomeValidationConfig()` | Controls score preconditioning and mixed total-protein guardrails. |
| `output` | `SignalomeOutputConfig` | `SignalomeOutputConfig()` | Controls network thresholding and network policy. |
| `performance` | `SignalomePerformanceConfig` | `SignalomePerformanceConfig()` | Controls scale guardrails for exact tree construction and full candidate scoring. |

Use presets for common behaviour:

```python
strict = SignalomeConfig.strict()
permissive = SignalomeConfig.permissive_missing_scores()
sampled = SignalomeConfig.sampled_candidate_scoring()
```

`SignalomeConfig.sampled_candidate_scoring()` changes candidate module-count
scoring only and keeps scale guards enabled.

## Scientific Parameters

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `substrate_support_cutoff` | `float` | `0.5` | `0.0 <= value <= 1.0` | Prediction support cutoff for kinase-supported substrates. |
| `assignment_policy` | `str` | `"cutoff_binary"` | `"cutoff_binary"`, `"weighted_top"` | Controls how sites are assigned to kinase-supported signalome modules. |

Example:

```python
scientific = SignalomeScientificConfig(
    substrate_support_cutoff=0.5,
    assignment_policy="cutoff_binary",
)
```

Use `"weighted_top"` when you want top-supported assignments to carry weighted
support rather than only cutoff membership.

## Clustering Parameters

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `module_count` | `int` or `None` | `None` | Integer `>= 1`, or `None` | Explicit module count. Use `None` for automatic module selection. |
| `module_selection_primary_correlation_threshold` | `float` | `0.5` | `0.0 <= value <= 1.0` | Primary correlation threshold for automatic module selection. |
| `module_selection_fallback_correlation_threshold` | `float` | `0.1` | `0.0 <= value <= 1.0` | Fallback threshold when the primary threshold does not select a usable module count. |
| `module_selection_max_clusters` | `int` | `10` | Integer `>= 1` | Largest candidate module count considered during automatic selection. |
| `candidate_scoring_policy` | `str` | `"full"` | `"full"`, `"sampled"` | Controls candidate module-count scoring cost. `"sampled"` reduces candidate scoring cost only. |
| `clustering_engine` | `str` | `"scipy_hierarchical"` | `"scipy_hierarchical"`, `"exact_python"` | Selects the clustering implementation. `"scipy_hierarchical"` is preferred for production runs. |

Automatic module selection example:

```python
clustering = SignalomeClusteringConfig(
    module_count=None,
    module_selection_primary_correlation_threshold=0.5,
    module_selection_fallback_correlation_threshold=0.1,
    module_selection_max_clusters=10,
    candidate_scoring_policy="full",
    clustering_engine="scipy_hierarchical",
)
```

Explicit module-count example:

```python
clustering = SignalomeClusteringConfig(
    module_count=6,
    candidate_scoring_policy="full",
    clustering_engine="scipy_hierarchical",
)
```

Sampled candidate scoring example:

```python
clustering = SignalomeClusteringConfig(
    candidate_scoring_policy="sampled",
    clustering_engine="scipy_hierarchical",
)
```

`candidate_scoring_policy="sampled"` reduces candidate module-count scoring
cost only. It does not remove exact tree construction, and it does not make the
full signalome workflow approximate.

## Validation Parameters

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `score_preconditioning_policy` | `str` | `"error_on_drop"` | `"error_on_drop"`, `"allow_and_report"` | Handles all-missing downstream score rows before signalome analysis. |
| `allow_mixed_total_protein_quantitative_meaning` | `bool` | `False` | `True`, `False` | Allows signalome analysis on datasets containing both total-corrected and uncorrected rows. Keep `False` unless that mixed interpretation is intentional. |

Strict validation example:

```python
validation = SignalomeValidationConfig(
    score_preconditioning_policy="error_on_drop",
    allow_mixed_total_protein_quantitative_meaning=False,
)
```

Allow-and-report example:

```python
validation = SignalomeValidationConfig(
    score_preconditioning_policy="allow_and_report",
    allow_mixed_total_protein_quantitative_meaning=False,
)
```

When `"allow_and_report"` is used, inspect
`result.score_preconditioning_diagnostics` after the run.

## Output Parameters

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `network_correlation_threshold` | `float` | `0.5` | `0.0 <= value <= 1.0` | Correlation threshold used by the kinase-network policy. |
| `network_policy` | `str` | `"signed"` | `"positive_only"`, `"absolute_threshold"`, `"signed"` | Controls which kinase correlations are represented as network edges. |

Example:

```python
output = SignalomeOutputConfig(
    network_correlation_threshold=0.5,
    network_policy="signed",
)
```

Network policy guide:

| Value | Behaviour | When to Use It |
| --- | --- | --- |
| `"signed"` | Keeps signed relationships meeting the configured threshold policy. | Good default when positive and negative relationships both matter. |
| `"positive_only"` | Keeps positive relationships only. | Useful for simpler activation-oriented network views. |
| `"absolute_threshold"` | Uses absolute correlation magnitude. | Useful when both strong positive and strong negative associations are relevant. |

## Performance Parameters

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `max_exact_tree_sites` | `int` | `2000` | Integer `>= 1` | Hard guard for exact tree construction. |
| `max_full_candidate_scoring_sites` | `int` | `2000` | Integer `>= 1` | Hard guard for full candidate module-count scoring. |

Example:

```python
performance = SignalomePerformanceConfig(
    max_exact_tree_sites=2000,
    max_full_candidate_scoring_sites=2000,
)
```

Raise these values only after considering runtime and memory for your dataset.

## Full Signalome Workflow Example

```python
from phospy import SignalomeWorkflow
from phospy.api import (
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeOutputConfig,
    SignalomePerformanceConfig,
    SignalomeScientificConfig,
    SignalomeValidationConfig,
    SignalomeWorkflowRequest,
)

config = SignalomeConfig(
    scientific=SignalomeScientificConfig(
        substrate_support_cutoff=0.5,
        assignment_policy="cutoff_binary",
    ),
    clustering=SignalomeClusteringConfig(
        module_count=None,
        module_selection_primary_correlation_threshold=0.5,
        module_selection_fallback_correlation_threshold=0.1,
        module_selection_max_clusters=10,
        candidate_scoring_policy="full",
        clustering_engine="scipy_hierarchical",
    ),
    validation=SignalomeValidationConfig(
        score_preconditioning_policy="error_on_drop",
        allow_mixed_total_protein_quantitative_meaning=False,
    ),
    output=SignalomeOutputConfig(
        network_correlation_threshold=0.5,
        network_policy="signed",
    ),
    performance=SignalomePerformanceConfig(
        max_exact_tree_sites=2000,
        max_full_candidate_scoring_sites=2000,
    ),
)

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=config,
    )
)

print(signalome_result.module_assignments.table.shape)
print(signalome_result.signalome_modules.table.shape)
print(signalome_result.kinase_network.edges.shape)
```

## Output Tables

Common result tables include:

```python
module_assignments = signalome_result.module_assignments.table
signalome_modules = signalome_result.signalome_modules.table
network_edges = signalome_result.kinase_network.edges
network_nodes = signalome_result.kinase_network.nodes
```

Sidecar tables help with interpretation and auditing:

```python
expanded_signalome = signalome_result.to_dataframe()
site_membership = signalome_result.site_membership_dataframe()
protein_site_context = signalome_result.protein_site_context_dataframe()
```

Runtime diagnostics are recorded in provenance:

```python
scale_guard = signalome_result.provenance.workflow_parameters["scale_guard"]
print(scale_guard["clustering_engine"])
print(scale_guard["tree_generation_mode"])
print(scale_guard["tree_generation_is_approximate"])
print(scale_guard["candidate_scoring_mode"])
print(scale_guard["candidate_scoring_is_approximate"])
print(scale_guard["max_exact_tree_sites"])
print(scale_guard["max_full_candidate_scoring_sites"])
```

Signalome module and network scores are derived summaries over upstream kinase
score profiles. They are not probabilities, calibrated confidence values, or
direct proof of causal regulation.
