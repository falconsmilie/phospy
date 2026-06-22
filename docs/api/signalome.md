# Signalome Workflow

`SignalomeWorkflow` interprets kinase workflow output into score-derived
signalome module assignments, module summaries, kinase score-profile association
tables, and site/protein context sidecars. Use it after `KinaseWorkflow` has
produced a `KinaseWorkflowResult`.

## When to Use This Workflow

Use this workflow when you want exploratory module and network-style summaries
from kinase scoring/prediction output.

Good fits:

- site-keyed kinase prediction and downstream score matrices
- datasets with explicit `protein_id` grouping metadata for interpreted sites
- module assignment and kinase score-profile association summaries for
  exploratory analysis

Signalome outputs are derived summaries. They are not probabilities, calibrated
confidence values, causal evidence, or experimental validation of signalling
relationships.

## Inputs

`SignalomeWorkflowRequest.kinase_result` must be a `KinaseWorkflowResult` from
the kinase workflow.

The upstream result must provide:

- a valid `AnalysisReadyPhosphoDataset` with `site_key` row identity
- `prediction_result.pred_mat`
- an authoritative downstream score matrix from `scoring_result`
- overlapping site keys and kinase columns between prediction and score
  matrices
- non-empty `dataset.site_metadata.protein_id` values for interpreted sites

Signalome does not reinterpret `display_id` as row identity and does not repair
missing protein grouping metadata.

## Localisation Prerequisite

Signalome should inherit site-level localisation checks from dataset building.
A typical upstream policy is:

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

With this policy, dataset build fails when localisation metadata is missing,
invalid, missing per row, or below threshold.

## Request Object

Use `SignalomeWorkflowRequest`.

Important fields:

| Field | Meaning |
| --- | --- |
| `kinase_result` | Upstream `KinaseWorkflowResult`. |
| `config` | `SignalomeConfig` grouped by scientific, clustering, validation, output, and performance intent. |

Constructing the request records intent only. `SignalomeWorkflow.run(...)`
validates the upstream kinase result, matrix alignment, site identity, protein
grouping metadata, and config before execution.

## Request Configuration

Use `SignalomeConfig`.

Config sections:

| Section | Class | Notes |
| --- | --- | --- |
| `scientific` | `SignalomeScientificConfig` | `substrate_support_cutoff`, `assignment_policy`. |
| `clustering` | `SignalomeClusteringConfig` | module count, automatic selection thresholds, candidate scoring policy, clustering engine. |
| `validation` | `SignalomeValidationConfig` | score preconditioning, localisation, mixed total-protein guardrails. |
| `output` | `SignalomeOutputConfig` | network threshold and network policy. |
| `performance` | `SignalomePerformanceConfig` | scale guardrails for exact tree construction and full candidate scoring. |

Useful presets:

```python
strict = SignalomeConfig.strict()
permissive = SignalomeConfig.permissive_missing_scores()
sampled = SignalomeConfig.sampled_candidate_scoring()
```

Important fields:

| Field | Default | Notes |
| --- | --- | --- |
| `scientific.substrate_support_cutoff` | `0.5` | Prediction support cutoff for kinase-supported substrates. |
| `scientific.assignment_policy` | `"cutoff_binary"` | Also supports `"weighted_top"`. |
| `clustering.module_count` | `None` | Use `None` for automatic module selection. |
| `clustering.candidate_scoring_policy` | `"full"` | `"sampled"` approximates candidate module-count scoring only. |
| `clustering.clustering_engine` | `"scipy_hierarchical"` | `"exact_python"` is also available. |
| `validation.score_preconditioning_policy` | `"error_on_drop"` | `"allow_and_report"` drops all-missing score rows and reports counts. |
| `output.network_policy` | `"signed"` | Also supports `"positive_only"` and `"absolute_threshold"`. |
| `performance.max_exact_tree_sites` | `2000` | Exact tree scale guardrail. |

## Running the Workflow

```python
from phospy import SignalomeWorkflow
from phospy.api import SignalomeWorkflowRequest

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(kinase_result=kinase_result)
)
```

Use an explicit config when you need to document a policy choice:

```python
from phospy.api import SignalomeConfig

request = SignalomeWorkflowRequest(
    kinase_result=kinase_result,
    config=SignalomeConfig.strict(),
)
```

## Result Object

`SignalomeWorkflow.run(...)` returns `SignalomeWorkflowResult`.

Important fields and helpers:

| Field or helper | Meaning |
| --- | --- |
| `dataset` | Input dataset from the kinase result. |
| `kinase_result` | Upstream kinase result. |
| `module_assignments` | `SignalomeAssignments`; use `.table` or `.to_pandas()`. |
| `signalome_modules` | `SignalomeModules`; use `.table` or `.to_pandas()`. |
| `kinase_network` | `KinaseNetwork`; use `.edges`, `.nodes`, and optional `.candidate_correlations`. |
| `module_selection_diagnostics` | Module-count selection diagnostics. |
| `score_preconditioning_diagnostics` | Score-row preconditioning diagnostics. |
| `alignment_diagnostics` | Dataset/score/prediction alignment diagnostics. |
| `expanded_signalome` / `to_dataframe()` | Optional flattened signalome table. |
| `site_membership` / `site_membership_dataframe()` | Optional site-membership sidecar. |
| `protein_site_context` / `protein_site_context_dataframe()` | Optional protein-site context sidecar. |
| `provenance` | Workflow provenance. |

Site-level public sidecars include `site_key` and `display_id` where
applicable. `site_key` remains the row identity.

## Interpreting the Result

Signalome runs on the shared intersection of `site_key` values across the
dataset, prediction matrix, and downstream score matrix. Repeated `display_id`
values can appear in outputs when distinct `site_key` rows preserve different
protein context.

Missing kinase correlations remain missing. A value of `0.0` means a finite
near-zero correlation was estimated.

Key output meanings:

| Output | Meaning |
| --- | --- |
| `module_assignments.top_kinase` | Top supported kinase candidate for a site by prediction score. Ties are reported through candidate and selection-policy columns. |
| `module_assignments.module_id` | Score-derived candidate kinase-supported module ID for the current dataset and config. |
| `module_assignments.module_top_kinase` | Top supported kinase candidate summarized across the candidate module. It is a label, not a causal mechanism claim. |
| `signalome_modules` | Module-by-kinase percentages derived from the configured assignment policy. Values summarize support shares within modules. |
| `kinase_network.nodes` | Kinases retained in the aligned prediction and downstream score matrices. `degree` counts retained correlation edges and `n_substrates` counts predicted substrates above the support cutoff. |
| `kinase_network.edges` | Correlation edges between kinase score profiles. `source_kinase` and `target_kinase` are deterministic table labels, not inferred direction. |
| `kinase_network.edges.correlation` | Edge weight derived from pairwise finite downstream score-profile correlations. `signed`, `positive_only`, and `absolute_threshold` policies control thresholding and whether sign is retained or stored as magnitude. |
| `kinase_network.candidate_correlations` | Candidate pairwise correlations before edge filtering, with status and valid-observation counts. |
| `expanded_signalome.regulated_module_ids` | Stable legacy field name for score-supported module IDs linked to a focal kinase. It does not mean experimental regulation was shown. |

Network edges are exploratory score-profile associations. Correlations do not
establish causality, direction, or experimental validation of signalling
relationships.

`candidate_scoring_policy="sampled"` approximates candidate module-count scoring
only. It does not make every signalome step approximate, and scale-guard
diagnostics are recorded in provenance.

## Provenance and Reproducibility

Workflow provenance records upstream kinase provenance, resolved config,
alignment diagnostics, score-preconditioning diagnostics, scale-guard decisions,
scientific policy records, table fingerprints, and signalome score semantics.
The semantics include the network threshold, policy, correlation basis, edge
directionality, and interpretation limits.

## Limitations

- Requires explicit `protein_id` grouping metadata for interpreted sites.
- Does not infer protein identity from display labels.
- Does not run kinase scoring or prediction itself.
- Module and network-style outputs are derived summaries, not causal proof or
  experimental evidence of signalling relationships.
- Mixed corrected/uncorrected total-protein quantitative meaning is rejected by
  default unless explicitly allowed.

## Minimal Example

```python
from phospy import SignalomeWorkflow
from phospy.api import (
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeWorkflowRequest,
)

config = SignalomeConfig(
    clustering=SignalomeClusteringConfig(
        module_count=None,
        clustering_engine="scipy_hierarchical",
    )
)

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=config,
    )
)

module_assignments = signalome_result.module_assignments.table
network_edges = signalome_result.kinase_network.edges
print(module_assignments.shape)
print(network_edges.head())
```
