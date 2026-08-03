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
- datasets with explicit `protein_group_id` grouping metadata for interpreted
  sites
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
- non-empty `dataset.site_metadata.protein_group_id` values for interpreted
  sites

Signalome requires `protein_group_id` because it groups retained phosphosites
into protein-level module and protein-site context summaries. `protein_group_id`
is Signalome-specific grouping metadata, not core protein identity. The legacy
`dataset.site_metadata.protein_id` column is accepted only as a migration alias
for older datasets and bundles. Core protein identity remains the dataset-level
`organism`, `protein_namespace`, and `protein_identifier` metadata carried by
the `site_key` row identity contract. Do not replace `protein_namespace` or
`protein_identifier` with `protein_group_id`.

Signalome does not reinterpret `display_id` as row identity, does not infer
`protein_group_id` from `gene_symbol` or `display_id`, and does not repair
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

For production workflow validation, use the signalome preset that carries the
same site-level threshold requirement into the request:

```python
from phospy.api import SignalomeConfig

config = SignalomeConfig.production()
```

The default `SignalomeConfig()` is the production mode and is the recommended
entry point. It requires site-level localisation evidence with the production
threshold and uses five paired finite observations for new network edges.

Historical exploratory behavior is still available, but it must be named
explicitly:

```python
from phospy.api import SignalomeConfig

config = SignalomeConfig.compatibility()
```

Reference-context compatibility remains conservative in both modes: unknown
context fails unless
`ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT` is set
explicitly.

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
| `output` | `SignalomeOutputConfig` | network threshold, network policy, and optional minimum paired finite observations. |
| `performance` | `SignalomePerformanceConfig` | scale guardrails for exact tree construction and full candidate scoring. |

Useful presets:

```python
strict = SignalomeConfig.strict()
production = SignalomeConfig.production()
compatibility = SignalomeConfig.compatibility()
permissive = SignalomeConfig.permissive_missing_scores()
sampled = SignalomeConfig.sampled_candidate_scoring()
```

Important fields:

| Field | Default | Notes |
| --- | --- | --- |
| `mode` | `"production"` | Recommended mode. Use `"exploratory_compatibility"` only through `SignalomeConfig.compatibility()` for legacy/exploratory behavior. |
| `scientific.substrate_support_cutoff` | `0.5` | Prediction support cutoff for kinase-supported substrates. |
| `scientific.assignment_policy` | `"cutoff_binary"` | Also supports `"weighted_top"`. |
| `clustering.module_count` | `None` | Use `None` for automatic module selection. |
| `clustering.candidate_scoring_policy` | `"full"` | `"sampled"` approximates candidate module-count scoring only. |
| `clustering.clustering_engine` | `"scipy_hierarchical"` | `"exact_python"` is also available. |
| `validation.score_preconditioning_policy` | `"error_on_drop"` | `"allow_and_report"` drops all-missing score rows and reports counts. |
| `validation.localisation_requirement` | `LocalisationRequirement.production_site_level()` | Workflow-level localisation requirement. Production requires present site-level localisation with minimum probability 0.75. |
| `output.network_policy` | `"signed"` | Also supports `"positive_only"` and `"absolute_threshold"`. |
| `output.network_min_paired_finite_observations` | `5` | Production requires at least five paired finite observations. Compatibility mode may lower this to the public floor of three. |
| `performance.max_exact_tree_sites` | `2000` | Exact tree scale guardrail. |

## Running the Workflow

```python
from phospy import SignalomeWorkflow
from phospy.api import SignalomeConfig, SignalomeWorkflowRequest

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=SignalomeConfig.production(),
    )
)
```

Use an explicit config when you need to document a policy choice:

```python
from phospy.api import SignalomeConfig

request = SignalomeWorkflowRequest(
    kinase_result=kinase_result,
    config=SignalomeConfig.compatibility(),
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
| `kinase_network` | `KinaseNetwork`; use `.edges`, `.nodes`, `.correlation_diagnostics`, and optional `.candidate_correlations`. |
| `module_selection_diagnostics` | Module-count selection diagnostics, including the automatic-selection `stability_report`. |
| `score_preconditioning_diagnostics` | Score-row preconditioning diagnostics. |
| `alignment_diagnostics` | Dataset/score/prediction alignment diagnostics. |
| `expanded_signalome` / `to_dataframe()` | Optional flattened signalome table. |
| `site_membership` / `site_membership_dataframe()` | Optional site-membership sidecar. |
| `protein_site_context` / `protein_site_context_dataframe()` | Optional protein-site context sidecar. |
| `provenance` | Workflow provenance. |

Site-level public sidecars include `site_key` and `display_id` where
applicable. Signalome protein-group sidecars and assignment tables use
`protein_group_id` for Signalome grouping identity. `site_key` remains the row
identity.

## Interpreting the Result

Signalome runs on the shared intersection of `site_key` values across the
dataset, prediction matrix, and downstream score matrix. Repeated `display_id`
values can appear in outputs when distinct `site_key` rows preserve different
protein context.

Missing kinase correlations remain missing. A value of `0.0` means a finite
near-zero correlation was estimated.

`output.network_min_paired_finite_observations` controls candidate edge
eligibility. For each kinase pair, PhosPy counts rows where both downstream
score profiles are finite. Pairs below the effective minimum are skipped before
edge thresholding. The public floor is three and the built-in default is five;
legacy bundle payloads that recorded threshold two remain readable as historical
results, but new replay/re-execution must migrate to at least three. Constant
score profiles, missing scores, non-finite scores, and finite correlations below
the configured network policy are also skipped and reported in
`kinase_network.correlation_diagnostics` and workflow provenance.
The optional `kinase_network.candidate_correlations` table lists candidate
pairs with `correlation_status`, `valid_observations`, and
`correlation_reason`.

Signalome clustering prepares its score matrix by dropping fully missing
kinase/dimension columns first, then median-imputing only partially missing
cells in the retained columns. Dropped all-missing labels, dropped-cell counts,
per-column imputation counts, retained labels, and the prepared-matrix
fingerprint are recorded in structured diagnostics and provenance. If no
dimension remains after dropping fully missing columns, execution fails before
tree construction.

When `clustering.module_count` is `None`, automatic module-count selection also
returns `module_selection_diagnostics.stability_report`. This typed report
separates deterministic reproducibility from descriptive stability. It records
the seeded perturbation policy, perturbation count, selected-count frequency,
pairwise coassignment similarity, threshold-grid sensitivity, status
(`stable`, `unstable`, or `not_computable`), and limitations. Frequencies and
similarity scores are descriptive diagnostics; they are not p-values,
confidence probabilities, or biological validation. Inputs with insufficient
sample structure report `not_computable` rather than a fabricated stability
score.

Key output meanings:

| Output | Meaning |
| --- | --- |
| `module_assignments.top_kinase` | Top supported kinase candidate for a site by prediction score. Ties are reported through candidate and selection-policy columns. |
| `module_assignments.module_id` | Score-derived candidate kinase-supported module ID for the current dataset and config. |
| `module_assignments.module_top_kinase` | Top supported kinase candidate summarized across the candidate module. It is a label, not a causal mechanism claim. |
| `signalome_modules` | Module-by-kinase percentages derived from the configured assignment policy. Values summarize support shares within modules. |
| `kinase_network.nodes` | Kinases retained in the aligned prediction and downstream score matrices. `degree` counts retained correlation edges and `n_substrates` counts predicted substrates above the support cutoff. |
| `kinase_network.edges` | Correlation edges between kinase score profiles. `source_kinase` and `target_kinase` are deterministic table labels, not inferred direction. Each accepted edge includes `valid_observations`. |
| `kinase_network.edges.correlation` | Edge weight derived from pairwise finite downstream score-profile correlations. `signed`, `positive_only`, and `absolute_threshold` policies control thresholding and whether sign is retained or stored as magnitude. |
| `kinase_network.candidate_correlations` | Candidate pairwise correlations before edge filtering, with status and paired finite observation counts. |
| `kinase_network.correlation_diagnostics` | Retained edge count plus skipped-edge counts for below-threshold, insufficient paired observations, constant profiles, missing scores, non-finite scores, and undefined correlations. |
| `expanded_signalome.regulated_module_ids` | Stable legacy field name for score-supported module IDs linked to a focal kinase. It does not mean experimental regulation was shown. |

Network edges are exploratory score-profile associations. Correlations do not
establish causality, direction, or experimental validation of signalling
relationships. They are not physical or causal interactions.

`candidate_scoring_policy="sampled"` approximates candidate module-count scoring
only. It does not make every signalome step approximate, and scale-guard
diagnostics are recorded in provenance.

## Provenance and Reproducibility

Workflow provenance records upstream kinase provenance, resolved config,
alignment diagnostics, score-preconditioning diagnostics, clustering preparation
diagnostics, scale-guard decisions, scientific policy records, table
fingerprints, Signalome mode, Signalome grouping identity, and signalome score
semantics. The semantics include the network threshold, requested/effective
paired-observation minimum, stable policy identifier, correlation basis, edge
directionality, skipped-edge diagnostics, prepared clustering matrix
fingerprint, and interpretation limits.

Signalome provenance uses the same causal site-row attrition contract as other
workflows. `row_attrition`, when present, records only stage-local site-row
removal proven from that stage's input and output site indexes. Compatibility
diagnostics such as `row_attrition_metrics`, alignment diagnostics, and
correlation skipped-edge counts remain separate diagnostics and are not treated
as site-row removal.

## Limitations

- Requires explicit `protein_group_id` grouping metadata for interpreted sites so
  retained sites can be summarized by protein in signalome module context.
- Accepts legacy `protein_id` only as a migration alias and rejects conflicting
  `protein_group_id`/`protein_id` values.
- Does not infer `protein_group_id` or core protein identity from `gene_symbol`
  or display labels.
- Does not run kinase scoring or prediction itself.
- Module and network-style outputs are derived summaries, not causal proof or
  experimental evidence of signalling relationships.
- Signalome network output does not report p-values, confidence intervals, or
  multiple-testing correction; adding inferential statistics requires a separate
  statistical policy decision.
- Mixed corrected/uncorrected total-protein quantitative meaning is rejected by
  default unless explicitly allowed.

## Minimal Example

```python
from phospy import SignalomeWorkflow
from phospy.api import (
    ReferenceContextCompatibilityPolicy,
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeValidationConfig,
    SignalomeWorkflowRequest,
)

config = SignalomeConfig(
    clustering=SignalomeClusteringConfig(
        module_count=None,
        clustering_engine="scipy_hierarchical",
    ),
    validation=SignalomeValidationConfig(
        reference_context_compatibility_policy=(
            ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
        )
    ),
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
