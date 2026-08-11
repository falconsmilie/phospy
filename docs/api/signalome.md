# Signalome analysis workflow

## Plain-language introduction

`SignalomeWorkflow` turns a `KinaseWorkflowResult` into module assignments,
module-by-kinase summaries, kinase score-profile association tables, and
site/protein context sidecars.

Use it after kinase scoring and prediction when you want exploratory summaries
of score-supported phosphosite modules and kinase relationships. The workflow
expects an upstream kinase result with aligned site and kinase score matrices.
It returns a `SignalomeWorkflowResult` with result tables, diagnostics,
provenance, and caveats.

Signalome outputs are derived summaries. They are not probabilities, causal
evidence, pathway activation proof, or experimental validation of signalling
relationships.

## Input and dataset requirements

`SignalomeWorkflowRequest.kinase_result` must be a `KinaseWorkflowResult` from
[`KinaseWorkflow`](kinase.md). The upstream result must provide:

- a valid `AnalysisReadyPhosphoDataset` with `site_key` row identity;
- `prediction_result.pred_mat`;
- an authoritative downstream score matrix from `scoring_result`;
- overlapping `site_key` rows and kinase columns between prediction and score
  matrices;
- required `site_sequence` context from the upstream dataset;
- non-empty `dataset.site_metadata.protein_group_id` values for interpreted
  sites.

`protein_group_id` is Signalome-specific grouping metadata. Core protein
identity still comes from `organism`, `protein_namespace`, and
`protein_identifier`. Signalome does not infer `protein_group_id` from
`gene_symbol` or `display_id`. Legacy `protein_id` is accepted only as a
migration alias when it does not conflict with `protein_group_id`.

Signalome should inherit localisation checks from dataset building. Production
mode also requires site-level localisation evidence with minimum probability
`0.75`; production mode fails when localisation is missing or below threshold.

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

Reference-context compatibility is conservative by default. Unknown context
fails unless `ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT` is
set deliberately in the upstream kinase request and in the signalome config.

## Minimal end-to-end example

```python
from dataclasses import replace

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow, SignalomeWorkflow
from phospy.advanced import (
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
    SignalomeConfig,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    IntensityScaleKind,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
    SignalomeWorkflowRequest,
)

phospho = pd.DataFrame(
    {
        "sample_a": [1.00, 0.70, 0.85, 0.92, 0.66],
        "sample_b": [1.10, 0.80, 0.88, 0.96, 0.69],
        "sample_c": [0.95, 0.75, 0.92, 0.90, 0.72],
        "sample_d": [1.05, 0.79, 0.90, 0.94, 0.68],
        "sample_e": [1.08, 0.81, 0.91, 0.97, 0.70],
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
            ("A" * 15) + "S" + ("A" * 15),
            ("A" * 15) + "S" + ("A" * 15),
            ("A" * 15) + "Y" + ("A" * 15),
            ("A" * 15) + "T" + ("A" * 15),
            ("A" * 15) + "Y" + ("A" * 15),
        ],
        "display_id": [
            "TSC2;S939;",
            "GSK3A;S21;",
            "MAPK14;Y182;",
            "AKT1;T308;",
            "SRC;Y416;",
        ],
        "organism": ["rat"] * 5,
        "protein_namespace": ["protein_id"] * 5,
        "protein_identifier": ["TSC2", "GSK3A", "MAPK14", "AKT1", "SRC"],
        "protein_group_id": ["TSC2", "GSK3A", "MAPK14", "AKT1", "SRC"],
        "localisation_confidence": [0.95, 0.94, 0.96, 0.93, 0.92],
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
        scoring_config=KinaseScoringConfig(
            reliability_profile=KinaseReliabilityProfile.CUSTOM,
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        activity_config=None,
    )
)

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

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=signalome_config,
    )
)

module_assignments = signalome_result.module_assignments.table
network_edges = signalome_result.kinase_network.edges
print(module_assignments.loc[:, ["site_key", "display_id", "module_id", "top_kinase"]])
print(network_edges.head())
```

## Request model

Use `SignalomeWorkflowRequest`.

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `kinase_result` | `KinaseWorkflowResult` | Required | Upstream kinase result to summarize. | Must include aligned prediction and score matrices and a dataset with `protein_group_id` for interpreted sites. |
| `config` | `SignalomeConfig` | Default: `SignalomeConfig()` | Signalome policy grouped into scientific, clustering, validation, output, and performance sections. | Default mode is production. Invalid production policies fail at config construction or workflow validation. |

`SignalomeConfig`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `scientific` | `SignalomeScientificConfig` | Default: `SignalomeScientificConfig()` | Substrate-support and assignment policy. | See table below. |
| `clustering` | `SignalomeClusteringConfig` | Default: `SignalomeClusteringConfig()` | Module selection and clustering policy. | See table below. |
| `validation` | `SignalomeValidationConfig` | Default: `SignalomeValidationConfig()` | Score, localisation, reference-context, and quantitative-meaning guardrails. | Production requires localisation evidence with minimum probability `0.75`. |
| `output` | `SignalomeOutputConfig` | Default: `SignalomeOutputConfig()` | Network threshold and paired-observation settings. | Production requires effective `network_min_paired_finite_observations >= 5`. |
| `performance` | `SignalomePerformanceConfig` | Default: `SignalomePerformanceConfig()` | Scale guardrails. | See table below. |
| `mode` | `"production" | "exploratory_compatibility"` | Default: `"production"` | Selects production or named compatibility behavior. | Use `SignalomeConfig.compatibility()` only for explicit exploratory compatibility runs. |

Useful constructors:

- `SignalomeConfig.production()` returns the recommended production policy.
- `SignalomeConfig.strict()` keeps strict score preconditioning.
- `SignalomeConfig.permissive_missing_scores()` allows and reports all-missing
  score-row drops.
- `SignalomeConfig.sampled_candidate_scoring()` samples candidate module-count
  scoring while keeping scale guardrails.
- `SignalomeConfig.compatibility()` opts into exploratory compatibility behavior.

`SignalomeScientificConfig`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `substrate_support_cutoff` | `float` | Default: `0.5` | Prediction support cutoff for assigning kinase-supported substrates. | Must be in `[0, 1]`. |
| `assignment_policy` | `"cutoff_binary" | "weighted_top"` | Default: `"cutoff_binary"` | How site-to-kinase support contributes to module summaries. | `"weighted_top"` uses top-kinase weights when available. |

`SignalomeClusteringConfig`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `module_count` | `int | None` | Default: `None` | Requested number of modules. | `None` enables automatic selection; an integer must be at least `1`. |
| `module_selection_primary_correlation_threshold` | `float` | Default: `0.5` | Primary threshold used by automatic module-count selection. | Must be in `[0, 1]`. |
| `module_selection_fallback_correlation_threshold` | `float` | Default: `0.1` | Fallback threshold for automatic module-count selection. | Must be in `[0, 1]`. |
| `module_selection_max_clusters` | `int` | Default: `10` | Maximum cluster count evaluated by automatic selection. | Must be at least `1`. |
| `candidate_scoring_policy` | `"full" | "sampled"` | Default: `"full"` | Candidate module-count scoring strategy. | `"sampled"` approximates candidate scoring only. |
| `clustering_engine` | `"scipy_hierarchical" | "exact_python"` | Default: `"scipy_hierarchical"` | Tree-building engine. | Both use the public signalome clustering contract. |

`SignalomeValidationConfig`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `score_preconditioning_policy` | `"error_on_drop" | "allow_and_report"` | Default: `"error_on_drop"` | Controls all-missing score-row handling. | `"allow_and_report"` drops all-missing score rows and reports counts. |
| `localisation_requirement` | `LocalisationRequirement` | Default: `LocalisationRequirement.production_site_level()` | Site-level localisation requirement. | Production requires present localisation with minimum probability `0.75`. |
| `allow_mixed_total_protein_quantitative_meaning` | `bool` | Default: `False` | Allows mixed total-protein quantitative meaning. | Keep `False` unless intentionally documented. |
| `reference_context_compatibility_policy` | `"require_known_match" | "allow_unknown_with_caveat"` | Default: `"require_known_match"` | Controls unknown dataset/reference context. | Allowing unknown context records a caveat. |

`LocalisationRequirement`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `require_present` | `bool` | Default from production preset: `True` | Requires localisation evidence to be present. | Missing localisation fails when `True`. |
| `minimum_probability` | `float | None` | Default from production preset: `0.75` | Minimum acceptable localisation probability. | If set, must be in `[0, 1]`. |

`SignalomeOutputConfig`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `network_correlation_threshold` | `float` | Default: `0.5` | Threshold for retaining kinase score-profile association edges. | Must be in `[0, 1]`. |
| `network_policy` | `"signed" | "positive_only" | "absolute_threshold"` | Default: `"signed"` | How correlation sign is handled. | Controls edge filtering and stored correlation values. |
| `network_min_paired_finite_observations` | `int | None` | Default: `5` | Minimum paired finite observations for a kinase pair. | Public floor is `3`; production requires the effective value to be at least `5`. `None` resolves to the default. |

`SignalomePerformanceConfig`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `max_exact_tree_sites` | `int` | Default: `2000` | Scale guardrail for exact tree construction. | Must meet the configured floor. |
| `max_full_candidate_scoring_sites` | `int` | Default: `2000` | Scale guardrail for full candidate scoring. | Must meet the configured floor. |

## Running the workflow

Call `SignalomeWorkflow().run(request)`.

```python
from phospy import SignalomeWorkflow
from phospy.api import SignalomeWorkflowRequest

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=SignalomeConfig.production(),
    )
)
```

For the same upstream kinase result, config, and package version, the production
signalome path is deterministic. If the upstream kinase result came from an
adaptive seeded prediction run, that upstream seed affects the signalome input.

The workflow can raise `WorkflowValidationError` before table construction when
the upstream result is missing required matrices, site keys do not align,
kinase columns do not overlap, `protein_group_id` is missing, localisation
requirements fail, all score rows would be dropped under a strict policy, or
network observation settings are below the supported floor.

## Response model and output formats

`SignalomeWorkflow.run(...)` returns `SignalomeWorkflowResult`.

Top-level result:

| Attribute or helper | Python type | Always present? | Meaning |
| --- | --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | Yes | Dataset from the upstream kinase result. |
| `kinase_result` | `KinaseWorkflowResult` | Yes | Upstream kinase result. |
| `module_assignments` | `SignalomeAssignments` | Yes | Site-level module assignment table. |
| `signalome_modules` | `SignalomeModules` | Yes | Module-by-kinase percentage table. |
| `kinase_network` | `KinaseNetwork` | Yes | Network-style edge/node/correlation-diagnostic tables. |
| `module_selection_diagnostics` | object | Yes | Module-count selection diagnostics. |
| `clustering_preparation_diagnostics` | object | Yes | Score-matrix preparation diagnostics. |
| `score_preconditioning_diagnostics` | object | Yes | Score preconditioning diagnostics. |
| `alignment_diagnostics` | object | Yes | Dataset/prediction/score alignment diagnostics. |
| `provenance` | object or mapping | Yes | Workflow provenance. |
| `caveats` | `tuple` | Yes | Structured caveats and warnings. |
| `expanded_signalome` | `pandas.DataFrame | None` | Optional | Flattened table for reporting. |
| `to_dataframe()` | `pandas.DataFrame | None` | Optional | Defensive snapshot of `expanded_signalome`. |
| `site_membership` | `pandas.DataFrame | None` | Optional | Site-membership sidecar. |
| `site_membership_dataframe()` | `pandas.DataFrame | None` | Optional | Defensive snapshot of the site-membership sidecar. |
| `protein_site_context` | `pandas.DataFrame | None` | Optional | Protein-site context sidecar. |
| `protein_site_context_dataframe()` | `pandas.DataFrame | None` | Optional | Defensive snapshot of protein-site context. |
| `input_dataset_preprocessing_report` | object or `None` | Optional | Preprocessing report carried from the upstream dataset. |
| `scientifically_equals(...)` | `bool` | Yes | Comparison helper for scientific/result equivalence checks. |

`module_assignments.table`:

| Column or index | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| index | Site identity. | `site_key` string labels | Yes |
| `site_key` | Same site identity as the index. | string | Yes |
| `display_id` | Human-readable site label. | string | Yes |
| `gene_symbol` | Display gene symbol. | string | Yes |
| `site` | Residue-position label. | string | Yes |
| `protein_group_id` | Signalome protein grouping identity. | string | Yes |
| `protein_accession` | Protein accession when available. | string or missing | Yes |
| `isoform_id` | Isoform identifier when available. | string or missing | Yes |
| `module_id` | Assigned module. | non-negative integer-compatible label | Yes |
| `top_kinase` | Top supported kinase candidate for the site. | string | Yes |
| `top_score` | Top kinase support score. | finite numeric except unsupported rows | Yes |
| `top_kinase_candidates` | Candidate top kinases. | sequence/list-like | Yes |
| `top_kinase_weights` | Candidate weights. | sequence/list-like pairs | Yes |
| `top_kinase_tie_count` | Number of tied top candidates. | non-negative integer-like | Yes |
| `top_kinase_is_ambiguous` | Whether top-kinase selection was ambiguous. | boolean-like | Yes |
| `top_kinase_selection_policy` | Selection policy used for site-level top kinase. | string | Yes |
| `module_top_kinase` | Top kinase summarized across the module. | string | Yes |
| `module_top_kinase_candidates` | Candidate module-top kinases. | sequence/list-like | Yes |
| `module_top_kinase_tie_count` | Number of tied module-top candidates. | non-negative integer-like | Yes |
| `module_top_kinase_is_ambiguous` | Whether module-top selection was ambiguous. | boolean-like | Yes |
| `module_top_kinase_selection_policy` | Selection policy used for module-top kinase. | string | Yes |

`signalome_modules.table`:

| Column or index | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| index | Module identity. | non-negative integer-compatible `module_id` labels | Yes |
| columns | Kinase labels. | strings | Yes |
| cells | Percent support share for the kinase within the module. | finite numeric values from `0` to `100` | Yes |

Rows summarize within-module kinase support. A row total is approximately `100`
when support is assigned, or `0` when no support is available.

`kinase_network.edges`:

| Column | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| `source_kinase` | First kinase label. | string | Yes |
| `target_kinase` | Second kinase label. | string | Yes |
| `correlation` | Pairwise score-profile association after network policy. | numeric in `[-1, 1]` | Yes |
| `valid_observations` | Paired finite site count used for the correlation. | non-negative integer-like | Yes |

Edges are association summaries. `source_kinase` and `target_kinase` are table
labels, not inferred direction.

`kinase_network.nodes`, when present:

| Column or index | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| index | Kinase label. | string | Yes |
| `degree` | Number of retained network edges for the kinase. | non-negative integer-like | Yes |
| `n_substrates` | Number of predicted substrates above support cutoff. | non-negative integer-like | Yes |

`kinase_network.candidate_correlations`, when present:

| Column | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| `source_kinase` | First kinase label. | string | Yes |
| `target_kinase` | Second kinase label. | string | Yes |
| `correlation` | Candidate pairwise correlation. | finite numeric in `[-1, 1]` only when status is `finite`; otherwise missing | Yes |
| `correlation_status` | Candidate correlation status. | `finite`, `constant_profile`, `insufficient_observations`, `missing_values`, `non_finite_values`, or `undefined` | Yes |
| `valid_observations` | Paired finite observation count. | non-negative integer-like | Yes |
| `correlation_reason` | User-facing reason for the status. | string | Yes |

`site_membership`, when present:

| Column | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| `site_key` | Site identity. | string | Yes |
| `display_id` | Human-readable site label. | string | Yes |
| `site_id` | Site identifier used by the sidecar. | string | Yes |
| `site` | Residue-position label. | string | Yes |
| `protein_group_id` | Signalome grouping identity. | string | Yes |
| `protein_accession` | Protein accession when available. | string or missing | Yes |
| `isoform_id` | Isoform identifier when available. | string or missing | Yes |
| `site_cluster` | Site cluster label. | string or integer-like | Yes |
| `protein_module_id` | Protein-level module label. | integer-like or missing | Yes |
| `included_in_module_table` | Whether the site contributes to module output. | boolean-like | Yes |
| `excluded_reason` | Reason for exclusion, if any. | string or missing | Yes |
| `gene_symbol` | Display gene symbol. | string | Yes |
| `top_kinase` | Top kinase label for the site. | string | Yes |
| `top_kinase_score` | Top kinase support score. | numeric or missing | Yes |
| `top_kinase_weight` | Top kinase weight. | numeric or missing | Yes |
| `n_supported_kinases` | Number of supported kinases for the site. | non-negative integer-like | Yes |

`protein_site_context`, when present:

| Column | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| `protein_group_id` | Signalome grouping identity. | string | Yes |
| `n_sites` | Number of sites for the protein group. | non-negative integer-like | Yes |
| `site_ids` | Site identifiers in the group. | sequence/list-like | Yes |
| `site_keys` | `site_key` values in the group. | sequence/list-like | Yes |
| `display_ids` | Display labels in the group. | sequence/list-like | Yes |
| `site_clusters` | Site cluster labels in the group. | sequence/list-like | Yes |
| `n_distinct_site_clusters` | Distinct site-cluster count. | non-negative integer-like | Yes |
| `protein_module_id` | Protein module label. | integer-like or missing | Yes |
| `multi_site_protein` | Whether more than one site is present. | boolean-like | Yes |
| `ambiguous_module_context` | Whether the protein has ambiguous module context. | boolean-like | Yes |
| `gene_symbol` | Display gene symbol. | string | Yes |
| `site` | Site label summary. | string or sequence-like | Yes |
| `protein_accession` | Protein accession summary. | string or sequence-like or missing | Yes |
| `isoform_id` | Isoform identifier summary. | string or sequence-like or missing | Yes |
| `top_kinases_by_site` | Top kinases by site. | mapping/list-like | Yes |
| `module_ids_by_site` | Module labels by site. | mapping/list-like | Yes |
| `site_key_to_display_id` | Mapping from `site_key` to display label. | mapping-like | Yes |

`expanded_signalome`, when present, is a flattened reporting table. It always
contains at least `site_key` and `display_id`; rows labelled as site rows by
`row_kind == "site"` have non-empty values for both fields. The remaining
columns depend on the materialized signalome output and should be treated as a
reporting view rather than the primary result contract.

## Interpreting the result

Signalome aligns the dataset, prediction matrix, and downstream score matrix by
`site_key`. Missing rows after alignment mean the site did not pass the upstream
or signalome eligibility path; absence from output does not prove biological
absence.

`module_id` labels are score-derived groups for the current dataset and config.
They are not stable biological pathway identifiers. `top_kinase` and
`module_top_kinase` are highest-supported labels under the selected policy; they
are not causal mechanism claims.

Network correlations summarize kinase score-profile associations across paired
finite observations. A missing correlation means the pair was not estimable or
did not pass filtering. A value near `0.0` means a finite near-zero association
was estimated. Edges do not establish causality, physical interaction, or
direction.

Automatic module-count selection diagnostics are descriptive. Stability
frequencies and similarities are not p-values, confidence probabilities, or
biological validation.

## Common problems

| Problem | What to check |
| --- | --- |
| Missing `protein_group_id` | Add non-empty `site_metadata.protein_group_id` before running kinase and signalome. Do not rely on gene symbols or display labels. |
| Localisation failure | Build the upstream dataset with localisation metadata and `DatasetLocalisationConfig`; production signalome fails when localisation is missing. |
| No overlapping rows | Confirm `prediction_result.pred_mat` and the authoritative scoring matrix share `site_key` rows. |
| No overlapping kinases | Confirm prediction and score matrices use matching kinase labels. |
| All-missing score rows | Use strict mode to fail, or `score_preconditioning_policy="allow_and_report"` only when dropping all-missing rows is acceptable and reported. |
| Empty network edges | Check `network_correlation_threshold`, `network_policy`, paired finite observations, constant profiles, and candidate-correlation diagnostics. |
| Too few retained sites for production | Production network settings require at least five paired finite observations for retained edges. Use more data or explicitly choose compatibility mode for exploratory demonstrations. |

## Related documentation

- [Kinase analysis](kinase.md)
- [Preparing a dataset](dataset-build-workflow.md)
- [Reference data](../reference_bundles.md)
- [Scientific interpretation and limitations](../scientific-interpretation.md)
