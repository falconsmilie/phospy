# Enrichment Workflow

## Plain-language introduction

`EnrichmentWorkflow` runs offline over-representation analysis (ORA) with
identifier sets that you provide. Use it when you already have a selected list
of genes, proteins, or phosphosites and want to ask whether that list overlaps
local enrichment sets more than expected under an explicit background universe.

The workflow expects one selected-identifier source, one homogeneous set
collection, and one explicit background universe. It returns an
`EnrichmentWorkflowResult` with one row per tested term, ORA p-values,
multiple-testing-adjusted p-values, overlap diagnostics, caveats, and run
provenance.

This workflow does not fetch GO, KEGG, Reactome, PTM-SEA, Enrichr, gseapy, or
other online resources. It does not implement ranked-list enrichment, GSEA,
ssGSEA, or PTM-SEA.

## Input and dataset requirements

Enrichment does not take an `AnalysisReadyPhosphoDataset` directly. It works on
identifiers that you pass either as a sequence or in a pandas DataFrame column.
If those identifiers came from a PhosPy dataset or workflow result, the upstream
workflow owns the dataset requirements.

Workflow-specific requirements:

- Choose one identifier level: gene/protein identifiers or PTM/phosphosite
  identifiers.
- Set `identifier_kind` to one of `"gene_symbol"`, `"protein_id"`,
  `"site_key"`, `"display_id"`, or `"phosphosite"`.
- Use a matching collection:
  - `GeneSetCollection` for `"gene_symbol"` or `"protein_id"`.
  - `PtmSetCollection` for `"site_key"`, `"display_id"`, or `"phosphosite"`.
- Provide exactly one selected-identifier source:
  - `selected_identifiers`, or
  - `input_table` plus `identifier_column`.
- Provide a non-empty `background_universe`. PhosPy does not infer it from a
  dataset, reference bundle, or set collection.
- Keep selected identifiers, background identifiers, and set members in the
  same namespace. PhosPy does not map between gene symbols, accessions,
  `display_id` labels, and `site_key` values.
- For PTM-level enrichment derived from PhosPy data, prefer analysis-ready
  `site_key` identifiers. Upstream dataset construction requires
  `site_sequence`, but enrichment itself does not inspect site sequences.
- ORA uses selected/background membership only. There is no replicate,
  condition, or intensity-scale requirement unless you provide optional
  provenance for an identifier set derived from a quantitative PhosPy result.

The main scientific assumption is that your background universe represents the
identifiers that could have been selected. A poor background can dominate the
result.

For shared dataset preparation, see
[Preparing a dataset](dataset-build-workflow.md).

## Minimal end-to-end example

```python
from phospy.api import (
    EnrichmentConfig,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
)

collection = GeneSetCollection(
    sets={
        "kinase_response": ("AKT1", "MAPK1", "MTOR"),
        "cell_cycle": ("CDK1", "CDK2", "MAPK1"),
    },
    identifier_kind="gene_symbol",
    term_names={
        "kinase_response": "Kinase response",
        "cell_cycle": "Cell cycle",
    },
    source_name="example in-memory gene sets",
    source_version="2026-06",
)

request = EnrichmentWorkflowRequest(
    identifier_column="gene_symbol",
    identifier_kind="gene_symbol",
    set_collection=collection,
    selected_identifiers=("AKT1", "MAPK1"),
    background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
    config=EnrichmentConfig(
        method="over_representation",
        multiple_testing_correction="benjamini_hochberg",
    ),
)

result = EnrichmentWorkflow().run(request)
print(result.table.loc[:, ["term_id", "input_overlap_count", "p_value"]])
```

`EnrichmentWorkflow` is public through `phospy.api`; it is not a top-level
`phospy` convenience export.

## Request model

Use `EnrichmentWorkflowRequest`.

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `identifier_column` | `str` | Required | Column to read when `input_table` is the selected-identifier source. Also records the identifier column name in request semantics. | Must be non-empty after trimming. The column must exist when `input_table` is supplied. |
| `identifier_kind` | `Literal["gene_symbol", "protein_id", "site_key", "display_id", "phosphosite"]` | Required | Declares the identifier namespace used by selected identifiers, the background universe, and the set collection. | Must match `set_collection.identifier_kind`. Gene-level kinds use gene sets; PTM-level kinds use PTM sets. |
| `set_collection` | `EnrichmentSetCollection` | Required | Local enrichment set collection supplied by the caller. `GeneSetCollection` and `PtmSetCollection` are the usual constructors. | Must be non-empty, homogeneous, and match `identifier_kind`. Set IDs must be unique. Set names and members must be non-empty strings. |
| `background_universe` | `Sequence[str]` | Required | Explicit universe for the hypergeometric ORA test. | Missing the argument raises `TypeError`. An empty sequence is rejected by `EnrichmentWorkflow.run(...)`. Values are trimmed and de-duplicated for execution. |
| `config` | `EnrichmentConfig` | Default: `EnrichmentConfig()` | Method, multiple-testing, background-policy, and set-size filtering options. | Must be an `EnrichmentConfig`. |
| `input_table` | `pandas.DataFrame \| None` | Default: `None` | Optional table containing selected identifiers in `identifier_column`. | Mutually exclusive with `selected_identifiers`. The current workflow path accepts a non-empty pandas DataFrame with the required column. |
| `selected_identifiers` | `Sequence[str] \| None` | Default: `None` | Optional explicit selected identifiers. | Mutually exclusive with `input_table`. Must be non-empty when used. Values are trimmed and de-duplicated for execution. |
| `selected_identifier_provenance` | `EnrichmentIdentifierSetProvenance \| None` | Default: `None` | Optional typed provenance for the selected identifier set. | If supplied, `identifier_count` must match the normalized selected identifier count. |
| `background_identifier_provenance` | `EnrichmentIdentifierSetProvenance \| None` | Default: `None` | Optional typed provenance for the background universe. | If supplied, `identifier_count` must match the normalized background identifier count. |

### `EnrichmentConfig`

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `method` | `Literal["over_representation"]` | Default: `"over_representation"` | Enrichment method. | ORA is the only supported method. |
| `multiple_testing_correction` | `Literal["benjamini_hochberg", "bonferroni", "holm", "benjamini_yekutieli", "none"]` | Default: `"benjamini_hochberg"` | Adjusts p-values across tested sets. | Correction is applied after background policies and optional set-size filters. `"none"` returns adjusted p-values equal to the raw p-values. |
| `min_set_size` | `int \| None` | Default: `None` | Drops sets with fewer retained members than this value. | `None` disables the lower bound. Integers must be `>= 1`. Size is measured after intersecting each set with `background_universe`. |
| `max_set_size` | `int \| None` | Default: `None` | Drops sets with more retained members than this value. | `None` disables the upper bound. Integers must be `>= 1`. If both bounds are set, `min_set_size <= max_set_size`. |
| `selected_outside_background_policy` | `Literal["error", "drop"]` | Default: `"error"` | Controls selected identifiers that are absent from the background. | `"error"` raises a validation error. `"drop"` runs ORA on the selected/background intersection and records dropped identifiers. |
| `set_member_outside_background_policy` | `Literal["error", "drop"]` | Default: `"drop"` | Controls set members that are absent from the background. | `"drop"` supports broad local collections tested against an experiment-specific background. `"error"` requires every set member to be inside the background. |
| `minimum_retained_foreground_fraction` | `float \| None` | Default: `None` | Optional guard on selected-identifier attrition when dropping outside-background identifiers. | `None` disables the guard. Values must be within `[0.0, 1.0]`. The workflow fails if the retained selected fraction is lower. |

`EnrichmentConfig.publishing()` returns
`EnrichmentConfig(min_set_size=5, max_set_size=500, ...)` unless you override
those bounds.

### Set collections

Use `GeneSetCollection` for gene/protein-level ORA and `PtmSetCollection` for
phosphosite/PTM-level ORA.

| Constructor | Main parameters | Defaults | Constraints |
| --- | --- | --- | --- |
| `GeneSetCollection` | `sets`, `identifier_kind`, `term_names`, `source_name`, `source_version`, `descriptions` | `identifier_kind="gene_symbol"`, `source_name="user"` | `identifier_kind` must be `"gene_symbol"` or `"protein_id"`. |
| `PtmSetCollection` | `sets`, `identifier_kind`, `term_names`, `source_name`, `source_version`, `descriptions` | `identifier_kind="site_key"`, `source_name="user"` | `identifier_kind` must be `"site_key"`, `"display_id"`, or `"phosphosite"`. |
| `EnrichmentSet` | `set_id`, `name`, `identifiers`, `identifier_kind`, `source_name`, `source_version`, `description` | Optional source and description fields default to `None`. | Use when you want to build an `EnrichmentSetCollection` from explicit set records. |

For mapping-style `sets`, keys are term IDs and values are identifier
sequences. `term_names` and `descriptions` may name only IDs present in `sets`.
Duplicate identifiers within one set are collapsed.

### Optional identifier-set provenance

`EnrichmentIdentifierSetProvenance` is optional. Use it when you want the result
to record where the selected or background identifiers came from.

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `source_type` | `EnrichmentIdentifierSetSourceType` | Required | Source category. | Accepted values are `"manual"`, `"raw_identifier_list"`, and `"phospy_derived_quantitative"`. |
| `source_label` | `str` | Required | Human-readable label for the identifier set source. | Must be non-empty after trimming. |
| `identifier_count` | `int` | Required | Number of normalized identifiers represented by this provenance record. | Must be non-negative and must match the normalized selected or background count used by the workflow. |
| `upstream_workflow_id` | `str \| None` | Default: `None` | Optional upstream workflow identifier. | Must be non-empty if supplied. |
| `upstream_result_id` | `str \| None` | Default: `None` | Optional upstream result identifier. | Must be non-empty if supplied. |
| `input_intensity_scale_evidence` | `InputIntensityScaleEvidence \| None` | Default: `None` | Required only for `"phospy_derived_quantitative"` provenance. | Must be omitted for `"manual"` and `"raw_identifier_list"`. |
| `derived_quantitative_provenance` | `EnrichmentDerivedQuantitativeSetProvenance \| None` | Default: `None` | Required only for `"phospy_derived_quantitative"` provenance. | Must be omitted for `"manual"` and `"raw_identifier_list"`. Its identifier namespace must match `identifier_kind`. |

Manual provenance can be supplied entirely through the stable public API:

```python
from phospy.api import (
    EnrichmentIdentifierSetProvenance,
    EnrichmentIdentifierSetSourceType,
)

selected_provenance = EnrichmentIdentifierSetProvenance(
    source_type=EnrichmentIdentifierSetSourceType.MANUAL,
    source_label="curated hit list",
    identifier_count=2,
)
```

## Running the workflow

Run enrichment with the public workflow object:

```python
from phospy.api import EnrichmentWorkflow

result = EnrichmentWorkflow().run(request)
```

`run(...)` receives one `EnrichmentWorkflowRequest` and returns one
`EnrichmentWorkflowResult`. The ORA calculation is deterministic; there is no
seed parameter. Run provenance records the offline/no-online-resource policy;
online resources are not used.

Before ORA runs, PhosPy validates the request. Common validation failures
include:

- the request is not an `EnrichmentWorkflowRequest`;
- `identifier_kind` is unsupported or does not match the set collection;
- both `input_table` and `selected_identifiers` are supplied, or neither is
  supplied;
- the selected source, background universe, or set collection is empty;
- selected identifiers are outside the background while
  `selected_outside_background_policy="error"`;
- all selected identifiers would be dropped by the background policy;
- retained selected-identifier fraction is below
  `minimum_retained_foreground_fraction`;
- set members are outside the background while
  `set_member_outside_background_policy="error"`;
- optional provenance counts do not match the normalized identifier counts.

These are reported as `WorkflowValidationError` during workflow execution.
Invalid `EnrichmentConfig` construction raises `ContractValidationError`.

## Response model and output formats

`EnrichmentWorkflow.run(...)` returns `EnrichmentWorkflowResult`.

| Attribute or helper | Python type | Always present? | Meaning |
| --- | --- | --- | --- |
| `identifier_kind` | `str` | Yes | Identifier namespace used for the run. |
| `set_collection` | `EnrichmentSetCollection` | Yes | The set collection represented by the result. |
| `config` | `EnrichmentConfig` | Yes | Resolved enrichment configuration. |
| `records` | `tuple[EnrichmentResultRecord, ...]` | Yes; may be empty | One record per tested term. |
| `table` | `pandas.DataFrame` | Yes; may be empty | Defensive snapshot of the result table. |
| `result_table` | `pandas.DataFrame` | Yes; may be empty | Alias for `table`. |
| `to_dataframe()` | `pandas.DataFrame` | Yes; may be empty | Returns the same result table shape as a defensive snapshot. |
| `unmatched_identifiers` | `tuple[str, ...]` | Yes; may be empty | Selected identifiers retained in the background that did not overlap any tested set. |
| `warnings` | `tuple[str, ...]` | Yes; may be empty | User-facing warning strings if present. |
| `caveats` | `tuple[ResultCaveat, ...]` | Yes; may be empty | Structured scientific caveats attached to the result. |
| `diagnostics` | `Mapping[str, object]` | Yes | Execution diagnostics, including foreground/background overlap and multiple-testing status. |
| `method_metadata` | `Mapping[str, object]` | Yes | Resolved method metadata, including ORA and multiple-testing settings. |
| `background_summary` | `Mapping[str, object]` | Yes | Background-universe and selected-source summary. |
| `set_collection_summary` | `Mapping[str, object]` | Yes | Set collection size, source, and attrition summary. |
| `selected_identifier_provenance` | `EnrichmentIdentifierSetProvenance \| None` | Optional | Selected-identifier provenance when supplied on the request. |
| `background_identifier_provenance` | `EnrichmentIdentifierSetProvenance \| None` | Optional | Background-universe provenance when supplied on the request. |
| `provenance` | `RunProvenance \| None` | Present for workflow-created results | Run provenance, including parameters, table fingerprints, offline/no-online-resource policy, and limitations. |

`table`, `result_table`, and `to_dataframe()` return defensive snapshots.
Changing the returned DataFrame does not mutate the result.

### Result table schema

The result table has a default row index. Each row is one tested term. Rows are
ordered by raw `p_value`, then by `term_id` for ties.

| Column or index | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| row index | Positional row number in the sorted result table. | pandas index | Yes |
| `term_id` | Set identifier. | non-empty string | Yes |
| `term_name` | Human-readable set name. | string | Yes |
| `collection_kind` | Set collection level. | `"gene_set"` or `"ptm_set"` | Yes |
| `identifier_kind` | Identifier namespace used by this row. | `"gene_symbol"`, `"protein_id"`, `"site_key"`, `"display_id"`, or `"phosphosite"` | Yes |
| `input_overlap_count` | Number of selected identifiers that overlap this tested set after background filtering. | integer count | Yes |
| `background_overlap_count` | Number of set members retained in the background universe. | integer count | Yes |
| `set_size` | Raw number of identifiers in the set before background intersection. | integer count | Yes |
| `overlap_identifiers` | Sorted identifiers in the selected/set overlap. | tuple of strings | Yes; may be empty |
| `p_value` | Hypergeometric ORA p-value for over-representation. | float in `[0.0, 1.0]` | Yes for workflow-created rows |
| `adjusted_p_value` | Multiple-testing-adjusted p-value. | float in `[0.0, 1.0]` | Yes for workflow-created rows |
| `correction_method` | Multiple-testing method used for this row. | configured correction method | Yes for workflow-created rows |
| `enrichment_ratio` | Observed overlap divided by expected overlap under the background. | non-negative float or `None` | Yes; `None` when the expected overlap is undefined |

If all sets are removed by configured set-size filters, `records` is empty and
`table` is an empty DataFrame with the same columns.

### Main diagnostics and summaries

`diagnostics["foreground_background"]` reports the selected/background
intersection:

| Key | Meaning |
| --- | --- |
| `identifier_kind` | Identifier namespace used for the run. |
| `foreground_size_before_intersection` | Normalized selected-identifier count before background intersection. |
| `background_size` | Normalized background-universe size. |
| `usable_foreground_size_after_background_intersection` | Selected identifiers retained in the background and used by ORA. |
| `retained_foreground_fraction` | Fraction of selected identifiers retained in the background. |
| `foreground_identifiers_missing_from_background_count` | Count of selected identifiers absent from the background. |
| `foreground_identifiers_missing_from_background` | Selected identifiers absent from the background. |
| `selected_outside_background_policy` | Resolved selected-identifier policy. |
| `set_member_outside_background_policy` | Resolved set-member policy. |
| `minimum_retained_foreground_fraction` | Configured retained-foreground guard, when any. |
| `tested_set_count` | Number of sets tested after optional set-size filters. |
| `dropped_set_count` | Number of sets dropped by optional set-size filters. |
| `set_identifiers_missing_from_background_count` | Count of distinct set members absent from the background. |
| `set_identifiers_missing_from_background` | Bounded, sorted preview of absent set members. |
| `set_identifiers_missing_from_background_truncated` | Whether the preview was shortened. |
| `set_identifiers_missing_from_background_preview_limit` | Preview limit. |

`diagnostics["multiple_testing_correction"]` records the correction method,
whether correction was applied, and the number of tested records. When
`min_set_size` or `max_set_size` is configured,
`diagnostics["set_size_filter"]` records the applied bounds, input/tested/drop
counts, drop reasons, and dropped-set summaries.

`background_summary` includes the explicit background source, universe size,
selected identifier count, selected source, retained foreground fraction, and
outside-background policy results. `set_collection_summary` includes collection
kind, identifier kind, set count, member counts, source metadata, empty-set
counts after background filtering, and dropped-set counts when set-size filters
are used.

## Interpreting the result

ORA asks whether the selected identifiers overlap a set more than expected if
the selected identifiers were drawn from the explicit background universe.

Low `p_value` or `adjusted_p_value` values mean the observed overlap is
unlikely under that ORA model. Higher `enrichment_ratio` values mean the
observed overlap is larger than the expected overlap under the same background.
An `enrichment_ratio` of `0.0` means no selected identifiers overlapped that
set. `None` means the ratio is undefined, usually because the selected set or
background-retained set has size zero.

The result does not prove pathway activation, regulation, causality, or kinase
activity. ORA is sensitive to the background universe, identifier namespace,
set curation, and the rule used to select foreground identifiers. Absence from
the output can also mean a set was filtered out or had no overlap in the chosen
background; it is not evidence of biological absence.

Gene-level and PTM-level ORA results are not interchangeable. Do not interpret
a gene-symbol set as phosphosite evidence, and do not collapse PTM-level
results to gene-level claims unless that transformation is part of your
analysis design.

For broader caveats, see
[Scientific interpretation and limitations](../scientific-interpretation.md).

## Common problems

| Problem | Most likely fix |
| --- | --- |
| Selected identifiers are outside the background. | Check that the selected list and background use the same namespace. Use `"drop"` only when selected/background intersection is intended, and review retained-foreground diagnostics. |
| The workflow says exactly one identifier source is required. | Provide either `selected_identifiers` or `input_table`, not both. |
| `input_table` is rejected. | Pass a non-empty pandas DataFrame and make sure `identifier_column` exists. |
| Gene and PTM semantics do not match. | Use `GeneSetCollection` with `"gene_symbol"` or `"protein_id"`; use `PtmSetCollection` with `"site_key"`, `"display_id"`, or `"phosphosite"`. |
| Empty results after filtering. | Review `min_set_size`, `max_set_size`, and `diagnostics["set_size_filter"]`. Set size is measured after background intersection. |
| Unexpectedly weak enrichment. | Revisit the background universe, selected-identifier rule, and set collection. ORA results are only as meaningful as those inputs. |
| Provenance count mismatch. | Set `identifier_count` to the normalized selected or background count after trimming blanks and removing duplicates. |

## Related documentation

- [Preparing a dataset](dataset-build-workflow.md)
- [Differential analysis](differential-analysis.md)
- [Scientific interpretation and limitations](../scientific-interpretation.md)
- [Reference data](../reference_bundles.md)
- [API guide](guide.md)
