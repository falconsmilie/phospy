# Enrichment

Enrichment asks whether a selected list of identifiers overlaps a collection
of biological sets more often than expected. `EnrichmentWorkflow` performs
offline over-representation analysis (ORA) against an explicit background
universe.

This walkthrough uses gene symbols and two small in-memory gene sets. The data
are illustrative and are not a biological benchmark. Run the Python blocks in
order in one fresh process.

## Define the Set Collection

A set collection maps stable term identifiers to their member identifiers. All
members use one namespace. Here that namespace is `gene_symbol`.

```python
from phospy.api import GeneSetCollection

set_collection = GeneSetCollection(
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
```

In a study, use a locally curated collection whose version and licence are
appropriate for the analysis. PhosPy does not fetch or bundle GO, KEGG,
Reactome, PTM-SEA, Enrichr, or gseapy resources.

## Choose the Selected Identifiers

The foreground is the list whose set overlap you want to test. These two gene
symbols might represent a study-defined hit list.

```python
selected_identifiers = ("AKT1", "MAPK1")
```

The background represents identifiers that could reasonably have been
selected in the experiment. It defines the reference population for the
hypergeometric test and must be supplied explicitly.

```python
background_universe = ("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2")
```

Selected identifiers, background identifiers, and set members must share the
same namespace. PhosPy does not translate gene symbols, accessions,
`display_id` labels, or `site_key` values.

## Configure the Test

ORA is the supported method. Benjamini-Hochberg correction adjusts the raw
*p* values across the terms tested in this run.

```python
from phospy.api import EnrichmentConfig

enrichment_config = EnrichmentConfig(
    method="over_representation",
    multiple_testing_correction="benjamini_hochberg",
)
```

## Construct the Request

`identifier_column` records the selected identifier field. It is also the
column name used when the alternative `input_table` route is selected.

```python
from phospy.api import EnrichmentWorkflowRequest

enrichment_request = EnrichmentWorkflowRequest(
    identifier_column="gene_symbol",
    identifier_kind="gene_symbol",
    set_collection=set_collection,
    selected_identifiers=selected_identifiers,
    background_universe=background_universe,
    config=enrichment_config,
)
```

## Run the Workflow

```python
from phospy.api import EnrichmentWorkflow

enrichment_workflow = EnrichmentWorkflow()
enrichment_result = enrichment_workflow.run(enrichment_request)
```

## Inspect and Interpret the Result

The result contains one row per tested term. The table accessor returns an
independent pandas snapshot.

```python
result_columns = [
    "term_id",
    "input_overlap_count",
    "background_overlap_count",
    "p_value",
    "adjusted_p_value",
]
result_view = enrichment_result.table.loc[:, result_columns]

print("Displayed result shape:", result_view.shape)
print(result_view)
```

The displayed shape is `(2, 5)`. `kinase_response` overlaps both selected
genes and has a raw *p* value of `0.3`; `cell_cycle` overlaps one and has a raw
*p* value of `0.9`. These values show how ORA behaves on a tiny universe. They
do not support a biological conclusion.

A smaller raw or adjusted *p* value means stronger evidence under the stated
foreground, background, and set collection. It is statistical evidence, not an
effect size. ORA does not prove pathway activation, regulation, kinase
activity, or causality.

## Choose the Identifier Level

Use one level consistently across the selected list, background, and
collection:

| Identifier Kind | Collection |
| --- | --- |
| `"gene_symbol"` or `"protein_id"` | `GeneSetCollection` |
| `"site_key"`, `"display_id"`, or `"phosphosite"` | `PtmSetCollection` |

For PhosPy-derived site-level analyses, prefer `site_key`. Provide selected
identifiers through either `selected_identifiers` or `input_table` with
`identifier_column`, never both.

PhosPy does not implement GSEA, ssGSEA, or PTM-SEA. ORA does not imply support
for those methods. Workflow-created provenance records that execution was
offline and used no online-resource lookup.

`EnrichmentWorkflow` is public through `phospy.api`; it is not exported from
the package root.

## Request

Create an `EnrichmentWorkflowRequest`.

| Parameter | Type | Required or Default | Description | Main Constraint |
| --- | --- | --- | --- | --- |
| `identifier_column` | `str` | Required | Column used when `input_table` supplies the selected identifiers. It also records the selected identifier field name. | Must be non-empty and present in `input_table` when that source is used. |
| `identifier_kind` | Identifier-kind string | Required | Namespace shared by selected identifiers, background, and sets. | Must match the collection level. |
| `set_collection` | `GeneSetCollection` or `PtmSetCollection` | Required | Caller-supplied local sets. | Must be non-empty and use the same identifier kind. |
| `background_universe` | `Sequence[str]` | Required | Universe used by the hypergeometric test. | Must be non-empty. Values are trimmed and de-duplicated. |
| `config` | `EnrichmentConfig` | `EnrichmentConfig()` | ORA, correction, attrition, and size-filter policy. | Only supported policies are accepted. |
| `input_table` | `pandas.DataFrame` or `None` | `None` | Optional table containing selected identifiers. | Mutually exclusive with `selected_identifiers`. |
| `selected_identifiers` | `Sequence[str]` or `None` | `None` | Optional explicit selected list. | Mutually exclusive with `input_table`; must be non-empty when used. |
| `selected_identifier_provenance` | `EnrichmentIdentifierSetProvenance` or `None` | `None` | Optional typed provenance for the selected set. | `identifier_count` must match the normalised selected count. |
| `background_identifier_provenance` | `EnrichmentIdentifierSetProvenance` or `None` | `None` | Optional typed provenance for the background. | `identifier_count` must match the normalised background count. |

### `EnrichmentConfig`

| Parameter | Default | Description |
| --- | --- | --- |
| `method` | `"over_representation"` | ORA is the only supported method. |
| `multiple_testing_correction` | `"benjamini_hochberg"` | Supports `"none"`, `"benjamini_hochberg"`, `"bonferroni"`, `"holm"`, and `"benjamini_yekutieli"`. |
| `min_set_size` | `None` | Optional lower set-size bound after background intersection. |
| `max_set_size` | `None` | Optional upper set-size bound after background intersection. |
| `selected_outside_background_policy` | `"error"` | Use `"drop"` only when selected/background intersection is intentional and will be reviewed. |
| `set_member_outside_background_policy` | `"drop"` | Use `"error"` to require every set member to occur in the background. |
| `minimum_retained_foreground_fraction` | `None` | Optional lower bound on the selected fraction retained after background intersection. |

`EnrichmentConfig.publishing()` supplies default set-size bounds of 5 and 500
unless you override them.

<details markdown="1">
<summary><strong>Set Collection Parameters</strong></summary>

`GeneSetCollection` and `PtmSetCollection` accept:

| Parameter | Description |
| --- | --- |
| `sets` | Mapping from term ID to identifiers, or explicit `EnrichmentSet` records. |
| `identifier_kind` | Namespace used by every member. Defaults to `"gene_symbol"` for gene sets and `"site_key"` for PTM sets. |
| `term_names` | Optional readable names keyed by term ID. |
| `source_name`, `source_version` | Optional source identity recorded with the collection. |
| `descriptions` | Optional term descriptions keyed by term ID. |

Set IDs must be unique and non-empty. Duplicate members within a set are
collapsed. `term_names` and `descriptions` may refer only to existing set IDs.

Optional `EnrichmentIdentifierSetProvenance` records where a selected or
background set came from, including an upstream workflow, selection rule,
identifier count, quantitative context, caveats, and upstream provenance
fingerprints.

</details>

## Execution Behaviour and Errors

The primary journey above shows the named workflow and request boundary. After
running it, the same request can also be passed through the concise public form:

```python
repeated_result = EnrichmentWorkflow().run(enrichment_request)
```

Execution is deterministic for the same normalised inputs and configuration.
Invalid local config values raise `ContractValidationError`; cross-input
problems raise `WorkflowValidationError` when the workflow runs.

Typical failures include an empty background, both selected-input routes being
set, mismatched identifier kinds, selected identifiers outside the background
under the `"error"` policy, or provenance counts that do not match normalised
identifiers.

## Response

`EnrichmentWorkflow.run(...)` returns an `EnrichmentWorkflowResult`.

| Attribute or Helper | Format | Meaning |
| --- | --- | --- |
| `table`, `result_table`, `to_dataframe()` | `pandas.DataFrame` | Independent snapshots of the sorted result table. |
| `records` | `tuple[EnrichmentResultRecord, ...]` | One typed record per tested term; may be empty. |
| `identifier_kind` | `str` | Namespace used for the run. |
| `set_collection`, `config` | Typed objects | Resolved collection and configuration. |
| `unmatched_identifiers` | `tuple[str, ...]` | Retained selected identifiers that overlap no tested set. |
| `warnings`, `caveats` | Tuples | User-facing warnings and structured scientific limits. |
| `diagnostics` | Mapping | Foreground/background, set filtering, and correction diagnostics. |
| `method_metadata` | Mapping | Resolved method name and method-specific execution metadata. |
| `background_summary`, `set_collection_summary` | Mappings | Compact input and attrition summaries. |
| `selected_identifier_provenance`, `background_identifier_provenance` | Typed provenance or `None` | Optional upstream identifier-set evidence. |
| `provenance` | `RunProvenance` or `None` | Workflow parameters, input fingerprints, and limitations. |

Provenance records the offline and no-online-resource policy for a workflow-created
result.

### Result Table Format

Each row is one tested term. Rows are ordered by raw `p_value`, then by
`term_id` when values tie.

| Column | Meaning |
| --- | --- |
| `term_id`, `term_name` | Set identity and readable name. |
| `collection_kind` | `"gene_set"` or `"ptm_set"`. |
| `identifier_kind` | Namespace used by the test. |
| `input_overlap_count` | Selected identifiers overlapping the term after background filtering. |
| `background_overlap_count` | Term members retained in the background. |
| `set_size` | Original set size before background intersection. |
| `overlap_identifiers` | Sorted selected identifiers in the overlap. |
| `p_value` | Hypergeometric ORA *p* value. |
| `adjusted_p_value` | Multiple-testing-adjusted *p* value. |
| `correction_method` | Correction applied to the tested terms. |
| `enrichment_ratio` | Observed overlap divided by expected overlap; `None` when undefined. |

When set-size filtering removes every set, the table is empty but keeps the same
columns. Returned DataFrames are independent snapshots.

<details markdown="1">
<summary><strong>Key Diagnostic Groups</strong></summary>

`diagnostics["foreground_background"]` reports normalised selected and
background counts, retained foreground fraction, outside-background identifiers,
resolved policies, and tested/dropped set counts.

`diagnostics["multiple_testing_correction"]` records the method and number of
tested terms. When size bounds are active, `diagnostics["set_size_filter"]`
records the bounds, dropped terms, and reasons.

</details>

## Interpret the Result

ORA asks whether the selected identifiers overlap a set more than expected if
they were drawn from the stated background universe.

A smaller `p_value` or `adjusted_p_value` indicates stronger evidence under that
model. A larger `enrichment_ratio` indicates more overlap than expected. These
values do not prove pathway activation, regulation, kinase activity, or
causality.

ORA does not imply GSEA, ssGSEA, or PTM-SEA support.

Results are sensitive to the background, namespace, set curation, and foreground
selection rule. A missing term may have been filtered or may have had no usable
overlap; it is not evidence of biological absence.

Gene-level and PTM-level results are not interchangeable. Do not turn a
gene-level overlap into a phosphosite-level claim without an explicit and
scientifically justified transformation.

## Common Issues

| Issue | What to Check |
| --- | --- |
| Selected identifiers fall outside the background. | Confirm that both use the same namespace. Choose `"drop"` only deliberately, then review the retained fraction. |
| The request has two selected sources. | Provide either `selected_identifiers` or `input_table`, not both. |
| A table input is rejected. | Pass a non-empty pandas DataFrame containing `identifier_column`. |
| Collection and identifiers do not match. | Use a gene collection for gene/protein kinds and a PTM collection for site kinds. |
| The result is empty. | Review set-size filters and `diagnostics["set_size_filter"]`. |
| Enrichment is unexpectedly weak. | Revisit the background universe, selection rule, and set curation. |

## Complete Runnable Example

```python
from phospy.api import GeneSetCollection

set_collection = GeneSetCollection(
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

selected_identifiers = ("AKT1", "MAPK1")

background_universe = ("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2")

from phospy.api import EnrichmentConfig

enrichment_config = EnrichmentConfig(
    method="over_representation",
    multiple_testing_correction="benjamini_hochberg",
)

from phospy.api import EnrichmentWorkflowRequest

enrichment_request = EnrichmentWorkflowRequest(
    identifier_column="gene_symbol",
    identifier_kind="gene_symbol",
    set_collection=set_collection,
    selected_identifiers=selected_identifiers,
    background_universe=background_universe,
    config=enrichment_config,
)

from phospy.api import EnrichmentWorkflow

enrichment_workflow = EnrichmentWorkflow()
enrichment_result = enrichment_workflow.run(enrichment_request)

result_columns = [
    "term_id",
    "input_overlap_count",
    "background_overlap_count",
    "p_value",
    "adjusted_p_value",
]
result_view = enrichment_result.table.loc[:, result_columns]

print("Displayed result shape:", result_view.shape)
print(result_view)
```

## Continue

- [Differential Analysis](differential-analysis.md)
- [Scientific Interpretation and Limitations](../scientific-interpretation.md)
- [Reference Data](../reference_bundles.md)
- [Public Python API](guide.md)
