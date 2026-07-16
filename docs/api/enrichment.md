# Enrichment Workflow

`EnrichmentWorkflow` runs offline over-representation analysis against
caller-supplied gene-set or PTM-set collections. It is intentionally explicit:
you provide the selected identifiers, the background universe, identifier
semantics, and the set collection.

## When to Use This Workflow

Use this workflow when you already have a set of selected identifiers and want a
simple offline ORA result over local, caller-supplied sets and a caller-supplied
background universe.

Good fits:

- gene-level ORA over `gene_symbol` or `protein_id`
- PTM/site-level ORA over `site_key`, `display_id`, or `phosphosite`
- small local collections created in Python or read from local files

This workflow does not fetch GO, KEGG, Reactome, PTM-SEA, PTMsigDB, Enrichr,
gseapy, clusterProfiler, or other online resources. It does not implement
ranked-list enrichment, GSEA, ssGSEA, or PTM-SEA. A possible future ranked-list
mode is deferred by
[ADR-0030](../adr/adr_0030_ranked_list_enrichment_prerequisites.md).

## Inputs

Provide exactly one identifier source:

- `selected_identifiers`, or
- `input_table` plus `identifier_column`

Also provide:

- an explicit `identifier_kind`
- a homogeneous `GeneSetCollection`, `PtmSetCollection`, or
  `EnrichmentSetCollection`
- an explicit non-empty `background_universe`
- an `EnrichmentConfig`

The selected identifiers, background universe, and set collection must use the
same identifier namespace. Gene-level and PTM-level collections are deliberately
separate. Selected identifiers that are not present in the background fail by
default. Intentional foreground intersection requires
`selected_outside_background_policy="drop"` in `EnrichmentConfig`.

## Request Object

Use `EnrichmentWorkflowRequest`.

Important fields:

| Field | Meaning |
| --- | --- |
| `identifier_column` | Column used when `input_table` is provided. |
| `identifier_kind` | Explicit identifier semantics, such as `"gene_symbol"` or `"site_key"`. |
| `set_collection` | Caller-supplied enrichment collection. |
| `background_universe` | Explicit universe for the hypergeometric test. |
| `config` | `EnrichmentConfig` for method and multiple-testing policy. |
| `input_table` | Optional table source; mutually exclusive with `selected_identifiers`. |
| `selected_identifiers` | Optional explicit identifier list; mutually exclusive with `input_table`. |
| `selected_identifier_provenance` | Optional typed provenance for the selected identifier set. |
| `background_identifier_provenance` | Optional typed provenance for the explicit background set. |

Construction stores the request payload and can coerce provenance mappings into
typed dataclasses. Workflow validation owns identifier-source selection,
collection/background compatibility, identifier-count checks, and provenance
compatibility. Workflow execution calculates enrichment statistics.

## Identifier-Set Provenance

Identifier-set provenance is optional for legacy and ordinary manual enrichment
requests. A request with no selected or background provenance keeps those fields
absent; PhosPy does not silently label omitted provenance as manual.

Use `EnrichmentIdentifierSetProvenance` when you want the result to distinguish
manual/raw identifier lists from identifier sets derived from PhosPy
quantitative workflows. The same typed model is used for both selected and
background sets.

`source_type` is one of:

- `EnrichmentIdentifierSetSourceType.MANUAL`
- `EnrichmentIdentifierSetSourceType.RAW_IDENTIFIER_LIST`
- `EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE`

For `MANUAL` and `RAW_IDENTIFIER_LIST`, intensity-scale evidence is not needed
and should not be supplied. For `PHOSPY_DERIVED_QUANTITATIVE`, typed
`InputIntensityScaleEvidence` is mandatory so downstream provenance preserves
whether the upstream quantitative scale came from an observed transformation or
from a user declaration.

`identifier_count` must match the workflow-normalized identifier count:
distinct, nonblank identifiers after enrichment input normalization and before
reference mapping, annotation-universe overlap, term filtering, statistical
testing, or multiple-testing correction. For the background role, provenance is
valid only when an explicit `background_universe` is supplied.

Manual selected list with optional manual provenance:

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

request = EnrichmentWorkflowRequest(
    identifier_column="gene_symbol",
    identifier_kind="gene_symbol",
    set_collection=collection,
    selected_identifiers=("AKT1", "MAPK1"),
    background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
    selected_identifier_provenance=selected_provenance,
)
```

PhosPy-derived quantitative selected set with observed transformation evidence:

```python
from phospy.api import (
    EnrichmentIdentifierSetProvenance,
    EnrichmentIdentifierSetSourceType,
)
from phospy.provenance.models import InputIntensityScaleEvidence

selected_provenance = EnrichmentIdentifierSetProvenance(
    source_type=EnrichmentIdentifierSetSourceType.PHOSPY_DERIVED_QUANTITATIVE,
    source_label="differential significant genes",
    identifier_count=2,
    upstream_workflow_id="differential-workflow",
    upstream_result_id="contrast-a-vs-b",
    input_intensity_scale_evidence=InputIntensityScaleEvidence(
        input_intensity_scale="log2",
        input_intensity_scale_evidence_level="observed_transformation",
        input_intensity_scale_source="phospy.science.transformations.transformers.log2",
        input_intensity_scale_source_detail="log2_transform",
    ),
)
```

If `input_intensity_scale_evidence_level="declared_by_user"`, the enrichment
result emits a role-specific caveat for the selected or background set. Observed
transformation evidence is recorded in the result and run provenance, but does
not emit the declared-scale caveat.

## Request Configuration

Use `EnrichmentConfig`.

| Field | Default | Supported values |
| --- | --- | --- |
| `method` | `"over_representation"` | `"over_representation"` |
| `multiple_testing_correction` | `"benjamini_hochberg"` | `"benjamini_hochberg"`, `"bonferroni"`, `"holm"`, `"benjamini_yekutieli"`, `"none"` |
| `min_set_size` | `None` | `None` or an integer `>= 1` |
| `max_set_size` | `None` | `None` or an integer `>= 1` |
| `selected_outside_background_policy` | `"error"` | `"error"` or `"drop"` |
| `set_member_outside_background_policy` | `"drop"` | `"error"` or `"drop"` |
| `minimum_retained_foreground_fraction` | `None` | `None` or a number in `[0.0, 1.0]` |

The default multiple-testing correction is `"benjamini_hochberg"`. Available
methods are `"benjamini_hochberg"`, `"bonferroni"`, `"holm"`,
`"benjamini_yekutieli"`, and `"none"`. Correction is applied across the sets
that are actually tested, after the configured universe policy has been applied
and after optional set-size filters drop any sets.

## Universe and Attrition Policy

The selected foreground is conservative by default. If any selected identifier
is outside `background_universe`, workflow validation raises
`WorkflowValidationError`. This catches common namespace and universe mistakes
before they become smaller analyses.

Use `selected_outside_background_policy="drop"` only when the intended analysis
is the selected/background intersection. When drop is enabled, PhosPy records
the outside-background selected identifiers, retained counts, and retained
foreground fraction in diagnostics and provenance. If all selected identifiers
would be dropped, validation fails.

Reference or set-member behavior is controlled separately with
`set_member_outside_background_policy`. The default `"drop"` supports broad
caller-supplied collections tested against an experiment-specific background.
Set `"error"` when every set member must already be inside the background.

Use `minimum_retained_foreground_fraction` to reject excessive foreground
attrition even when selected dropping is enabled. For example, `0.8` requires at
least 80% of normalized selected identifiers to remain in the background.

Related collection classes:

- `EnrichmentSet`
- `EnrichmentSetCollection`
- `GeneSetCollection`
- `PtmSetCollection`

## Set-Size Filters

Use `min_set_size` and `max_set_size` when you want to test only sets within a
specific size range. These filters are optional. When both are `None`, PhosPy
tests the same sets as before.

Set size is measured after each set is intersected with the explicit
`background_universe`. This matters when a local set contains identifiers that
are not in your background. For example, a set with 20 raw members but only 8
members in the background has a filtered set size of 8.

Sets outside the configured range are excluded before ORA p-values are
calculated. Multiple-testing correction then uses only the tested sets.

```python
from phospy.api import EnrichmentConfig, EnrichmentWorkflow, EnrichmentWorkflowRequest

request = EnrichmentWorkflowRequest(
    identifier_column="gene_symbol",
    identifier_kind="gene_symbol",
    set_collection=collection,
    selected_identifiers=("AKT1", "MAPK1"),
    background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
    config=EnrichmentConfig(
        min_set_size=2,
        max_set_size=500,
    ),
)

result = EnrichmentWorkflow().run(request)
print(result.diagnostics["set_size_filter"]["dropped_set_count"])
```

## Running the Workflow

`EnrichmentWorkflow` is exported from `phospy.api`, not from the top-level
`phospy` package.

```python
from phospy.api import EnrichmentWorkflow

result = EnrichmentWorkflow().run(request)
```

## Result Object

`EnrichmentWorkflow.run(...)` returns `EnrichmentWorkflowResult`.

Important fields and helpers:

| Field or helper | Meaning |
| --- | --- |
| `table` | Defensive snapshot of the enrichment table. |
| `result_table` | Alias for the same table snapshot. |
| `to_dataframe()` | Defensive enrichment table snapshot. |
| `records` | Tuple of `EnrichmentResultRecord` values. |
| `unmatched_identifiers` | Selected identifiers not matched into tested sets. |
| `warnings` | User-facing warnings. |
| `diagnostics` | Execution diagnostics. |
| `method_metadata` | Method metadata. |
| `background_summary` | Background-universe summary. |
| `set_collection_summary` | Collection summary. |
| `selected_identifier_provenance` | Typed selected identifier-set provenance, when supplied. |
| `background_identifier_provenance` | Typed background identifier-set provenance, when supplied. |
| `provenance` | Run provenance populated by workflow execution. |

The result table includes one row per tested term with overlap counts, overlap
identifiers, p-values, adjusted p-values, correction method, and enrichment
ratio.

`diagnostics["foreground_background"]` reports how the selected identifiers,
background, and set collection overlap:

| Field | Meaning |
| --- | --- |
| `identifier_kind` | Identifier kind used for the run. |
| `foreground_size_before_intersection` | Number of selected identifiers before intersecting with the background. |
| `background_size` | Number of identifiers in the explicit background. |
| `usable_foreground_size_after_background_intersection` | Number of selected identifiers present in the background and used by ORA. |
| `retained_foreground_fraction` | Fraction of normalized selected identifiers retained in the background. |
| `foreground_identifiers_missing_from_background_count` | Count of selected identifiers absent from the background. |
| `foreground_identifiers_missing_from_background` | The selected identifiers absent from the background. |
| `selected_outside_background_policy` | Resolved selected-identifier outside-background policy. |
| `set_member_outside_background_policy` | Resolved set-member outside-background policy. |
| `minimum_retained_foreground_fraction` | Configured retained-foreground threshold, when any. |
| `tested_set_count` | Number of sets tested after optional set-size filters. |
| `dropped_set_count` | Number of sets dropped by optional set-size filters. |
| `set_identifiers_missing_from_background_count` | Count of distinct set identifiers absent from the background. |
| `set_identifiers_missing_from_background` | Bounded preview of set identifiers absent from the background. |
| `set_identifiers_missing_from_background_truncated` | Whether the set-identifier preview was shortened. |

When set-size filters are configured, `diagnostics["set_size_filter"]` reports:

- `applied_after_background_intersection`
- `min_set_size` and `max_set_size`
- `input_set_count`, `tested_set_count`, and `dropped_set_count`
- `dropped_set_reason_counts`
- `dropped_sets` with each set ID, reason, raw size, and background-overlap size

The same run also records `tested_set_count`, `dropped_set_count`,
`dropped_set_ids`, and `dropped_set_reason_counts` in
`set_collection_summary`.

## Interpreting the Result

ORA tests whether the selected identifiers overlap a set more than expected
under the explicit background universe. Results depend strongly on the
background you supply. PhosPy does not infer that universe from a dataset,
reference bundle, or set collection.

The background matters because it defines the identifiers that could have been
selected. If a foreground identifier is absent from the background, PhosPy
raises a validation error by default instead of inventing a mapping, expanding
the universe, or silently shrinking the foreground. When the caller explicitly
sets `selected_outside_background_policy="drop"`, the dropped identifiers and
retained fraction are recorded in diagnostics and provenance.

The enrichment ratio is a descriptive overlap summary for the selected
identifiers, set members, and background used in the run. Adjusted p-values
describe statistical evidence under the ORA model and selected multiple-testing
correction; they are not proof that a pathway is active, regulated, causal, or
mechanistically responsible for the observed phosphoproteomics pattern.

Gene-level and site-level results are not interchangeable. A gene-symbol set is
not reinterpreted as a PTM set, and a PTM set is not collapsed to gene symbols.
Do not mix gene-level and site-level interpretation in the same result.

## Provenance and Reproducibility

Workflow provenance records method, identifier column and kind, collection kind,
analysis level, explicit background size, selected identifier count, selected
identifier source, set-collection source metadata when provided,
multiple-testing correction, table fingerprints, offline/no-online-resource
policy, and limitations. It also records `universe_policy`, including the
resolved selected and set-member outside-background policies, selected
outside-background identifiers, retained foreground fraction, and deterministic
set-member outside-background preview metadata.

When identifier-set provenance is supplied, `result.provenance.workflow_parameters`
also stores compact `selected_identifier_provenance` and
`background_identifier_provenance` payloads. These payloads include source type,
source label, normalized identifier count, upstream workflow/result IDs when
provided, and the nested input-intensity-scale evidence payload. They do not
copy full identifier lists; identifier lists remain represented by the existing
input-table fingerprints.

## Limitations

- Offline ORA only.
- No bundled curated enrichment resources for this feature.
- No online service calls.
- No ranked-list enrichment, GSEA, ssGSEA, PTM-SEA, or pathway activity
  inference.
- ORA does not prove pathway activation, regulation, or biological causality.
- Background choice is caller-owned and scientifically important.

## Minimal Example

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
        min_set_size=2,
    ),
)

result = EnrichmentWorkflow().run(request)
print(result.table.loc[:, ["term_id", "input_overlap_count", "p_value"]])
```
