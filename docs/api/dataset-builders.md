---
search:
  exclude: true
---

# Dataset Builders

`AnalysisReadyDatasetBuilder` is the supported way to construct an
analysis-ready dataset from ordinary input tables.

Most users should follow [Build an Analysis-Ready Dataset](dataset-build-workflow.md).
This page summarises the stricter identity and diagnostic contracts for readers
who need them.

## Use the Public Builder

The builder accepts input tables, applies the requested preprocessing, validates
the resulting scientific state, and returns an
`AnalysisReadyPhosphoDataset`.

```python
from phospy import AnalysisReadyDatasetBuilder

builder = AnalysisReadyDatasetBuilder()
dataset = builder.run(dataset_request)
```

The direct `AnalysisReadyPhosphoDataset(...)` constructor raises immediately.

## Understand Row Identity

Analysis-ready rows use `site_key` as their unique identity.

- `dataset.phospho.index` and `dataset.site_metadata.index` contain `site_key`
  values.
- `dataset.site_metadata["site_key"]` exactly matches the metadata index.
- `display_id` is a human-readable label, such as `GENE;SITE;`.
- A `display_id` may repeat when the corresponding `site_key` values differ.
- Site-level workflow outputs include both identifiers where applicable.

Builder input may use a display-indexed shape when enough protein context is
available to derive a unique `site_key`. The final dataset does not silently
fall back to display labels as row identity.

Kinase references may remain display-ID keyed at the reference boundary. The
Kinase workflow maps those reference labels to dataset `site_key` values before
scoring.

## Resolve Duplicate Sites Deliberately

Rows that resolve to the same `site_key` are scientifically ambiguous and fail
by default. PhosPy does not make them unique by renaming them.

The supported non-error policies are deliberate scientific choices:

- `max_mean_signal`
- `first`
- `aggregate_mean`
- `aggregate_median`

After using one of these policies, inspect:

- `dataset.preprocessing_report.duplicate_site_resolution`
- `dataset.preprocessing_report.metadata_conflicts`

Peptide-evidence input follows the policies documented in the
[dataset-build workflow](dataset-build-workflow.md#use-peptide-evidence-input).

## Use Trusted Construction Only When Appropriate

`AnalysisReadyPhosphoDataset.from_trusted_tables(...)` is an advanced route for
callers who already own prepared, `site_key`-indexed tables and their evidence.
It is not a shortcut around validation.

Trusted construction requires:

- aligned phosphosite and metadata tables;
- non-empty `site_sequence` values;
- `organism`, `protein_namespace`, `protein_identifier`, `gene_symbol`, and
  `site` metadata;
- encoded `site_key` indexes; and
- `TrustedDatasetConstructionAssertions` covering identity, intensity scale,
  quantitative meaning, localisation, sequence, reference context, and table
  structure.

Localisation evidence must record its source, policy, and threshold. Sequence
presence is not localisation evidence. Supplied provenance must match the
actual table fingerprints.

These checks establish that the caller supplied a structurally coherent claim.
They cannot prove that user-asserted biological provenance is correct.

## Inspect Builder Diagnostics

Use the documented output objects:

- the returned `AnalysisReadyPhosphoDataset`;
- `dataset.preprocessing_report`;
- `dataset.provenance`; and
- the result fields documented by each workflow.

Do not import internal processing-state classes from implementation modules.
Those classes are not public extension points.

For the formal identity decision, see
[ADR-0024: Protein-Scoped Phosphosite Row Identity](../adr/adr_0024_protein_scoped_phosphosite_row_identity.md).
