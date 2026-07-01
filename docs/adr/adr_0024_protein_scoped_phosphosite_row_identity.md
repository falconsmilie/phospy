# ADR: Protein-Scoped Phosphosite Row Identity

## Document Control

- **ADR ID:** ADR-0024
- **Title:** Protein-Scoped Phosphosite Row Identity
- **Status:** Accepted
- **Date:** 2026-05-27
- **Decision Type:** Architecture Decision Record

Supersedes ADR-0023 for analysis-ready phosphosite row identity scope and
amends ADR-0021 and ADR-0003 where they discuss analysis-ready row identity.

## Status

Accepted.

## Context

ADR-0023 documented the previous boundary where analysis-ready phosphosite
rows were indexed by display-oriented identifiers such as `GENE;SITE;`.

That boundary does not preserve protein-scoped scientific identity when the
same display token appears across distinct protein contexts. The repository now
needs one explicit analysis-ready row identity that is protein-scoped, unique,
and stable for scientific workflows, while preserving human-readable display
labels.

## Decision

PhosPy adopts a two-key model for analysis-ready phosphosite rows:

- `site_key` is the analysis-ready row identity.
- `display_id` is the human-readable display label, typically `GENE;SITE;`.

`site_key` must be unique in analysis-ready datasets.

`display_id` is not the analysis-ready row identity and may repeat.

`AnalysisReadyPhosphoDataset.phospho.index` is `site_key`.

`AnalysisReadyPhosphoDataset.site_metadata.index` is `site_key`.

`AnalysisReadyPhosphoDataset.site_metadata["display_id"]` must be present.

Auditable protein context metadata is required to construct and validate
`site_key`. At the advanced/trusted direct analysis-ready boundary,
`site_metadata` must include non-empty:

- `site_key`
- `display_id`
- `organism`
- `protein_namespace`
- `protein_identifier`
- `gene_symbol`
- `site`
- `site_sequence`

Direct analysis-ready datasets must not silently fall back to display-site row
identity. Ordinary user construction should use `AnalysisReadyDatasetBuilder`.

Builder input may remain user-friendly and accept legacy display-indexed input
only when enough protein context is available to deterministically derive
`site_key`.

Peptide-evidence protein_accession is row-identity metadata. It must be
preserved as protein_accession or explicit protein_namespace/protein_identifier
metadata. It must not be rewritten into protein_id, which remains available for
grouping semantics.

Workflows operate on `site_key`. Site-level workflow outputs that materialize
row identity include both `site_key` and `display_id`.

## Boundary Clarifications

Once `site_key` becomes row identity, repeated `display_id` values are allowed.
`display_id` duplicates are valid only under unique `site_key` rows.

This decision applies to analysis-ready dataset identity boundaries and
downstream workflow contracts that consume those datasets.

Display-indexed input is a builder compatibility input only. It is not valid
direct analysis-ready identity. Direct `AnalysisReadyPhosphoDataset`
construction is advanced/trusted use and requires `phospho.index`,
`site_metadata.index`, and `site_metadata["site_key"]` to all use the same
unique encoded `site_key` values.

Kinase reference resources may continue to use display IDs at the reference
boundary. The kinase workflow must match those display IDs through an explicit
reference-mapping layer from dataset `display_id` metadata to internal
`site_key` rows. Reference display IDs remain reference/display identifiers and
are not converted into analysis-ready row identity.

## Explicitly Out of Scope for This Migration

This migration does not include:

- automatic protein remapping
- peptide-evidence modelling changes
- automatic duplicate aggregation beyond existing explicit duplicate policies
- opaque hashes as the only public row key
- `pandas.MultiIndex` as the public row index

## Consequences

Analysis-ready identity becomes protein-scoped and unambiguous, while
human-readable display labels remain available for interpretation and reporting.

Direct dataset constructors now require enough identity context to produce
`site_key` and must fail explicitly when that context is missing. They also
validate that `site_key` matches the metadata-derived protein-scoped key.

Builder pathways can preserve backward-friendly ingestion shape only when they
can derive the required `site_key` without ambiguity.
