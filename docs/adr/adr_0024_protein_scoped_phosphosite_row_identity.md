# ADR: Protein-Scoped Phosphosite Row Identity

## Document Control

- **ADR ID:** ADR-0024
- **Title:** Protein-Scoped Phosphosite Row Identity
- **Status:** Accepted
- **Date:** 2026-05-27
- **Decision Type:** Architecture Decision Record

Supersedes ADR-0023 for analysis-ready phosphosite row identity scope and
amends ADR-0021 and ADR-0003 where they discuss analysis-ready row identity.

Update note (2026-07-26, sealed runtime construction boundary): direct
`AnalysisReadyPhosphoDataset(...)` construction now raises immediately. The
successful construction paths are
`AnalysisReadyDatasetBuilder.run(DatasetBuildRequest(...))` for ordinary
construction and `AnalysisReadyPhosphoDataset.from_trusted_tables(...)` for
advanced trusted reconstruction with complete
`TrustedDatasetConstructionAssertions`. Supplied trusted provenance must match
the actual represented-table fingerprints. Trusted assertions are audit
evidence supplied by the caller, not proof of biological correctness.

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

Analysis-ready datasets are single-organism datasets. The dataset has one
resolved `Organism`, and every row-level `organism` value plus every decoded
`site_key.organism` value must agree with it. Uniform row metadata may infer an
omitted dataset-level `organism`; an explicit dataset-level `organism` must
match every row. Mixed human/mouse/rat rows are invalid for
`AnalysisReadyPhosphoDataset`; if mixed-species support is added later, it must
use a separate explicit type and contract.

`AnalysisReadyPhosphoDataset.phospho.index` is `site_key`.

`AnalysisReadyPhosphoDataset.site_metadata.index` is `site_key`.

`AnalysisReadyPhosphoDataset.site_metadata["display_id"]` must be present.

Auditable protein context metadata is required to construct and validate
`site_key`. At the advanced/trusted factory reconstruction boundary,
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
`site_key`. Builder and trusted-factory construction paths normalize supported
organism aliases and case variants to the shared `Organism` enum before storing
analysis-ready state; arbitrary organism strings are not valid internal
organism state.

Update note (2026-07-17, provenance coherence): the same construction-boundary
organism rule applies to run provenance. If a direct, trusted, derived, builder,
or restored dataset path supplies `RunProvenance.reference_context` or selected
reference provenance, those organism values must resolve to the dataset
`Organism`. Workflows consume this sealed state and must not repair mismatches.

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
analysis-ready identity. Direct `AnalysisReadyPhosphoDataset(...)`
construction raises immediately. Advanced trusted reconstruction uses
`AnalysisReadyPhosphoDataset.from_trusted_tables(...)` and requires
`phospho.index`, `site_metadata.index`, and `site_metadata["site_key"]` to all
use the same unique encoded `site_key` values.

Dataset construction owns organism coherence. Workflows must not repair,
reinterpret, or independently normalize `dataset.organism`,
`site_metadata["organism"]`, decoded `site_key.organism`, or provenance
reference-context organism; they consume the resolved single-organism dataset
or fail earlier at the dataset boundary.

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

Trusted dataset reconstruction now requires enough identity context to produce
`site_key` and must fail explicitly when that context is missing. It also
validates that `site_key` matches the metadata-derived protein-scoped key and
that dataset, row, and decoded site-key organisms resolve to one `Organism`.

Builder pathways can preserve backward-friendly ingestion shape only when they
can derive the required `site_key` without ambiguity.
