# ADR: Phosphosite Identity and Protein Context Policy

## Document Control

- **ADR ID:** ADR-0021
- **Title:** Phosphosite Identity and Protein Context Policy
- **Status:** Accepted
- **Date:** 2026-05-12
- **Decision Type:** Architecture Decision Record

## Context

PhosPy historically centered site identity on display-oriented `GENE;SITE;`
tokens. That format is convenient and remains useful for matrix indexing, but
it is not sufficient scientific identity for accession-aware phosphoproteomics.

The same display token can correspond to different proteins, isoforms, source
records, or namespaces. If those contexts are silently collapsed, downstream
kinase/signalome/sequence-aware interpretation can become irreproducible or
scientifically ambiguous.

## Decision

PhosPy introduces an explicit phosphosite identity model that separates:

- **Display identity:** standardized `GENE;SITE;` token
- **Scientific identity context:** protein/accession/provenance fields

Identity is represented in a dedicated site-domain module and validated at
dataset/preprocessing boundaries before workflow execution.

### Identity Model

The identity model includes:

- `display_id`
- `gene_symbol`
- `residue`
- `position`
- optional context fields:
  - `organism`
  - `protein_id`
  - `protein_accession`
  - `isoform_id`
  - `source_namespace`
  - `source_site_id`

### Default Site-Token Strictness

By default, phosphosite `site` values must be strict phosphoproteomics tokens:
`S`, `T`, or `Y` followed by a positive integer (for example `S123`, `T45`,
`Y999`).

Malformed or non-phosphorylatable tokens (for example `FOO`, `A123`, `S0`,
blank, or null values) are rejected before workflow execution.

Opaque site values are only valid behind an explicit waiver
(`allow_opaque_site_values=True`), and that waiver must be provenance-visible
in any workflow path that enables it.

### Compatibility

Standardized `GENE;SITE;` parsing remains supported and is still required for
display IDs at strict dataset boundaries. Display identity remains an indexing
surface, not complete scientific identity.

### Ambiguity and Collision Policy

When multiple rows map to the same display site ID, PhosPy now rejects rows
that carry conflicting scientific identity contexts (for example conflicting
`protein_id`/`protein_accession` values) before duplicate-site collapse.

Semantically identical repeated identities remain valid.

## Workflow Requirements

- **Signalome:** remains strict about explicit protein identity metadata and
  consumes validated phosphosite identities.
- **Kinase workflows:** consume validated phosphosite identities but do not
  require protein/accession identity for all lanes yet.
- **Differential/simple matrix workflows:** continue to accept gene-site-only
  datasets when scientifically acceptable.

## Responsibility Audit

Ownership boundaries are:

- **Site identity domain module:** parsing + identity construction + collision
  validation.
- **Dataset/table schemas:** structural shape and column contracts.
- **Workflow validators:** workflow-specific identity requirements (for example
  strict signalome protein identity expectations).
- **Executors/interpreters:** consume validated identity; they do not infer or
  synthesize identity from fallback heuristics.

## Consequences

### Positive

- Scientific ambiguity in same-display-ID rows is explicit and testable.
- Identity expectations become workflow-aware and auditable.
- Compatibility with legacy display IDs is preserved where acceptable.

### Tradeoffs

- Some duplicate-site datasets that previously collapsed silently now fail fast
  until identity ambiguity is resolved upstream.
- Users may need to provide richer identity metadata for strict workflows.
