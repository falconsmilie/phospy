# ADR: Supported Phosphosite Display-Site Identity Scope

## Document Control

- **ADR ID:** ADR-0023
- **Title:** Supported Phosphosite Display-Site Identity Scope
- **Status:** Superseded
- **Date:** 2026-05-26
- **Decision Type:** Architecture Decision Record

Superseded by ADR-0024 for analysis-ready phosphosite row identity.

## Status

Superseded.

## Context

PhosPy uses normalised phosphosite display identifiers such as `GENE;SITE;` as
the analysis-ready row index surface for `AnalysisReadyPhosphoDataset`.

Phosphosite rows may also carry protein and mapping context metadata, including
`protein_id`, `protein_accession`, `protein_group`, `isoform_id`, `source`,
`source_namespace`, and `site_id`.

Those fields are useful for downstream interpretation and reporting, but they
do not currently define analysis-ready dataset row identity.

Without an explicit decision record, contributors may incorrectly treat
protein-, isoform-, source-, or peptide-evidence-level context as the current
row identity boundary, and may place duplicate-site handling in downstream
workflow validators instead of the dataset-construction boundary.

## Decision

PhosPy currently uses normalised phosphosite display identifiers, for example
`GENE;SITE;`, as the analysis-ready dataset row identity.

`AnalysisReadyPhosphoDataset` supports one row per normalised display-site
identifier.

Duplicate display-site identifiers are rejected during dataset construction.
This rejection applies even when duplicate rows have identical protein, isoform,
or source metadata values.

Protein, isoform, source, and mapping-context fields may be stored as metadata,
but they do not currently define row identity.

## Supported behaviour

PhosPy currently supports all of the following:

- one analysis-ready row per normalised display-site identifier
- display-site identifiers such as `GENE;SITE;`
- protein/isoform/source fields stored as metadata only
- duplicate display-site rejection during dataset construction
- duplicate display-site rejection even when protein metadata is identical
- downstream workflows receiving only datasets that already passed this
  identity rule

## Unsupported behaviour

PhosPy does not currently support:

- protein-scoped row identity
- isoform-scoped row identity
- source-scoped row identity
- peptide-evidence-scoped row identity
- automatic duplicate phosphosite aggregation
- automatic duplicate phosphosite renaming
- signalome-specific duplicate repair
- downstream workflow validators as the primary place for duplicate phosphosite
  detection

## Consequences

Users must resolve duplicate display-site rows before dataset construction.

Ambiguous protein or isoform mappings must be handled upstream.

PhosPy intentionally fails early at dataset construction rather than silently
collapsing scientifically distinct rows.

Signalome workflows may require protein metadata, but that metadata does not
currently make dataset row identity protein-scoped.

Future protein- or isoform-scoped row identity support would require a
deliberate dataset-model change, not a small downstream validator change.

## Migration path

A future migration to protein- or isoform-scoped row identity would require
explicit architectural work and governance. At minimum, that migration must
include:

- a first-class scientific site identity key
- a preserved display label such as `GENE;SITE;`
- dataset-builder changes
- dataset table/model changes
- validation changes
- workflow identity-contract changes
- documentation updates
- migration guidance for existing users

This migration is out of scope for the current decision and is not implemented
here.
