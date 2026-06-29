# ADR: Supported Phosphosite Display-Site Identity Scope

## Document Control

- **ADR ID:** ADR-0023
- **Title:** Supported Phosphosite Display-Site Identity Scope
- **Status:** Superseded
- **Date:** 2026-05-26
- **Decision Type:** Architecture Decision Record

Superseded by ADR-0024 for analysis-ready phosphosite row identity. This ADR is
retained only as historical context for the pre-`site_key` boundary.

## Status

Superseded.

## Context

This ADR recorded an earlier display-label-based boundary. ADR-0024 supersedes
that boundary because display labels do not preserve protein-scoped scientific
identity when the same `GENE;SITE;` label appears in multiple protein or isoform
contexts.

## Current Decision After Supersession

The supported analysis-ready identity model is now:

- `site_key` is the unique analysis-ready row identity.
- `display_id` is the human-readable label, typically `GENE;SITE;`.
- `AnalysisReadyPhosphoDataset.phospho.index` is `site_key`.
- `AnalysisReadyPhosphoDataset.site_metadata.index` is `site_key`.
- `display_id` may repeat when `site_key` values differ.
- Direct `AnalysisReadyPhosphoDataset` construction is advanced/trusted use,
  requires `site_key`, and must not silently fall back to display labels.
- Direct construction also requires auditable protein context metadata:
  `organism`, `protein_namespace`, `protein_identifier`, `gene_symbol`, `site`,
  and `site_sequence`.
- Builder ingestion may accept legacy display-indexed input only when enough
  protein context exists to derive `site_key`.
- Workflows operate on `site_key`; site-level outputs that materialize identity
  include both `site_key` and `display_id`.
- Kinase references may use display IDs only through an explicit mapping layer
  to dataset `site_key` rows.

## Consequences After Supersession

This ADR must not be used as implementation guidance for current
analysis-ready row identity. Dataset/model/validation boundaries enforce
`site_key`; workflow validators compose those guarantees and do not repair
display-indexed datasets.

Duplicate `display_id` values are not invalid by themselves. They are valid
when distinct `site_key` values preserve the protein-scoped identity. Duplicate
rows still fail when they resolve to the same `site_key` without a configured
preprocessing policy that resolves them before final dataset construction.
