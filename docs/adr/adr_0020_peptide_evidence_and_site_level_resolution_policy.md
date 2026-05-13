# ADR: Peptide Evidence and Site-Level Resolution Policy

## Document Control

- **ADR ID:** ADR-0020
- **Title:** Peptide Evidence and Site-Level Resolution Policy
- **Status:** Accepted
- **Date:** 2026-05-12
- **Decision Type:** Architecture Decision Record

## Context

PhosPy workflows downstream of dataset building consume site-level rows.
However, some upstream pipelines provide peptide-level evidence where one
observation can map to one or more phosphosites. Treating that ambiguity as
already resolved silently changes scientific meaning and prevents reproducible
audit of row-level decisions.

## Decision

Dataset building now supports two explicit input modes:

1. `site_level_resolved`
2. `peptide_evidence`

When `peptide_evidence` is used, `multi_site_policy` is required and must be
one of:

- `reject`
- `exclude_from_sequence_scoring`
- `keep_joint`
- `split`

Policy mapping reuses the existing `phospy.science.evidence` models and multi-site
resolution logic (`PeptideEvidenceTable`, `MultiSiteHandlingConfig`,
`SiteEvidenceMapping`) instead of introducing parallel ambiguity models.

## Peptide vs Site Ambiguity

- **Peptide-level ambiguity:** one peptide evidence row references multiple site
  tokens (for example `S10,T12`).
- **Site-level data contract:** downstream lanes require explicit site rows with
  explicit provenance of how ambiguity was resolved.

Ambiguity is not hidden by implicit defaults inside normalisation, kinase
scoring, or signalome scoring.

## Where Resolution Happens

Ambiguity resolution is owned by evidence/preprocessing boundary logic during
dataset interpretation, before `AnalysisReadyPhosphoDataset` construction.

`dataset.preprocessing_report` and dataset run provenance record the decision
and counts.

## Provenance Requirements

The dataset report/provenance includes:

- peptide observations received
- unique site IDs produced
- ambiguous observations
- excluded observations
- split observations (when applicable)
- selected multi-site policy

## Consequences

### Positive

- Input boundary is explicit about whether site-level resolution is already done.
- Multi-site handling is reproducible and audit-ready.
- Downstream workflows can trust dataset rows without owning peptide collapse.

### Tradeoffs

- Call sites that provide peptide evidence must now provide
  `peptide_evidence_sample_intensity_columns` and `multi_site_policy`.
- `keep_joint` rows preserve ambiguous site tokens and must be interpreted as
  ambiguous by consumers.

## Responsibility Audit

This ADR explicitly keeps ownership boundaries as:

- evidence resolution: `phospy.science.evidence` + dataset preprocessing/builder
- dataset model: no peptide collapse
- kinase workflow: no peptide ambiguity resolution
- signalome workflow: no peptide ambiguity resolution
