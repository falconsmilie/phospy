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

### Explicit Aggregation Semantics

Peptide-to-site aggregation is scientifically explicit and owned by
`phospy.science.evidence.dataset_resolution`:

- **Aggregation policy name:** `mapping_weighted_mean`
- **Formula:** `site_intensity = mean(per_peptide_intensity * mapping_weight)`
- **Mapping-weight representation:** `site_mapping.mapping_weight`
- **Weight normalisation contract:** mapping weights must sum to `1.0` per
  `peptide_row_id`
- **Derived default when absent:** equal weight per mapped site
- **Duplicate peptide handling:** retain all peptide rows as independent
  observations (no sequence-level de-duplication)
- **Mixed ambiguous/unambiguous handling:** both contribute through the same
  weighted-mean aggregation once mapping is resolved

This replaces implicit split-plus-mean behaviour with named, provenance-visible
semantics.

### Protein Identity Metadata

Peptide-evidence protein_accession is row-identity metadata. It must be
preserved as protein_accession or explicit protein_namespace/protein_identifier
metadata. It must not be rewritten into protein_id, which remains available for
grouping semantics.

If peptide-evidence rows collapse to one resolved site, the resolved site must
not aggregate multiple distinct accessions silently. Conflicting accessions for
one resolved site are rejected until a future ADR defines explicit aggregation
semantics.

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

## Sequence Context Policy

Peptide-evidence resolution may normalise sequence strings by trimming
surrounding whitespace and uppercasing letters, but it must not repair
biological sequence context. In particular, the resolver must never rewrite the
centre residue of a provided `site_sequence` to match a resolved site token.

When multiple peptide-evidence rows collapse to one resolved site, supplied
`site_sequence` context is resolved as a set, not by row order. All non-null,
non-blank supplied values for that resolved site are normalized with the same
rules. They must normalize to exactly one unique valid value, and there must be
no invalid supplied non-blank values in the same group.

Conflicting valid contexts are rejected with `PhosPyInputError`. PhosPy must not
guess by choosing the first row, the most common value, the lexicographically
smallest value, or any other implicit precedence. Mixed evidence containing one
valid supplied value and one or more invalid supplied values is also rejected,
because reducing that group to the valid value would hide source-evidence
disagreement.

Split multi-site peptide-context derivation is a fallback for absent or
unusable split-site context when existing split-domain rules can derive the
resolved-site window from `peptide_sequence` plus `site_string`. It is not a
precedence rule over conflicting valid supplied contexts. When all supplied
contexts for a split target are invalid, derivation may succeed only if it
produces exactly one distinct normalized sequence; multiple derived candidates
or no derived candidate fail deterministically.

For unambiguous centred sequence windows, centre-residue mismatch is a hard
input error. Users must remove the peptide-evidence `site_sequence` value to
allow trusted reference derivation, or correct the upstream evidence. A future
repair path would require an explicit ADR and machine-readable provenance for
the repair decision.

## Provenance Requirements

The dataset report/provenance includes:

- peptide observations received
- unique analysis-ready `site_key` rows produced
- ambiguous observations
- excluded observations
- split observations (when applicable)
- selected multi-site policy
- aggregation policy + formula
- mapping-weight source + normalisation contract
- duplicate peptide policy + duplicate peptide row count
- mixed-ambiguity handling policy

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
- Explicit mapping-weight contracts reject malformed mappings where per-peptide
  weights do not sum to `1.0`.

## Responsibility Audit

This ADR explicitly keeps ownership boundaries as:

- evidence resolution: `phospy.science.evidence` + dataset preprocessing/builder
- dataset model: no peptide collapse
- kinase workflow: no peptide ambiguity resolution
- signalome workflow: no peptide ambiguity resolution
- dataset builder: invokes evidence resolution and records provenance, but does
  not duplicate aggregation policy logic
- downstream workflows: consume resolved site-level data after provenance is
  attached
