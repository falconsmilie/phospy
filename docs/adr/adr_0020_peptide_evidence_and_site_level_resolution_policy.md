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

## Amendment: Differential Uncertainty Scope

**Date:** 2026-07-29

ADR-0020 covers dataset-construction resolution from peptide evidence to
analysis-ready site-level intensity rows. It remains the preferred
PhosPy-origin lane for peptide-to-site differential analysis: resolve peptide
evidence at the sample-intensity level first, then fit the core site-level
differential model.

Post-hoc peptide-level differential estimate combination is governed separately
by [ADR-0041: Peptide-to-Site Differential Uncertainty Policy](adr_0041_peptide_to_site_differential_uncertainty_policy.md).
The supported post-hoc route requires typed standard errors, original
finite-degree-of-freedom uncertainty, source experiment/run identifiers,
dependence policy, and mapping policy. Same-experiment peptide estimates are
rejected unless a future dependence-aware method is explicitly added.

Scientific rationale:

- same-experiment peptide-level differential statistics are not independent
  external studies and cannot be treated as independent production
  meta-analysis input;
- finite-degree-of-freedom t statistics must not be used directly as z
  statistics;
- a single-estimate pass-through is not meta-analysis;
- multi-site allocation remains explicit provenance, not hidden statistical
  repair.

### Explicit Aggregation Semantics

The semantics in this section are dataset-construction signal-allocation and
site-resolution semantics. They are not post-hoc differential statistic
meta-analysis semantics.

Peptide-to-site aggregation is scientifically explicit and owned by
`phospy.science.evidence.dataset_resolution`:

- **Mapping-weight source policy:**
  `explicit_mapping_weight_when_supplied_else_equal_fraction_per_resolved_site`
- **Mapping-weight normalisation policy:**
  `sum_to_one_per_peptide_evidence_row`
- **Signal allocation policy:** `multiply_peptide_signal_by_mapping_fraction`
- **Site summarisation policy:** `arithmetic_mean_of_allocated_signals`
- **Duplicate evidence policy:**
  `retain_duplicate_peptide_evidence_rows_as_separate_observations`
- **Mixed ambiguity policy:**
  `combine_ambiguous_and_unambiguous_allocated_signals_in_site_mean`
- **Localisation aggregation policy:**
  `arithmetic_mean_of_finite_reported_localisation_values`
- **Legacy aggregation-policy alias:**
  `legacy_alias_for_arithmetic_mean_of_allocated_signals`
- **Mapping-fraction representation:** `site_mapping.mapping_weight`
- **Weight normalisation contract:** mapping weights must sum to `1.0` per
  `peptide_row_id`
- **Derived default when absent:** equal weight per mapped site

The current supported mathematics is:

```text
a[p,s,j] = w[p,s] * x[p,j]
y[s,j] = mean over retained evidence rows p mapped to s of a[p,s,j]
```

where:

- `p` is a retained peptide evidence row.
- `s` is a resolved site.
- `j` is a sample.
- `x[p,j]` is the peptide-row signal.
- `w[p,s]` is the mapping allocation fraction for that peptide row and site.
- `a[p,s,j]` is the allocated peptide-row signal.
- `y[s,j]` is the site-level signal.

`w[p,s]` is an allocation fraction, not a statistical inverse-variance weight
or localisation-confidence weight. The final arithmetic mean is taken over
allocated evidence-row signals. This is not equivalent in general to either
`sum(w * x)` or `sum(w * x) / sum(w)`.

This policy does not generally conserve total signal after site-level
summarisation because the allocated signals are averaged per resolved site.
Duplicate retained rows affect the arithmetic mean because they are treated as
separate evidence observations. Ambiguous and unambiguous evidence rows mapped
to the same site are summarised by the same arithmetic mean over allocated
signals.

This ADR documents the currently supported policy. It does not claim that this
policy is universally optimal; alternatives such as summation, conventional
normalised weighted means, robust medians, best-peptide methods, or hierarchical
models require separate scientific review and a future ADR.

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
- mapping-weight source policy and actual source observed for the run
- mapping-weight normalisation policy
- signal allocation policy
- site summarisation policy
- legacy aggregation-policy alias + formula
- duplicate evidence policy + duplicate peptide row count
- mixed-ambiguity handling policy
- localisation aggregation policy

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
