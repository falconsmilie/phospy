# ADR: Localisation Confidence and Site-Level Eligibility Policy

## Document Control

- **ADR ID:** ADR-0018
- **Title:** Localisation Confidence and Site-Level Eligibility Policy
- **Status:** Accepted
- **Date:** 2026-05-11
- **Decision Type:** Architecture Decision Record

## Context

Site-level phosphoproteomics interpretation depends on confidence that the
reported modified residue is correctly localised on the peptide. Ambiguous
localisation can materially change biological interpretation at both kinase and
signalome levels. If low-confidence or unknown localisation is allowed to enter
site-level pipelines silently, downstream analyses can look precise while
including uncertain site assignments.

PhosPy needs an explicit, auditable policy for localisation confidence at
dataset construction time so workflow components consume already eligible data.

Update note (2026-07-16, trusted construction evidence): Trusted
analysis-ready construction must not infer localisation confidence from
`site_sequence` presence or from any generic "trusted" switch. The preferred
advanced factory requires explicit localisation evidence with source, policy,
and threshold, or an explicit localisation waiver. Compatibility direct
construction remains advanced/internal and records missing assertion metadata
when callers bypass that preferred factory.

Update note (2026-07-17, no inferred localisation): The seven-dimension
trusted construction assertion schema keeps localisation separate from sequence
and aligned-structure evidence. `site_sequence` presence remains mandatory at
the dataset boundary but is never accepted as localisation evidence. The
localisation dimension must record source, policy, and threshold, or an
explicit localisation waiver.

## Decision

PhosPy enforces localisation eligibility at the dataset preprocessing/validation
boundary using an explicit dataset policy:

- `require_threshold`
- `allow_missing_with_waiver`
- `ignore`

Default analysis-ready behavior is conservative:

- `mode="require_threshold"`
- `min_confidence=0.75`
- `confidence_column="localisation_confidence"`

### Enforcement Boundary

Localisation confidence rules are enforced in:

- preprocessing configuration validation
- preprocessing-stage execution
- dataset/site-metadata validation and evidence-model range checks

Localisation filtering or correction is **not** introduced in kinase or
signalome executors.

### Waiver Rules

`allow_missing_with_waiver` is acceptable only when an explicit waiver reason is
provided. Waiver mode can retain missing or below-threshold localisation
entries, but invalid confidence values remain rejected.

### Provenance and Reporting

When localisation validation is waived, PhosPy records this decision in:

- preprocessing trace diagnostics (mode, threshold, waiver reason, counts)
- row-audit/report rows for retained-with-waiver sites
- preprocessing operations/provenance payload (`localisation_mode`,
  `localisation_min_confidence`, `localisation_confidence_column`,
  `localisation_waiver_reason`)

For trusted `AnalysisReadyPhosphoDataset.from_trusted_tables(...)`
construction, localisation is recorded in
`TrustedDatasetConstructionAssertions.localisation` as either:

- typed evidence: source, policy, and threshold
- explicit waiver: waiver reason, with optional policy/details

This assertion payload is immutable and linked into construction provenance by
fingerprint. Sequence presence alone is not localisation evidence.

## Consequences

### Positive

- Prevents silent inclusion of unknown/low-confidence localisation in
  site-level dataset construction.
- Makes relaxation of localisation rules explicit and reviewable.
- Preserves clear separation of responsibilities: dataset preprocessing validates
  eligibility; workflows consume prepared datasets.

### Tradeoffs

- Some legacy datasets now require additional metadata or explicit waiver
  configuration.
- Users must choose policy intentionally for exploratory vs. analysis-ready
  builds.

## Scope Boundaries

This ADR does not:

- move localisation eligibility into kinase workflow execution
- move localisation eligibility into signalome clustering execution
- apply hidden default filtering without provenance

## References

Cox, J., & Mann, M. (2008). MaxQuant enables high peptide identification rates,
individualized p.p.b.-range mass accuracies and proteome-wide protein
quantification. *Nature Biotechnology, 26*(12), 1367-1372.
https://doi.org/10.1038/nbt.1511

Olsen, J. V., Blagoev, B., Gnad, F., Macek, B., Kumar, C., Mortensen, P., &
Mann, M. (2006). Global, in vivo, and site-specific phosphorylation dynamics in
signaling networks. *Cell, 127*(3), 635-648.
https://doi.org/10.1016/j.cell.2006.09.026

Sharma, K., D'Souza, R. C. J., Tyanova, S., Schaab, C., Wisniewski, J. R., Cox,
J., & Mann, M. (2014). Ultradeep human phosphoproteome reveals a distinct
regulatory nature of Tyr and Ser/Thr-based signaling. *Cell Reports, 8*(5),
1583-1594. https://doi.org/10.1016/j.celrep.2014.07.036
