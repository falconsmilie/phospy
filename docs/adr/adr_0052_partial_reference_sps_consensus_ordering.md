# ADR-0052: Partial-Reference SPS Consensus Ordering

## Status

- **ADR ID:** ADR-0052
- **Title:** Partial-Reference SPS Consensus Ordering
- **Status:** Accepted
- **Date:** 2026-09-18
- **Decision Type:** Scientific Algorithm and Reproducibility
- **Refines:** ADR-0051

## Context

SPS discovery permits a site to contribute when it has valid condition-level
evidence in at least `minimum_datasets_per_site` references. This retains sites
that are absent or unrankable in some independent reference datasets. The
existing Fisher-style consensus score was calculated from the available ranks,
but scores calculated from different numbers of references were compared
directly. That left missing-reference patterns able to change selection without
an explicit cross-count scientific policy.

## Decision

Partial-reference contribution remains supported. Let (R_i) be the references
that provide valid evidence for site (i), (k_i=|R_i|), and (n_r) the number
of candidate sites ranked in reference (r). If \(\bar{a}_{ir}\) is the
tie-aware midrank of site (i) in reference (r), with smaller stability
magnitudes ranked first, define

\[
q_{ir} = \frac{n_r - \bar{a}_{ir} + 1/2}{n_r}, \qquad
T_i = -2\sum_{r\in R_i}\log(q_{ir}),
\]

and retain the fixture-backed PhosR-compatible within-stratum consensus
statistic

\[
C_i = \Pr\{\chi^2_{2(k_i-1)} \ge T_i\}.
\]

Sites are ordered lexicographically by

\[
(-k_i,\;-C_i,\;\text{site_key}_i).
\]

Thus, more independent contributing references take precedence. The consensus
score compares stability only among sites with the same evidence count, and
ascending `site_key` resolves an exact tie. The rationale is
conservative corroboration: absence of reference evidence is not evidence of
stability and cannot improve a site's position relative to a site supported by
more references. `minimum_datasets_per_site` remains the explicit eligibility
floor.

The policy is named
`contributing_dataset_count_descending_then_consensus_stability_descending` in
provenance. Algorithm version `2.0.0` owns this ordering. Version `1.0.0`
serialized results retain their historical score-first ordering when loaded;
they are not silently reinterpreted.

## Consequences

- A very stable site supported by fewer references can rank below a less stable
  site supported by more references.
- Reordering rows, samples, or reference inputs cannot change the result.
- When every site is supported by every reference, the evidence-count key is
  constant. The numeric consensus scores and complete-reference ordering are
  therefore unchanged, preserving the pinned PhosR comparison.
- Partial-reference behavior remains a PhosPy policy outside the external PhosR
  parity claim and is protected by synthetic validation.

## Related Records

- [ADR-0051: SPS Reference Quantitative State and Discovery Separation](adr_0051_sps_reference_quantitative_state_and_discovery_separation.md)
- [Native SPS/RUV-Style Discovery and Correction](../api/sps-ruv.md)
- [Comparison With PhosR](../parity.md)
