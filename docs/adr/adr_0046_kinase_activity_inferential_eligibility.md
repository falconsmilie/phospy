# ADR-0046: Kinase Activity Inferential Eligibility

## Status

- **ADR ID:** ADR-0046
- **Title:** Kinase Activity Inferential Eligibility
- **Status:** Accepted
- **Date:** 2026-08-05
- **Decision Type:** Scientific Workflow Contract

## Context

KSEA-style normal-approximation p-values assume that the tested substrate set is
fixed independently of the quantitative profile used in the z-score. In the
default kinase workflow, substrate membership can be selected by ranking or
thresholding prediction scores derived from the same phosphosite matrix later
tested by KSEA. Combined profile/motif scoring and leave-one-out profile scoring
also consume the tested matrix before membership selection.

Emitting ordinary p-values after this adaptive selection would present a null
model that was not actually used. A warning alone is insufficient because finite
p/q values are easily reused downstream as if they were ordinary inferential
evidence.

## Decision

KSEA-style activity requires typed substrate-membership provenance before
ordinary p/q values can be emitted. The activity science domain derives and
enforces the final inferential-eligibility decision from those facts. Workflow
code supplies facts only; it does not supply the final scientific conclusion.
Executors must not infer eligibility from raw matrices, method names, or
descriptive labels.
Source category and supporting facts must form one coherent state before the
eligibility decision is derived. Contradictory combinations are rejected rather
than downgraded or caveated.

The membership-selection record includes:

- source category;
- selection method and version;
- score source;
- threshold/top-k policy;
- source reference fingerprint(s);
- selection quantitative matrix fingerprint when quantitative data were
  consumed while deriving or selecting membership;
- tested quantitative matrix fingerprint for the exact KSEA background matrix;
- whether selection consumed the tested matrix;
- selected kinase universe;
- selected substrate universe; and
- the derived inferential decision, including availability, reason/status,
  missing evidence, and policy version.

The activity-domain KSEA membership policy is fail-closed:

| Membership source category | Ordinary KSEA p/q policy |
| --- | --- |
| `profile_derived` | Ineligible. |
| `fused_profile_motif` | Ineligible. |
| `prediction_selected` with data-adaptive or quantitative selection evidence | Ineligible. |
| `unknown` or incomplete provenance | Ineligible. |
| leave-one-out profile-derived membership | Ineligible. |
| `fixed_external_reference` | Eligible only with source-reference fingerprints, tested-matrix fingerprint, selected universes, and explicit fixed-reference independence policy evidence. |
| `sequence_only_motif` | Eligible only with source-reference fingerprints, tested-matrix fingerprint, selected universes, `kinase_library_motif_scores` score source, and explicit sequence-only independence policy evidence. |

Missing source fingerprints, missing tested-matrix fingerprint evidence, or
incomplete independence-policy evidence produce descriptive KSEA output only.
They must not silently become eligible.

The source-specific coherence gate rejects records such as
`fixed_external_reference` combined with data-adaptive membership,
`consumed_tested_matrix=True`, a selection quantitative-matrix fingerprint,
known profile-derived or fused score sources, sequence-only motif score sources,
or source-specific method/policy tokens from another category. Built-in score
sources have one interpretation: `profile_scores` and
`rank_weighted_fusion_scores` are profile-derived,
`combined_profile_motif_scores` is fused profile/motif-derived, and
`kinase_library_motif_scores` is sequence-only motif-derived. Arbitrary
externally supplied method strings do not establish independence; fixed
external eligibility requires the explicit fixed-reference independence-policy
token plus supported version and complete evidence.

For ineligible membership, KSEA computes descriptive z-scores and records
missing p/q values. `p_value_matrix` and `q_value_matrix` are unavailable, and
statistics tables use missing p/q cells plus structured inferential-status
columns. Benjamini-Hochberg adjustment is run only over eligible finite
p-values.

Every constructed KSEA result carries an explicit membership-selection record.
Missing or legacy membership provenance is represented by an explicit missing,
ineligible record for descriptive-only output. Finite p/q matrices or finite
statistics-table p/q cells without eligible membership provenance are rejected
during result construction and bundle reconstruction.

A future nested resampling or sample-splitting KSEA method may emit valid
inferential values for adaptive membership only if substrate selection is
repeated inside the null procedure or otherwise separated from the tested data.
That would be a new method/policy contract, not an interpretation change to the
ordinary normal approximation.

Serialized membership payloads are reconstructed from underlying facts, then
the decision is recomputed. Serialized eligibility, reason, status, and nested
decision fields are accepted only when they match the recomputed decision.
Contradictory or tampered payloads are rejected, including relabelling an
adaptive record as `fixed_external_reference`, retaining adaptive facts under a
fixed-external label, or adding favourable independence tokens to an adaptive
record. Legacy payloads that lack the new tested-matrix or independence
evidence are not upgraded to ordinary inference eligibility, and missing
fingerprints are not fabricated during migration.

The ordinary normal-approximation assumptions remain scientific assumptions
even after membership eligibility is established; eligibility only establishes
that the ordinary KSEA p/q output is allowed under this policy gate.

## Consequences

Positive consequences:

- Default `phosr_rank_weighted` profile-derived KSEA no longer emits ordinary
  p/q values under an invalid null.
- Combined profile/motif and leave-one-out profile-derived membership are
  explicitly marked descriptive for ordinary KSEA inference.
- Motif-only or fixed-reference membership can still use the documented
  approximation when the membership provenance is independent.
- Provenance and bundles preserve the complete selection chain and eligibility
  state for audit and reload.

Negative consequences:

- Existing downstream consumers that expected finite KSEA p/q values from the
  default adaptive workflow must switch to z-score interpretation or provide an
  independent membership source.

Neutral consequences:

- KSEA z-score calculation and computability statuses are unchanged.
- ssGSEA permutation p-values and differential/enrichment p-value policies are
  unchanged.

## Implementation Notes

The typed membership contract lives in
`src/phospy/science/activities/membership.py`. The KSEA method gates p/q
emission in `src/phospy/science/activities/methods/ksea_zscore.py`.

Kinase workflow code constructs membership provenance in
`src/phospy/workflows/kinase/membership.py` from the resolved score source,
prediction threshold/top-k policy, reference fingerprints, the quantitative
matrix actually consumed by selection when applicable, the exact
`request.ksea_background_phospho_matrix` tested by KSEA, and selected
kinase/substrate universes.

Kinase result caveats include a structured adaptive-membership caveat when KSEA
p/q values are unavailable. Kinase bundle manifests persist
`outputs.activity.membership_selection` so reload preserves eligibility state.

The KSEA scientific policy record is versioned as
`ksea_zscore_activity_v1` policy version `4`. The membership-selection policy
version is `3`, and the KSEA membership inferential policy version is `3`.

## Related Records

- [ADR-0033: Result Caveats and Scientific Warning Ownership](adr_0033_result_caveats.md)
- [ADR-0034: Quantitative State, Motif Scoring Semantics, and Reference Context](adr_0034_quantitative_state_motif_semantics_and_reference_context.md)
- [ADR-0035: Provenance Immutability and Stable Serialization](adr_0035_provenance_immutability_and_stable_serialization.md)
- [Kinase Workflow](../api/kinase.md)
- [Workflow Contracts](../workflow_contracts.md)
