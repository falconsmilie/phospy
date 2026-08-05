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
ordinary p/q values can be emitted. The activity science domain owns the final
inferential-eligibility decision. Workflow code supplies typed provenance, and
executors must not infer eligibility from raw matrices, method names, or
descriptive labels.

The membership-selection record includes:

- source category;
- selection method and version;
- score source;
- threshold/top-k policy;
- source reference fingerprint(s);
- quantitative dataset fingerprint when the tested matrix was consumed;
- whether selection consumed the tested matrix;
- selected kinase universe;
- selected substrate universe; and
- inferential eligibility plus reason/status.

Membership selected from a fixed external reference or a sequence-only motif
source can be eligible when it is demonstrably independent of the tested
quantitative matrix. Membership selected from profile-derived prediction,
combined profile/motif scoring, prediction-selected adaptive support, or
leave-one-out profile scoring is inferentially ineligible for ordinary KSEA
p/q values because selection consumed the tested matrix.

For ineligible membership, KSEA computes descriptive z-scores and records
missing p/q values. `p_value_matrix` and `q_value_matrix` are unavailable, and
statistics tables use missing p/q cells plus structured inferential-status
columns. Benjamini-Hochberg adjustment is run only over eligible finite
p-values.

A future nested resampling or sample-splitting KSEA method may emit valid
inferential values for adaptive membership only if substrate selection is
repeated inside the null procedure or otherwise separated from the tested data.
That would be a new method/policy contract, not an interpretation change to the
ordinary normal approximation.

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
prediction threshold/top-k policy, reference fingerprints, activity matrix
fingerprint when consumed, and selected kinase/substrate universes.

Kinase result caveats include a structured adaptive-membership caveat when KSEA
p/q values are unavailable. Kinase bundle manifests persist
`outputs.activity.membership_selection` so reload preserves eligibility state.

The KSEA scientific policy record is versioned as
`ksea_zscore_activity_v1` policy version `2` to record the p/q eligibility gate.

## Related Records

- [ADR-0033: Result Caveats and Scientific Warning Ownership](adr_0033_result_caveats.md)
- [ADR-0034: Quantitative State, Motif Scoring Semantics, and Reference Context](adr_0034_quantitative_state_motif_semantics_and_reference_context.md)
- [ADR-0035: Provenance Immutability and Stable Serialization](adr_0035_provenance_immutability_and_stable_serialization.md)
- [Kinase Workflow](../api/kinase.md)
- [Workflow Contracts](../workflow_contracts.md)
