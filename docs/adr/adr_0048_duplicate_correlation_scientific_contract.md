# ADR-0048: Duplicate-Correlation Scientific Contract

## Status

- **ADR ID:** ADR-0048
- **Title:** Duplicate-Correlation Scientific Contract
- **Status:** Accepted
- **Date:** 2026-08-18

## Context

PhosPy already supports `reject` and `fixed_block` at the public paired-design
request boundary. `fixed_block` is a valid fixed-effects model: block identity
is represented by ordinary nuisance coefficients in the design matrix.

`duplicate_correlation` is also a valid model, but with different assumptions.
It treats block identity as correlation-group metadata rather than fixed design
coefficients. This ADR originally defined the internal contract before public
exposure; the completed implementation now exposes the policy through the
existing differential request/configuration contract.

## Decision

For `duplicate_correlation`, the fixed-effects design contains condition terms
and any already-supported non-block covariates. Block dummy variables must not
be included in the same model as a block-correlation structure.

Observations from different blocks are independent. Observations in the same
block use compound-symmetry correlation:

- diagonal entries equal `1`;
- off-diagonal entries equal the consensus correlation;
- singleton blocks remain valid observations but contribute no correlated pair.

One consensus correlation is estimated for the full design and reused for all
contrasts derived from that fit.

The first estimator implementation must use feature-wise REML correlation
estimates, Fisher `atanh` transformation, a `15%` trimmed mean from each tail,
and inverse `tanh` transformation. The trim value is fixed at `0.15`; it is not
a public tuning parameter.

For the estimator to run, the full non-block fixed-effects design must leave at
least two residual degrees of freedom, i.e. the analysed sample count must be at
least two larger than the design rank. Feature-level REML estimates with fewer
than two observed residual degrees of freedom are not eligible for consensus
aggregation.

The estimator uses the exact authoritative matrix entering differential model
fitting after approved workflow preprocessing. It must not read a hidden
alternate snapshot. Duplicate-correlation workflow provenance records the
`TableFingerprint` for that fitted matrix so the estimator input and final
differential model input can be audited as the same table. If imputed values are
present in that matrix, their participation is recorded explicitly in
provenance; imputed cells must not be silently reinterpreted as original
observations, and a pre-imputation matrix must not be silently substituted.

Failure outcomes are explicit and never silently choose another paired-design
policy, fixed block, ordinary least squares, or correlation zero:

- no repeated blocks;
- missing or empty block identities;
- insufficient observations relative to design rank;
- rank-deficient fixed-effects design;
- no feature with an estimable correlation;
- numerical non-convergence;
- invalid or non-positive-definite covariance;
- unsupported observation weighting;
- all eligible feature estimates being non-finite.

Authority boundaries:

- workflow validation owns public request admissibility and block metadata
  completeness;
- the experimental-design interpreter owns condition/covariate design assembly
  and excludes block columns for `duplicate_correlation`;
- the shared linear-model decomposition owns rank, conditioning, residual
  degrees of freedom, and contrast estimability;
- the duplicate-correlation estimator owns feature-wise REML estimates,
  Fisher-scale aggregation, consensus result, and covariance validity;
- the GLS fitter owns use of the consensus compound-symmetry covariance;
- empirical-Bayes moderation remains owned by the existing engine and is not
  rewritten for this policy.

A feature-level REML failure only means that feature does not contribute to the
consensus. It is not final differential row attrition by itself: the same
feature may still be fitted later using a successfully estimated consensus
correlation.

Internal frozen typed models define feature-level status, consensus result,
failure reasons, block-structure summary, and duplicate-correlation workflow
provenance. Internal estimator results may retain feature-wise estimates for
testing and aggregation. Public-sized workflow provenance uses a separate
feature-free consensus summary and must retain summaries and counts, not
thousands of feature-wise values.

## Consequences

- `fixed_block` and `duplicate_correlation` are clearly distinct valid models.
- Public request behaviour supports explicit `duplicate_correlation` selection
  without changing the default `reject` policy or auto-selecting from metadata.
- Scientific decisions stay in the differential science domain; workflow
  validation remains in the validation/workflow domains; dataset validation
  remains private.
- Provenance uses typed immutable contracts instead of arbitrary dictionaries.
- The first estimator method and `0.15` trim value are fixed by contract rather
  than accepted as user-facing tuning parameters.

## Alternatives Considered

1. Add `duplicate_correlation` to public config before the estimator, GLS,
   provenance, validation, and fixture parity were complete. Rejected because
   that would have exposed a selectable but unfinished policy.
2. Treat duplicate-correlation as a fixed-block variant. Rejected because the
   statistical assumptions and design matrix differ.
3. Fall back to OLS, fixed block, or correlation zero when estimation fails.
   Rejected because silent fallback changes the scientific model.
4. Make trim user-configurable now. Rejected because the first implementation
   uses the fixed published `15%` trimmed-mean contract.

## Implementation Notes

- Public paired-design constants include `reject`, `fixed_block`, and
  `duplicate_correlation` in `src/phospy/science/configs/differential.py` and
  are re-exported through the supported configuration facades.
- Internal typed contracts live in
  `src/phospy/science/differential/models/duplicate_correlation.py`.
- Public-sized workflow provenance uses the existing `TableFingerprint`
  provenance model for the authoritative fitting matrix and serializes through
  narrow `to_payload()` methods that omit retained feature-wise estimates.
- Linear-model decomposition remains in
  `src/phospy/science/differential/linear_model.py`.
- Empirical-Bayes moderation remains in
  `src/phospy/science/differential/empirical_bayes.py`.
- Related ADRs: [ADR-0019](adr_0019_experimental_design_and_contrast_contract.md)
  and [ADR-0044](adr_0044_differential_replicate_reliability_policy.md).
