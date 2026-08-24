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

PhosPy intentionally uses the implemented limma-aligned residual-space
variance-component REML formulation. It does not define the estimator as a
direct one-dimensional profile-likelihood optimiser over correlation, and this
ADR does not claim algebraic equivalence to that alternative formulation.
PhosPy's scientific authority for this implementation is the residual-space
estimator in `src/phospy/science/differential/duplicate_correlation.py` together
with the version-pinned limma fixtures that verify feature-level estimates and
the consensus correlation.

For one feature, finite observations are selected first. Let `y_g` be the
observed feature vector, `X_g` the corresponding non-block fixed-effects design,
and `B_g` the corresponding block-indicator matrix. The rank of `X_g` is
recomputed after this feature-specific subsetting. A residual-space basis `Q2_g`
is formed for the orthogonal complement of `col(X_g)`, representing removal of
the fixed-effects design. The estimator then works with `z_g = Q2_g.T @ y_g`
and `W_g = Q2_g.T @ B_g`. After SVD of `W_g`, squared rotated residual-space
coordinates are fitted against a two-column component design containing an
intercept and the block eigenvalues. The fitted quantities are the residual
component and block component; the raw feature correlation is calculated after
that fit as:

```text
rho_raw = block_component / (residual_component + block_component)
```

Initial component coefficients come from least squares on the residual-space
component design. When the residual-space information is adequate, the
implementation refines them with the Gamma mean-variance iteration used by the
current estimator. Non-finite components, non-finite objective values, and
non-finite raw correlations become typed
`non_finite_objective_or_estimate` feature outcomes. Finite component fits that
do not converge become typed `optimisation_failed` feature outcomes.

Feature-specific missingness is part of the estimator contract. Removing
non-finite observations can change `X_g`, the observed block assignments, design
rank, residual degrees of freedom, the maximum observed repeated-block size, and
therefore the admissible lower correlation boundary for that feature. Features
that lose fixed-effect estimability, have two or fewer observed residual
degrees of freedom, have no repeated observations after subsetting, or have zero
or unusable residual variation receive typed failure statuses and do not
contribute to the consensus.

Successful feature estimates include both interior estimates and boundary
converged estimates. Their clamped correlations are Fisher transformed with
`atanh`, sorted, trimmed by `floor(n * 0.15)` from each tail, averaged on the
Fisher scale, and transformed back with `tanh`. Boundary-converged feature
estimates are valid consensus contributors.

The feature-level correlation clamp is applied after `rho_raw` is calculated
from the fitted components. For a feature whose maximum observed repeated-block
size is `m >= 2`, the compound-symmetry positive-definite lower boundary is
`-1 / (m - 1)`. PhosPy clamps feature estimates to:

```text
lower = -1 / (m - 1) + 0.01
upper = 0.99
```

The `0.01` lower-bound offset keeps lower-bound feature estimates strictly
inside the positive-definite compound-symmetry interval. The `0.99` upper cap
keeps estimates away from the unit-correlation singularity. Together these
clamps avoid singular or effectively singular covariance matrices, keep Fisher
aggregation finite and deterministic, preserve the limma-aligned numerical
behaviour verified by fixtures, and make boundary handling reproducible.

If the clamped estimate is within the boundary-detection tolerance of either
clamp, the feature status is `boundary_converged` and the boundary is recorded
as `lower` or `upper`; otherwise the feature status is `estimated`. Boundary
counts and the workflow-level positive-definite interval are retained in typed
diagnostics.

After Fisher aggregation, the consensus correlation is validated against the
full workflow block structure before GLS. This check uses the full maximum
repeated-block size, not any one feature's missingness-reduced block structure.
For full maximum repeated-block size `M >= 2`, the consensus lower bound is
`-1 / (M - 1) + 1e-10`; the upper bound is `0.99`. A consensus outside that
interval produces the typed
`invalid_or_non_positive_definite_covariance` outcome. The GLS fitter then
validates the supplied consensus again by requiring it to be finite, inside
`(-1, 1)`, strictly above the compound-symmetry lower boundary for the supplied
blocks, and Cholesky-factorable as a block covariance matrix.

For the estimator to run, the full non-block fixed-effects design must leave
more than two residual degrees of freedom, i.e. the analysed sample count must be
more than two larger than the design rank. Feature-level REML estimates with two
or fewer observed residual degrees of freedom are not eligible for consensus
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

The limma duplicate-correlation fixtures have explicit parity scopes. Fixtures
A-C validate the complete supported public path through feature-wise REML,
consensus aggregation, compound-symmetry GLS, contrasts, empirical Bayes, and
final moderated statistics. Fixture D validates estimator and GLS parity for
controlled feature-level missingness and failure cases, including raw
feature-specific GLS residual degrees of freedom and final feature-fit
statuses. Fixture D does not claim public empirical-Bayes parity: its actual
missing values are outside the public `AnalysisReadyPhosphoDataset` input
contract, and the internal duplicate-correlation executor also fails closed
before moderation if final GLS returns any non-fit feature status. PhosPy does
not expose a partial feature-specific moderated result after final GLS
failures. This fixture-D boundary does not weaken the complete-case public
parity claim made for fixtures A-C.

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
5. Replace the parity-verified residual-space variance-component REML estimator
   with a direct bounded scalar optimiser over correlation. Rejected for this
   corrective round because the implemented residual-space estimator already
   matches the version-pinned limma feature and consensus fixtures.

## Implementation Notes

- Public paired-design constants include `reject`, `fixed_block`, and
  `duplicate_correlation` in `src/phospy/science/configs/differential.py` and
  are re-exported through the supported configuration facades.
- Internal typed contracts live in
  `src/phospy/science/differential/models/duplicate_correlation.py`.
- Fixture parity scopes are declared in
  `tests/support/duplicate_correlation_parity_scopes.py` and checked against the
  committed fixture manifest by the release-gated fixture integrity tests.
- Public-sized workflow provenance uses the existing `TableFingerprint`
  provenance model for the authoritative fitting matrix and serializes through
  narrow `to_payload()` methods that omit retained feature-wise estimates.
- Linear-model decomposition remains in
  `src/phospy/science/differential/linear_model.py`.
- Empirical-Bayes moderation remains in
  `src/phospy/science/differential/empirical_bayes.py`.
- Related ADRs: [ADR-0019](adr_0019_experimental_design_and_contrast_contract.md)
  and [ADR-0044](adr_0044_differential_replicate_reliability_policy.md).
