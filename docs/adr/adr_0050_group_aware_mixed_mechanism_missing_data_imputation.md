# ADR: Group-Aware Mixed-Mechanism Missing-Data Imputation

## Document Control

- **ADR ID:** ADR-0050
- **Title:** Group-Aware Mixed-Mechanism Missing-Data Imputation
- **Status:** Accepted
- **Date:** 2026-09-15
- **Decision Type:** Architecture Decision Record

## Context

Whole-matrix missing-data policies force every eligible missing value through
one numerical mechanism. Experimental group structure can instead support an
explicit, auditable routing decision before imputation. Missingness patterns
are evidence for that routing decision; they are not proof that values are
missing at random (MAR) or missing not at random (MNAR).

Because one route uses MinProb, every route must consume the same explicitly
established log2 matrix. Running KNN on linear values and MinProb on log2 values
within one policy would make the routing result quantitatively incoherent.

## Decision

PhosPy exposes `missing_data.policy="impute_group_aware"` as an opt-in
policy. Callers must name the sample-metadata group column explicitly; group
membership is never inferred from sample names.

Routing decisions are made from the original missingness pattern before any
numerical imputation:

- partial, sufficiently supported groups are assigned to KNN
- strongly asymmetric, fully missing groups with adequate reference support
  are assigned to MinProb
- ambiguous cases are not forcibly imputed

KNN and MinProb operate only on the target cells explicitly assigned to them.
Synthetic values produced by either mechanism must never become routing
evidence, donor evidence, reference evidence, or numerical input for the other
mechanism.

The policy requires log2 input, MinProb parameters (`q`, `width`, and `seed`),
KNN parameters (`k` and `distance="nan_euclidean"`), explicit group-routing
thresholds, and `no_overlap_policy="error"`. The standalone KNN column-mean
fallback is not available because it would add a third numerical mechanism to
the mixed routing contract. Group-aware routing owns row eligibility, so the
standalone `max_missing_fraction_per_row` control is not accepted.

An omitted `no_overlap_policy` resolves to `"error"`; callers may also supply
`"error"` explicitly. No other value is supported for this policy.

The missing-data stage implements the complete policy. It captures the original
missingness profile, resolves groups from aligned `dataset.sample_metadata`,
routes cells, and removes unsupported rows before numerical imputation. KNN and
MinProb each consume the same original retained matrix independently; their
disjoint results are merged only after both mechanisms finish. The stage then
validates target/imputation-mask invariants, preservation of originally observed
values, and a complete output matrix.

## Consequences

Automatic preprocessing order places a configured log2 transform before
group-aware missing-data processing. If input is explicitly declared already
log2, no additional transform is required.

Existing standalone `forbid`, row-median, MinProb, and KNN policies retain
their current validation and numerical semantics, including the standalone KNN
fallback documented by ADR-0037.

The group-aware KNN route always uses the conservative no-overlap behavior: a
target without a valid overlapping donor fails with actionable site/cell
context. It never uses the standalone KNN column-mean fallback.

The policy consumes `dataset.sample_metadata` in addition to the phospho and
site-metadata tables; standalone missing-data policies do not acquire that
dependency. Stage provenance records this plan-dependent dependency, the random
seed, routing target-mask evidence, and that both mechanisms consumed the
original retained matrix. Because MinProb is part of the configured policy,
group-aware execution is `SEEDED_STOCHASTIC` even when an input routes no cells
to MinProb. The dataset-level observation mask remains binary (`True` for
originally observed and `False` for imputed), while mechanism-specific routing
and mask details remain in audit diagnostics.

The policy does not claim PhosR numerical parity.

## Validation

Contract tests cover required routing, KNN, and MinProb parameters;
threshold domains; the log2-only scale rule; rejection of standalone row
eligibility and column-mean fallback controls; seeded-stochastic planning; and
automatic stage ordering. Numerical integration tests prove mixed-mechanism
execution, filtering before both mechanisms, preservation of observed values,
complete output, observation-mask semantics, seeded reproducibility, and that
neither mechanism consumes the other's synthetic outputs.
