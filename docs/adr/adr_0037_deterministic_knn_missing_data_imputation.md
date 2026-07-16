# ADR: Deterministic KNN Missing-Data Imputation

## Document Control

- **ADR ID:** ADR-0037
- **Title:** Deterministic KNN Missing-Data Imputation
- **Status:** Accepted
- **Date:** 2026-07-16
- **Decision Type:** Architecture Decision Record

## Context

PhosPy exposes missing-data handling through dataset preprocessing. The KNN
policy had to remain semantically stable while becoming practical enough for
release-gated preprocessing workloads. Replacing it with scikit-learn without a
formal semantic parity decision would change donor selection, tie behavior, or
fallback semantics in ways that are difficult for users to audit.

## Decision

`missing_data.policy="impute_knn"` remains a custom deterministic preprocessing
implementation owned by
`src/phospy/science/datasets/preprocessing/stages/missing_data/knn.py`.

The implementation must document and preserve these semantics:

- rows above `max_missing_fraction_per_row` are dropped before donor search
- only `distance="nan_euclidean"` is supported
- donors for a missing cell must be observed in that cell's column
- distance uses target/donor shared observed columns and scales by
  `n_columns / shared_observed_column_count`
- exact donor ties are ordered by `(str(row_id), original_position)`
- selected donor values are averaged without distance weighting
- when no donor has any shared observed column, the retained-column mean is the
  deterministic fallback

The implementation must be vectorized or chunked and guarded. Current public
guardrails are:

- retained site rows `<= 50,000`
- sample columns `<= 64`
- estimated distance-feature operations `<= 2,000,000,000`
- target-by-donor distance chunks sized to about `96 MiB`
- release-gate KNN peak-memory budget `< 384 MiB`

## Consequences

KNN imputation is suitable for sparse retained missingness at large site counts
and fails fast for broad missingness that would require impractical all-pairs
distance work. The error must report retained rows, rows with missing values,
sample columns, estimated work, and the configured budget.

The policy does not silently fall back to row-median, MinProb, or scikit-learn
behavior. Users must choose a different missing-data policy explicitly when KNN
is outside the guardrail.

Workflow orchestration must not implement KNN imputation details. Workflows may
consume preprocessing-owned imputation observation metadata, but KNN donor,
distance, tie, and fallback semantics stay in the preprocessing missing-data
domain.

## Validation

Required coverage:

- parity with small reference-loop fixtures
- deterministic donor-tie fixtures
- all-missing retained rows
- sparse-overlap and no-overlap fallback cases
- guardrail rejection for impractical requests
- 10k, 25k, and 50k-site release-gate benchmarks
- peak-memory regression coverage
- workflow-orchestration ownership checks
