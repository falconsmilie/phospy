# ADR: Deterministic KNN Missing-Data Imputation

## Document Control

- **ADR ID:** ADR-0037
- **Title:** Deterministic KNN Missing-Data Imputation
- **Status:** Accepted
- **Date:** 2026-07-16
- **Last Amended:** 2026-08-07
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
- no-overlap behavior is explicit and versioned through
  `missing_data.no_overlap_policy`:
  - `column_mean_with_caveat` (policy version `1`) uses the retained-column
    mean as the deterministic fallback when no eligible donor has any shared
    observed column
  - `error` rejects the retained missing cell with a `PhosPyInputError`

The implementation must be vectorized or chunked and guarded. Current public
guardrails are:

- retained site rows `<= 50,000`
- sample columns `<= 64`
- estimated distance-feature operations `<= 2,000,000,000`
- target-by-donor distance chunks sized to about `48 MiB`
- release-gate KNN peak-memory budget `< 384 MiB`

KNN execution must return typed imputation-mechanism provenance owned by the
preprocessing missing-data stage. Every imputed cell must belong to exactly one
mechanism:

- nearest-neighbour-derived imputation mask
- retained-column-mean fallback imputation mask

The outcome and serialized missing-data diagnostics must include exact
mechanism counts, affected row IDs, affected column IDs, stable mechanism mask
hashes, configured no-overlap policy and policy version, and rows whose imputed
cells were entirely column-mean fallback-imputed. If column-mean fallback is
used, diagnostics must carry a typed caveat code rather than only free text.

## Consequences

KNN imputation is suitable for sparse retained missingness at large site counts
and for moderate retained missing-target workloads that remain within the
distance-work guardrail. It fails fast for broad missingness that would require
impractical all-pairs distance work. The error must report retained rows, rows
with missing values, sample columns, estimated work, and the configured budget.

The policy does not silently fall back to row-median, MinProb, or scikit-learn
behavior. Users must choose a different missing-data policy explicitly when KNN
is outside the guardrail.

Workflow orchestration must not implement KNN imputation details. Workflows may
consume preprocessing-owned imputation observation metadata, but KNN donor,
distance, tie, and fallback semantics stay in the preprocessing missing-data
domain.

Differential or downstream workflows must not reconstruct fallback cells from
numeric imputed values. They may consume dataset-owned KNN imputation facts from
preprocessing provenance.

## Validation

Required coverage:

- parity with small reference-loop fixtures
- deterministic donor-tie fixtures
- all-missing retained rows
- sparse-overlap and no-overlap fallback cases
- explicit `error` no-overlap policy rejection
- mixed nearest-neighbour and column-mean fallback cells in one retained row
- deterministic mechanism masks and hashes under row reordering
- bundle serialization round-trips for KNN mechanism diagnostics
- row audit and report caveat coverage for column-mean fallback
- guardrail rejection for impractical requests
- 10k, 25k, and 50k-site release-gate benchmarks for sparse 12-sample and
  moderate 24-sample retained missing-target tiers
- output equivalence across distance chunk sizes
- peak-memory regression coverage
- workflow-orchestration ownership checks
