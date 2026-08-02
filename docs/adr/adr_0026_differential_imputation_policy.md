# ADR-0026: Differential Imputation Policy

## Status

- **ADR ID:** ADR-0026
- **Title:** Differential Imputation Policy
- **Status:** Accepted
- **Date:** 2026-06-14

## Context

Analysis-ready datasets may be complete because upstream preprocessing imputed
missing phosphosite intensities. A complete imputed matrix is not the same as a
fully observed matrix. Treating imputed cells as ordinary observations by
default would make residual degrees of freedom, fitted effects, and p-values
look more certain than the measurement evidence supports.

Dataset preprocessing now preserves imputation facts: an observed-cell mask and
per-feature `imputed_cell_count`, `observed_cell_count`, and
`imputed_fraction`. Preprocessing also records the imputation input scale and
whether imputation happened before or after intensity transformation. These are
dataset-owned scientific facts; differential analysis needs an explicit policy
for any use of imputed datasets.

## Decision

`DifferentialAnalysisConfig.imputed_value_policy` defaults to `"reject"`.
Imputed datasets are rejected unless callers explicitly request a supported
non-default policy.

Imputation science remains owned by preprocessing. Missing-data configuration
and preprocessing plan interpretation decide the imputation method, required or
selected input scale, operation order, and observation-mask preservation.
Differential analysis does not reinterpret or rerun imputation; it only
validates downstream eligibility from the dataset-owned imputation metadata.

The supported non-default policy is
`"withhold_imputed_features"`. It requires dataset-owned imputation observation
metadata. During interpretation, PhosPy computes feature-level imputation
metadata over the actual analysed sample subset and assigns one status per
feature:

- `tested`
- `withheld_high_imputation`
- `withheld_insufficient_observed_values`

Rows are withheld when their imputed fraction is above
`imputed_value_max_fraction`, or when any requested contrast condition has fewer
originally observed samples than `minimum_condition_replicates`.

Only `tested` rows are passed to the statistical executor. Withheld rows are
reattached to public contrast tables with missing `logFC`, `t`, `P.Value`, and
`adj.P.Val`.

Benjamini-Hochberg correction is computed only over tested rows. Withheld and
otherwise non-testable rows are not included in the denominator.

`observed_only` fitting is not supported in this ADR. PhosPy does not subtract
imputed cells from residual degrees of freedom while still fitting imputed
values and does not claim feature-specific residual degrees of freedom for this
policy.

## Consequences

- **Positive**
  - Default behavior remains scientifically conservative.
  - Non-default imputed differential analysis is explicit and auditable.
  - Result tables expose imputation counts, fractions, policy, threshold, and
    status.
  - Adjusted p-value behavior is deterministic and test-pinned.
- **Negative**
  - Projects that want observed-only modelling still need a future
    implementation with per-feature design fitting and residual degrees of
    freedom.
- **Neutral**
  - The withhold policy may still fit a tested feature on an analysis-ready row
    containing imputed values when the imputed fraction is within the configured
    limit and observed support is sufficient. That limitation is recorded in
    provenance and documentation.

## Implementation Notes

- Public config:
  `src/phospy/contracts/configs/differential.py`
- Eligibility validation:
  `src/phospy/validation/workflows/differential.py`
- Imputation metadata alignment and status assignment:
  `src/phospy/workflows/differential/interpreter.py`
- Tested-row statistical execution and public result reattachment:
  `src/phospy/workflows/differential/executor.py`
- Result/provenance/table models:
  `src/phospy/science/differential/models/results.py`,
  `src/phospy/science/differential/models/provenance.py`, and
  `src/phospy/science/differential/models/tables.py`
- Provenance construction:
  `src/phospy/workflows/differential/provenance.py`
- Tests:
  `tests/unit/test_differential_imputation_policy.py`
