# PhosPy Release Notes

## Version 1.7.6

Release date: 2026-09-26.

These notes describe only the changes since Version 1.7.5.

## Release Overview

PhosPy 1.7.6 makes two important workflow boundaries easier to discover and
reduces the preparation cost of native RUV-style correction. It does not add or
remove public symbols, change configuration defaults, change a supported
scientific estimator, or broaden the package's external parity claims.

## Differential Input Boundaries

Dataset imputation and differential fitting remain separate decisions. A
dataset produced with `impute_row_median`, `impute_minprob`, `impute_knn`, or
`impute_group_aware` can carry imputed cells even though its numeric matrix is
complete. `DifferentialAnalysisWorkflow` rejects such a dataset by default.

The existing explicit alternative,
`imputed_value_policy="withhold_imputed_features"`, first withholds a feature
when its imputed fraction exceeds `imputed_value_max_fraction` or when a
condition used by any requested contrast has fewer than
`minimum_condition_replicates` originally observed values. Remaining tested
features are fitted on the workflow-approved matrix. If the configured positive
threshold retains imputed cells, those values participate in fitting; this is
not observed-only fitting. The dataset guide, quickstart, differential guide,
and public docstrings now state this boundary before users select imputation.

Suspicious declared-log2 input also has two independent acknowledgements.
`DatasetBuildRequest.allow_suspicious_declared_input_intensity_scale=True`
permits dataset construction and preserves the warning in provenance. It does
not authorize scale-sensitive differential output. Differential analysis still
requires the separate existing
`DifferentialAnalysisConfig.allow_suspicious_declared_input_scale=True`
override. Builder warnings and differential validation errors now name that
scope and recovery path directly; the underlying default rejection policy is
unchanged.

## RUV-Style Preparation Internals

Governed observation-mask materialisation and temporary row-median completion
now use array-oriented implementations in both RUV-III and native SPS/RUV-style
correction paths; mask materialisation is factored into a shared internal
helper. The scientific and numerical contracts are unchanged: exact
scalar-reference and complete-executor tests cover prepared and corrected
matrices, missingness classification, restored positions, ordered coordinates,
statuses, diagnostics, warnings, provenance, and fingerprints.

Checked-in same-machine reports record materially lower preparation times for
the governed missingness cases at small, representative, and stress tiers. The
reports also record the temporary dense-workspace memory tradeoff and a separate
complete-input fast-path comparison. These measurements are local benchmark
evidence, not a general runtime guarantee; kernel timing is outside their
preparation measurement.

The release adds a focused benchmark command:

```bash
make benchmark-ruv-iii-missingness
```

See [RUV-III missingness-preparation evidence](https://github.com/falconsmilie/phospy/blob/main/benchmarks/evidence/ruv-iii-missingness-preparation-2026-09-22.md)
for the environment, exact measurements, source hashes, and scope. The README
also now includes the project repository banner.

## Compatibility and Scientific Scope

No compatibility shim or migration is required. Existing public request and
configuration names remain in place, and default differential eligibility and
RUV-style scientific behaviour remain unchanged. The RUV preparation work does
not establish new PhosR or `ruv` parity; the bounded claims documented in
[Parity](parity.md) and [Scientific Coverage](scientific-coverage.md) still
apply.

Next: [Quickstart](quickstart.md) or [API Guide](api/guide.md).
