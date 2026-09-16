# PhosPy Release Notes

## Version 1.7.4

Release date: 2026-09-16.

These notes describe the changes since Version 1.7.3.

## Release Overview

PhosPy 1.7.4 adds an opt-in, PhosPy-native group-aware mixed-mechanism
missing-data policy:

```python
policy = "impute_group_aware"
```

The policy uses exact, aligned sample-group metadata and the original observed
missingness pattern to route eligible targets. It does not detect or prove MAR
or MNAR, and it does not claim numerical parity with PhosR `scImpute` or
`ptImpute`.

## Kinase Scientific-Policy Versions

The current implementation owns these policy and schema versions:

| Policy | Implemented version |
| --- | ---: |
| KSEA activity policy | 5 |
| Membership-selection policy | 4 |
| Inferential policy | 4 |
| Membership payload schema | 2 |
| Membership-independence policy | 2 |

These versions govern the KSEA scientific contract recorded in provenance and
bundles: membership evidence, whether substrate membership was selected
independently of the tested matrix, whether ordinary KSEA p/q output is
eligible, and compatibility for persisted membership and provenance payloads.
They are compatibility and interpretation contract identifiers, not empirical
proof of scientific validity.

## Group-Aware Imputation

- Complete groups require no imputation.
- Sufficiently supported partially observed groups route missing targets to
  KNN.
- Fully missing groups route to MinProb only when another group for the same
  site provides sufficient reference observation, expressing an explicit
  left-censored modelling assumption.
- Rows containing unsupported or ambiguous patterns are removed rather than
  forced through an imputation mechanism.

The policy requires aligned `sample_metadata`, a valid `group_column`,
established log2 quantitative input, explicit routing thresholds, KNN settings
(`k` and `distance="nan_euclidean"`), and MinProb settings (`q`, `width`, and
`seed`). Group-aware KNN uses strict `no_overlap_policy="error"`; the standalone
KNN column-mean fallback is not available in this mixed policy.

KNN and MinProb consume independent copies of the same original retained
matrix. Their target masks cannot overlap, and synthetic values from one
mechanism never become routing evidence, donor evidence, or numerical input for
the other. Originally observed values and the existing binary observation-mask
semantics are preserved.

See [Dataset Build Workflow](api/dataset-build-workflow.md#group-aware-knn-minprob)
for the complete user-facing configuration and interpretation contract.

## Diagnostics and Compatibility

- Typed routing, mechanism provenance, row audits, and diagnostics record the
  group-aware decision path.
- Diagnostics schema v2 preserves exact dataset-facing sample and group labels
  through serialization and bundle reconstruction.
- Historical schema-v1 parsing, normalization, round trips, and bundle
  reconstruction retain their previous semantics. Schema-v1 payloads are not
  migrated automatically, and schema v1 cannot represent
  `impute_group_aware` diagnostics.
- Existing standalone `forbid`, `impute_row_median`, `impute_knn`, and
  `impute_minprob` policies remain available with unchanged behavior.

## Scientific Scope

Observed missingness patterns are modelling evidence used to select an explicit
route; they do not establish the missingness mechanism for any value. The
MinProb route makes a left-censored assumption only for eligible asymmetric
fully missing groups. The partial-observation KNN route is a PhosPy mechanism,
not a reimplementation of PhosR `scImpute`.

The broader PhosPy scientific boundaries remain in force: differential analysis
is limited to tested design and contrast envelopes; bundled runtime references
are rat-only; enrichment is offline ORA over caller-supplied collections; native
SPS/RUV-style correction is not PhosR-equivalent RUV/SPS/RUV-III parity;
`duplicate_correlation` is a narrow paired-design GLS route rather than a
mixed-effects framework; and kinase/signalome outputs should be interpreted
through the documented workflow assumptions and caveats.

Next: [Quickstart](quickstart.md) or [API Guide](api/guide.md).
