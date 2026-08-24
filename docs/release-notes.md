# PhosPy Release Notes

## Version 1.7.1

Release date: 2026-08-24.

These notes describe the changes since Version 1.7.0.

## Release Overview

PhosPy 1.7.1 is a paired-design differential-analysis and release-documentation
release. It adds the explicit `duplicate_correlation` paired-design policy for
blocked differential designs, documents the final estimator and GLS contract,
names the exact version-pinned R/limma `duplicateCorrelation` fixture scope, and
corrects release documentation so MkDocs remains a standalone
documentation-maintenance command rather than a package release gate.

No breaking public API change is introduced relative to 1.7.0. The default
paired-design policy remains `reject`; PhosPy still does not infer or
automatically select a paired-design model from block metadata.

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

## Compatibility and Migration

- `DifferentialAnalysisConfig.paired_design_policy` now supports
  `"duplicate_correlation"` in addition to the existing `"reject"` and
  `"fixed_block"` values.
- `fixed_block` remains a valid supported design. It models block identity as
  ordinary fixed nuisance coefficients and does not estimate within-block
  correlation.
- `duplicate_correlation` is opt-in. It uses block identity as a correlation
  group, not as fixed block coefficients, and it is not selected automatically.
- Unsupported duplicate-correlation requests and fit failures fail closed. There
  is no fallback to `fixed_block`, ordinary least squares, or correlation zero.
- MkDocs remains available only through the standalone documentation path such
  as `make docs-build`. It is not part of `make release-check`, package
  building, wheel/source-distribution verification, or installed-distribution
  smoke checks.

## Major Additions

- Explicit `paired_design_policy="duplicate_correlation"` support for blocked
  differential designs through the public `DifferentialAnalysisRequest` and
  `DifferentialAnalysisConfig` contract. The supported constant is
  `PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION` from `phospy.advanced` and the
  configuration facades.
- Feature-wise duplicate-correlation estimation using PhosPy's implemented
  residual-space variance-component REML formulation. Eligible feature
  correlations are combined into one robust consensus on the Fisher `atanh`
  scale with the fixed 15% trim policy.
- Final duplicate-correlation fitting by compound-symmetry GLS. Contrasts and
  empirical-Bayes moderation reuse the existing differential pipeline after the
  GLS fit.
- Typed result provenance for duplicate correlation, including selected and
  normalized policy, covariance treatment and structure, estimator policy,
  matrix/design/block fingerprints, block summaries, consensus correlation,
  attempted/eligible/failed feature counts, typed failure summaries,
  convergence and boundary summaries, imputation participation, design rank,
  and GLS fit status.
- Version-pinned R/limma duplicate-correlation fixtures. Fixtures A-C validate
  the complete supported public path through final moderated output; fixture D
  validates estimator and GLS behavior for controlled feature-level
  missingness/failure cases only.
- Installed wheel and source-distribution smoke coverage now exercises the
  public duplicate-correlation workflow outside the source checkout.

## Fixes and Hardening

- The duplicate-correlation estimator contract now documents the residual-space
  REML formulation and correlation clamp policy. Feature-level estimates use
  the observed repeated-block size after feature-specific subsetting, clamp the
  lower bound to `-1 / (m - 1) + 0.01`, cap the upper bound at `0.99`, and
  validate the final consensus against the full workflow block structure before
  GLS.
- Validation rejects missing block IDs, designs without repeated blocks,
  fixed-block columns in duplicate-correlation designs, unsupported precision
  weights, rank-deficient non-block designs, insufficient residual degrees of
  freedom, and non-estimable contrasts before unsupported models can run.
- The workflow distinguishes feature-level REML failures from final GLS
  failures. A feature that fails to contribute to the consensus may still be
  fitted by GLS if a valid consensus is obtained from other features; a final
  GLS failure stops the duplicate-correlation workflow.
- Release process documentation now consistently states that successful package
  release checks do not validate rendered documentation. Documentation builds
  are maintained separately.

## Scientific Scope

`duplicate_correlation` is a narrow paired-design differential model. It uses
one shared consensus within-block correlation and compound-symmetry covariance
for the final GLS fit. It is not a general mixed-effects framework, does not
fit random slopes or multiple random effects, does not support arbitrary
longitudinal covariance, and does not combine fixed block coefficients with a
block-correlation structure.

The committed duplicate-correlation limma fixtures are implementation evidence
for their declared fixture scopes. They are not independent biological
validation and do not imply general limma equivalence.

The broader PhosPy scientific boundaries from 1.7.0 remain in force:
differential analysis is limited to tested design and contrast envelopes;
bundled runtime references are rat-only; enrichment is offline ORA over
caller-supplied collections; native SPS/RUV-style correction is not
PhosR-equivalent RUV/SPS/RUV-III parity; and kinase/signalome outputs should be
interpreted through the documented workflow assumptions and caveats.

Next: [Quickstart](quickstart.md) or [API Guide](api/guide.md).
