# PhosPy Release Notes

## Version 1.7.3

Release date: 2026-09-10.

These notes describe the changes since Version 1.7.2.

## Release Overview

PhosPy 1.7.3 adds opt-in quantification-depth-aware empirical-Bayes
moderation for differential analysis. This mode can model feature-specific
empirical-Bayes prior variance against per-site quantification depth, using
PSM count or peptide count metadata supplied for the tested sites.

Existing empirical-Bayes behavior remains the default. `EmpiricalBayesConfig()`
still uses the global prior, and `EmpiricalBayesConfig(trend=True)` still uses
the existing mean-intensity variance trend. Quantification-depth moderation is
selected only when the caller explicitly requests the quantification-depth trend
covariate and declares the depth kind.

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

- No migration is required for existing empirical-Bayes users. The global and
  mean-intensity trend modes keep their existing configuration semantics.
- Quantification-depth-aware moderation is selected with
  `EmpiricalBayesConfig(trend=True, trend_covariate="quantification_depth",
  quantification_depth_kind="psm_count")` or
  `quantification_depth_kind="peptide_count"`.
- The authoritative metadata source is
  `site_metadata["quantification_depth"]`. Valid finite integer counts greater
  than or equal to one are required for every feature that enters differential
  testing and empirical-Bayes moderation. Features withheld before testing may
  still appear in re-expanded contrast tables with missing
  quantification-depth trend diagnostics.
- The depth trend is fit against `log2(quantification_depth)`. Result
  diagnostics expose both the raw depth and the transformed trend covariate.
- MaxQuant and FragPipe/PTMProphet importers can populate quantification-depth
  metadata only when callers explicitly map a source column and declare the
  depth kind. The importers do not infer PSM counts from row counts, spectrum
  identifiers, PSM identifiers, or arbitrary numeric columns.

## Major Additions

- Added quantification-depth-aware empirical-Bayes moderation across ordinary
  differential analysis, robust empirical Bayes, duplicate-correlation
  analysis, and the experimental protein-aware differential lane.
- Added separate `QuantificationDepthTrendDiagnostics` so depth trend
  diagnostics are reported independently from mean-intensity trend diagnostics.
- Added empirical-Bayes provenance that records the selected trend covariate,
  depth kind, depth transformation, empirical-Bayes method, robust mode, and
  winsor settings as applicable.
- Added conservative site-level quantification-depth handling for MaxQuant and
  FragPipe/PTMProphet imports. Unambiguous identical mapped counts are
  preserved; ambiguous, split, or conflicting site evidence leaves depth
  unavailable for downstream validation.
- Added a checked-in R/DEqMS `spectraCounteBayes` reference fixture and
  release-gated tests for DEqMS-inspired quantification-depth-aware empirical
  Bayes.

## Fixes and Hardening

- Invalid or non-finite empirical-Bayes trend covariates are rejected before
  fitting.
- `site_metadata["quantification_depth"]` is the authority for
  quantification-depth trend moderation.
- Ordinary differential regression coverage protects the default global prior,
  the existing mean-intensity trend, duplicate-correlation execution, importer
  opt-in behavior, provenance payloads, and re-expanded result diagnostics.

## Scientific Scope

The new mode is quantification-depth-aware empirical Bayes inspired by DEqMS.
It retains PhosPy's empirical-Bayes trend architecture and should not be
interpreted as exact numerical or API compatibility with DEqMS.

Quantification depth is feature metadata for the sites entering differential
testing and moderation. PhosPy does not derive this metadata from unrelated
evidence structure, and missing or invalid depth for a tested site fails closed
instead of silently falling back to another trend.

The broader PhosPy scientific boundaries remain in force: differential analysis
is limited to tested design and contrast envelopes; bundled runtime references
are rat-only; enrichment is offline ORA over caller-supplied collections; native
SPS/RUV-style correction is not PhosR-equivalent RUV/SPS/RUV-III parity;
`duplicate_correlation` is a narrow paired-design GLS route rather than a
mixed-effects framework; and kinase/signalome outputs should be interpreted
through the documented workflow assumptions and caveats.

Next: [Quickstart](quickstart.md) or [API Guide](api/guide.md).
