# ADR-0027: Target Future Native PhosR-Style SPS/RUV-III Correction

## Status

- **ADR ID:** ADR-0027
- **Title:** Target Future Native PhosR-Style SPS/RUV-III Correction
- **Status:** Accepted
- **Date:** 2026-06-14
- **Decision Type:** Architecture and Scientific Roadmap

## Context

PhosPy already exposes limited batch-related functionality. Dataset
preprocessing can run `linear_residualize_batch`, which is fixed-effect
residualisation of batch terms while preserving condition effects by design.
Dataset construction can also report `ruv_readiness` metadata when requested.

Those features do not implement PhosR-style stable-phosphosite (SPS) selection,
RUV model fitting, or RUV-III correction. They also do not provide output
parity with PhosR's `RUVphospho` lane. PhosR's documented RUV path uses stable
phosphosites as controls, a design matrix, an unwanted-factor count, and
temporary imputation behavior around RUV-III normalisation.

Users coming from PhosR therefore need a direct answer: current PhosPy does not
support native RUV/SPS/RUV-III correction, but this is a scientifically
important preprocessing direction that should be planned explicitly.

## Decision

PhosPy will add native PhosR-style SPS/RUV-III correction as a future
preprocessing/normalisation implementation epic.

This is a future architecture commitment, not a current feature-support claim.
Until executable implementation, public contracts, documentation, validation,
and tests exist, PhosPy must continue to state that native RUV/SPS/RUV-III
correction is not supported.

Future correction logic belongs before downstream scientific workflows consume
the quantitative matrix. It must not be owned by differential analysis, kinase
analysis, enrichment, or signalome execution.

## Current Behaviour

Current PhosPy behavior is:

- no native RUV/SPS/RUV-III correction support
- no PhosR-equivalent SPS/RUV-III batch-correction lane
- no bundled SPS/control reference set for correction
- no RUV-III numerical correction kernel
- no batch-correction workflow for RUV/SPS/RUV-III
- no PhosR parity fixtures for RUV/SPS/RUV-III corrected outputs

`linear_residualize_batch` remains the only executable preprocessing method in
the `batch_correction` config group. It is fixed-effect residualisation, not
SPS selection, not RUV factor estimation, and not RUV-III correction.

Differential fixed-effect batch covariates are ordinary model terms. They are
not preprocessing correction and are not equivalent to PhosR-style RUV/SPS
correction.

`ruv_readiness` and similarly named diagnostics are report-only readiness
signals. They must not be described as correction, partial correction, or
support for RUV/SPS/RUV-III.

## Future Architecture Direction

Future SPS/RUV-III support should be implemented as preprocessing and
normalisation infrastructure that produces a corrected analysis-ready matrix
with explicit provenance and eligibility metadata.

Likely future ownership areas are:

```text
src/phospy/science/batch_correction/
src/phospy/science/normalisation/
src/phospy/science/references/
src/phospy/workflows/batch_correction/
src/phospy/validation/workflows/batch_correction/
src/phospy/provenance/
tests/
```

These paths are future-facing suggestions only. This ADR does not create them
and does not require them to exist today.

Expected future responsibilities:

- science-layer correction contracts and numerical kernels should own
  SPS/RUV-III model semantics.
- normalisation/preprocessing orchestration should decide when correction runs
  relative to missing-data handling, scaling, total-protein correction, and
  dataset construction.
- reference modules should own user-supplied or bundled SPS/control-set
  contracts and provenance checks.
- workflow modules should own a validator -> interpreter -> executor path only
  after public request/result contracts are defined.
- validation modules should reject incomplete, confounded, or scientifically
  unsafe correction requests before execution.
- provenance modules should record enough information to reproduce correction
  decisions and downstream eligibility.
- tests should include method-level validation, workflow contract checks,
  before/after diagnostics, imputation-mask behavior, and PhosR comparison
  fixtures where practical.

## Scientific Constraints

Future support must define explicit contracts for:

1. SPS/control phosphosite selection.
2. User-supplied or bundled control sets, including source provenance,
   redistribution status, organism compatibility, identifier semantics, and
   versioning.
3. Batch labels and replicate/design structure.
4. RUV/SPS model parameters, including unwanted-factor count where applicable.
5. Temporary imputation rules if the correction method requires a complete
   matrix.
6. Restoration, removal, or flagging of originally missing/imputed values after
   correction.
7. Before/after diagnostics.
8. Downstream workflow eligibility rules.
9. Provenance sufficient to reproduce correction decisions.
10. Parity or validation fixtures against PhosR where practical.

Future implementation must also align with PhosPy's imputation policy. If
RUV/SPS/RUV-III correction requires temporary imputation to create a complete
matrix, those imputed values must not silently become fully observed biological
measurements downstream. The implementation must do at least one of the
following:

1. restore missingness after correction;
2. carry an observation/imputation mask forward;
3. flag imputed positions in downstream outputs; or
4. withhold features from downstream statistical testing where imputation makes
   inference unsafe.

Temporary imputation for correction mechanics is not equivalent to observed
data.

## Consequences

Positive consequences:

- PhosPy has a clear long-term direction for users who need PhosR-style
  SPS/RUV-III correction.
- Current documentation can be honest that residualisation and readiness
  diagnostics are not correction support.
- Future implementation work has a preprocessing/normalisation ownership
  boundary before downstream workflows.

Negative consequences:

- PhosR users still need external PhosR or other tooling for SPS/RUV-III
  correction until a future implementation lands.
- The future epic is scientifically larger than adding a single numerical
  routine because it requires control-set provenance, design contracts,
  missingness handling, diagnostics, workflow eligibility, and validation
  evidence.

Neutral consequences:

- `linear_residualize_batch` continues to exist as limited fixed-effect
  residualisation.
- `ruv_readiness` can continue to report metadata readiness, but only with
  wording that keeps it separate from correction support.

## Future Implementation Work

Future tickets should be split into separate implementation steps:

1. Define batch-correction workflow contracts.
2. Define SPS/control-set input contracts and provenance.
3. Define replicate/design mapping contracts.
4. Define temporary-imputation and missingness-restoration policy.
5. Implement a native RUV-III numerical kernel.
6. Add a validator -> interpreter -> executor batch-correction workflow.
7. Add before/after diagnostics.
8. Add PhosR parity fixtures where practical.
9. Document supported and unsupported correction modes.

## Non-Goals

This ADR does not:

- implement RUV, SPS selection, control-set selection, or RUV-III correction
- add bundled SPS/control references
- add a batch-correction workflow
- add numerical correction kernels
- change differential analysis behavior
- change existing residualisation behavior
- claim PhosR parity
- move correction ownership into differential, kinase, enrichment, or
  signalome workflows
- add public API options for a future feature

## References

Bioconductor. (2025). *PhosR: A set of methods and tools for comprehensive
analysis of phosphoproteomics data* (Version 1.20.0) [R package manual].
https://bioconductor.org/packages/release/bioc/manuals/PhosR/man/PhosR.pdf

Gagnon-Bartsch, J. (2019). *ruv: Detect and remove unwanted variation using
negative controls* (Version 0.9.7.1) [R package].
https://CRAN.R-project.org/package=ruv

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R.,
& Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.
