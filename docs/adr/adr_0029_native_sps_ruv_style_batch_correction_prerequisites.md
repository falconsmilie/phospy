# ADR-0029: Native SPS/RUV-Style Batch Correction Prerequisites

## Status

- **ADR ID:** ADR-0029
- **Title:** Native SPS/RUV-Style Batch Correction Prerequisites
- **Status:** Accepted
- **Date:** 2026-06-23
- **Decision Type:** Scientific Architecture and Roadmap
- **Supersedes:** ADR-0027

## Context

PhosPy currently has two batch-related preprocessing surfaces:

- `linear_residualize_batch`, an opt-in fixed-effect residualisation method
  under dataset preprocessing.
- `ruv_readiness`, report-only metadata that helps audit whether future
  SPS/RUV-style preprocessing inputs might be present.

Neither surface selects stable or control phosphosites, estimates unwanted
variation factors, or applies RUV-style correction. Differential batch
covariates are also ordinary model terms, not preprocessing correction.

Native SPS/RUV-style correction would change the quantitative phosphosite
matrix consumed by downstream workflows. It therefore belongs on the
preprocessing side of the strict `AnalysisReadyPhosphoDataset` boundary, with
explicit validation, interpretation, execution, provenance, and tests before it
can become a supported feature.

## Problem Statement

Users may reasonably want batch-effect correction based on stable or control
phosphosites. A small flag such as `use_ruv=True` would be unsafe because the
method depends on control-site representation, batch and replicate structure,
missing-value handling, imputation policy, stage order, and downstream
eligibility.

PhosPy must define these contracts before implementing native SPS/RUV-style
correction. Until those prerequisites exist, PhosPy should keep stating that
RUV/SPS/RUV-III correction is not supported.

## Decision

PhosPy defers native SPS/RUV-style batch correction.

No public configuration flag, workflow, bundled control set, or numerical
correction kernel should be added until the prerequisites in this ADR are in
place. Future implementation, if added, must:

1. remain a preprocessing or normalisation concern before downstream workflows
   consume the matrix;
2. preserve the strict `AnalysisReadyPhosphoDataset` boundary;
3. use a validator -> interpreter -> executor workflow shape where a workflow
   is needed;
4. expose component and stage entrypoints through `run()`;
5. keep validation logic in validation-focused modules; and
6. record enough provenance to reproduce the correction decision.

This ADR defines prerequisites only. It does not implement RUV correction.

## Risks of Casual Batch-Correction Flags

Casual flags are risky because they can hide scientific choices that should be
audited. In particular, they can:

- treat arbitrary rows as stable controls without validating site identity,
  organism, source, or selection rules;
- run on confounded batch and condition structures where correction could remove
  biological signal;
- convert temporary imputation into apparent observed evidence;
- leave downstream workflows unable to tell which cells, sites, or samples were
  corrected or imputed;
- blur fixed-effect residualisation, diagnostic readiness reporting, and real
  RUV-style correction; and
- create support claims before contract tests and documentation can defend
  them.

## Relationship to Current Residualisation and Readiness Reporting

`linear_residualize_batch` remains the only executable batch-related
preprocessing method. It is fixed-effect residualisation that preserves
condition effects by including condition terms in the residualisation design.
It is not SPS selection, unwanted-factor estimation, RUV/SPS/RUV-III
correction, ComBat, limma `removeBatchEffect` parity, or mixed-effects
modelling.

`ruv_readiness` may continue to report metadata readiness signals. It must stay
report-only. It does not select controls, apply correction, modify the matrix,
or make a dataset partially corrected.

Differential workflow batch covariates remain model terms inside differential
analysis. They must not be documented as preprocessing correction.

## Required SPS and Control-Site Representation

Future native support requires an explicit representation for stable or control
phosphosites before any correction can run.

That representation must define:

1. how control status maps to `site_key` rows;
2. whether controls are caller-supplied, selected by a PhosPy preprocessing
   step, or packaged with PhosPy;
3. organism, identifier namespace, source, version, license, and redistribution
   status for any packaged control set;
4. selection method, thresholds, and exclusion reasons when controls are
   selected from the input data;
5. how ambiguous, duplicate, missing, or incompatible control annotations are
   rejected;
6. whether controls are binary labels, weighted controls, grouped controls, or a
   richer typed object; and
7. how control eligibility is recorded for downstream audit.

No future implementation should fetch control sites from online resources
behind the user's back.

## Missing-Value and Imputation Requirements

The analysis-ready dataset boundary requires a complete matrix, but a complete
matrix can contain values that were imputed upstream. Future SPS/RUV-style
correction must not treat imputed values as fully observed measurements unless
that policy is explicit and validated.

If a correction method needs a complete matrix internally, PhosPy must define:

- whether temporary imputation is allowed;
- which imputation method, parameters, and random seed are used;
- which rows or cells were originally missing;
- whether originally missing values are restored, masked, flagged, or withheld
  after correction;
- how row and sample eligibility changes after imputation;
- how the policy interacts with downstream differential imputation handling; and
- when missingness or imputation makes correction unsafe and should be rejected.

Temporary imputation for correction mechanics is not observed biological
evidence.

## Provenance Requirements

Future correction must record provenance for:

- requested method and resolved parameters;
- preprocessing stage order;
- control-site source, version, selection method, and selected `site_key` rows;
- batch, replicate, and design metadata used for correction;
- missing-value and imputation policy, including observation masks;
- input and output matrix fingerprints;
- before/after diagnostics and warnings;
- rejected or withheld rows and reasons; and
- PhosPy version and relevant dependency versions when available.

Provenance should live in the same provenance-oriented ownership areas used by
other workflows and preprocessing stages. Result tables alone are not enough.

## Validation Requirements

Future validation must happen before execution and must live in validation
modules or focused validation collaborators.

Validation must reject at least:

- missing, duplicate, or incompatible `site_key` control mappings;
- control sites with incompatible organism or identifier metadata;
- too few eligible controls for the requested method;
- confounded or rank-deficient batch, condition, replicate, or design
  structures;
- sample metadata that is missing, duplicated, or misaligned with the matrix;
- missing values without an explicit supported policy;
- imputation policies that cannot preserve observation provenance;
- unsupported stage orders;
- attempts to run correction after downstream workflows have consumed the
  matrix; and
- requests that would weaken the analysis-ready dataset boundary.

Validators should produce clear user-facing errors rather than silently dropping
controls, samples, or sites.

## Future Test Strategy

Any future implementation should add focused tests before claiming support:

- configuration and request validation tests for accepted and rejected cases;
- interpreter tests for resolved design, controls, missingness, and stage order;
- numerical executor tests with small deterministic fixtures;
- provenance tests for method parameters, controls, masks, diagnostics, and
  matrix fingerprints;
- dataset-builder integration tests that prove the output still satisfies
  `AnalysisReadyPhosphoDataset`;
- downstream workflow tests showing corrected datasets are consumed only through
  existing strict boundaries;
- documentation checks that `ruv_readiness` and `linear_residualize_batch`
  remain distinct from RUV-style correction; and
- external comparison fixtures only when their scope is narrow and documented
  as limited method evidence, not package-level equivalence.

## Future Ownership

Likely future ownership areas are:

```text
src/phospy/contracts/configs/
src/phospy/science/batch_correction/
src/phospy/science/datasets/preprocessing/
src/phospy/workflows/batch_correction/
src/phospy/validation/workflows/batch_correction/
src/phospy/provenance/
tests/
```

These paths are future-facing. This ADR does not require them to exist today.

Numerical method semantics should live in science-layer batch-correction code.
Workflow orchestration, if needed, should keep the validator -> interpreter ->
executor pattern. Validation modules should own scientific and contract
rejections. Differential, kinase, enrichment, and signalome workflows should
consume already corrected analysis-ready data only; they should not own
correction logic.

## Consequences

Positive consequences:

- Current docs can be precise about what is and is not supported.
- Future work has clear prerequisites instead of a single unsafe feature flag.
- The analysis-ready boundary and validation architecture remain intact.

Negative consequences:

- Users who need SPS/RUV-style correction must use external, explicit tooling
  before building a PhosPy analysis-ready dataset for now.
- Future implementation will require more than a numerical routine because it
  must include representation, validation, imputation policy, provenance, and
  tests.

Neutral consequences:

- `linear_residualize_batch` continues as limited fixed-effect residualisation.
- `ruv_readiness` continues as report-only readiness metadata.

## Non-Goals

This ADR does not:

- implement RUV, SPS selection, control-site selection, or RUV-III correction;
- add public API flags for RUV-style correction;
- add bundled stable or control phosphosite references;
- add hidden online fetching;
- add large dependencies;
- claim equivalence to PhosR;
- weaken preprocessing or dataset validation constraints; or
- move correction ownership into differential, kinase, enrichment, or signalome
  workflows.

## Related Records

- [ADR-0003: Analysis-Ready Dataset and Preprocessing Boundary for PhosPy](adr_0003-analysis_ready_dataset_and_preprocessing_boundary.md)
- [ADR-0007: Validation Domain Architecture for PhosPy](adr_0007_validation_domain_architecture.md)
- [ADR-0025: Competitive Phosphoproteomics Workflow Coverage Roadmap](adr_0025_competitive_phosphoproteomics_workflow_coverage.md)
- [ADR-0026: Differential Imputation Policy](adr_0026_differential_imputation_policy.md)
- [ADR-0027: Target Future Native PhosR-Style SPS/RUV-III Correction](adr_0027_target_future_native_phosr_style_sps_ruv_iii_correction.md)
- [Workflow Contracts](../workflow_contracts.md)
- [Scientific Coverage](../scientific-coverage.md)
