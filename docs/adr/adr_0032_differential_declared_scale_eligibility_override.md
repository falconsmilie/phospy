# ADR-0032: Differential Declared Scale Eligibility Override

## Status

- **ADR ID:** ADR-0032
- **Title:** Differential Declared Scale Eligibility Override
- **Status:** Accepted
- **Date:** 2026-07-01

## Context

Differential analysis reports `logFC`, so it requires established log2
phosphosite intensities. Dataset construction can establish log2 state by
applying a supported PhosPy transformation, by preserving a caller declaration,
or by restoring trusted bundle provenance.

Declared input scale can be suspicious even when it is structurally established.
The dataset builder records suspicious declared-scale diagnostics in
`IntensityScaleEstablishmentProvenance.diagnostic_warnings`. Directly
constructed or bundle-restored datasets can therefore carry enough provenance
for downstream workflow eligibility to make a conservative decision.

## Decision

`DifferentialAnalysisConfig.allow_suspicious_declared_input_scale` defaults to
`False`.

Differential eligibility rejects an established log2 dataset when all of the
following are true:

- `dataset.intensity_scale_state.establishment_provenance` exists
- the provenance mode is `declared`
- the provenance records one or more diagnostic warnings
- the differential override is not set

Declared log2 datasets without diagnostic warnings remain valid. PhosPy
transformed log2 state remains valid. Bundle-restored trusted provenance remains
valid unless it records declaration mode and warnings.

When the override is set, differential policy provenance records
`statistical_testing.allow_suspicious_declared_input_scale=True`.

## Consequences

- **Positive**
  - Differential analysis is conservative by default for suspicious declared
    log2 input.
  - The override is explicit in public configuration and machine-readable
    result policy provenance.
  - Builder diagnostics remain owned by dataset construction; differential only
    consumes recorded provenance.
- **Negative**
  - Trusted suspicious declarations require an additional explicit workflow
    option.
- **Neutral**
  - This policy does not alter the analysis-ready dataset boundary or site
    sequence requirements.

## Implementation Notes

- Public config:
  `src/phospy/contracts/configs/differential.py`
- Eligibility validation:
  `src/phospy/validation/workflows/differential.py`
- Workflow request validation:
  `src/phospy/workflows/differential/validator.py`
- Policy provenance:
  `src/phospy/science/differential/models/provenance.py` and
  `src/phospy/workflows/differential/provenance.py`
- Tests:
  `tests/unit/test_differential_analysis.py` and
  `tests/integration/workflows/differential/test_differential_result_provenance.py`
