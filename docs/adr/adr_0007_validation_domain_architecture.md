# ADR-0007: Validation Domain Architecture for PhosPy

## Status

- **ADR ID:** ADR-0007
- **Title:** Validation Domain Architecture for PhosPy
- **Status:** Accepted
- **Date:** 2026-05-13

## Context

PhosPy has three validation layers:

- shared structural primitives
- domain-specific validators
- workflow boundary validators

Without strict ownership, generic helpers become dumping grounds, scientific
rules drift into structural modules, and workflow validators take on execution
behavior.

Ownership is now maintained as an executable map in
`docs/validation-ownership.md` with invariant owner, enforcement point,
exclusions, and associated tests.

## Decision

Validation ownership is explicit and enforced by module boundaries:

1. `validation/common` owns structural primitives only (DataFrame type/shape,
   required columns, uniqueness, finite/missing checks, and generic alignment).
2. Domain modules own scientific and domain semantics (for example phosphosite
   identity, reference compatibility, replicate policy, and localisation
   eligibility).
3. Workflow validators compose structural and domain validators at request
   boundaries, then pass validated objects downstream.
4. Dataset construction owns dataset invariants. Validation-domain adapters may
   delegate through model construction, but do not duplicate model-boundary
   invariant logic.
5. Reference/resource model construction owns narrow structural invariants for
   the value object being constructed. Validation modules may re-export or
   compose those model-boundary validators, but the value object must not import
   `phospy.validation` to validate itself.
6. Workflow validators do not execute scoring/prediction/clustering science and
   do not perform data-transformation side effects.
7. Public presets/config objects are still required to pass the same validator
   boundaries as manually constructed configs.
8. Dataset validation modules under `phospy.validation.datasets` are internal
   support for builders, importers, preprocessing, and workflow boundaries.
   They are not promoted through `phospy.api` and are not a supported user
   validation route.
9. Public request DTOs remain passive command payloads. Their constructors store
   payloads only; request type compatibility, contextual compatibility, and
   scientific validation belong at builder and workflow boundaries.
10. Public contract/config scalar coercion helpers live under
    `phospy.contracts.configs` when contract construction needs them.
    `phospy.contracts` must not import `phospy.validation`.
11. Generic public contract value objects that do validate during construction
    enforce only context-free local invariants and raise
    `ContractValidationError`. Dataset and reference constructors keep their
    boundary-specific validation exceptions. Workflow validators own
    workflow-context and scientific validation and raise workflow validation
    errors.

The ownership map in `docs/validation-ownership.md` is part of ADR governance,
not optional commentary.

## Consequences

- **Positive**
  - Validation responsibilities are reviewable and auditable.
  - Domain rules remain near scientific ownership.
  - Workflow validators stay predictable and easier to test.
- **Negative**
  - Convenience refactors that move rules into unrelated modules should be
    rejected, even if they seem shorter.
  - Some existing helpers may need ownership cleanup when those files are
    touched.
- **Neutral**
  - Validation remains an internal architecture concern; this ADR clarifies
    governance rather than changing public API surfaces.
  - Constructor validation timing is predictable: passive request DTO
    construction never validates scientific compatibility, while validated
    value-object construction reports local contract failures with non-workflow
    exceptions.
  - Users validate data by constructing datasets with supported builders or by
    running workflows. Internal validators remain directly importable only for
    package implementation and focused tests.

## Alternatives Considered

1. Keep a broad "reusable validation" guideline without owner mapping.
   Rejected because it did not prevent drift.
2. Centralize all validation in workflow validators.
   Rejected because workflow validators would become oversized and domain rules
   would lose locality.
3. Move scientific validation into executors/interpreters.
   Rejected because boundary validation must fail early before execution.

## Implementation Notes

- Ownership registry: `docs/validation-ownership.md`.
- Shared structural primitives: `src/phospy/validation/common/dataframes.py`.
- Phosphosite-specific identifier/coherence validation owner:
  `src/phospy/science/sites/validation.py`.
- Differential design/contrast validation ownership:
  `src/phospy/validation/workflows/differential.py` plus
  `src/phospy/workflows/differential/validator.py`.
- Technical-replicate policy owner:
  `src/phospy/workflows/differential/replicates.py`.
- Intensity-scale establishment and dataset coherence owners:
  `src/phospy/validation/transformations/state.py` and
  `src/phospy/science/datasets/models.py`.
- Reference bundle and Kinase Library resource construction invariants:
  `src/phospy/science/references/models.py` and
  `src/phospy/science/references/kinase_library.py`; validation-package routes
  re-export these validators for internal compatibility.
- Workflow-boundary owners include
  `src/phospy/workflows/kinase/validator.py` and
  `src/phospy/workflows/signalome/validator.py`.
- Public API boundary checks live in
  `tests/unit/api/test_validation_not_public_api.py` and
  `tests/architecture/test_validation_boundaries.py`.
- Package dependency checks live in
  `tests/architecture/test_package_dependency_dag.py`.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R.,
& Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356. https://doi.org/10.1093/bioinformatics/btz306

YangLab. (n.d.). *PhosR* (Version release) [Computer software]. GitHub.
https://github.com/PYangLab/PhosR
