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

Update note (2026-07-17, validation ownership under package DAG): Validation
remains private and independently owned. Science and table modules no longer
import concrete `phospy.validation` implementations; they either own local
domain invariants or consume protocols wired by API/workflow adapters.

Update note (2026-07-29, phosphosite identity and sequence-context ownership):
Phosphosite identity contracts and workflow sequence-context semantics are
science-owned. `phospy.science.sites.identity_contracts` owns the concrete
`PhosphositeIdentityContract` class and shared identity enforcement.
`phospy.science.sites.sequence_context` owns the concrete
`SequenceContextContract` class, centered sequence-context enforcement,
residue/window policy, and known sequence-source checks. Validation-package
routes may compose or identity-preservingly re-export those objects, but must
not define copies.

Update note (2026-07-30, executable ownership governance): The ownership map is
now checked by architecture tests for current module paths, concrete
symbol-level duplicate definitions across owner trees, identity-preserving
compatibility re-exports, workflow-validator composition wording, and
import-graph diagnostics that report the AST line of the import.

## Decision

Validation ownership is explicit and enforced by module boundaries:

1. `phospy.frames` owns neutral structural dataframe primitives (DataFrame
   type/shape, required columns, uniqueness, finite/missing checks, and generic
   alignment). `validation/common` remains a compatibility route and must not
   regain ownership of neutral frame behavior.
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
12. Validation modules may depend on consumer-owned protocols and science-owned
    value models, but concrete workflow/build/preprocessing collaborators must
    not depend back on validation implementations.
13. Validators do not own contextual interpretation. They may enforce
    construction-time and boundary invariants for public config objects, but
    interpreter stages own translation from public DTOs to resolved execution
    models.
14. Generic scalar configuration validation primitives live in the lowest
    applicable shared layer, such as `phospy.policies`. Domain-specific config
    helpers may live with their domain owner and be re-exported by higher
    layers, but duplicated helper implementations are not allowed.

The ownership map in `docs/validation-ownership.md` is part of ADR governance,
not optional commentary. Updates to ADR-0007 and that map must be kept in sync
when ownership or compatibility routes change.

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

- Ownership registry: `docs/validation-ownership.md`, enforced by
  `tests/architecture/test_validation_ownership_governance.py`.
- Shared structural primitives: `src/phospy/frames/validation.py`; legacy
  validation routes under `src/phospy/validation/common/` are
  identity-preserving compatibility wrappers.
- Phosphosite-specific identifier/coherence validation owners:
  `src/phospy/science/sites/validation.py`,
  `src/phospy/science/sites/metadata_validation.py`,
  `src/phospy/science/sites/identity_contracts.py`, and
  `src/phospy/science/sites/sequence_context.py`.
  `src/phospy/validation/identity_contracts.py` and selected
  `src/phospy/validation/datasets/site_metadata.py` names are compatibility
  routes and must preserve object identity with those science-owned
  implementations.
- Differential design/contrast validation ownership:
  `src/phospy/validation/workflows/differential.py` plus
  `src/phospy/workflows/differential/validator.py`.
- Technical-replicate policy owner:
  `src/phospy/workflows/differential/replicates.py`.
- Intensity-scale establishment and dataset construction coherence owners:
  `src/phospy/validation/transformations/state.py`,
  `src/phospy/science/datasets/construction/analysis_ready.py`,
  `src/phospy/science/datasets/construction/validation.py`,
  `src/phospy/science/datasets/construction/trusted_assertions.py`, and
  `src/phospy/science/datasets/construction/fingerprints.py`.
- Reference manifest parsing, bundled resource integrity, redistribution
  release-gate policy, and bundle semantics owners:
  `src/phospy/science/references/validation/manifest_schema.py`,
  `src/phospy/science/references/validation/resource_integrity.py`,
  `src/phospy/science/references/validation/redistribution_policy.py`, and
  `src/phospy/science/references/validation/bundle_semantics.py`. The
  `src/phospy/science/references/validation/__init__.py` route is a
  compatibility facade only.
- Reference bundle and Kinase Library resource construction invariants:
  passive models live in `src/phospy/science/references/models.py` and
  `src/phospy/science/references/kinase_library_models.py`; validator service
  routes live in `src/phospy/science/references/validation/bundle.py` and
  `src/phospy/science/references/validation/kinase_library.py`, with
  validation-package routes re-exporting these validators for internal
  compatibility.
- Workflow-boundary owners include
  `src/phospy/workflows/kinase/validator.py` and
  `src/phospy/workflows/signalome/validator.py`.
- Public API boundary checks live in
  `tests/unit/api/test_validation_not_public_api.py` and
  `tests/architecture/test_validation_boundaries.py`.
- Package dependency checks live in
  `tests/architecture/test_package_dependency_dag.py`.
- Import-graph records must retain AST line numbers for diagnostics; package
  dependency failures should not collapse all imports to line 1.
- Duplicate-definition checks for validation/scientific ownership use concrete
  symbols and expected owner files, not package-name-only ownership claims.
- Compatibility import routes are allowed only when they preserve object
  identity with the owner symbol.
- Config ownership and duplicate-definition checks live in
  `tests/architecture/test_config_ownership.py`.
- Adapter wiring checks live in `tests/unit/test_protocol_adapter_wiring.py`.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R.,
& Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356. https://doi.org/10.1093/bioinformatics/btz306

YangLab. (n.d.). *PhosR* (Version release) [Computer software]. GitHub.
https://github.com/PYangLab/PhosR
