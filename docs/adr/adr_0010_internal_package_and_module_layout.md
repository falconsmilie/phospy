# ADR: Internal Package and Module Layout for PhosPy

## Document Control

- **ADR ID:** ADR-0010
- **Title:** Internal Package and Module Layout for PhosPy
- **Status:** Accepted
- **Date:** 2026-05-02
- **Decision Type:** Architecture Decision Record

## Abstract

This ADR defines module-splitting governance and compatibility-shim governance
for internal PhosPy packages. The architecture uses split component modules and
keeps compatibility layers limited to explicit serialization/bundle concerns.

## Status

Accepted.

This ADR supersedes earlier high-level package-layout guidance.

## Context and Problem Statement

Internal workflow code has been split into focused modules (validators,
interpreters, executors, provenance builders, result assembly, exporters) while
keeping public workflow classes stable. Bundle-format replay still requires
explicit compatibility adapters in selected IO paths.

Without governance:

- file-size-driven splits create incoherent module boundaries
- compatibility adapters can accumulate business logic
- compatibility names can leak into accidental public API

## Decision Drivers

1. Keep internal modules aligned to responsibility, not file length.
2. Preserve workflow public contracts while allowing internal movement.
3. Prevent compatibility modules from becoming new logic owners.

## Decision

### Module Splitting Rules

1. Modules are split by responsibility, not arbitrary file length.
2. Workflow public classes remain orchestration shells.
3. Validators, interpreters, executors, provenance builders, result models, and
   exporters remain separate responsibilities.
4. New behavior is implemented in owning modules.

### Compatibility-Shim Governance

1. Workflow packages must not expose compatibility re-export modules.
2. Compatibility shims must not own logic.
3. Compatibility shims must not become public API unless exported through
   `phospy.api`.
4. Compatibility shims require a clear reason to exist.
5. New code must not add behaviour to compatibility modules.

### Policy-Enum Ownership

1. Policy enums must be defined in the domain that owns the behavior.
2. Shared policy-enum parsing infrastructure lives under `phospy.policies`.
3. Root-level dumping-ground modules (for example `phospy.policy_models`) are
   not allowed.

### Scientific-Policy Record Ownership

1. Shared scientific policy record models belong in
   `phospy.provenance.scientific_policy_models`.
2. Domain-specific scientific policy builders/records must live in the owning
   domain module.
3. Root-level scientific policy dumping-ground modules (for example
   `phospy.scientific_policies`) are not allowed.
4. Compatibility re-export shims for removed root scientific policy modules are
   not allowed.

### Acceptable Compatibility Shim Examples

- bundle-format compatibility adapters where the compatibility concern is
  explicit (for example payload normalization/parsing modules under
  `io/bundles/.../compatibility/`)

### Unacceptable Compatibility Shim Examples

- modules that combine validation, execution, and provenance responsibilities
- generic `helpers.py`, `utils.py`, or `compatibility.py` modules containing
  new domain logic
- public aliases added only to avoid updating internal imports
- `workflows/*/components.py` compatibility re-export modules

## Consequences

### Positive Consequences

- Responsibility boundaries stay inspectable.
- Internal refactors stay safer without public-surface drift.
- Compatibility concerns stay explicit and limited.

### Negative Consequences

- Teams must update internal imports instead of parking logic in shims.
- Shim modules need periodic review and removal when no longer justified.

## Affected Modules

- `src/phospy/workflows/kinase/validator.py`
- `src/phospy/workflows/kinase/interpreter.py`
- `src/phospy/workflows/kinase/executor.py`
- `src/phospy/workflows/kinase/provenance.py`
- `src/phospy/workflows/signalome/validator.py`
- `src/phospy/workflows/signalome/interpreter.py`
- `src/phospy/workflows/signalome/executor.py`
- `src/phospy/workflows/signalome/provenance.py`
- `src/phospy/io/bundles/_signalome/compatibility/`
- `src/phospy/validation/references/compatibility.py`

## Scope Boundaries

This ADR governs internal module boundaries and compatibility shim usage. It
does not change public API export policy (ADR-0001), validation ownership
(ADR-0007), or test/performance governance (ADR-0014).

## Validation and Review Criteria

Future changes must satisfy all of the following:

1. Is the split responsibility-driven?
2. Does the public workflow class stay orchestration-only?
3. Is a compatibility shim thin, internal, and justified?
4. Is any new behavior added only in owning modules?
5. Is any shim exposed publicly without explicit `phospy.api` export? If yes,
   reject.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi,
R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub.
https://github.com/PYangLab/PhosR
