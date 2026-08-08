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
keeps compatibility layers limited to explicit serialization/bundle concerns
and documented identity-preserving import routes.

## Status

Accepted.

This ADR supersedes earlier high-level package-layout guidance.

Update note (2026-07-17, actual package DAG): The enforced package graph is now
the complete top-level `phospy` graph, not a curated subset. Architecture tests
parse every source file and compare actual package edges with an explicit
allowed-edge table.

Update note (2026-08-05, scientific model decomposition): oversized
transformation, activity, reference-manifest, and batch-correction model/helper
modules are split by change reason. Historical aggregate routes remain only as
identity-preserving imports; owner modules define the single class objects.

Update note (2026-08-08, contracts facade dependency rule):
`phospy.contracts` is a transport/public facade. It may re-export exact
science-owned objects from designated public science modules, but it must not
import private science modules, executors, construction services, internal
views, or validation implementation modules.

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
4. Keep the package dependency graph acyclic and enforceable in CI.

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

### Configuration Ownership

1. Each configuration dataclass, enum, or validation helper has exactly one
   implementation owner.
2. Algorithm and scientific policy types used directly by numerical code live
   in the owning science domain. Higher layers may re-export those exact
   objects, but must not copy the source or wrap them in subclasses.
3. Public transport/request configuration lives under `phospy.contracts.configs`
   when it describes API shape, optional user input, or workflow composition.
   `phospy.science` must not import those public DTOs.
4. Interpreters translate public config DTOs into distinctly named resolved
   execution models where defaults, references, seeds, eligible sets, or
   contextual choices are interpreted.
5. Shared scalar validation helpers live in `phospy.policies` or another lower
   shared layer. Generic `common.py` dumping grounds for unrelated policies are
   not acceptable.
6. Re-export shims are allowed only when they preserve object identity and
   respect the package DAG.

### Scientific-Policy Record Ownership

1. Shared scientific policy record models belong in
   `phospy.provenance.scientific_policy_models`.
2. Domain-specific scientific policy builders/records must live in the owning
   domain module.
3. Root-level scientific policy dumping-ground modules (for example
   `phospy.scientific_policies`) are not allowed.
4. Compatibility re-export shims for removed root scientific policy modules are
   not allowed.

### Package-Ownership Boundaries

1. `phospy.science` owns scientific domain code:
   domain models, phosphosite/sequence logic, references/evidence handling,
   prediction/activity/differential/signalome computation, dataset-domain
   preprocessing, and transformation-domain models.
2. `phospy.workflows` owns orchestration only.
3. `phospy.validation` owns validation infrastructure only.
4. `phospy.provenance` owns shared provenance infrastructure only.
5. `phospy.tables` owns identity-preserving compatibility imports only. It
   must not own scientific table logic or define duplicate table/schema
   classes.
6. `phospy.frames` owns pandas frame ownership helpers and generic table-schema
   infrastructure.
7. `phospy.data` owns packaged static resources only.
8. `phospy.contracts` owns stable transport facade modules: passive request
   DTOs, public result containers, and identity-preserving re-export routes for
   science-owned domain values.
9. `phospy.api` and `phospy.advanced` own aggregate public namespaces only.

### Enforceable Package Dependency DAG

Package imports must remain acyclic. The architecture tests parse all AST import
statements, including local imports, relative imports, `TYPE_CHECKING` imports,
and static `importlib.import_module(...)` / `__import__(...)` calls, so cycles
cannot be hidden by deferring imports.

The enforced package rules are:

1. `phospy.errors` is a leaf. It must not import any other `phospy` package.
2. `phospy.contracts` must not import `phospy.validation` or
   `phospy.workflows`.
3. `phospy.science` must not import `phospy.contracts`, `phospy.io`,
   `phospy.tables`, `phospy.validation`, or `phospy.workflows`.
4. `phospy.contracts` may import `phospy.science` only through designated
   public domain modules: public science config policy modules; science-owned
   model/result/table-schema modules; and narrow reference, evidence,
   transformation, and result-caveat domain modules listed in the architecture
   test. It must not import private science modules, executors, construction
   services, internal views, or validation implementation modules.
5. Concrete local readers, reference source loaders, and nested workflow runners
   are injected by API/workflow orchestration adapters.
6. Compatibility modules may re-export moved names, but they must not reintroduce
   reverse imports or own new behavior.
7. API adapters may import private validation modules for composition, but those
   validators must not be exported from public namespaces.

The current top-level orientation is:

- `api -> contracts, errors, io, science, validation, workflows`
- `contracts -> errors, frames, policies, provenance, science`
- `io -> contracts, errors, provenance, science, validation`
- `science -> errors, frames, policies, provenance`
- `tables -> frames, science`
- `validation -> contracts, errors, frames, provenance, science`
- `workflows -> contracts, errors, provenance, science, validation`

The CI graph check lives in
`tests/architecture/test_package_dependency_dag.py`. Config ownership checks
also live in `tests/architecture/test_config_ownership.py`. New exemptions
require an ADR update and a named test change; unexplained exemptions are not
allowed.

### Acceptable Compatibility Shim Examples

- bundle-format compatibility adapters where the compatibility concern is
  explicit (for example payload normalization/parsing modules under
  `io/bundles/.../compatibility/`)
- `phospy.tables.*` identity-preserving re-export modules that keep supported
  legacy scientific table import routes available while owned
  implementations live under `phospy.science.tables` and generic schema
  infrastructure lives under `phospy.frames`

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
- `src/phospy/io/bundles/_signalome/{config,diagnostics,tables}.py`
- `src/phospy/validation/references/compatibility.py`
- `src/phospy/science/references/models.py`
- `src/phospy/science/references/manifest.py`
- `src/phospy/science/references/manifest_common.py`
- `src/phospy/science/references/manifest_files.py`
- `src/phospy/science/references/manifest_model.py`
- `src/phospy/science/references/manifest_policy.py`
- `src/phospy/science/references/redistribution.py`
- `src/phospy/science/transformations/policy.py`
- `src/phospy/science/transformations/provenance.py`
- `src/phospy/science/transformations/scale_state.py`
- `src/phospy/science/transformations/scale_values.py`
- `src/phospy/science/activities/method_models.py`
- `src/phospy/science/activities/inputs.py`
- `src/phospy/science/activities/diagnostics.py`
- `src/phospy/science/activities/results.py`
- `src/phospy/science/activities/result_validation.py`
- `src/phospy/science/datasets/preprocessing/batch_correction.py`
- `src/phospy/science/datasets/preprocessing/batch_correction_models.py`
- `src/phospy/science/datasets/preprocessing/batch_correction_engine.py`
- `src/phospy/science/datasets/preprocessing/batch_correction_provenance.py`
- `src/phospy/science/datasets/preprocessing/batch_correction_provenance_payloads.py`
- `src/phospy/science/datasets/preprocessing/batch_correction_provenance_validation.py`
- `src/phospy/science/references/kinase_library.py`
- `src/phospy/io/readers/dataset_inputs.py`
- `src/phospy/workflows/batch_correction/preprocessing_adapter.py`
- `tests/architecture/test_package_dependency_dag.py`

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
6. Does the AST package dependency graph remain acyclic?
7. Does `science` avoid concrete I/O and workflow imports by using injected
   protocols?

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi,
R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub.
https://github.com/PYangLab/PhosR
