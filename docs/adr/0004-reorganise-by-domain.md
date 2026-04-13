# ADR 0004: Reorganise `src/phospy/` by Domain Capability and Process

- **Status:** Proposed
- **Date:** 2026-04-10

## Context

The root application package, `src/phospy/`, has grown in a largely flat and historically shaped way. This made early
iteration easier, but it is now creating the same class of problem that previously existed in the validation layer
before that package was reorganised by validation type.

At present, the root package mixes multiple kinds of concerns at the same level, including:

- dataset models and loading
- preprocessing logic
- prediction execution
- activity analysis
- signalome analysis
- reference resolution
- public workflow entry points
- shared I/O concerns

This shape creates several issues.

First, discoverability is weakening. A contributor cannot quickly answer basic questions such as where preprocessing
ends and workflow orchestration begins, or whether a given module represents stable product behaviour or only wiring.

Second, orchestration and domain logic are too close together. This increases the likelihood that top-level convenience
flows accumulate their own processing behaviour instead of delegating to one stable implementation path. The recent
concern around `SimpleKinaseWorkflow` is an example of this structural risk.

Third, mixed-responsibility root modules are becoming difficult to review. Files such as `workflow.py`, `dataset.py`,
`activities.py`, and related modules carry both domain behaviour and public composition concerns, which increases
coupling and makes change boundaries unclear.

Fourth, the current layout does not reflect the actual product shape of PhosPy. PhosPy is not a generic collection of
helpers and services. It is a phosphoproteomics application with clear capability areas, including data ingestion,
preprocessing, prediction, downstream activity analysis, signalome analysis, and biological reference handling.

A clearer package structure is needed so that:

- top-level packages reflect stable domain capabilities
- processes inside each domain have obvious homes
- orchestration becomes thinner and easier to reason about
- future growth does not recreate root-level dumping-ground modules

## Decision

`src/phospy/` will be reorganised by **domain capability first**, and then by **process within that domain**.

The target structure will introduce explicit top-level packages for the main product capabilities, such as:

    src/phospy/
        __init__.py

        api/
        datasets/
        preprocessing/
        prediction/
        activities/
        signalomes/
        references/
        io/
        validation/
        errors/
        internal/

The responsibility of each package will be as follows.

### `api/`

This package will contain the supported public entry points of the application. Its role is thin orchestration and
public API stability, not ownership of core domain behaviour.

Examples include:

- stable user-facing workflows
- public façades over lower-level processing
- compatibility wrappers where needed during migration

Business logic should not accumulate here.

### `datasets/`

This package will contain dataset-shaped models and dataset construction concerns.

Examples include:

- `PhosphoDataset`
- `AnalysisReadyPhosphoDataset`
- dataset loading and building
- dataset result containers

This package will not own preprocessing strategy.

### `preprocessing/`

This package will contain the transformation of raw or semi-structured phosphoproteomic inputs into analysis-ready
forms.

Examples include:

- full-mode preprocessing
- phospho-only preprocessing
- filtering
- correction
- normalisation
- imputation
- site matrix construction
- stable preprocessing result objects

This package will own preprocessing modes directly rather than allowing workflow façades to assemble preprocessing ad
hoc.

### `prediction/`

This package will contain prediction-related logic and execution.

Examples include:

- prediction engines
- sampling behaviour
- scoring
- prediction results
- predMat-related execution paths

This package will not own unrelated workflow or reference orchestration.

### `activities/`

This package will contain kinase activity analysis logic built on top of prediction outputs and compatible
phospho-derived inputs.

Examples include:

- activity scoring
- enrichment-style calculations
- result models

### `signalomes/`

This package will contain signalome construction and related downstream analysis.

Examples include:

- signalome analysis
- map generation
- network-style outputs
- result models

### `references/`

This package will contain biological reference handling as its own explicit domain.

Examples include:

- bundled reference assets
- species/reference resolution
- substrate maps
- motif resources
- site sequence resources

This package exists to prevent reference concerns from being scattered across prediction, workflow, and preprocessing
modules.

### `io/`

This package will contain shared input and output concerns that are genuinely about structured data access rather than a
specific scientific domain.

Examples include:

- table reading
- mapping file loading
- path-oriented I/O utilities

### `validation/`

This package will continue to hold validation concerns, organised by validation type as established in the validation
refactor.

### `errors/`

This package will contain the application’s error hierarchy, grouped in a way that supports clearer boundaries between
validation, domain, and other application errors.

### `internal/`

This package will contain truly internal support code that does not belong to a domain package and is not part of the
supported public surface.

This package must remain narrow. It must not become a replacement dumping ground for miscellaneous logic.

## Decision Rules

The following rules will govern the reorganisation.

1. **Top-level packages must represent domain capabilities, not implementation styles.**
   Packages such as `services/`, `helpers/`, `common/`, or `utils/` will not be introduced as general-purpose buckets.

2. **Processes should live inside the relevant domain package.**
   For example, filtering or imputation related to phosphoproteomic preparation belongs under `preprocessing/`, not in a
   flat root module.

3. **Public orchestration must remain thin.**
   The `api/` package may coordinate domain operations, but it must not become a second implementation layer.

4. **Reference handling is a first-class domain.**
   Reference bundle and species-resolution logic will not remain scattered across unrelated modules.

5. **Mixed-responsibility files must be split while moving.**
   This change is not a file relocation exercise. Where a root module currently mixes multiple concerns, those concerns
   will be separated before or during the move.

6. **Stable public imports may be preserved temporarily where useful.**
   Compatibility shims are allowed during migration when they reduce churn, but they must not permanently hide poor
   internal structure.

## Rationale

This decision is intended to make the package structure reflect the actual product architecture of PhosPy.

A domain-led package layout improves discoverability because contributors can navigate by capability rather than by
historical file names. It also reduces drift by making it harder for orchestration modules to absorb domain behaviour.
This is particularly important for a PhosR port, where the Python implementation should become clearer and more explicit
than the original workflow sprawl, not inherit or reproduce it.

Separating `references/` as its own package is especially important. Reference resolution, bundled assets, and substrate
or sequence resources are used across multiple parts of the application, but they are not merely implementation details
of prediction. Treating them as a distinct domain makes that relationship explicit and reduces coupling.

Creating an `api/` package also makes the supported user-facing surface easier to define. This helps distinguish between
stable entry points and internal composition code.

## Alternatives Considered

### Keep the Current Flat Root Structure

This was rejected because the current layout is already showing the same symptoms previously seen in the validation
layer: broad modules, unclear ownership, and growing review complexity.

### Reorganise by Technical Layer

An alternative would be to structure the application using packages such as `models/`, `services/`, `managers/`, or
`utils/`.

This was rejected because those package names do not describe PhosPy’s product capabilities. In practice, they tend to
become dumping grounds and make it harder, not easier, to understand where functionality belongs.

### Reorganise by Workflows Only

Another alternative would be to keep the root structure mostly flat and group modules around specific workflows.

This was rejected because workflows are not the only stable shape of the system. Core domain capabilities such as
preprocessing, prediction, activities, and reference handling are broader and more stable than any single workflow path.
Organising by workflows would also risk repeating the same over-coupling problem in a different form.

## Consequences

### Positive Consequences

- The package structure will better reflect the actual capabilities of the application.
- Contributors will be able to navigate by domain and process more quickly.
- Orchestration code will have clearer boundaries.
- Reference handling will become easier to reason about and extend.
- Future growth will be less likely to recreate root-level dumping grounds.
- The public API surface will be easier to define and protect.

### Negative Consequences

- The refactor will touch a large number of imports.
- Some files will need to be split rather than simply moved, which increases the size of the change.
- Temporary compatibility shims may be needed during migration.
- The reorganisation may create short-term churn in tests and documentation.

## Scope

This ADR covers the structural reorganisation of the root application package, `src/phospy/`, around domain capability
and process.

This ADR does not define the exact implementation of every moved module. It sets the target architecture and package
responsibilities that subsequent tickets must follow.

## Out of Scope

The following are out of scope for this ADR:

- redesigning the scientific behaviour of preprocessing, prediction, or downstream analyses
- changing supported public behaviour unless required by a specific follow-up ticket
- eliminating all compatibility imports in one step
- introducing new public workflows purely as part of the package move

## Implementation Notes

Implementation should proceed incrementally rather than as one uncontrolled rename.

A sensible order is:

1. define the target package responsibilities
2. move the worst mixed-responsibility root modules first
3. push core behaviour down into domain packages
4. preserve stable public import paths where they meaningfully reduce churn
5. remove temporary shims once the new structure is stable

Particular care should be taken with current orchestration-heavy modules so that the refactor improves boundaries rather
than simply moving existing coupling into new folders.

## Outcome Sought

The intended outcome is that `src/phospy/` becomes readable as a product architecture:

- here is how data enters
- here is how data is prepared
- here is how prediction works
- here is how activity analysis works
- here is how signalome analysis works
- here is how references are resolved
- here is how validation works
- here is the thin public API

That is the structural direction PhosPy should follow going forward.