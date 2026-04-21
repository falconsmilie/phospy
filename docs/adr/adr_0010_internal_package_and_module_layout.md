# ADR: Internal Package and Module Layout for PhosPy

## Document Control

- **ADR ID:** ADR-010
- **Title:** Internal Package and Module Layout for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines how the PhosPy codebase should be physically organised into packages and modules. The package is being developed as a maintainable Python port of PhosR. To support that goal, the internal layout must reflect the architectural boundaries already established for the public API, workflows, dataset boundary, references, result models, transformation handling, validation, builders, and exceptions.

The decision is to organise the codebase around clear internal domains rather than around ad hoc technical helpers or historical implementation seams. The physical module layout should reinforce the product shape of one dataset model and two workflows, while keeping shared concerns such as validation, transformations, references, builders, and errors in dedicated internal domains.

## Status

Accepted.

This ADR defines the package and module layout that supports the earlier ADR set and prevents structural drift during the rewrite.

## Context and Problem Statement

The previous ADRs define a clear architectural direction:

- one public dataset model
- two public workflows
- workflow stages built around validator, interpreter, and executor
- strict dataset boundary
- explicit reference resolution
- disciplined result models
- typed transformation state
- private validation domain
- flexible builder path below the strict dataset boundary
- explicit exception taxonomy

However, architecture decisions only remain stable if the physical code layout supports them. A poor module structure causes otherwise good design to degrade over time. In practice, this leads to:

- responsibilities bleeding across modules
- helper-heavy packages with unclear ownership
- duplicated logic placed where it happens to fit rather than where it belongs
- public and internal concepts being mixed together
- difficulty locating related code
- pressure to create wrapper layers just to cross awkward package seams

PhosPy therefore needs an explicit decision on internal package and module layout.

## Decision Drivers

The decision is driven by the following considerations:

1. **Architectural alignment.** The code layout should reflect the earlier ADR decisions rather than undermine them.
2. **Maintainability.** Related code should live together in obvious domains.
3. **Discoverability.** Developers should be able to find dataset, workflow, validation, reference, and builder logic easily.
4. **Boundary clarity.** Public API modules and internal implementation modules should be visibly distinct.
5. **Extensibility.** New supported input paths, organisms, or workflow internals should fit naturally into the layout.
6. **Simplicity.** The structure should stay boring and readable rather than becoming a maze of technical helper packages.

## Proposed Decision

PhosPy will use a domain-oriented internal package layout with a small public API surface and clearly separated internal implementation domains.

The top-level internal layout should reflect the main architectural responsibilities:

- public API surface
- datasets and builders
- workflows
- references
- transformations
- validation
- errors
- domain services for prediction, activity, and signalome logic

The physical module structure should reinforce the rule that the package is centred on one analysis-ready dataset and two primary workflows.

## Core Design Principle

The internal package layout should follow **domain and boundary ownership**, not **incidental implementation convenience**.

A module should exist because it owns a coherent responsibility, not because code happened to be added there first.

## Proposed Top-Level Package Direction

A likely healthy top-level structure is conceptually similar to:

```text
phospy/
  api/
  datasets/
  workflows/
  references/
  transformations/
  validation/
  errors/
  prediction/
  activities/
  signalomes/
  io/
```

The exact names may vary slightly, but this is the intended shape.

## Public API Package

The `api/` package should define the canonical public API ownership surface.

It should contain the public models and entry points that are defined and
organised as package-owned API contracts.

Normal user-facing imports should stay anchored at top-level `phospy`, which
acts as a minimal convenience entrypoint for the main product objects only.

A likely shape is:

```text
api/
  datasets.py
  requests.py
  results.py
  workflows.py
  enums.py
```

This package should stay intentionally small.

It should not become a mirror of the entire internal package structure.

User-handleable exceptions remain implemented in `errors/`, but supported
public exception imports should be centralised under `phospy.api`.

## Datasets Package

The `datasets/` package should own the analysis-ready dataset boundary and the builder path beneath it.

A likely shape is:

```text
datasets/
  models.py
  builders/
  services/
```

### Responsibilities

- `models.py` owns `AnalysisReadyPhosphoDataset`
- `builders/` owns the public builder route and its specialised collaborators
- `services/` may hold dataset-shaping or preprocessing helpers that belong specifically to the dataset boundary

Builder collaborators should remain under the dataset builder area rather than being spread across the codebase. This keeps the builder story consistent and prevents ingestion logic from fragmenting into multiple unrelated packages.

The dataset package should not own workflow execution logic.

## Workflows Package

The `workflows/` package should own the two public workflow families and their internal stage components.

A likely shape is:

```text
workflows/
  kinase/
    contracts.py
    validator.py
    interpreter.py
    executor.py
  signalome/
    contracts.py
    validator.py
    interpreter.py
    executor.py
```

### Responsibilities

- workflow request/result contract support that is internal to workflow staging
- validator, interpreter, and executor components
- workflow-specific orchestration logic only

`contracts.py` should primarily hold internal stage DTOs and related data contracts. Workflow-specific protocol definitions should live there only if they are tightly scoped to that workflow domain and improve locality. In general, behavioural protocols should not be mixed in casually with DTO contracts. If protocol definitions start to accumulate or serve a broader role, they should move to a more explicit internal module such as `interfaces.py` within the same workflow domain.

The workflow package should not absorb shared validation, reference loading, or transformation establishment logic that belongs elsewhere.

## References Package

The `references/` package should own reference models, bundled-resource access, and provider logic.

A likely shape is:

```text
references/
  models.py
  providers/
  resources.py
```

### Responsibilities

- `ReferenceBundle`
- bundled reference resource access
- provider implementations for built-in organisms
- future provider-oriented extension points for additional organisms or external APIs

Reference resolution should remain separate from workflow execution and should not be hidden in generic helpers.

## Transformations Package

The `transformations/` package should own transformation-state types and transformer logic.

A likely shape is:

```text
transformations/
  models.py
  transformers/
```

### Responsibilities

- `TransformationState`
- transformer interfaces and implementations
- transformation-specific internal helpers

This keeps transformation concerns distinct from dataset models and workflows while still serving the dataset-building boundary.

## Validation Package

The `validation/` package should exist as its own internal domain, in line with ADR-007.

A likely shape is:

```text
validation/
  common/
  datasets/
  references/
  transformations/
  workflows/
```

### Responsibilities

- reusable lower-level validation components
- dataset validation building blocks
- reference validation building blocks
- transformation validation building blocks
- workflow-level validator composition support

Validation components should not depend on each other across subdomains. Composition should happen at higher levels.

## Errors Package

The `errors/` package should own the explicit exception taxonomy defined by ADR-009.

A likely shape is:

```text
errors/
  base.py
  validation.py
  input.py
  build.py
  references.py
  transformations.py
  workflows.py
```

### Responsibilities

- package-level base exception
- domain-specific base exceptions
- no unrelated helper logic

This keeps failure semantics centralised and avoids exception definitions being scattered across unrelated modules.

## Prediction Package

The `prediction/` package should own kinase scoring and prediction domain logic.

A likely shape is:

```text
prediction/
  models.py
  scoring/
  prediction/
```

### Responsibilities

- scoring-stage domain logic
- prediction-stage domain logic
- prediction-related internal DTOs and models

This package should not own public workflow orchestration.

## Activities Package

The `activities/` package should own kinase activity analysis logic.

A likely shape is:

```text
activities/
  models.py
  services.py
```

This package remains smaller than the prediction domain, but it should still remain clearly separate rather than being folded into generic workflow helpers.

If `activities/` grows meaningfully, it should evolve in a way that stays consistent with other scientific domain packages rather than remaining a catch-all `services.py` pattern indefinitely.

## Signalomes Package

The `signalomes/` package should own downstream signalome-specific analysis logic and domain models.

A likely shape is:

```text
signalomes/
  models.py
  services.py
```

This package should not become a second workflow package. Workflow orchestration remains in `workflows/signalome/`, while the signalome domain logic itself stays here.

If `signalomes/` grows meaningfully, it should evolve consistently with the `activities/` package and the broader domain-layout rules rather than diverging into an unrelated structure.

## IO Package

The `io/` package should own file and output concerns that are not part of the core workflow/result contract.

A likely shape is:

```text
io/
  readers/
  publishing.py
```

Reader logic should remain under `io/readers/` internally. That is the cleanest physical home for file-reading concerns.

However, file-path builder entry points should hide those reader components from the public usage story. Users should experience file-path handling through the builder route, not through direct exposure to internal reader modules.

This package should remain secondary to the main scientific domains.

It should not become the place where the builder or workflows quietly re-home their responsibilities.

## Public vs Internal Boundary Direction

Public product-facing code should be shallow and easy to understand.

Internal implementation packages may be richer, but their structure should not leak upward into the public import story.

In practice:

- `api/` stays small and curated as the canonical API-definition namespace
- top-level `phospy` stays a small convenience entrypoint
- internal packages stay domain-oriented
- callers should not need to understand the internal package layout to use the package correctly

## Module Granularity Direction

Granularity should favour clarity and ease of composition.

The codebase should avoid both extremes:

### Too coarse

- giant modules that mix many responsibilities
- workflow files that absorb half the package
- helper dumps with unclear ownership

### Too fragmented

- one-file-per-tiny-concept when the split adds no clarity
- deep nesting that forces constant cross-import churn
- micro-packages created only for theoretical neatness

The rule should be:

- split when ownership becomes clearer
- do not split merely to appear clean

## Naming Direction

Names should stay concise and honest.

This aligns with the project preference for avoiding overly long type and module names.

Examples of the intended direction:

- `models.py`
- `services.py`
- `validator.py`
- `interpreter.py`
- `executor.py`
- `providers/`
- `transformers/`

Avoid long compound names where a shorter honest name works.

## Dependency Direction

The package layout should reinforce sane dependency flow.

A healthy direction is:

- `api/` depends on stable public-facing models and entry points
- `workflows/` depends on datasets, references, transformations, validation, prediction, activities, and signalomes
- `datasets/builders/` depends on validation, transformations, references where needed, and dataset models
- validation components do not depend on each other across subdomains
- domain packages such as `prediction/`, `activities/`, and `signalomes/` do not depend on workflow packages

This helps prevent cycles and keeps orchestration above domain logic rather than buried inside it.

## Consequences

### Positive consequences

- The code layout starts to match the architecture rather than fight it.
- Responsibilities become easier to locate and maintain.
- Public and internal concepts remain more clearly separated.
- New work has obvious homes instead of creating structural drift.
- The rewrite gains a stable physical structure to target.

### Negative consequences

- Existing code may need meaningful movement across modules and packages.
- Some current implementation seams may disappear entirely.
- The project must stay disciplined so the layout does not re-accumulate helper-heavy clutter.

### Neutral consequences

- Exact filenames and some subpackage details may still evolve.
- Internal modules may still use additional DTOs or helper components where justified, as long as they respect the broader layout rules.

## Rejected Alternatives

### Alternative 1: Keep the current structure and gradually tidy it up opportunistically

This option was rejected because a drifting structure tends to preserve past mistakes and makes the rewrite less coherent.

### Alternative 2: Organise everything around technical layers only

This option was rejected because the package is shaped around domain workflows and scientific boundaries, not around abstract technical layering alone.

### Alternative 3: Collapse everything into a few large modules for speed

This option was rejected because it would undermine the clear domain boundaries established by the earlier ADRs.

### Alternative 4: Over-fragment the code into many tiny packages immediately

This option was rejected because it would create navigation overhead and complexity without necessarily improving ownership.

## Resolved Decisions

The following decisions are now resolved for this ADR.

1. Builder collaborators should remain under the dataset-builder area rather than being spread across the codebase.
2. User-handleable exceptions should be publicly centralised through `phospy.api`.
3. Workflow `contracts.py` modules should primarily hold internal stage DTOs; tightly scoped workflow-local protocols may live there when justified, but broader behavioural protocols should move to a more explicit internal module.
4. Reader logic should remain under `io/readers/` internally, while builder file-path entry points hide that structure from the public usage story.
5. `activities/` and `signalomes/` should evolve consistently if they grow.

## Implementation Guidance

A likely healthy first-pass rewrite target is:

- stabilise `api/`
- create clear `datasets/`, `workflows/`, `validation/`, `references/`, `transformations/`, and `errors/` domains
- keep `prediction/`, `activities/`, and `signalomes/` as the primary scientific domain packages
- let builders sit below the dataset boundary rather than leaking into workflows
- keep module names concise and boring

Reviewers should reject changes that place code according to convenience alone when a clear domain owner already exists.

## Scope Boundaries

This ADR defines internal package and module layout only.

It does not define:

- the full rewrite plan or migration sequence
- the exact final public API contents beyond earlier ADRs
- the full file-IO strategy
- packaging and release mechanics
- deployment structure

Those concerns should be addressed separately.

## Validation and Review Criteria

Future code and review work should check proposed changes against the following questions:

1. Does this module have a clear domain owner?
2. Does this placement reinforce or weaken the earlier ADR boundaries?
3. Does this keep public and internal concerns clearly separated?
4. Does this reduce or increase helper sprawl?
5. Does this make the package easier or harder to navigate?

If the answers are weak or negative, the design should be reconsidered.

## Relationship to Earlier ADRs

This ADR complements the earlier architecture decisions.

- ADR-001 defines the intended public API contract.
- ADR-002 defines the internal workflow architecture.
- ADR-003 defines the dataset and preprocessing boundary.
- ADR-004 defines the reference resolution strategy and `ReferenceBundle` contract.
- ADR-005 defines result-model design.
- ADR-006 defines the transformation-state and transformer contract.
- ADR-007 defines the validation-domain architecture.
- ADR-008 defines the analysis-ready dataset builder architecture.
- ADR-009 defines the exception and error taxonomy.
- ADR-010 defines the internal package and module layout that supports all of the above.

Together, these ADRs establish:

- one public dataset model
- two public workflows
- one strict dataset boundary
- one flexible builder path below that boundary
- one private validation domain
- one stronger transformation-state contract
- one explicit exception taxonomy
- one domain-oriented physical code layout

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR*. GitHub repository. https://github.com/PYangLab/PhosR
