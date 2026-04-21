# ADR: Rewrite Roadmap and Fresh-Start Plan for PhosPy

## Document Control

- **ADR ID:** ADR-012
- **Title:** Rewrite Roadmap and Fresh-Start Plan for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines how the PhosPy rewrite should proceed. The package is being developed as a maintainable Python port of PhosR. To support that goal, the rewrite must follow the architecture established by the ADR set rather than attempt to preserve the structure of the current application.

The decision is to treat the rewrite as a fresh start. The existing application is a source of scientific logic, behavioural examples, and useful implementation ideas, but it is not a structural baseline to migrate forward. The roadmap should focus on building the new architecture directly, validating scientific behaviour as it is reintroduced, and resisting compatibility-driven carryover of the current design.

## Status

Accepted.

This ADR defines the rewrite execution strategy that follows the architectural decisions already made for the public API, workflows, dataset boundary, references, results, transformation state, validation, builders, exceptions, and package layout.

## Context and Problem Statement

The ADR set now defines a clear target architecture for PhosPy:

- one public analysis-ready dataset model
- two public workflows
- workflow stages built around validator, interpreter, and executor
- strict dataset boundary
- flexible builder path below that boundary
- explicit reference-resolution path
- disciplined result models
- typed transformation state
- private validation domain
- explicit exception taxonomy
- domain-oriented internal package layout

At the same time, the existing application reflects an older structural state that the project no longer intends to preserve. That codebase may still contain useful scientific logic and useful examples, but it also contains architectural patterns the project has now explicitly decided against.

If the rewrite is treated as a migration from the current application, several risks follow:

- historical structure will be preserved by habit rather than by design
- wrappers, aliases, and pass-through layers will survive because they already exist
- old internal seams will distort the new package layout
- implementation effort will drift into adaptation work instead of building the intended architecture directly
- the rewrite will become a refactor swamp rather than a clean restart

PhosPy therefore needs an explicit decision about how to execute the rewrite.

## Decision Drivers

The decision is driven by the following considerations:

1. **Architectural clarity.** The ADR set already defines the intended target architecture.
2. **Freedom to simplify.** The rewrite should not preserve current application structure by default.
3. **Scientific fidelity.** Useful scientific logic may still be reused where correct.
4. **Maintainability.** The new codebase should be organised around the fresh architecture from the start.
5. **Execution discipline.** The team needs a clear rule for what is reused and what is discarded.
6. **Avoiding pseudo-migration.** The project should not drift into carrying old architecture forward under the label of a rewrite.

## Proposed Decision

PhosPy will be rewritten as a fresh-start implementation.

This is not a migration of the existing application structure.

The existing application may be used as:

- a source of scientific logic
- a source of validation rules
- a source of behavioural examples
- a source of test ideas
- a source of implementation techniques where they fit the new design

It must not be treated as:

- the structural template for the rewrite
- a compatibility baseline for internal architecture
- a package layout to preserve
- a workflow/result shape to retain automatically

## Core Design Principle

Reuse **behaviour and scientific logic** where appropriate, but rebuild **structure and boundaries** from the ADR-defined architecture.

The rewrite should ask:

- is this logic scientifically useful?
- does it fit the new architecture cleanly?

It should not ask:

- how do we keep this old class alive?
- how do we adapt this old layer instead of replacing it?

## Fresh-Start Rule

The fresh-start rule is:

- begin from the intended public API and internal architecture
- build the new package layout directly in a separate clean package area
- port logic selectively into the new boundaries
- do not preserve old structure merely because it already exists

This means the rewrite is design-led, not legacy-led.

## Repository and Codebase Stance

The rewrite should begin in a separate clean package area rather than trying to reshape the existing application structure in place.

The old application may remain temporarily as reference material only while useful logic is being reintroduced. However, it should be removed as soon as that reference value is no longer needed.

The goal is not long-term coexistence. The goal is to minimise the period during which the old structure can continue influencing new design decisions.

## What May Be Reused

The following are acceptable sources of reuse from the existing application:

- scientific calculations that are still correct
- domain-specific transformation or scoring logic
- reference-resource handling ideas that fit the new reference architecture
- sequence derivation logic that fits the builder boundary
- validation rules that fit the validation domain
- examples and tests that help confirm intended behaviour

Reused logic should be adapted into the new boundaries rather than wrapped in compatibility layers.

## What Should Not Be Carried Forward Automatically

The following should not be preserved by default:

- existing internal package layout
- wrapper-heavy orchestration
- duplicated result accessors
- compatibility aliases
- helper sprawl
- ad hoc validation placement
- old builder or ingestion entry points that conflict with the new public builder story
- legacy method naming patterns that conflict with the `run(...)` convention

Any carryover of these patterns requires a positive architectural justification, not simple historical inertia.

## Rewrite Outcome Target

The rewrite should converge directly on the architecture already defined by ADRs 001 through 011.

In practical terms, the target is:

- public API first
- dataset boundary below it
- builder path below the dataset boundary
- workflows built around validator, interpreter, and executor
- reference resolution in the interpreter path
- domain-oriented package layout
- private validation domain
- explicit exception taxonomy

The implementation sequence should reinforce this structure from the start.

## Roadmap Strategy

The rewrite should proceed in architecture-aligned stages.

### Stage 1: Establish the package skeleton

Build the new package structure first.

This includes the main domains defined by the ADR set, such as:

- public API
- datasets
- workflows
- references
- transformations
- validation
- errors
- prediction
- activities
- signalomes
- io

The purpose of this stage is to make the physical structure match the intended architecture before deeper logic is introduced.

### Stage 2: Stabilise the public contract

Build the public surface next:

- public enums and errors
- `AnalysisReadyPhosphoDataset`
- request/result DTOs
- builder public contract
- `KinaseWorkflow`
- `SignalomeWorkflow`

The purpose of this stage is to lock the external product story before deeper internals are added.

### Stage 3: Reintroduce builder and dataset flow

Implement the path from messy inputs to analysis-ready dataset:

- public builder class
- `DatasetBuildRequest`
- specialised internal builder collaborators
- sequence derivation support
- transformation-state establishment
- validation composition

This stage should deliver a working dataset boundary before the full workflows are completed.

### Stage 4: Reintroduce kinase workflow

Implement the end-to-end kinase path:

- validator
- interpreter
- executor
- scoring logic
- prediction logic
- activity analysis
- result assembly

Scientific logic from the existing application may be reused here where correct, but only within the new architecture.

### Stage 5: Reintroduce signalome workflow

Implement the downstream signalome path:

- signalome validator
- signalome interpreter
- signalome executor
- signalome domain services
- signalome result assembly

This stage should assume the upstream kinase result contract already exists.

### Stage 6: Harden and simplify

After the main paths work, tighten the implementation:

- remove temporary scaffolding
- remove old reference material that is no longer needed
- collapse any accidental wrapper drift
- improve error messaging
- align tests to the new contracts
- confirm package layout matches ADR-010

This stage exists to keep the rewrite from accumulating transitional clutter.

## Testing Strategy During Rewrite

The rewrite should use tests to confirm scientific and contract behaviour, not to preserve legacy structure.

A healthy testing direction is:

- test the new public contracts directly
- test validation and exception behaviour explicitly
- test builder output against the dataset contract
- test workflow results against the new result model design
- reuse existing tests or examples only when they validate correct behaviour under the new architecture

Tests should not force preservation of old internal seams.

Parity should be considered separately as its own architectural/testing question rather than assumed as part of this roadmap. If parity requirements are needed, they should be defined in a dedicated ADR.

## Completion Rule

The rewrite should be considered complete enough to remove any remaining temporary bridging code or reference scaffolding once the proposed functionality exists.

For PhosPy, that means:

- public dataset creation
- public kinase workflow
- public signalome workflow

Once those are present in the new architecture and working to the expected standard, temporary old-structure scaffolding should be removed rather than preserved out of caution.

## Migration Stance

The project stance is explicitly:

- no migration of the existing application structure
- no compatibility-preservation obligation for old internals
- no requirement to adapt the old package layout into the new one

This is a fresh-start rewrite.

If a piece of old logic is useful, it should be ported into the new design as though it were source material, not infrastructure being carried forward.

## Documentation Strategy During Rewrite

The ADR set should act as the decision baseline.

During implementation:

- architecture docs should describe the new system only
- the docs should not become mixed narratives about old and new structure
- temporary implementation notes may exist, but they should not blur the rewrite direction

The documentation should present PhosPy as the new design, not as a migration story.

## Consequences

### Positive consequences

- The rewrite stays aligned with the architecture instead of being pulled backward by the old structure.
- Scientific logic can still be reused where valuable.
- The team gets a clear rule for what to keep and what to discard.
- Implementation effort is focused on the target design rather than compatibility glue.
- The resulting package is more likely to be coherent and maintainable.

### Negative consequences

- More code may be rewritten or re-homed than in a migration-style approach.
- Some useful old abstractions may be discarded even if they previously worked.
- The team must stay disciplined and avoid slipping into legacy-led decisions.

### Neutral consequences

- Old code may still remain useful as reference material during implementation.
- Some names or isolated implementation ideas may survive if they fit the new architecture cleanly.

## Rejected Alternatives

### Alternative 1: Incremental migration from the current application structure

This option was rejected because it would preserve too much of the existing architecture by inertia and weaken the value of the ADR-driven redesign.

### Alternative 2: Hybrid rewrite that preserves old internals while changing only the public API

This option was rejected because the current internal structure is part of the problem, not a neutral substrate.

### Alternative 3: Big-bang logic copy without architectural staging

This option was rejected because it would risk reproducing old structure accidentally and make implementation harder to reason about.

### Alternative 4: Preserve old tests and structures as the main rewrite guide

This option was rejected because tests and code from the current application should inform behaviour, not dictate architecture.

## Resolved Decisions

The following decisions are now resolved for this ADR.

1. The rewrite should happen in a separate clean package area.
2. The old application should be removed as soon as it is no longer needed as reference material.
3. The preferred sequencing for scientific logic reuse is to establish the package skeleton first.
4. Parity should be addressed in its own ADR rather than being assumed implicitly in the rewrite roadmap.
5. The rewrite should be considered complete enough to remove remaining temporary scaffolding once it delivers public dataset creation, kinase workflow, and signalome workflow functionality.

## Implementation Guidance

A likely healthy execution rule set is:

- build the new package structure first
- implement the new public API and internal boundaries directly
- port old scientific logic selectively into the new boundaries
- delete or ignore structural carryover that does not fit
- keep the rewrite narrative honest: this is a fresh-start implementation

Reviewers should reject work that preserves old structure without a clear ADR-aligned reason.

## Scope Boundaries

This ADR defines the rewrite execution strategy only.

It does not define:

- the full implementation details of each package
- the exact scientific algorithms to port
- release planning or versioning policy
- repository management details beyond the rewrite stance
- deployment sequencing

Those concerns should be addressed separately.

## Validation and Review Criteria

Future code and review work should check proposed changes against the following questions:

1. Does this implementation follow the ADR-defined target architecture directly?
2. Is this reusing old logic, or accidentally preserving old structure?
3. Does this move the rewrite forward, or recreate a migration path by habit?
4. Does this keep the public and internal contracts consistent with the ADR set?
5. Does this reduce or increase compatibility-driven clutter?

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
- ADR-008 defines the internal analysis-ready dataset builder architecture.
- ADR-009 defines the exception and error taxonomy.
- ADR-010 defines the internal package and module layout.
- ADR-011 defines the public builder API contract.
- ADR-012 defines how the rewrite should be executed as a fresh-start implementation.

Together, these ADRs establish:

- one public dataset model
- two public workflows
- one strict dataset boundary
- one coherent builder story
- one private validation domain
- one explicit failure taxonomy
- one domain-oriented code layout
- one fresh-start rewrite strategy

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR*. GitHub repository. https://github.com/PYangLab/PhosR

