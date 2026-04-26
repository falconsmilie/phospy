# ADR: Validation Domain Architecture for PhosPy

## Document Control

- **ADR ID:** ADR-007
- **Title:** Validation Domain Architecture for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines how validation should be structured in PhosPy. The package is being developed as a maintainable Python port of PhosR. To support that goal, validation must be explicit, reusable, private, and easy to compose without leaking into the public API as a separate product surface.

The decision is to establish validation as its own internal domain. Shared lower-level validation concerns should live in that domain, while workflow-level validators should act as composers that call shared validation components and then add workflow-specific rules. Validation methods should follow the project-wide `run(...)` convention.

## Status

Accepted.

This ADR defines the validation architecture that supports the public API, dataset boundary, reference handling, intensity-scale and processing-state contract, and workflow architecture established by earlier ADRs.

## Context and Problem Statement

Earlier ADRs have already made validation central to the architecture:

- the dataset boundary depends on strong dataset invariants
- reference handling depends on structural and compatibility validation
- intensity scale and processing state depend on validated typed contracts
- workflows are built around validator, interpreter, and executor stages

At the same time, the project direction is clear that validation should remain private. It should not become a standalone public API surface that users are expected to wire manually.

Without a dedicated validation domain, several problems tend to appear:

- repeated checks spread across workflows, builders, and services
- drift between similar validation rules implemented in different places
- workflow validators becoming bloated with low-level data checks
- public models carrying too much validation detail directly
- helper classes with overlapping responsibilities and unclear ownership

PhosPy therefore needs an explicit decision about where validation lives, how it is composed, and how it interacts with the workflow architecture.

## Decision Drivers

The decision is driven by the following considerations:

1. **Reusability.** Basic validation rules should be written once and reused.
2. **Clarity.** Validation responsibility should be easy to locate and reason about.
3. **Maintainability.** Shared rules should not be copied into multiple workflows or services.
4. **Workflow simplicity.** Workflow validators should compose shared validation rather than reimplement everything themselves.
5. **Privacy of implementation.** Validation should not become a separate public product surface.
6. **PhosR alignment.** The package should feel like a workflow-driven scientific tool, not a toolkit of free-floating validation helpers.

## Proposed Decision

Validation in PhosPy will be implemented as its own internal domain.

This validation domain will contain reusable lower-level validation components for shared concerns such as:

- dataset structure
- matrix and metadata alignment
- metadata column requirements
- reference bundle validity
- organism compatibility
- intensity-scale validity
- processing-state validity
- general numeric and structural constraints

Workflow-level validators will remain part of the workflow architecture, but they will act primarily as composers. They will call the shared validation-domain components and then apply workflow-specific rules. Their job is to validate successfully or fail clearly, not to reshape the request into a new form.

Validation remains private. Users should not be expected to work directly with validation classes as part of the normal public API.

## Core Design Principle

Validation should be **centralised and composable**, not **duplicated and scattered**.

The goal is not to create one giant validator. The goal is to create a validation domain where shared rules can be defined once and assembled cleanly where needed.

## Validation Domain Scope

The validation domain should own reusable validation concerns that are broader than any one workflow stage.

Examples include:

- checking required DataFrame properties
- checking index uniqueness
- checking matrix and metadata alignment
- checking required metadata columns
- checking non-empty structured resources
- checking enum-like support constraints
- checking intensity-scale-state presence and type
- checking processing-state presence and type
- checking organism compatibility rules

These are not workflow-specific business rules. They are shared validation building blocks.

## Workflow Validator Scope

Workflow validators remain important, but their role is narrower.

A workflow validator should:

- accept the public workflow request DTO
- call relevant shared validators from the validation domain
- apply workflow-specific constraints that only make sense in that workflow
- return the original request unchanged for pipeline continuity, or otherwise complete successfully without reshaping it into a new DTO

A workflow validator should not:

- reimplement basic dataset, reference, or intensity-scale checks already owned by the validation domain
- become a dumping ground for low-level helper logic
- expand into a second general-purpose validation framework
- create dedicated validated DTOs as part of its normal job

## Example Composition Direction

A workflow validator should read more like this:

```python
class KinaseWorkflowValidator:
    def __init__(self, dataset_validator, reference_validator, transform_validator) -> None:
        self._dataset_validator = dataset_validator
        self._reference_validator = reference_validator
        self._transform_validator = transform_validator

    def run(self, request: KinaseWorkflowRequest) -> KinaseWorkflowRequest:
        self._dataset_validator.run(request.dataset)
        self._transform_validator.run(request.dataset.intensity_scale_state)
        self._reference_validator.run(request.reference, request.dataset.organism)
        self._validate_workflow_specific_rules(request)
        return request
```

The exact class names may differ, but the shape should remain.

## Method Naming Convention

All validators should expose a `run(...)` method.

This keeps validation aligned with the broader workflow architecture convention already established elsewhere in the project.

Examples:

- `self._dataset_validator.run(dataset)`
- `self._reference_validator.run(bundle)`
- `self._workflow_validator.run(request)`

Verb-specific method names such as `validate(...)` are discouraged for these components.

## Validation Result Direction

Validation components should validate by succeeding or failing clearly.

The default direction is:

- validation components check conditions
- they either complete successfully or raise a validation exception
- they do not create new DTOs as part of normal validation work

Where pipeline continuity benefits from returning a value, the validator should return the original input unchanged rather than manufacturing a separate validated DTO.

This keeps the validator focused on validation rather than type reshaping. Meaningful reshaping belongs to later stages such as interpretation, not to validation itself.

## Error Strategy

Validation should fail clearly and specifically.

The validation domain should prefer domain-appropriate validation errors over accidental low-level exceptions such as:

- raw `KeyError`
- raw `IndexError`
- ambiguous `ValueError` with little context

Error messages should identify:

- what object or field was invalid
- what requirement failed
- enough context for the caller to fix the issue

A shared base validation exception type should exist for the validation domain so validation failures are explicit and can be handled consistently.

## Dataset Validation Direction

Dataset validation should be shared, central, and private.

It should cover concerns such as:

- phospho matrix presence and numeric shape
- metadata alignment
- required site metadata columns
- optional sample metadata alignment
- optional total matrix alignment
- intensity-scale-state presence and validity
- processing-state presence and validity

`AnalysisReadyPhosphoDataset` may still perform local invariant checks at construction time.

However, dataset-model construction should not become the place where broad shared validation logic is orchestrated. The preferred direction is:

- dataset construction keeps local invariant checks that are truly about the model itself
- builders and preprocessing paths use validation-domain components for broader shared validation before final construction
- workflows may rely on that established boundary and compose further workflow-specific checks where needed

This keeps the dataset model smaller, reduces coupling between the model and the wider validation domain, and avoids repeating deep validation orchestration in multiple places.

## Reference Validation Direction

Reference validation should live in the validation domain.

It should cover concerns such as:

- required `ReferenceBundle` fields
- non-empty resources
- organism support
- dataset/reference organism compatibility
- preset-resolution compatibility checks where appropriate

Workflow validators should compose these checks rather than duplicate them.

## Transformation Validation Direction

Intensity-scale-state validation should live in the validation domain.

It should cover concerns such as:

- presence of `IntensityScaleState`
- supported intensity scale kind
- consistency with dataset expectations
- rejection of unsupported or loosely declared intensity-scale claims

This keeps intensity-scale validation out of workflow executors and out of scattered dataset helpers.

## Builder and Preprocessing Interaction

Builders and preprocessing services may call validation-domain components directly where needed.

This is appropriate because they sit at the ingestion and dataset-construction boundary.

However, this should not dissolve the distinction between:

- reusable shared validation
- workflow-specific validation composition

The validation domain provides the building blocks. Builders and workflow validators each compose them for their own boundary.

## Public API Boundary

Validation remains outside the public product surface.

This means:

- users are not expected to instantiate validators directly
- validators are not a first-class public workflow alternative
- public API documentation should focus on datasets, workflows, requests, and results rather than teaching validation internals

Validation is critical to the product, but it is still internal architecture.

## Internal Package Direction

A likely package direction is to give validation its own internal domain package.

A likely shape is conceptually similar to:

```text
validation/
  datasets/
  references/
  transformations/
  common/
  workflows/
```

This keeps shared validation organised by concern without making workflow validators swallow everything.

The exact package layout may vary, but validation should be recognisable as its own domain.

## Consequences

### Positive consequences

- Shared validation rules are centralised and reusable.
- Workflow validators stay smaller and more readable.
- Validation ownership becomes clearer.
- Builder and workflow boundaries both benefit from the same lower-level rules.
- Error quality improves because validation becomes intentional rather than accidental.

### Negative consequences

- The architecture introduces another explicit internal domain that must be maintained.
- Some current validation logic will need to be extracted and reorganised.
- The project must stay disciplined to avoid both over-fragmentation and validator bloat.

### Neutral consequences

- Some model-level invariant checks may still remain on public models where appropriate.
- Builders, datasets, and workflows may all depend on the same validation domain without making it public.

## Rejected Alternatives

### Alternative 1: Keep validation embedded inside each workflow and builder

This option was rejected because it duplicates rules, encourages drift, and makes workflows harder to maintain.

### Alternative 2: Put all validation directly on public models only

This option was rejected because many reusable validation concerns cross model boundaries and should not be trapped inside one class.

### Alternative 3: Expose validators as a public API surface

This option was rejected because the package is intended to present a clean workflow-oriented product surface, not a toolbox of standalone validators.

### Alternative 4: Create one monolithic validator for the whole package

This option was rejected because it would centralise validation in the wrong way and would become a maintenance bottleneck.

## Resolved Decisions

The following decisions are now resolved for this ADR.

1. The validation domain should have a shared base validation exception type.
2. `AnalysisReadyPhosphoDataset` should keep local invariant checks at construction time, while builders and preprocessing paths use validation-domain components for broader shared validation.
3. Common validators should remain granular so they are easy to compose and reuse.
4. Validators should not create dedicated validated DTOs as part of their normal job.
5. Validation-domain components should not depend on each other across subdomains. Composition should happen at higher levels.

## Implementation Guidance

A likely healthy split is:

- validation domain for reusable lower-level rules
- workflow validators as composers
- builder and preprocessing services reusing shared validation where appropriate
- public models retaining only the invariant checks that are truly local to the model

Reviewers should reject changes that duplicate shared validation logic in multiple workflows or that attempt to make validators part of the public product story.

## Scope Boundaries

This ADR defines validation-domain architecture only.

It does not define:

- the full exception hierarchy
- the exact dataset contract
- the exact workflow request and result contracts
- exporter or visualisation behaviour
- migration strategy from current code

Those concerns should be handled separately.

## Validation and Review Criteria

Future code and review work should check proposed changes against the following questions:

1. Does this validation logic belong in the shared validation domain or in a workflow-specific validator?
2. Does this change reduce or increase duplication?
3. Does this keep validators private rather than public-facing?
4. Does this keep workflow validators acting as composers rather than low-level utility bins?
5. Does this improve or weaken validation clarity and error quality?

If the answers are weak or negative, the design should be reconsidered.

## Relationship to Earlier ADRs

This ADR complements the earlier architecture decisions.

- ADR-001 defines the intended public API contract.
- ADR-002 defines the internal workflow architecture.
- ADR-003 defines the dataset and preprocessing boundary.
- ADR-004 defines the reference resolution strategy and `ReferenceBundle` contract.
- ADR-005 defines result-model design.
- ADR-006 defines the intensity-scale and dataset-processing-state contract.
- ADR-007 defines how validation is organised and composed across the package.

Together, these ADRs establish:

- one public dataset model
- two public workflows
- one consistent internal workflow pattern
- one explicit reference-resolution path
- one disciplined result-model approach
- one stronger intensity-scale and processing-state contract
- one private validation domain with workflow-level composition

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR*. GitHub repository. https://github.com/PYangLab/PhosR
