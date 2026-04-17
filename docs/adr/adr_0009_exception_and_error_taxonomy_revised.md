# ADR: Exception and Error Taxonomy for PhosPy

## Document Control

- **ADR ID:** ADR-009
- **Title:** Exception and Error Taxonomy for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Authors:** OpenAI ChatGPT with project direction from the PhosPy maintainer
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines how failures should be represented and communicated in PhosPy. The package is being developed as a maintainable Python port of PhosR. To support that goal, exception handling must be explicit, consistent, and aligned with the architecture decisions already made for validation, dataset building, reference resolution, transformation handling, and workflow execution.

The decision is to establish a clear internal exception taxonomy with a small number of meaningful base exception types. Validation failures, input failures, builder/preprocessing failures, reference-resolution failures, transformation failures, and workflow-execution failures should each have explicit domain-appropriate exception types. Errors should fail early, fail clearly, and avoid leaking low-level implementation exceptions as the normal public failure story.

## Status

Accepted.

This ADR defines the exception and error taxonomy that supports the public API, validation domain, builder architecture, reference handling, transformation contract, and workflow architecture established by earlier ADRs.

## Context and Problem Statement

Several earlier ADRs already depend on explicit, high-quality failure behaviour:

- the validation domain should raise clear validation failures
- the builder should fail quickly on unsupported input conventions
- reference resolution should fail explicitly for unsupported or mismatched organism scenarios
- transformation state should be established through supported PhosPy paths rather than inferred or loosely declared
- workflows should remain simple and should not become dumping grounds for ad hoc failure translation

Without a defined exception taxonomy, several problems appear quickly:

- callers receive inconsistent exception types for similar failures
- low-level exceptions such as `KeyError` or pandas-specific failures leak through as the normal API story
- builders, validators, and workflows each invent their own failure patterns
- public API users get poor messages and unclear recovery guidance
- error handling becomes hard to test and hard to document

PhosPy therefore needs an explicit decision on how failures are grouped, where they originate, and how much of the taxonomy should be visible at the public boundary.

## Decision Drivers

The decision is driven by the following considerations:

1. **Clarity.** Similar failures should fail in similar ways.
2. **Maintainability.** Exception handling should not be improvised in each module.
3. **User experience.** Failures should be understandable and actionable.
4. **Architectural alignment.** Exception types should reinforce the boundaries already established in earlier ADRs.
5. **Testability.** Failure modes should be explicit enough to assert cleanly in tests.
6. **Privacy of implementation.** Internal library and framework exceptions should not define the public failure story by default.

## Proposed Decision

PhosPy will use a small, explicit exception taxonomy with one package-level base exception and a limited number of meaningful domain-level base exceptions.

A likely direction is:

- `PhosPyError` as the package-level base exception
- `PhosPyValidationError` for validation-domain failures
- `PhosPyInputError` for input-reading and input-shape failures at the ingestion boundary
- `PhosPyBuildError` for builder and preprocessing failures that are not merely input failures
- `PhosPyReferenceError` for reference-resolution and reference-compatibility failures
- `PhosPyTransformationError` for transformation-state and transformer failures
- `PhosPyWorkflowError` for workflow-stage orchestration and execution failures

These types may have narrower subtypes later where useful, but the initial taxonomy should remain intentionally small.

## Core Design Principle

Failures should be **domain-meaningful** and **boundary-aware**.

The exception type should help answer:

- what kind of failure occurred?
- at what boundary did it occur?
- what should the caller do next?

The goal is not to create a huge class hierarchy. The goal is to make failures predictable and useful.

## Package-Level Base Exception

A single package-level base exception should exist:

- `PhosPyError`

This creates a clean common parent for package-defined failures and gives callers a simple catch point when they need one.

All PhosPy-defined domain exceptions should derive from this base type.

## Validation Exceptions

Validation-domain failures should derive from:

- `PhosPyValidationError`

This aligns with ADR-007, which already established that the validation domain should have a shared base validation exception type.

Validation failures include cases such as:

- missing required fields
- invalid DataFrame structure
- invalid metadata alignment
- unsupported enum-like values
- missing `site_sequence` at the dataset boundary
- invalid transformation-state shape or type

Validation components should prefer raising this explicit type rather than leaking raw `KeyError`, `TypeError`, or generic `ValueError` as the intended failure mode.

## Input Exceptions

Input-boundary failures should derive from:

- `PhosPyInputError`

These failures include cases such as:

- unreadable file paths
- unsupported file types at the ingestion boundary
- input file structure that cannot be interpreted by the supported builder path
- malformed raw input payloads before dataset-building can meaningfully proceed

This keeps input failures separate from builder failures. An input that cannot even be read or interpreted into the supported ingestion shape is not the same thing as a dataset build failure.

## Builder and Preprocessing Exceptions

Builder and preprocessing failures should derive from:

- `PhosPyBuildError`

These failures include cases such as:

- failed site-sequence derivation at the builder boundary
- inability to complete dataset construction even though input was successfully read
- builder-path failure after ingestion has already produced supported raw structures
- preprocessing-stage failure while shaping a dataset toward the analysis-ready contract

This keeps build-boundary failures distinct from pure validation failures and distinct from input-reading failures.

## Reference Exceptions

Reference-related failures should derive from:

- `PhosPyReferenceError`

These failures include cases such as:

- unsupported organism for bundled reference resolution
- missing dataset organism when `AUTO` is used
- dataset/reference organism mismatch
- malformed `ReferenceBundle`
- unresolved preset-to-bundle mapping

This aligns with ADR-004 and keeps reference failures clearly separate from general workflow failures.

## Transformation Exceptions

Transformation-related failures should derive from:

- `PhosPyTransformationError`

These failures include cases such as:

- unsupported transformation kind
- invalid transformation state
- transformer failure to establish supported transformation state
- inconsistent transformation assumptions during dataset build

This aligns with ADR-006 and keeps transformation concerns explicitly identifiable.

## Workflow Exceptions

Workflow execution failures that are not better described as validation, input, build, reference, or transformation failures should derive from:

- `PhosPyWorkflowError`

These are failures such as:

- inconsistent internal workflow-stage assumptions
- orchestration-level execution failures after request validation and interpretation have already succeeded
- impossible internal states that indicate workflow failure rather than user-input invalidity

`PhosPyWorkflowError` is sufficient for impossible internal states at this stage. A separate dedicated impossible-state exception type is not currently needed.

Workflow exceptions should not become the default wrapper for everything. If a failure is more precisely a validation, input, build, reference, or transformation problem, it should keep that more specific type.

## Exception Translation Policy

Low-level exceptions from libraries and infrastructure should not normally leak through unchanged as the intended public failure story.

Examples include:

- pandas structural exceptions
- raw `KeyError`
- raw `IndexError`
- low-context `ValueError`
- file-system exceptions with poor domain framing

PhosPy components should translate such failures into the appropriate domain exception where the boundary meaning is clear.

Translation should be purposeful and boundary-aware. The application should not be designed around debugging-first exception flow. The intended failure contract should be explicit and meaningful for callers.

However, translation should still preserve useful debugging context. The preferred pattern is:

- raise the domain-appropriate PhosPy exception
- preserve the original exception as the cause where appropriate

This keeps both usability and traceability.

## Message Quality Policy

Exception messages should be explicit, meaningful, consistent, and actionable.

An exception is the path by which an error reaches the end user, so message formatting is part of the architectural contract.

A good failure message should identify:

- what object, field, or boundary failed
- what rule or expectation was violated
- enough context for the caller to correct the issue

Examples of strong message style:

- dataset is missing required `site_sequence` in `site_metadata`
- `ReferencePreset.AUTO` requires `dataset.organism`, but none was provided
- unsupported builder input convention for phospho input; expected one of the supported column naming patterns
- transformation state could not be established from the provided dataset input
- input file could not be read as a supported phosphoproteomics table

Messages should avoid vague wording such as:

- invalid input
- unsupported data
- failed to process

unless more context is also provided.

These message-formatting rules should remain documented as part of this ADR.

## Boundary-Aware Failure Direction

### Dataset boundary

Failures at the strict dataset boundary should normally be validation failures unless they originate from broader builder orchestration.

### Input boundary

Failures while reading or initially interpreting supported file-path or raw input payloads should normally be input failures.

### Builder boundary

Failures while turning successfully ingested user input into a strict dataset should normally be build failures, with nested validation or transformation causes where relevant.

### Reference boundary

Failures while resolving presets or validating bundle compatibility should normally be reference failures.

### Transformation boundary

Failures while establishing supported transformation state should normally be transformation failures.

### Workflow boundary

Failures after a request has passed validation and interpretation should normally be workflow failures, unless a more specific domain exception remains the truer description.

## Public API Visibility

Any exception type the end user is expected to handle should be publicly exported.

This means public export should follow practical handling expectations rather than an arbitrary internal/public split.

Exceptions are part of the public failure contract, but they are not the main product entry point into the package.

Users should still primarily interact with:

- datasets
- builders
- workflows
- results

## Testing Direction

Tests should assert domain-meaningful exceptions rather than generic failures wherever reasonable.

Examples:

- invalid dataset structure raises `PhosPyValidationError`
- unreadable or unsupported ingestion input raises `PhosPyInputError`
- failed dataset shaping after successful ingestion raises `PhosPyBuildError`
- organism mismatch raises `PhosPyReferenceError`
- unsupported transformation path raises `PhosPyTransformationError`

This strengthens failure semantics and reduces accidental drift.

## Internal Package Direction

A likely internal direction is to define exception types in a small dedicated module or package area rather than scattering them across many domains.

A likely shape is conceptually similar to:

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

The exact module layout may vary, but the taxonomy should remain recognisable and intentionally small.

## Consequences

### Positive consequences

- Failure behaviour becomes more consistent across the package.
- Public and internal boundaries gain clearer error semantics.
- Tests become easier to write and maintain.
- Builders, validators, and workflows can fail in domain-appropriate ways without improvising.
- User-facing error messages improve.

### Negative consequences

- Some current code will need explicit failure translation rather than relying on whatever low-level exception occurs.
- The project must stay disciplined so the taxonomy does not grow into an unnecessary hierarchy.
- Developers will need to decide deliberately which boundary a failure belongs to.

### Neutral consequences

- Narrower subtypes may still be introduced later if a real need appears.
- Some internal low-level exceptions may still surface as causes during debugging, but they should not define the intended failure contract.

## Rejected Alternatives

### Alternative 1: Let each module raise whatever native exception is convenient

This option was rejected because it produces inconsistent behaviour and poor public failure semantics.

### Alternative 2: Use one generic package exception for all failures

This option was rejected because it collapses too many different failure meanings into one coarse type and weakens caller handling.

### Alternative 3: Create a deep and highly granular exception hierarchy immediately

This option was rejected because it would over-design the taxonomy before real usage pressure justifies it.

### Alternative 4: Keep exceptions completely internal and undocumented

This option was rejected because failure behaviour is part of the practical public contract, even if it is not the main product entry point.

## Resolved Decisions

The following decisions are now resolved for this ADR.

1. Any exception type the end user is expected to handle should be publicly exported.
2. File-reading and raw input-interpretation failures at the ingestion boundary should be represented as `PhosPyInputError`, not `PhosPyBuildError`.
3. `PhosPyWorkflowError` is sufficient for impossible internal states for now.
4. Internal components should translate low-level exceptions into domain-appropriate PhosPy exceptions at meaningful boundaries rather than letting debugging-oriented failure behaviour define the application contract.
5. Error message formatting rules should remain documented as part of this ADR.

## Implementation Guidance

A likely healthy split is:

- one package-level base exception
- one shared validation-domain base exception
- a small number of boundary-aware domain exceptions
- local translation of low-level exceptions at meaningful boundaries
- preservation of original causes where appropriate

Reviewers should reject changes that default back to generic exception behaviour or that create new exception types without a clear boundary or semantic reason.

## Scope Boundaries

This ADR defines exception and error taxonomy only.

It does not define:

- the full validation-domain implementation
- the full public builder API
- the exact workflow request and result contracts
- logging or telemetry strategy
- migration strategy from current code

Those concerns should be addressed separately.

## Validation and Review Criteria

Future code and review work should check proposed changes against the following questions:

1. Does this failure use the most appropriate domain exception type?
2. Does this translate a low-level exception at the right boundary?
3. Does the message clearly explain what failed and why?
4. Does this keep the taxonomy small and meaningful?
5. Does this improve or weaken caller handling and testability?

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
- ADR-009 defines how failures are classified and communicated across the package.

Together, these ADRs establish:

- one public dataset model
- two public workflows
- one strict dataset boundary
- one flexible builder path below that boundary
- one private validation domain
- one stronger transformation-state contract
- one explicit, boundary-aware failure taxonomy

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR*. GitHub repository. https://github.com/PYangLab/PhosR

