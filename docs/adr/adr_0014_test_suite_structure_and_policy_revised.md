# ADR: Test Suite Structure and Policy for PhosPy

## Document Control

- **ADR ID:** ADR-014
- **Title:** Test Suite Structure and Policy for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Authors:** OpenAI ChatGPT with project direction from the PhosPy maintainer
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines how the test suite should be structured and what each major test class is responsible for in PhosPy. The package is being developed as a maintainable Python port of PhosR. To support that goal, the test suite must reinforce the fresh-start architecture, distinguish clearly between ordinary behavioural testing and scientific parity testing, and keep failure signals easy to interpret.

The decision is to structure the test suite around three primary categories:

- unit tests
- integration tests
- parity tests

Each category should have a clear purpose, scope, and comparison policy. Unit tests should validate isolated components and local rules. Integration tests should validate multi-component flows within the new architecture. Parity tests should validate scientifically meaningful outputs against explicit reference expectations, primarily derived from PhosR outputs.

## Status

Accepted.

This ADR defines the test-suite policy that supports the fresh-start rewrite strategy in ADR-012 and the scientific parity strategy in ADR-013.

## Context and Problem Statement

The rewrite now has a strong architectural foundation:

- one strict dataset boundary
- one coherent builder story
- two public workflows
- validator, interpreter, executor workflow staging
- private validation domain
- explicit exception taxonomy
- scientific parity defined separately from architectural conformance

Without a clear test-suite structure, these good architectural decisions can still degrade during implementation. A poorly structured test suite tends to produce several problems:

- unit, integration, and parity concerns get mixed together
- failures become harder to interpret
- parity checks start acting as hidden contract tests or vice versa
- scientific comparisons become scattered and inconsistent
- contributors are unsure where new tests belong
- the suite becomes harder to maintain as the rewrite grows

PhosPy therefore needs an explicit decision about how the test suite is organised and what each layer of testing is supposed to prove.

## Decision Drivers

The decision is driven by the following considerations:

1. **Clarity of intent.** Test categories should answer different questions clearly.
2. **Architectural alignment.** The test suite should support the fresh-start ADR set rather than recreate legacy pressure.
3. **Maintainability.** Contributors should know where to place tests and how to write them.
4. **Scientific confidence.** Parity testing must be present without swallowing the rest of the suite.
5. **Failure readability.** Test failures should make it obvious what kind of expectation was violated.
6. **Scalability.** The suite should remain organised as the package grows.

## Proposed Decision

PhosPy will structure its test suite around three primary categories:

- unit
- integration
- parity

These categories should be physically separated in the test tree and conceptually separated in contributor expectations.

The package should not treat parity as just another form of integration testing, and it should not use unit tests to carry scientific-reference comparisons that belong in parity.

## Core Design Principle

Each test category should answer a different question.

- **Unit tests** ask whether a focused component behaves correctly in isolation.
- **Integration tests** ask whether multiple components work together correctly inside the new architecture.
- **Parity tests** ask whether the implementation produces the expected scientific outputs for meaningful fixtures.

If a test does not have a clear answer to one of those questions, it probably belongs in the wrong place.

## Proposed Test Tree Direction

A likely healthy top-level test structure is:

```text
tests/
  support/
  unit/
  integration/
  parity/
```

Additional internal organisation may exist beneath these directories, but this separation should remain the primary visible structure.

`tests/support/` is the preferred home for shared test helpers that are genuinely reusable across test categories.

## Unit Test Policy

Unit tests should focus on small, isolated responsibilities.

Examples include:

- model invariant checks
- validator behaviour for narrow rules
- builder collaborator behaviour in isolation
- reference provider behaviour for a focused case
- transformer behaviour for a focused case
- exception translation at a narrow boundary
- result-model contract checks at object level

Unit tests should:

- be fast
- minimise external dependencies
- avoid unnecessary fixture sprawl
- prefer focused assertions over broad end-to-end scenarios

Unit tests should not be used to validate full workflow parity or multi-stage scientific behaviour.

## Integration Test Policy

Integration tests should validate multi-component flows inside the new architecture.

Examples include:

- builder flow from supported input to `AnalysisReadyPhosphoDataset`
- end-to-end kinase workflow across validator, interpreter, executor, and domain services
- end-to-end signalome workflow using upstream kinase outputs
- reference-resolution flow for supported presets
- exception behaviour across meaningful package boundaries

Integration tests should prove that the new architecture works together as intended.

They should not be overloaded with the full responsibility of parity against external scientific references.

Fixture sharing between integration and parity suites is acceptable where it improves consistency and keeps both suites running against the same underlying data.

## Parity Test Policy

Parity tests should validate scientifically meaningful outputs against explicit reference expectations.

This ADR follows the parity strategy established in ADR-013.

Parity tests should:

- live in a dedicated parity directory
- use explicit fixtures
- use explicit reference-output provenance
- apply documented comparison rules
- focus on meaningful scientific outputs rather than incidental implementation details

The primary reference source for initial parity fixtures should be PhosR outputs.

The initial parity suite should aim to include at least the scientific output areas already present in the old application.

Parity fixtures do not need extra repository-level versioning or annotation machinery beyond being kept current and correct.

## Separation Rules

The suite should follow these separation rules.

### Unit tests must not

- depend on large scientific reference fixture sets unless the component itself is inherently fixture-driven
- attempt to prove full workflow parity
- encode broad end-to-end behaviour that belongs in integration tests

### Integration tests must not

- silently take on the entire parity burden
- act as a substitute for focused unit tests
- rely on unclear comparison rules for scientific outputs

### Parity tests must not

- become a hidden migration test suite for legacy structure
- assert internal class or module relationships
- preserve old wrapper or alias behaviour that the ADR set has rejected

## Comparison Policy by Test Category

### Unit tests

Prefer exact, narrow assertions wherever possible.

### Integration tests

Prefer exact behavioural assertions for contract and orchestration expectations, with focused tolerance-aware checks only where numerically necessary.

### Parity tests

Use explicit comparison rules per fixture or output family, including documented tolerances where exact comparison is not appropriate.

This keeps scientific comparison policy visible rather than hidden inside generic helper assertions.

## Fixture Policy

Fixture management should align with the three-part structure.

### Unit fixtures

Should stay small, local, and easy to understand.

### Integration fixtures

May be broader, but should still be targeted toward architectural flows rather than giant scientific archives.

### Parity fixtures

Should be curated carefully and include:

- explicit input data
- explicit expected output reference
- explicit provenance
- explicit comparison rules

Parity fixtures should be treated as high-value assets, not casual test data.

## Test Naming and Contributor Policy

Contributors should be able to answer two questions before adding a test:

1. What kind of test is this?
2. What architectural question is it answering?

A healthy contributor rule is:

- if the test targets one narrow component, it is probably a unit test
- if it validates a new-architecture flow across multiple components, it is probably an integration test
- if it compares scientific outputs to a reference baseline, it is probably a parity test

If that classification is unclear, the test design should be reconsidered.

For public API contract checks that are still narrow enough to be unit-like, the preferred direction is:

- keep them under `tests/unit/`
- group them under a clear `api/` subdirectory if helpful
- use names that make the contract focus explicit, such as:
  - `test_public_contract_dataset.py`
  - `test_public_contract_workflows.py`
  - `test_public_contract_results.py`

The key rule is that the filename should say `public_contract` when that is the primary concern.

## CI and Execution Direction

The suite should be runnable as a whole, but the categories should remain separable.

A healthy direction is:

- unit tests are expected to run frequently and quickly
- integration tests are expected to run regularly as part of architectural confidence
- parity tests are expected to run as a distinct scientific-confidence layer

Parity tests should not run by default. They should be run only when explicitly requested by the user at pytest invocation time.

The exact CI policy may be defined elsewhere, but the suite structure should support separate execution and reporting.

## Failure Interpretation Policy

Test failures should be easy to classify.

### Unit failure

Usually means a local rule or component behaviour is wrong.

### Integration failure

Usually means the new architecture is not composing correctly.

### Parity failure

Usually means either:

- the scientific output differs from the reference expectation
- the reference expectation is wrong or outdated
- the comparison rule needs review

This classification helps developers understand what to investigate first.

## Relationship to Exception and Validation Policy

The test suite should explicitly cover:

- public and internal exception taxonomy behaviour where appropriate
- validation-domain failure behaviour
- builder failure behaviour
- workflow failure boundaries

However, these checks should still be placed in the appropriate suite category rather than treated as a separate fourth top-level test class.

## Minimum Coverage Direction

At a minimum, the suite should ensure:

- public API contracts are covered
- dataset builder flow is covered
- kinase workflow flow is covered
- signalome workflow flow is covered
- parity coverage reaches at least the level already attained in the old application

This does not mean every area must have equal test volume. It means no core public capability should be left without meaningful test coverage.

## Consequences

### Positive consequences

- The test suite becomes easier to navigate and maintain.
- Failures become easier to interpret.
- Scientific parity remains visible without overwhelming ordinary development tests.
- Contributors have clearer guidance for where tests belong.
- The suite reinforces the fresh-start architecture rather than legacy structure.

### Negative consequences

- Contributors will need to think more deliberately about test placement.
- Some existing or future tests may need to be moved to preserve the intended separation.
- Parity fixtures will require more discipline than casual test data.

### Neutral consequences

- Additional sub-structure may still exist under each major test directory.
- Some tests may still require judgement calls when they sit near the boundary between categories.

## Rejected Alternatives

### Alternative 1: Keep all tests in one undifferentiated directory

This option was rejected because it weakens the distinction between contract, integration, and parity concerns.

### Alternative 2: Treat parity as a subtype of integration testing only

This option was rejected because parity answers a distinct scientific-confidence question and needs a clearer identity.

### Alternative 3: Add many more top-level test categories immediately

This option was rejected because it would overcomplicate the suite before the rewrite has stabilised.

### Alternative 4: Use only unit and parity tests

This option was rejected because integration tests play an important role in validating the architecture as assembled.

## Resolved Decisions

The following decisions are now resolved for this ADR.

1. Shared test helpers should live in a common support module under `tests/`.
2. Fixture sharing between integration and parity suites is acceptable.
3. Parity fixtures do not need additional repository-level versioning or metadata annotation machinery beyond being kept current.
4. Narrow public API contract tests should stay in the unit suite and should use `public_contract` in their naming.
5. Parity tests should not run by default and should be requested explicitly at pytest invocation time.

## Implementation Guidance

A likely healthy direction is:

- create `tests/support`, `tests/unit`, `tests/integration`, and `tests/parity`
- keep unit tests focused and fast
- use integration tests for assembled architectural flows
- treat parity fixtures and comparison helpers as deliberate scientific assets
- document contributor expectations clearly in repository documentation
- ensure parity execution is opt-in rather than default

Reviewers should reject tests that blur category boundaries without a strong reason.

## Scope Boundaries

This ADR defines test suite structure and policy only.

It does not define:

- the full CI execution policy
- the exact fixture repository layout beyond the high-level separation policy
- the detailed scientific parity rules beyond ADR-013
- performance benchmarking policy
- release gating thresholds

Those concerns should be addressed separately.

## Validation and Review Criteria

Future code and review work should check proposed changes against the following questions:

1. Is this test in the right category?
2. Does it answer the right kind of question for that category?
3. Does it keep scientific parity concerns separate from ordinary behavioural checks?
4. Does it make future maintenance easier or harder?
5. Does it reinforce or weaken the fresh-start architecture?

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
- ADR-012 defines the fresh-start rewrite roadmap.
- ADR-013 defines scientific parity strategy and parity-testing policy.
- ADR-014 defines how the overall test suite should be structured around unit, integration, and parity concerns.

Together, these ADRs establish:

- one public dataset model
- two public workflows
- one coherent builder story
- one fresh-start rewrite strategy
- one scientific-parity policy
- one explicit three-part test suite structure

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR*. GitHub repository. https://github.com/PYangLab/PhosR

