# ADR: Scientific Parity Strategy and Parity-Testing Policy for PhosPy

## Document Control

- **ADR ID:** ADR-013
- **Title:** Scientific Parity Strategy and Parity-Testing Policy for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines what parity means for PhosPy and how parity should be tested during the fresh-start rewrite. The package is being developed as a maintainable Python port of PhosR. To support that goal, parity must be framed around scientific behaviour and meaningful workflow outputs rather than around preservation of the existing application structure.

The decision is that PhosPy should pursue scientific parity, not structural or legacy-implementation parity. The primary parity target is the intended PhosR-aligned behaviour at stable workflow boundaries and scientifically meaningful outputs. The existing application may be used as a secondary reference source where useful, but it is not the authoritative definition of parity. Parity tests should be explicit, fixture-driven, tolerance-aware where necessary, and clearly separated from ordinary contract and unit tests.

## Status

Accepted.

This ADR defines the parity strategy that supports the fresh-start rewrite defined in ADR-012 and the broader architectural decisions already established for the public API, workflows, dataset boundary, references, transformation state, validation, builders, and exceptions.

## Context and Problem Statement

The rewrite is explicitly a fresh-start implementation rather than a migration of the current application structure. That decision creates an immediate follow-on question:

- what does the new implementation need to match?

For a PhosR port, the wrong answer is:

- every internal class, layer, and historical behaviour of the current application

That would recreate migration pressure and pull the rewrite back toward the old structure.

At the same time, the rewrite still needs a disciplined way to confirm that reintroduced scientific logic is behaving correctly. Without an explicit parity strategy, the team risks one of two bad extremes:

- treating the existing application as the hidden authoritative baseline for everything
- avoiding parity altogether and relying only on contract tests, leaving scientific confidence weaker than it should be

PhosPy therefore needs an explicit decision about what parity means, what outputs matter, what reference sources matter, and how parity tests should fit alongside ordinary tests.

## Decision Drivers

The decision is driven by the following considerations:

1. **Scientific fidelity.** PhosPy is intended to be a PhosR port, so scientifically meaningful behaviour matters.
2. **Fresh-start integrity.** Parity must not become a back door for preserving legacy structure.
3. **Test clarity.** Parity tests should be distinguishable from unit, contract, and validation tests.
4. **Practicality.** Numerical workflows often require tolerance-aware comparison rather than brittle exact equality.
5. **Maintainability.** The parity policy must stay focused on meaningful outputs, not every incidental intermediate detail.
6. **Confidence during rewrite.** Selective parity checks can reduce risk while scientific logic is reintroduced.

## Proposed Decision

PhosPy will pursue **scientific parity** rather than **structural parity**.

The primary parity target is the intended PhosR-aligned scientific behaviour at stable workflow and stage-output boundaries.

The existing application may be used as a secondary reference source where useful, but it is not the authoritative definition of parity and it must not dictate the architecture of the rewrite.

Parity tests should be:

- explicit
- fixture-driven
- targeted at meaningful scientific outputs
- tolerance-aware where numerical comparison requires it
- separate from ordinary unit, contract, and validation tests

## Core Design Principle

Parity should compare **scientific outcomes** and **stable workflow outputs**, not **legacy structure**.

The rewrite should ask:

- does this new implementation produce the expected scientific result?

It should not ask:

- does this preserve the old internal class graph?
- does this keep every historical convenience behaviour?
- does this match legacy implementation detail that the ADRs have already rejected?

## Definition of Parity

For PhosPy, parity means that the new implementation reproduces the intended scientific behaviour of the target analysis flow closely enough to be considered faithful to the PhosR-style workflow.

Parity does not mean exact preservation of:

- internal package layout
- helper boundaries
- wrapper classes
- legacy result aliases
- historical orchestration layers
- accidental implementation details of the previous application

## Primary Parity Target

The primary parity target should be PhosR-aligned scientific behaviour.

In practical terms, this means parity should be judged against outputs that matter scientifically, such as:

- analysis-ready dataset construction outcomes where the contract is stable
- kinase scoring outputs
- kinase prediction outputs
- kinase activity outputs
- signalome outputs

Parity should be judged at stable boundaries, not at arbitrary low-level implementation points.

## Secondary Reference Sources

The existing application may still be used as a secondary reference source where useful.

Appropriate uses include:

- extracting expected outputs for selected scientific fixtures
- confirming that logic previously believed correct still behaves as expected
- identifying existing test scenarios worth porting into the new parity suite

However, the existing application must not be treated as the unquestioned authority for parity. It is reference material, not the design baseline.

## Recommended Parity Levels

A healthy parity policy should distinguish several levels.

### Level 1: Public Contract Tests

These are not parity tests. They confirm that the new public API behaves as designed.

Examples include:

- builder returns `AnalysisReadyPhosphoDataset`
- workflows accept the correct request DTOs
- result models have the intended shape
- exceptions and validation behaviour follow the new contract

### Level 2: Scientific Parity Tests

These are the core parity tests.

They compare scientifically meaningful outputs for selected fixtures.

Examples include:

- expected scoring matrices
- expected prediction matrices
- expected activity outputs
- expected signalome outputs

### Level 3: Internal-Stage Checks Where Justified

These are allowed only where a stage output is both stable and scientifically meaningful.

Examples may include:

- a scoring-stage output that is central to downstream interpretation
- a reference-resolution outcome that must match expected supported behaviour

These should not expand into exhaustive testing of every internal intermediate.

## Parity Scope Policy

Parity should focus on outputs that satisfy all of the following:

1. scientifically meaningful
2. stable enough to compare repeatedly
3. important enough that drift would matter
4. tied to the public or stable internal story rather than incidental implementation detail

If an output does not meet those criteria, it probably does not belong in the parity suite.

## Dataset Parity Direction

For dataset construction, parity should focus on the strict analysis-ready boundary rather than on every ingestion step.

That means parity should compare outcomes such as:

- required fields present
- metadata alignment
- derived `site_sequence`
- established transformation state
- final analysis-ready structure

Parity should not be defined as matching every intermediate preprocessing decision from the old application if those decisions are not part of the stable dataset contract.

## Workflow Parity Direction

Workflow parity should focus on meaningful stage and final outputs.

Examples include:

- kinase scoring outputs
- prediction outputs such as `pred_mat`
- activity outputs
- signalome outputs

Where a workflow output depends on floating-point or tolerance-sensitive logic, parity should use explicit comparison rules rather than brittle exact equality.

## Numerical Comparison Policy

Parity tests must acknowledge that scientific and numerical workflows often require tolerance-aware comparison.

The policy should therefore be:

- use exact equality only where the output is truly deterministic and exact comparison is appropriate
- use explicit tolerances where floating-point behaviour, ordering, or platform-level differences make exact equality too brittle
- document the comparison method used for each parity fixture or parity test family

Tolerances must be deliberate, not vague. A tolerance should reflect meaningful scientific and numerical expectations rather than simply making tests pass.

Where a stable tolerance policy is not yet known, the rewrite should begin with a conservative documented comparison rule per output family and refine it as parity fixtures accumulate. Tolerance choice should remain explicit and reviewable rather than hidden in ad hoc assertions.

## Fixture Policy

Parity tests should be fixture-driven.

A fixture should include:

- clear input data
- clear expected output reference
- clear comparison rule
- clear provenance of how the expected output was obtained

Fixture sets should stay selective and meaningful. The goal is not to build a giant archive of every possible example.

## Legacy Science Coverage Inventory Policy (2026-04-20)

Legacy science tests are acceptable donors, but not active parity authorities.
When a donor scenario is still scientifically relevant, the blocking check must
be rewrite-owned under `tests/{unit,integration,parity}` with rewrite-owned
fixture paths under `tests/fixtures/rewrite_parity/**` or
`tests/fixtures/public_workflow_reference/**`.

Historical archive trees should not be the normal fixture source for active
rewrite parity.

Rewrite-side visibility of full legacy-science coverage should be maintained
via an explicit inventory and check that:

- classifies each legacy-science area as `PORTED`, `INTENTIONALLY_RETIRED`,
  `OPEN_GAP`, or `CONTRACT_CHANGED`
- distinguishes scoped parity passes from full legacy-science coverage
- prevents blanket "no open gaps" language when inventory entries remain open

## Reference Output Policy

Expected parity outputs should primarily come from PhosR outputs where available and appropriate.

Selected outputs from the existing application may still be used as a secondary reference source when they help confirm behaviour that already exists and is worth preserving scientifically.

The reference source for a parity fixture should always be explicit.

A parity fixture should not be treated as self-explanatory. Its origin matters.

## Separation From Ordinary Tests

Parity tests should be separate from ordinary unit and contract tests.

This separation is important because parity answers a different question.

Unit and contract tests ask:

- does this component satisfy the new architecture and API contract?

Parity tests ask:

- does this implementation still produce the expected scientific outcome for this meaningful scenario?

Keeping them separate prevents both kinds of tests from becoming muddled.

A separate parity test directory is the preferred direction. The overall testing story should therefore distinguish at least:

- unit tests
- integration tests
- parity tests

## Failure Interpretation Policy

When a parity test fails, the first question should be:

- is the new output scientifically wrong, or is the reference expectation no longer the correct target?

Parity failures should not automatically be treated as bugs in the new implementation. They may also reveal:

- incorrect or outdated reference outputs
- tolerance rules that need refinement
- scientifically improved behaviour that differs from the old application
- hidden assumptions in the previous implementation

Parity failures therefore require review, not blind rollback toward the legacy result.

## Relationship to Fresh-Start Rewrite

Parity must support the fresh-start rewrite, not undermine it.

This means:

- parity is used to confirm scientific behaviour
- parity is not used to justify preservation of old architecture
- parity is not used as a reason to reintroduce wrappers, aliases, or structural compromise

The architecture remains ADR-led. Parity provides scientific confidence within that architecture.

## Consequences

### Positive Consequences

- The rewrite gains a disciplined way to confirm scientific fidelity.
- Parity stays focused on what matters rather than preserving legacy structure.
- Contract tests and parity tests remain conceptually distinct.
- Reference outputs become more explicit and reviewable.
- Numerical comparisons gain a deliberate policy rather than ad hoc assertion style.

### Negative Consequences

- Parity fixtures and expected outputs must be curated carefully.
- The team will need to review parity failures thoughtfully rather than treating them as purely mechanical.
- Some scientifically meaningful outputs may still be hard to compare exactly.

### Neutral Consequences

- The old application may still contribute selected reference outputs without being treated as the architecture baseline.
- Not every stage or dataset behaviour will necessarily receive parity coverage.

## Rejected Alternatives

### Alternative 1: Treat the Existing Application as the Authoritative Parity Baseline for Everything

This option was rejected because it would turn parity into a hidden migration strategy and would weaken the fresh-start rewrite stance.

### Alternative 2: Avoid Parity Entirely and Rely Only on New Contract Tests

This option was rejected because contract tests alone are not enough to establish confidence in scientific fidelity for a PhosR port.

### Alternative 3: Require Parity for Every Intermediate Internal Detail

This option was rejected because it would over-constrain the rewrite and drag legacy implementation detail into the new architecture.

### Alternative 4: Use Vague Visual or Informal Comparison of Outputs Only

This option was rejected because parity needs explicit and testable comparison rules.

## Resolved Decisions

The following decisions are now resolved for this ADR.

1. The initial parity suite should aim to include at least the scientific output areas already present in the old application.
2. The primary reference source for initial parity fixtures should be PhosR outputs.
3. Parity tests should live in a separate test directory, distinct from unit and integration tests.
4. The minimum parity coverage required before the old application is removed entirely as reference material is the level already attained at that point.
5. Numerical comparison tolerances should remain explicit per output family and should be refined deliberately rather than guessed implicitly.

## Implementation Guidance

A likely healthy direction is:

- keep contract tests and parity tests separate
- begin with a small, meaningful parity fixture set
- choose explicit PhosR-derived reference outputs first
- compare stable, scientifically meaningful outputs first
- expand parity coverage at least to the level already attained in the old application
- refine tolerance rules deliberately as real fixtures are introduced

Reviewers should reject parity proposals that effectively preserve legacy structure instead of validating scientific behaviour.

## Scope Boundaries

This ADR defines scientific parity strategy and parity-testing policy only.

It does not define:

- the detailed implementation of the scientific algorithms
- the exact file layout of the test suite beyond the high-level separation policy
- release or acceptance thresholds for versions
- the detailed migration/removal timing of every old module beyond ADR-012

Those concerns should be addressed separately.

## Validation and Review Criteria

Future code and review work should check proposed changes against the following questions:

1. Does this parity check validate a scientifically meaningful outcome?
2. Does this avoid binding the rewrite to legacy structure?
3. Is the reference source explicit and justified?
4. Is the comparison rule explicit and appropriate for the output type?
5. Does this strengthen confidence without bloating the rewrite with incidental parity work?

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
- ADR-013 defines what scientific parity means and how parity testing should be applied during the rewrite.

Together, these ADRs establish:

- one public dataset model
- two public workflows
- one strict dataset boundary
- one coherent builder story
- one fresh-start rewrite strategy
- one scientific-parity policy that supports fidelity without forcing legacy structure forward

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR*. GitHub repository. https://github.com/PYangLab/PhosR
