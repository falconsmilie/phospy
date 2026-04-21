# ADR: Transformation-State and Transformer Contract for PhosPy

## Document Control

- **ADR ID:** ADR-006
- **Title:** Transformation-State and Transformer Contract for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines how transformation state should be represented and established in PhosPy. The package is being developed as a maintainable Python port of PhosR. To support that goal, the dataset boundary must carry stronger guarantees than a simple free-text intensity label.

The decision is to represent transformation state through a typed contract established by a transformer component. `AnalysisReadyPhosphoDataset` should carry validated transformation state rather than relying on an informal string such as `"log2"` or `"linear"` alone. Transformation logic belongs to preprocessing and builder paths, not to workflows.

## Status

Accepted.

This ADR defines the transformation-state contract that supports the dataset boundary established in ADR-003 and the public API direction established in ADR-001.

## Context and Problem Statement

ADR-003 established that a simple declared `intensity_scale` label is not sufficient as the long-term design for PhosPy. The dataset boundary should carry stronger guarantees about transformation state through a transformer-oriented design.

This matters because transformation state is not just descriptive metadata. It affects how input values are interpreted, how downstream logic reasons about expected data, and whether a dataset is truly analysis-ready for PhosR-style workflows.

If transformation state is represented only as a free-text label, several problems arise:

- the label may be inaccurate or stale
- users or builders may set it without actually establishing the corresponding state
- workflows may assume guarantees that were never enforced
- validation becomes weaker because the state is declarative rather than proven

PhosPy therefore needs an explicit contract for how transformation state is created, validated, and carried into the public dataset model.

## Decision Drivers

The decision is driven by the following considerations:

1. **Scientific correctness.** Transformation state affects the meaning of quantitative input data.
2. **Boundary honesty.** The analysis-ready dataset should carry real guarantees, not just descriptive labels.
3. **Maintainability.** Transformation logic should be centralised rather than reinterpreted in multiple places.
4. **Extensibility.** Different supported transformation paths should fit a common contract.
5. **Workflow simplicity.** Workflows should consume established transformation state rather than infer or repair it.
6. **Type clarity.** A dedicated transformation model is more explicit than a bare string field.

## Proposed Decision

Transformation state in PhosPy will be represented through a typed transformation contract established by a transformer component.

The public dataset boundary should carry validated transformation state rather than relying on a free-text label alone.

A transformer is responsible for establishing that state during preprocessing or dataset-building paths.

Workflows consume the established transformation state but do not create it, infer it heuristically, or reinterpret it.

## Core Design Principle

Transformation state must be **established**, not merely **declared**.

A dataset should not be considered analysis-ready simply because it carries a text label claiming a scale or transform. Instead, the dataset should carry a transformation state object that records the supported transformation contract applied by PhosPy.

## Proposed Public Direction

The dataset contract should move away from a bare public `intensity_scale` string as the primary source of truth.

The preferred long-term direction is:

- `AnalysisReadyPhosphoDataset` carries a typed transformation-state field
- any user-facing convenience label such as `log2` or `linear` is derived from that state rather than acting as the authoritative representation

This preserves usability while keeping the real guarantee explicit.

## Proposed Transformation-State Model

A transformation-state model should represent the established state of quantitative intensities.

A likely public shape is conceptually similar to:

- transformation kind
- any narrow metadata required to interpret the state safely

The authoritative public type should be named `TransformationState`.

Examples of transformation kinds that may be supported initially include:

- linear
- log2

The initial contract should stay intentionally small.

Separate phospho and total matrices should share one transformation-state contract. Distinct transformation states are not part of the intended design at this stage.

## Transformer Contract

A transformer is the component that establishes transformation state.

A likely interface direction is:

- transformer accepts quantitative input state
- transformer returns transformed data plus validated transformation state
- transformer exposes `run(...)` in line with the project naming convention

Conceptually:

```python
class TransformerInterface(Protocol):
    def run(self, data: pd.DataFrame) -> TransformationResult: ...
```

The class and DTO names may change, but the contract direction should remain.

## Transformer Responsibilities

A transformer may:

- apply a supported transformation path
- return the resulting transformation state object
- reject unsupported or invalid transformation scenarios

A transformer may not:

- hide transformation assumptions behind an unverified label
- offload transformation-state establishment to later workflow stages
- silently guess unsupported transformations from ambiguous input
- treat externally declared transformation state as equivalent to a PhosPy-established state

Transformation state should always be established within PhosPy through the supported transformer path.

## Dataset Responsibilities

`AnalysisReadyPhosphoDataset` should carry the resulting transformation state as part of its validated boundary.

The dataset may expose a derived convenience label for readability, but the authoritative contract should be the typed transformation state.

Any user-facing label should use honest naming and should not imply that a loose scale string is the true contract.

The dataset itself should not contain transformation logic beyond validating that a valid transformation state object is present.

## Preprocessing and Builder Responsibilities

Preprocessing and builder paths are responsible for invoking the appropriate transformer.

This means:

- raw or semi-structured input may be normalised and transformed before final dataset construction
- the final public dataset should already contain established transformation state
- builders should not allow a dataset to cross the public boundary with only an informal free-text transform claim

This supports the user goal of flexible input handling without weakening the analysis-ready dataset contract.

Transformation should not be exposed as a public builder choice. It should follow preprocessing policy and keep the public surface smaller.

## Workflow Responsibilities

Workflows may assume that transformation state has already been established on the incoming dataset.

Workflows may still validate workflow-specific compatibility constraints if needed, but they should not:

- transform intensities as part of normal workflow execution
- infer missing transformation state
- reinterpret a convenience label as though it were a validated guarantee

Transformation belongs to preprocessing and dataset-building, not workflow execution.

## Validation Responsibility

Validation for transformation state belongs in the validation domain.

This includes reusable concerns such as:

- whether the transformation state object is present
- whether it is supported
- whether it is internally consistent with the dataset contract

Workflow-level validators may compose those shared checks with workflow-specific rules, but they should not replace the underlying validation domain.

## Supported Initial Scope

The initial supported transformation scope should remain intentionally small.

A sensible starting point is:

- linear
- log2

Additional transformation types should be future work unless there is a strong scientific reason to introduce them immediately.

## Convenience Label Policy

A convenience label may still exist for usability, but it must not be the authoritative source of truth.

The correct direction is:

- typed transformation state is authoritative
- any label such as `log2` is derived from that state

The convenience label should use honest naming. `intensity_scale` should not remain the primary naming idea for this concept.

## Consequences

### Positive consequences

- The dataset boundary becomes more honest and scientifically meaningful.
- Transformation assumptions are centralised and testable.
- Workflows become simpler because they consume established state.
- The package gains a clear seam for future transformation support.
- User-facing convenience can still exist without weakening the contract.

### Negative consequences

- The design is more explicit and slightly more structured than a simple string field.
- Builders and preprocessing paths must do more real work before dataset construction.
- Some current code that treats scale as a loose label will need to be rewritten.

### Neutral consequences

- Internal transformers may still vary as long as they satisfy the same contract.
- The public API may still expose a human-readable derived transform view if it is not treated as the source of truth.

## Rejected Alternatives

### Alternative 1: Keep `intensity_scale` as a free-text label only

This option was rejected because it provides weak guarantees and allows the dataset boundary to claim a state that may never have been established.

### Alternative 2: Let workflows interpret or repair transformation state

This option was rejected because it leaks preprocessing concerns into workflow execution and weakens the meaning of the analysis-ready dataset.

### Alternative 3: Infer transformation state heuristically from input values

This option was rejected because heuristic inference is brittle, hard to explain, and easy to get wrong scientifically.

### Alternative 4: Accept externally declared transformation state as equivalent to a PhosPy-established state

This option was rejected because transformation state in PhosPy should always be established through the supported transformer path.

## Resolved Decisions

The following decisions are now resolved for this ADR.

1. The authoritative public type name should be `TransformationState`.
2. Transformation state should always be established within PhosPy through the supported transformer path.
3. Separate phospho and total matrices should share one transformation-state contract.
4. Any derived user-facing label should be renamed more honestly rather than retaining `intensity_scale` as the primary naming idea.
5. Transformation should not be exposed as a public builder choice. It should follow preprocessing policy and keep the public surface smaller.

## Implementation Guidance

A likely healthy split is:

- `TransformationState` as a typed value object or DTO
- transformer interface for establishing that state
- preprocessing or builder services that invoke the transformer
- validation-domain support for reusable transformation checks
- datasets that validate the presence of the resulting state

Reviewers should reject changes that reduce transformation state back to a bare free-text field or that move transformation establishment into workflows.

## Scope Boundaries

This ADR defines the transformation-state and transformer contract only.

It does not define:

- the full preprocessing pipeline
- the full dataset contract beyond the transformation aspect
- workflow result design
- exporter or visualisation APIs
- migration strategy from current code

Those concerns should be addressed separately.

## Validation and Review Criteria

Future code and review work should check proposed changes against the following questions:

1. Does this change strengthen or weaken the honesty of the dataset boundary?
2. Does this keep transformation establishment in preprocessing and builders?
3. Does this preserve a typed transformation contract rather than a loose label?
4. Does this avoid heuristic guessing of transformation state?
5. Does this make future transformer extension easier without bloating workflow APIs?

If the answers are weak or negative, the design should be reconsidered.

## Relationship to Earlier ADRs

This ADR complements the earlier architecture decisions.

- ADR-001 defines the intended public API contract.
- ADR-002 defines the internal workflow architecture.
- ADR-003 defines the dataset and preprocessing boundary.
- ADR-004 defines the reference resolution strategy and `ReferenceBundle` contract.
- ADR-005 defines result-model design.
- ADR-006 defines how transformation state is established and carried through the dataset boundary.

Together, these ADRs establish:

- one public dataset model
- two public workflows
- one consistent internal workflow pattern
- one explicit reference-resolution path
- one disciplined result-model approach
- one stronger transformation-state contract at the dataset boundary

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR*. GitHub repository. https://github.com/PYangLab/PhosR

