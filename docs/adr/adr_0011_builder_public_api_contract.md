# ADR: Builder Public API Contract for PhosPy

## Document Control

- **ADR ID:** ADR-0011
- **Title:** Builder Public API Contract for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines the intended public API contract for building `AnalysisReadyPhosphoDataset` in PhosPy. The package is being developed as a maintainable Python port of PhosR. To support that goal, the public builder contract must give users a simple and realistic route from messy phosphoproteomics inputs to the strict dataset boundary, without exposing the full complexity of the internal builder architecture.

The decision is to provide one consistent public builder story backed by a small family of specialised internal builders or collaborators. The public builder contract should be class-based, accept either file-path inputs or already-loaded `DataFrame` inputs, keep convention detection simple and boring, fail quickly on unsupported input shapes, and return only `AnalysisReadyPhosphoDataset`.

## Status

Accepted.

This ADR defines the public builder contract that sits above the internal
builder architecture established in ADR-0008 and below the dataset boundary
defined in ADR-0003.

Update note (2026-05-11): the builder contract is derive-or-fail for
`site_metadata.site_sequence` at the analysis-ready boundary. Ingestion may
omit sequences, but final dataset construction must not proceed unless every
row has a valid non-empty `site_sequence`.

Update note (2026-05-11, provenance addition): successful builds attach a
structured preprocessing sequence-resolution summary
(`SiteSequenceResolutionReport`) so callers can audit sequence origin,
conflicts, applied conflict policy, unresolved counts, and final
sequence-complete site totals.

Update note (2026-05-14, boundary clarification): builder output guarantees
required `site_sequence` presence at the analysis-ready dataset boundary. It
does not automatically imply sequence-aware centred-context suitability for
every workflow lane. Sequence-aware workflow validators own strict centred
context checks (odd length, central residue match to site token, strict
character policy unless an explicit relaxation is configured).

## Context and Problem Statement

Earlier ADRs established two important truths:

- the public workflow boundary must stay strict and accept only `AnalysisReadyPhosphoDataset`
- users should not be forced to manually normalise ugly phosphoproteomics inputs before using the package

ADR-0008 established that the package should provide a builder-oriented path
below the strict dataset boundary. However, that ADR mainly defined the
architectural role of the builder layer, not its public contract.

Without an explicit public builder contract, several problems are likely:

- users will not know the intended route from raw data to analysis-ready data
- internal builder flexibility may leak into the public API in inconsistent ways
- the package may accumulate multiple competing ingestion entry points
- workflows may remain clean on paper but dataset construction may still feel ad hoc to users

PhosPy therefore needs a clear decision about what the end user is supposed to call when they want to build an `AnalysisReadyPhosphoDataset`.

## Decision Drivers

The decision is driven by the following considerations:

1. **Usability.** Users need a clear and simple public route from raw input to an analysis-ready dataset.
2. **Boundary discipline.** The builder public API should be flexible, but not vague or magical.
3. **Consistency.** The package should present one coherent builder story rather than many loosely related ingestion helpers.
4. **Maintainability.** Internal builder specialisation should not force a fragmented public surface.
5. **PhosR alignment.** The product should feel like a scientific workflow package, not a general ingestion framework.
6. **Failure quality.** Unsupported input should fail early and clearly rather than being guessed into acceptance.

## Decision

PhosPy will expose one consistent public builder route for constructing `AnalysisReadyPhosphoDataset`.

Internally, this route may delegate to a small family of specialised builders or collaborators, but the external story should remain coherent and simple.

The public builder contract should:

- be exposed primarily as a class rather than mixing classes and functions
- accept either file-path input or already-loaded `DataFrame` input
- support a small set of known input conventions
- avoid aggressive auto-recognition or heuristic guessing
- keep transformation handling implicit under preprocessing policy rather than exposing it as a prominent public choice
- keep the supported site-matrix lane intentionally narrow so the public route still returns a missing-value-free `AnalysisReadyPhosphoDataset`
- fail quickly on unsupported input conventions with clear and actionable messages
- return only `AnalysisReadyPhosphoDataset`

## Core Design Principle

The builder public API should be **simple for users** and **strict in outcome**.

Users should not need to understand the internal builder family. They should only need to understand how to provide supported inputs and receive a validated analysis-ready dataset.

## Public Builder Story

The public story should be that PhosPy offers one recommended route for creating `AnalysisReadyPhosphoDataset` from real-world phosphoproteomics data.

That route should be more prominent than direct manual dataset construction for typical messy industry inputs.

Direct construction of `AnalysisReadyPhosphoDataset` remains valid for callers who already have fully prepared data, but it is not the recommended story for ordinary ingestion.

## Public Builder Shape

The public builder contract should be represented through a single clear class-based entry point.

A likely direction would look conceptually like:

```python
builder = AnalysisReadyDatasetBuilder(...)
dataset = builder.run(request)
```

The exact builder class name may still be refined, but the class-based single-route principle should remain.

The design should not mix a primary class-oriented story with a second competing function-oriented story.

## Public Builder Input Contract

The public builder should accept a single structured request model rather than many loosely related arguments.

The public request DTO should be named:

- `DatasetBuildRequest`

A likely request shape can describe:

- phospho input
- optional total input
- optional site metadata input
- optional sample metadata input
- narrow input overrides where genuinely necessary
- limited public input-format hints where required for supported operation

The builder should not present a long scalar-heavy call signature.

## Supported Input Forms

The builder public contract should support two broad input forms:

### File-Path Input

Users may provide file paths for supported ingestion routes.

This keeps the builder practical for common real-world use.

Initial support should stay simple and focus on straightforward tabular formats rather than broad reader coverage from the start.

### `DataFrame` Input

Users may provide already-loaded pandas `DataFrame` objects.

This keeps the builder practical for notebooks, programmatic pipelines, and advanced callers who already control loading.

## Convention Recognition Policy

Convention recognition should remain intentionally limited.

The builder should support a small set of known and documented conventions for common phosphoproteomics input shapes.

It should not try to be clever about ambiguous input.

The public contract should favour:

- documented supported conventions
- narrow explicit overrides where needed
- early failure when the input does not match supported shapes

This avoids opening the package up to incorrect silent acceptance of bad data.

## Override and Hinting Policy

The public builder may support narrow overrides or hints where they are genuinely necessary to disambiguate otherwise supported input.

However, hinting should not go overboard. It must remain limited enough that the builder does not turn into a highly configurable ingestion framework.

A healthy direction is:

- simple default path first
- narrow hint or override path when needed
- fail clearly when the input still falls outside supported bounds

## Return Contract

The public builder returns exactly one primary public product:

- `AnalysisReadyPhosphoDataset`

It should not return:

- a semi-ready intermediate dataset
- a separate diagnostics object in the initial public contract
- multiple alternative result shapes

This keeps the builder aligned with the strict dataset boundary established elsewhere.

## Failure Contract

The builder should fail quickly, clearly, and meaningfully.

Typical failure categories should include:

- input failures when files or raw inputs cannot be read or interpreted into supported shapes
- build failures when ingestion succeeds but the builder cannot complete construction of an analysis-ready dataset
- validation failures when shaped data violates the strict dataset contract
- transformation failures when supported transformation state cannot be established

The public builder should not try to rescue unsupported inputs through magic adaptation.

## Relationship to Internal Builder Architecture

Internally, the public builder route may delegate to specialised builders or collaborators.

Examples of internal responsibilities may include:

- file/input readers
- input normalisation
- metadata shaping
- sequence derivation
- transformation handling
- validation composition
- final dataset construction

However, the public contract should not expose this internal structure directly.

The end user should see one consistent builder path, not the full internal orchestration graph.

## Relationship to Workflow Public API

The builder public API exists to get the user to the strict workflow boundary.

The intended public flow remains:

1. build `AnalysisReadyPhosphoDataset`
2. run `KinaseWorkflow`
3. run `SignalomeWorkflow`

The builder must therefore stay aligned with the dataset contract and must not create alternate public dataset forms that confuse the workflow story.

## Public API Surface Direction

The builder contract should be part of the public API.

A healthy direction is:

- builder defined in the stable `phospy.api` ownership namespace and
  re-exported from top-level `phospy` as the primary user-facing import route
- internal builder collaborators remain internal
- public docs present builder usage as the standard route for messy inputs

## Consequences

### Positive Consequences

- Users get a clear and practical route from raw input to `AnalysisReadyPhosphoDataset`.
- The builder public surface stays smaller than the internal builder architecture.
- The strict dataset boundary remains intact.
- Failure behaviour becomes easier to explain and document.
- The overall product story becomes more coherent.

### Negative Consequences

- The public builder contract must be designed carefully so it remains stable.
- Some user requests will still fail if they fall outside supported input conventions.
- Internal builder flexibility must be hidden without becoming confusing to maintainers.

### Neutral Consequences

- Direct dataset construction remains available for already-prepared callers.
- Internal builder collaborators may evolve as long as the public builder story remains consistent.

## Rejected Alternatives

### Alternative 1: Expose Many Separate Public Ingestion Helpers

This option was rejected because it fragments the public story and makes the package harder to understand.

### Alternative 2: Force All Users to Build `AnalysisReadyPhosphoDataset` Manually

This option was rejected because it is unrealistic for messy industry input formats.

### Alternative 3: Make the Public Builder Highly Magical and Aggressively Auto-Detect Everything

This option was rejected because it increases the risk of accepting bad data under incorrect assumptions.

### Alternative 4: Return Multiple Public Result Shapes From the Builder

This option was rejected because it weakens the strict dataset boundary and complicates the product story.

## Resolved Decisions

The following decisions are now resolved for this ADR.

1. The public builder should be exposed primarily as a class.
2. The public request DTO should be named `DatasetBuildRequest`.
3. Explicit input-format hinting should remain narrow and should not become highly configurable.
4. The builder should be part of the public API.
5. Initial file-path support should stay simple rather than trying to cover broad reader scenarios from the start.

## Implementation Guidance

A likely healthy direction is:

- one public builder class
- one structured public request DTO
- internal delegation to specialised builder collaborators
- narrow supported hints or overrides only where necessary
- clear exceptions when the builder cannot proceed

Reviewers should reject changes that fragment the public builder story or that reintroduce ingestion flexibility into the workflow layer.

## Scope Boundaries

This ADR defines the public builder API contract only.

It does not define:

- the full internal builder architecture beyond ADR-0008
- the full dataset contract beyond ADR-0003
- the full IO subsystem design beyond what is needed for the builder contract
- reference resolution strategy
- migration strategy from current code

Those concerns should be addressed separately.

## Validation and Review Criteria

Future code and review work should check future changes against the following questions:

1. Does this keep the public builder story coherent and singular?
2. Does this keep the builder flexible enough for real-world inputs without becoming magical?
3. Does this preserve one strict dataset output contract?
4. Does this keep internal builder specialisation hidden from normal users?
5. Does this make the public path from raw input to workflow execution clearer or more confusing?

If the answers are weak or negative, the design should be reconsidered.

## Relationship to Earlier ADRs

This ADR complements the earlier architecture decisions.

- ADR-0001 defines the intended public API contract.
- ADR-0002 defines the internal workflow architecture.
- ADR-0003 defines the dataset and preprocessing boundary.
- ADR-0004 defines the reference resolution strategy and `ReferenceBundle`
  contract.
- ADR-0005 defines result-model design.
- ADR-0006 defines the transformation-state and transformer contract.
- ADR-0007 defines the validation-domain architecture.
- ADR-0008 defines the internal analysis-ready dataset builder architecture.
- ADR-0009 defines the exception and error taxonomy.
- ADR-0010 defines the internal package and module layout.
- ADR-0011 defines the public builder API contract that exposes the builder
  story coherently.

Together, these ADRs establish:

- one public dataset model
- three public workflows
- one strict dataset boundary
- one flexible builder path below that boundary
- one coherent public builder story above that boundary
- one private validation domain
- one explicit failure taxonomy

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub. https://github.com/PYangLab/PhosR
