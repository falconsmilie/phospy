# ADR: Analysis-Ready Dataset Builder Architecture for PhosPy

## Document Control

- **ADR ID:** ADR-0008
- **Title:** Analysis-Ready Dataset Builder Architecture for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines how PhosPy should turn messy phosphoproteomics inputs into `AnalysisReadyPhosphoDataset`. The package is being developed as a maintainable Python port of PhosR. To support that goal, the public workflow boundary must stay strict, while the ingestion and builder boundary remains flexible enough to handle real-world input variation.

The decision is to introduce a dedicated analysis-ready dataset builder path that accepts flexible user input, composes preprocessing services, validation-domain components, sequence derivation, and transformation handling, and returns a validated `AnalysisReadyPhosphoDataset`. The builder should be the main public route for converting industry-style inputs into the strict dataset contract.

## Status

Accepted.

This ADR defines the builder and ingestion architecture that supports the
dataset boundary established in ADR-0003, the transformation-state contract
established in ADR-0006, and the validation-domain architecture established in
ADR-0007.

Update note (2026-05-11): `site_sequence` may be omitted at ingestion, but it
is mandatory at the `AnalysisReadyPhosphoDataset` boundary. The builder owns
the derive-or-fail transition before final dataset construction.

Update note (2026-06-29): The builder is the documented supported construction
path for ordinary users. Builder-created datasets record construction
provenance that identifies the construction method, table identities, and
processing-state establishment. Direct dataset construction remains
advanced/trusted use for callers who already own fully prepared analysis-ready
tables.

Update note (2026-07-15, workflow-derived quantitative ownership): The builder
remains the ordinary public route for creating source analysis-ready datasets.
Workflow-derived quantitative datasets, such as technical-replicate aggregated
matrices, are not builder outputs and must not reuse builder preprocessing
reports or source builder provenance. They are internal derived dataset objects
with fresh derived-data provenance and explicit parent lineage.

## Context and Problem Statement

Earlier ADRs established a clear direction:

- workflows should accept only `AnalysisReadyPhosphoDataset`
- datasets should enforce a strict analysis-ready boundary
- `site_sequence` is mandatory in the final dataset contract
- transformation state must be established through PhosPy
- validation remains private and belongs to its own internal domain

At the same time, real phosphoproteomics inputs are often inconsistent and awkward. Column names vary, metadata is incomplete, site descriptions are inconsistent, and different sources structure phospho and total data differently. A strict dataset boundary is the right workflow contract, but it cannot be the first thing users have to construct by hand in every case.

PhosPy therefore needs a builder architecture that absorbs input messiness before the dataset boundary while keeping the public workflow surface simple.

## Decision Drivers

The decision is driven by the following considerations:

1. **User practicality.** Users should not need to manually normalise awkward input files before they can use the package.
2. **Strict workflow boundary.** Flexibility must stop before the dataset enters workflow execution.
3. **PhosR alignment.** The package should feel like a workflow-oriented scientific tool, not a low-level file-munging framework.
4. **Maintainability.** Ingestion complexity should be centralised rather than spread across workflows and models.
5. **Validation discipline.** Flexible ingestion must still converge on one strict validated dataset model.
6. **Extensibility.** New input conventions and data sources should be addable without reshaping workflow contracts.

## Decision

PhosPy will provide a dedicated builder-oriented path for producing `AnalysisReadyPhosphoDataset` from flexible user inputs.

This builder path is responsible for:

- accepting user-friendly raw or semi-structured inputs
- normalising column naming and input conventions
- shaping site and sample metadata
- preserving valid `site_sequence` values and deriving missing values before
  final dataset construction when supported
- invoking transformation handling through the supported transformer path
- composing shared validation-domain components
- returning a validated, missing-value-free `AnalysisReadyPhosphoDataset` in the supported public lane

The builder path should be the main public route from messy inputs to the strict dataset boundary.

## Core Design Principle

Input flexibility belongs at the builder boundary, not at the workflow boundary.

The builder exists to absorb variation. The dataset exists to represent a stable validated state. Workflows exist to operate on that state.

Each boundary must stay honest.

## Public Builder Direction

The public direction should favour a builder path over forcing users to instantiate `AnalysisReadyPhosphoDataset` manually from raw industry inputs.

Direct dataset construction may still exist and remain valid for callers who already have fully analysis-ready data, but the recommended public path should be the builder.

The public story should remain consistent, but internally the builder layer should be implemented as a small family of specialised builders or builder collaborators rather than one oversized builder object. This supports single responsibility without fragmenting the public experience.

## Builder Responsibilities

The builder path is responsible for converting flexible inputs into the final dataset contract.

This includes:

- accepting phospho input
- accepting optional total input
- accepting optional site and sample metadata
- normalising supported input column naming conventions
- shaping site metadata into the required public fields
- deriving `site_sequence` before final dataset construction when possible, and
  failing clearly when it cannot be resolved
- establishing transformation state through the supported transformer path
- composing shared validation-domain components
- constructing the final `AnalysisReadyPhosphoDataset`

The builder should not pass partially shaped raw structures into workflows.

## Input Acceptance Direction

The builder should accept either:

- already-loaded `DataFrame` inputs
- file-path inputs

This keeps the builder practical for both programmatic and file-driven usage.

File-reading support should exist only to the extent that it supports the builder's primary job of constructing the analysis-ready dataset. It should not turn the builder into a general-purpose file-ingestion framework.

## Input Flexibility Policy

The builder should intentionally tolerate reasonable input variation.

Examples include:

- varying column names for gene symbols, site identifiers, or sample identifiers
- different but recognisable site metadata conventions
- optional input separation between phospho and metadata tables
- presence or absence of total abundance input

This flexibility is a builder concern, not a workflow concern.

The builder should absorb these differences and converge on the standard public dataset contract.

However, the builder should stay simple and boring. It should not try to perform aggressive magic recognition of ambiguous conventions, because that increases the risk of bad data being accepted under incorrect assumptions.

## Column Naming Strategy

The builder architecture should support flexible column naming without forcing users to rename everything manually.

The preferred direction is a hybrid strategy:

- default support for a small set of known/common naming conventions
- internal mapping and normalisation logic within the builder path
- narrow explicit user overrides only where genuinely necessary

The public design should avoid making users specify many column-name arguments in normal cases.

Workflows should never accept repeated public arguments such as:

- `gene_col`
- `site_col`
- `sequence_col`

Those concerns belong entirely to the builder boundary.

## Unsupported Input Strategy

Unsupported input conventions should fail quickly.

The builder should prefer early, explicit failure with good messaging over trying to adapt unsupported structures through guesswork.

Failure messages should explain:

- what input shape or convention was not supported
- what was expected instead
- what narrow override or preparation step, if any, would allow the build to proceed

## Builder Output Contract

The builder returns exactly one primary public product:

- `AnalysisReadyPhosphoDataset`

That output must already satisfy the dataset boundary defined in ADR-0003.

In particular, the final dataset must already contain:

- required site metadata fields
- required `site_sequence` (validated as non-empty strings)
- established transformation state
- validated alignment across its components

The builder should not expose a second semi-ready dataset form as part of the public contract.

Diagnostics or separate build-report outputs are not a current concern and are outside the initial contract.

## Relationship to Trusted Reconstruction

Direct construction of `AnalysisReadyPhosphoDataset(...)` is sealed and raises
immediately. Trusted advanced/internal callers who already possess fully
analysis-ready data must use
`AnalysisReadyPhosphoDataset.from_trusted_tables(...)` with complete
`TrustedDatasetConstructionAssertions`. Under ADR-0024, fully analysis-ready
means `site_key` indexes plus the required auditable protein context metadata,
not display-indexed `GENE;SITE;` rows.

However:

- the builder should be the recommended public path for real-world ingestion
- trusted reconstruction should not become the expected path for messy industry
  inputs
- trusted reconstruction validates structure but cannot prove biological
  correctness of user-asserted provenance or assertions

This preserves the usefulness of the strict model without forcing all callers through a low-level manual preparation burden.

## Sequence Derivation Direction

Because `site_sequence` enrichment may still be useful, the builder path remains
the correct place to derive it when possible from supported resources.

This derivation must happen before final dataset construction.

If derivation cannot be completed from the available supported inputs and resources, the builder should fail clearly rather than constructing an incomplete analysis-ready dataset.

## Transformation Direction

The builder path must establish transformation state through the supported
transformer path defined by ADR-0006.

Transformation should not be exposed as a prominent public ingestion choice. It should follow preprocessing policy and keep the public surface small.

The builder must not treat a loose incoming transform label as equivalent to established PhosPy transformation state.

## Validation Composition Direction

The builder path should compose validation-domain components rather than reimplementing validation ad hoc.

This includes shared checks for:

- matrix structure
- metadata alignment
- required metadata fields
- derived sequence validity
- transformation-state validity

The builder may still contain orchestration-specific checks related to ingestion flow, but it should not duplicate the reusable rules that belong in the validation domain.

## Service Composition Direction

The builder path is expected to coordinate several internal services.

A likely healthy composition includes:

- input normalisation services
- metadata shaping services
- sequence derivation services
- transformer component
- validation-domain components
- final dataset construction

The exact service names may vary, but the architecture should stay boring and explicit.

## Public Surface Philosophy

The builder should not become a second application framework.

It should provide a simple, stable entry path that hides ingestion mess without hiding the final dataset contract.

The builder public surface should therefore remain smaller than the total number of internal services it coordinates.

## Consequences

### Positive Consequences

- Users get a realistic path from awkward real-world inputs to analysis-ready datasets.
- Workflow contracts stay strict and simple.
- Ingestion complexity is centralised where it belongs.
- Column-name and metadata variation stop leaking upward into workflows.
- The dataset boundary remains honest while still being practical.

### Negative Consequences

- The builder layer becomes an important internal coordination point that must be maintained carefully.
- The package must support a meaningful amount of input normalisation logic.
- Some ingestion edge cases may still require explicit user guidance or failure rather than silent adaptation.

### Neutral Consequences

- Direct dataset construction remains available for already-prepared data.
- Internal ingestion services may evolve without changing the workflow contract, provided they still converge on the same dataset boundary.

## Rejected Alternatives

### Alternative 1: Require Users to Construct `AnalysisReadyPhosphoDataset` Manually From All Raw Inputs

This option was rejected because it puts too much ingestion burden on users and is unrealistic for typical phosphoproteomics inputs.

### Alternative 2: Let Workflows Accept Raw Inputs Directly and Perform Ingestion Internally

This option was rejected because it weakens workflow contracts and reintroduces ingestion complexity into the workflow layer.

### Alternative 3: Support Multiple Semi-Ready Dataset Forms as Public Outputs

This option was rejected because it would weaken the meaning of the analysis-ready boundary and increase the public API burden.

### Alternative 4: Expose Extensive Builder Configuration for Every Ingestion Detail

This option was rejected because it would turn the public builder into a low-level configuration framework rather than a practical user entry point.

## Resolved Decisions

The following decisions are now resolved for this ADR.

1. The public builder story should be consistent, but internally it should use a small family of specialised builders or builder collaborators.
2. The builder should accept either file-path input or already-loaded `DataFrame` input.
3. The builder should remain simple and boring rather than aggressively auto-recognising ambiguous conventions.
4. The builder should return only `AnalysisReadyPhosphoDataset` in the initial contract.
5. Unsupported input conventions should fail quickly with clear and actionable messaging.

## Implementation Guidance

A likely healthy split is:

- public builder path as the recommended ingestion entry point
- internal specialised builders or builder collaborators
- internal services for input normalisation, shaping, sequence derivation, and transformation
- validation-domain reuse for shared checks
- final construction of `AnalysisReadyPhosphoDataset` only after the builder has established the required guarantees

Reviewers should reject changes that move ingestion flexibility into workflows or that force users to manually reproduce builder responsibilities in normal cases.

## Scope Boundaries

This ADR defines the analysis-ready dataset builder architecture only.

It does not define:

- the exact dataset contract beyond what ADR-0003 already establishes
- the full file-IO strategy for the package
- reference resolution strategy
- workflow result design
- migration strategy from current code

Those concerns should be addressed separately.

## Validation and Review Criteria

Future code and review work should check future changes against the following questions:

1. Does this keep ingestion flexibility below the dataset boundary?
2. Does this reduce or increase manual user burden for messy input formats?
3. Does this preserve one strict final dataset contract?
4. Does this reuse shared validation rather than duplicate it?
5. Does this keep workflows free from input-normalisation concerns?

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
- ADR-0008 defines how messy inputs are converted into
  `AnalysisReadyPhosphoDataset`.
- ADR-0018 defines the phosphosite identity/localisation policy that the
  builder must satisfy at the analysis-ready boundary.

Together, these ADRs establish:

- one public dataset model
- three public workflows
- one strict dataset boundary
- one flexible builder path below that boundary
- one consistent internal workflow pattern
- one private validation domain

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub. https://github.com/PYangLab/PhosR
