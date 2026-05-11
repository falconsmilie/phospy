# ADR: Analysis-Ready Dataset and Preprocessing Boundary for PhosPy

## Document Control

- **ADR ID:** ADR-0003
- **Title:** Analysis-Ready Dataset and Preprocessing Boundary for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines the intended boundary between preprocessing and workflow execution in PhosPy. The package is being developed as a maintainable Python port of PhosR. To support that goal, workflows must operate on a single, well-defined, analysis-ready dataset model rather than on loosely structured matrices and metadata fragments.

The decision is to make `AnalysisReadyPhosphoDataset` the only public dataset model used as workflow input. Preprocessing is responsible for transforming user-supplied phosphoproteomics data into that validated form. Public workflows should assume they receive a valid analysis-ready dataset and should not repeat preprocessing concerns except where a workflow introduces genuinely new constraints.

## Status

Accepted.

This ADR defines the dataset and preprocessing boundary that supports the public API and internal workflow architecture established by earlier ADRs.

Update note (2026-05-11): `site_sequence` may be omitted at ingestion, but it
is mandatory at the `AnalysisReadyPhosphoDataset` boundary. Builder
preprocessing owns the derive-or-fail transition before final dataset
construction.

Update note (2026-05-11, contract correction): final dataset construction must
reject any row with missing, blank, or invalid `site_sequence` values before
workflow execution.

Update note (2026-05-11, provenance correction): preprocessing now emits a
structured `SiteSequenceResolutionReport` that audits sequence origin and loss
(`provided_by_input`, `resolved_from_reference`, `resolved_from_fasta`,
`unresolved`, `conflicts`, `conflict_policy`,
`final_sequence_complete_sites`). Sequence-loss decisions are preprocessing
responsibility and must not be deferred to workflow scoring.

## Context and Problem Statement

PhosPy is intended to expose one public dataset model and three primary
workflows (differential, kinase, signalome). All three depend on a stable and
well-understood input boundary. Without a strict dataset contract, workflows
become polluted with repeated checks, fallback logic, column-name plumbing, and
interpretation of partially prepared data.

Historically, phosphoproteomics pipelines often allow flexible input structures during import and preprocessing. That flexibility is useful at the ingestion boundary, but it becomes harmful if it leaks into workflow execution. When workflows are required to interpret raw or semi-prepared inputs, the public API becomes harder to use, internal validation becomes repetitive, and the software becomes harder to maintain.

PhosPy therefore needs an explicit decision on what "analysis-ready" means and where that state is established.

## Decision Drivers

The decision is driven by the following considerations:

1. **Workflow clarity.** Workflows should accept a stable dataset model rather than infer structure from loosely coupled input fragments.
2. **PhosR alignment.** The package should model the conceptual transition from preprocessing into downstream analysis clearly.
3. **Maintainability.** Repeated input-shaping logic across workflows increases drift and complexity.
4. **Validation quality.** Public dataset invariants should be checked once at the correct boundary.
5. **Usability.** Users should know when their data is considered analysis-ready and what that implies.
6. **Extensibility.** Additional preprocessing paths should still converge on the same dataset contract.

## Decision

`AnalysisReadyPhosphoDataset` will be the only public dataset model accepted by public workflows.

Preprocessing is responsible for producing this model.

The dataset model represents a validated, analysis-ready phosphoproteomics state that is suitable for kinase workflow execution and, through that path, signalome workflow execution.

Public workflows must not accept raw phospho tables, raw metadata tables, or loose collections of data structures as alternative first-class inputs.

## Dataset Contract

The public dataset model contains:

- `phospho`
- `site_metadata`
- optional `sample_metadata`
- optional `total`
- optional `organism`
- `intensity_scale`

### Field Meanings

#### `phospho`

A numeric phosphosite-by-sample matrix representing the phosphoproteomics data used for downstream analysis.

#### `site_metadata`

A site-level metadata table indexed exactly like `phospho` and containing the standard metadata required for downstream workflows.

#### `sample_metadata`

An optional sample-level metadata table indexed exactly like the sample columns of `phospho`.

#### `total`

An optional protein total abundance matrix for workflows or analyses that require or benefit from paired total measurements.

#### `organism`

An optional explicit organism indicator that may assist reference resolution and validation.

#### `intensity_scale`

A simple declared label is not sufficient as the long-term design. The dataset boundary should carry stronger guarantees about transformation state through a typed transformation concept produced by a transformer component. In practice, this means the dataset should represent not only a nominal scale label, but a validated transformation state established during preprocessing.

## Required Dataset Invariants

The dataset is expected to enforce the following invariants at construction time.

### Core Matrix Invariants

- `phospho` must be a non-empty numeric `DataFrame`
- `phospho.index` must be unique
- `phospho.columns` must be unique

### Site Metadata Invariants

- `site_metadata` must be a non-empty `DataFrame`
- `site_metadata.index` must exactly match `phospho.index`
- required columns must be present:
  - `gene_symbol`
  - `site`
- `site_sequence` is required and must contain non-empty plausible sequence
  strings

### Sample Metadata Invariants

If `sample_metadata` is provided:

- it must be a `DataFrame`
- its index must exactly match `phospho.columns`

### Total Matrix Invariants

If `total` is provided:

- it must be a numeric `DataFrame`
- `total.columns` must exactly match `phospho.columns`
- `total.index` must be unique

### Transformation-State Invariant

- transformation state must be established by a supported transformer path
- the dataset must expose a validated transformation state rather than relying on an informal free-text label alone
- any remaining public scale label, if retained, must be derived from that stronger transformation state

## Meaning of "Analysis-Ready"

For the purposes of PhosPy, a dataset is analysis-ready when the following are true:

1. The phospho matrix and required metadata have been structurally validated.
2. Standard site metadata fields are normalised to the required public column names.
3. Matrix and metadata alignment invariants have been established.
4. Any optional total matrix has been aligned to the phospho samples.
5. The dataset is suitable for workflow validation without requiring raw-data interpretation.

This does not necessarily mean that every possible scientific normalisation or filtering decision has already occurred. It means the dataset has crossed the public boundary from raw or semi-structured input into a stable, validated analysis object.

For the primary supported lane, `dataset.site_metadata.site_sequence` is
required at this boundary. Ingestion may remain permissive, but unresolved
sequence context must not pass through this boundary.

## Preprocessing Responsibility

Preprocessing is responsible for converting user-provided raw or semi-structured data into `AnalysisReadyPhosphoDataset`.

This includes responsibilities such as:

- column normalisation at the ingestion boundary
- validation of required preprocessing fields
- site-level metadata shaping
- sample metadata shaping
- alignment of phospho and optional total data
- scale declaration or normalisation where required by preprocessing policy
- rejection of malformed inputs before workflow execution

Preprocessing may be implemented through builders or service classes, but those details are internal architecture concerns.

## Workflow Responsibility

Once a workflow receives `AnalysisReadyPhosphoDataset`, it may assume that the dataset contract has already been established.

Workflows may still validate constraints that are specific to the workflow itself, such as:

- whether required reference inputs are compatible with the dataset
- whether a requested activity grouping column exists in `sample_metadata`
- whether a downstream signalome analysis requires outputs from an earlier stage

Workflows should not be responsible for:

- interpreting arbitrary user column names
- reconstructing required site metadata fields from raw inputs
- aligning phospho and total matrices from scratch
- deciding whether the dataset is fundamentally analysis-ready

## Standardised Metadata Columns

To prevent column-name plumbing from leaking into workflow APIs, `site_metadata` must use standard public column names.

The required columns are:

- `gene_symbol`
- `site`
- required `site_sequence`

Additional metadata columns are allowed, but workflows should depend only on documented public fields unless they introduce a separate and explicit requirement.

This decision means public workflow requests should not accept repeated arguments such as:

- `gene_col`
- `site_col`
- `sequence_col`

Those concerns belong to preprocessing, not workflow execution.

`site_sequence` may be optional at ingestion, but preprocessing must derive it
or fail before final dataset construction.

## Construction and Validation Strategy

The design is that `AnalysisReadyPhosphoDataset` validates its own public
invariants at construction time.

This keeps the public dataset boundary honest and makes invalid datasets harder to construct accidentally.

Public construction should support a flexible builder path so users are not forced to pre-normalise awkward industry input formats themselves. In particular, preprocessing and builder services should absorb variability in column naming and related ingestion quirks before final dataset construction.

Preprocessing services or builders may perform additional shaping before construction, but they should not rely on workflows to complete the basic validity story.

## Optionality Rules

The following optionality rules apply.

### Required

- `phospho`
- `site_metadata`

### Optional

- `sample_metadata`
- `total`
- `organism`

### Rationale

This keeps the minimum public dataset contract small enough for kinase analysis while allowing richer downstream use where available.

## Implications for Signalome Workflow

The signalome workflow is downstream of the kinase workflow and should not introduce a second competing dataset boundary.

Signalome analysis should rely on the same `AnalysisReadyPhosphoDataset` that entered kinase analysis, carried through the kinase workflow result.

This keeps the public analysis flow coherent:

1. prepare dataset
2. run kinase workflow
3. run signalome workflow

## Consequences

### Positive Consequences

- Workflows have a clean and stable input contract.
- Preprocessing concerns stop leaking into workflow APIs.
- Column-name plumbing is eliminated from the public workflow surface.
- Validation becomes easier to reason about and easier to test.
- The package better reflects the conceptual transition from preprocessing to analysis.

### Negative Consequences

- Raw input flexibility is pushed earlier in the stack and no longer tolerated at workflow boundaries.
- Preprocessing must take responsibility for more of the public data-shaping story.
- Some current internal code that relies on semi-prepared data may need to be rewritten.

### Neutral Consequences

- Different preprocessing implementations may still exist internally as long as they converge on the same dataset model.
- Additional metadata fields may still be carried through the dataset without becoming required public contract fields.

## Rejected Alternatives

### Alternative 1: Allow Workflows to Accept Raw Matrices and Metadata Directly

This option was rejected because it weakens the public contract, encourages duplicated validation logic, and makes workflows responsible for input shaping that belongs earlier.

### Alternative 2: Support Multiple Public Dataset Shapes

This option was rejected because it would make the workflow boundary ambiguous and increase the documentation and validation burden.

### Alternative 3: Treat Analysis Readiness as a Soft Convention Rather Than a Strict Model

This option was rejected because the package needs a reliable and testable boundary, not an informal expectation.

## Resolved Decisions and Remaining Recommendation

The following decisions are now resolved for this ADR.

1. `site_sequence` is mandatory in the final `AnalysisReadyPhosphoDataset` contract.
2. A preprocessing path may accept missing `site_sequence` at ingestion, but it must derive it or fail before final dataset construction.
3. A simple declared `intensity_scale` label is not sufficient on its own; the dataset boundary should carry stronger guarantees about transformation state through a transformer-oriented design.
4. Public dataset construction should support a flexible builder path so users do not have to manually normalise difficult phosphoproteomics input formats or column naming schemes.
5. Dataset validation should remain private and should not become a standalone public surface.

The remaining design recommendation is about provenance.

### Recommendation on Provenance

Preprocessing provenance should remain outside the core public dataset contract.

The core dataset model should stay focused on validated analysis-ready state. Provenance is valuable, but it is better treated as optional accompanying metadata rather than a required part of the fundamental dataset identity.

This keeps the dataset model smaller and reduces the risk that operational or ingestion-history concerns dominate the core public contract.

If provenance is needed later, preferred options include:

- an optional metadata attachment on workflow results
- a separate provenance record produced by builders or preprocessing services
- internal audit metadata that does not reshape the core dataset model

## Implementation Guidance

The dataset model should remain a data container with validation, not a workflow service.

A likely healthy split is:

- `AnalysisReadyPhosphoDataset` as the validated public model
- preprocessing or builder services as internal shaping logic
- a dedicated validation domain for reusable lower-level validation concerns
- workflow-level validators as composers that call shared validation components and add workflow-specific rules

Reviewers should reject designs that attempt to move ingestion flexibility back into workflow requests.

## Scope Boundaries

This ADR defines the intended dataset and preprocessing boundary only.

It does not define:

- the exact internal structure of preprocessing services
- reference resolution strategy
- workflow result modelling
- export or visualisation APIs
- migration strategy from current code

Those concerns should be addressed separately.

## Validation and Review Criteria

Future code and review work should check future changes against the following questions:

1. Does this change preserve a single public analysis-ready dataset model?
2. Does this move raw-input interpretation into the correct boundary, or does it leak into workflows?
3. Does this strengthen or weaken dataset invariants?
4. Does this reduce or reintroduce column-name plumbing in public APIs?
5. Does this make workflow inputs clearer and harder to misuse?

If the answers are weak or negative, the design should be reconsidered.

## Relationship to Earlier ADRs

This ADR complements the earlier architecture decisions.

- ADR-0001 defines the intended public API contract.
- ADR-0002 defines the internal workflow architecture.
- ADR-0003 defines the dataset and preprocessing boundary that those workflows
  depend on.
- ADR-0018 defines phosphosite identity/localisation policy at the same
  analysis-ready boundary.

Together, these ADRs establish:

- one public dataset model
- three public workflows
- one consistent internal workflow pattern

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub. https://github.com/PYangLab/PhosR
