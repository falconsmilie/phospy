# ADR: Public API Contract for PhosPy

## Document Control

- **ADR ID:** ADR-001
- **Title:** Public API Contract for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Authors:** OpenAI ChatGPT with project direction from the PhosPy maintainer
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines the intended public API contract for PhosPy. The package is being developed as a maintainable Python port of PhosR. To support that goal, the public surface must remain small, workflow-oriented, and aligned with the conceptual model of phosphoproteomics analysis rather than shaped by historical implementation details.

The decision is to expose one analysis-ready dataset model and two primary workflows: a kinase workflow and a signalome workflow. Each workflow accepts a single typed request object and returns a single typed result object. Supporting configuration models are included only where they clarify the workflow contract.

## Status

Accepted.

This ADR establishes the target public API contract for the ongoing architecture reset and rewrite effort.

## Context and Problem Statement

PhosPy has accumulated internal structures that are more complex than necessary for the product it is intended to be. The project goal is not to build a generic phosphoproteomics framework or a highly abstract orchestration platform. The goal is to build a usable, scientifically faithful, and maintainable Python port of PhosR.

As the codebase evolved, design effort began to drift toward wrapper-heavy orchestration, duplicated result accessors, compatibility-oriented layering, and abstractions that made the software harder to understand and maintain. This produced friction in both implementation and review. It also created a risk that the public API would expand around internal architecture rather than around the scientific workflows that users actually care about.

The project therefore requires a clear decision about the public contract before internal implementation continues.

## Decision Drivers

The decision is driven by the following considerations:

1. **Conceptual fidelity to PhosR.** The public API should reflect the scientific workflow expected from a PhosR-style analysis package.
2. **Maintainability.** A small, explicit public API is easier to document, test, preserve, and reason about.
3. **Usability.** Users should be able to understand the package in terms of an analysis-ready dataset and a small number of meaningful workflows.
4. **Correctness.** Workflow contracts should make invalid usage difficult and valid usage obvious.
5. **Freedom to simplify internals.** The public API should not force preservation of historical internal structure.
6. **No backwards-compatibility burden.** The current effort is not constrained by the need to preserve older public shapes or compatibility aliases.

## Decision

PhosPy will expose a deliberately small public API centred on one dataset model and two primary workflows.

### Public API scope

The intended public surface is:

- `AnalysisReadyPhosphoDataset`
- `Organism`
- `ReferencePreset`
- `ReferenceBundle`
- `KinaseScoringConfig`
- `KinasePredictionConfig`
- `KinaseActivityConfig`
- `KinaseWorkflowRequest`
- `KinaseWorkflowResult`
- `KinaseWorkflow`
- `SignalomeConfig`
- `SignalomeWorkflowRequest`
- `SignalomeWorkflowResult`
- `SignalomeWorkflow`
- `PhosPyValidationError`

No additional public workflows, compatibility facades, or result-wrapper layers are part of the intended contract by default.

### Dataset contract

`AnalysisReadyPhosphoDataset` is the only public dataset model.

It represents validated, analysis-ready phosphoproteomics input data. It is the primary data boundary before workflow execution.

The dataset contract includes:

- `phospho`
- `site_metadata`
- optional `sample_metadata`
- optional `total`
- optional `organism`
- validated transformation state established through the supported transformer path

The dataset is expected to enforce core invariants at construction time.

Required site metadata columns are standardised and must include:

- `gene_symbol`
- `site`
- optional `site_sequence` (validated when present)

`site_sequence` is optional in the final public dataset contract. A preprocessing
path may derive it before final dataset construction, but workflows should not
assume it is present unless that specific workflow lane requires it.

This decision removes the need for repeated public column-name arguments across workflows.

### Kinase workflow contract

`KinaseWorkflow` is the primary public entry point for PhosR-style kinase analysis.

It accepts exactly one request object:

- `KinaseWorkflowRequest`

It returns exactly one result object:

- `KinaseWorkflowResult`

The request contains:

- dataset
- reference selection
- scoring configuration
- prediction configuration
- activity configuration

The result contains:

- dataset
- resolved references
- scoring result
- prediction result
- optional activity result

Nested stage results are exposed directly. The top-level result object should not mirror nested fields through duplicated convenience accessors unless a shortcut is clearly essential.

### Signalome workflow contract

`SignalomeWorkflow` is the primary public entry point for PhosR-style signalome analysis.

It accepts exactly one request object:

- `SignalomeWorkflowRequest`

It returns exactly one result object:

- `SignalomeWorkflowResult`

The request depends on the kinase workflow result so the public pipeline remains coherent, staged, and difficult to misuse.

The request contains:

- kinase workflow result
- signalome configuration

The result contains:

- dataset
- kinase workflow result
- module assignments
- signalome modules
- kinase network
- optional expanded signalome

## Dataset Construction and Validation Direction

Although `AnalysisReadyPhosphoDataset` is part of the public API contract, users should not be forced to normalise difficult phosphoproteomics inputs manually before using it.

Public construction should therefore support a flexible builder-oriented path that can absorb input variability, especially around column naming and related ingestion differences, before producing the final validated dataset.

Dataset validation remains a private concern. It is not intended to become a standalone public surface.

Transformation state must have stronger guarantees than a simple free-text scale label. The public contract should reflect a validated transformation state established through the supported transformer path, even if a derived user-facing scale label is retained for convenience.

These dataset-related decisions should be read together with ADR-003, which defines the analysis-ready dataset and preprocessing boundary in more detail.

## Rationale

This public contract matches the intended product shape closely.

First, it aligns the package with the actual scientific workflow rather than exposing implementation fragments. A user prepares an analysis-ready dataset, runs kinase analysis, and optionally runs signalome analysis downstream.

Second, it keeps the public surface small enough to be documented and maintained properly. Scientific packages tend to become difficult to use when the public API expands around helpers, wrappers, and historical aliases.

Third, it supports internal simplification. Because the contract is defined in terms of stable workflow concepts rather than internal service layers, implementation can be reorganised without disturbing the product surface.

Fourth, it encourages honest result modelling. Users should see the real nested structure of workflow results rather than being encouraged to depend on convenience mirrors that duplicate the same data through multiple access paths.

## Consequences

### Positive consequences

- The public API becomes easier to learn and document.
- The package aligns more clearly with its identity as a PhosR port.
- Internal implementation can be simplified aggressively without preserving historical architecture.
- Request and result models become stable anchors for validation and testing.
- The workflow-oriented user experience becomes clearer.

### Negative consequences

- Existing internal abstractions that were previously reflected in public shapes will be removed or hidden.
- Some convenience access patterns may be lost in favour of a smaller and more honest public model.
- Implementation work is required to make current code conform to the new contract.
- A strict public API boundary may require some users to adapt if they have been relying on unofficial or incidental access paths.

### Neutral consequences

- Internal service boundaries, validators, interpreters, executors, and providers remain implementation details unless intentionally promoted later.
- Supporting utilities may still exist internally, but they are not part of the product contract by default.

## Rejected Alternatives

### Alternative 1: Preserve the current public surface and refactor around it

This option was rejected because the current shape contains signs of architectural drift, including duplicated result accessors, wrapper-heavy orchestration, and compatibility-driven complexity. Preserving that surface would lock the rewrite into historical design compromises.

### Alternative 2: Expose multiple public workflows for different internal stages

This option was rejected because it overfits the public API to internal implementation details. The intended product is better represented by one dataset model and two main workflows.

### Alternative 3: Expose a broad utility-style API with many helper entry points

This option was rejected because it weakens the PhosR workflow model, increases documentation burden, and encourages fragmented usage patterns.

### Alternative 4: Build a highly abstract plugin-style public API from the start

This option was rejected because it would prioritise extensibility over clarity and would likely recreate the same problems of over-smart architecture that the current effort is trying to remove.

## Implementation Guidance

The public API defined in this ADR should be implemented through simple internal workflow stages:

- validator
- interpreter
- executor

These are internal implementation concerns, not public contract concerns.

Public workflows should expose a single `run(request)` method that returns a single result object.

Public result objects should be primarily data containers. They should not become service objects, export engines, or compatibility façades.

Public dataset construction should validate invariants early and consistently.

## Scope Boundaries

This ADR defines the intended public API contract only.

It does not define:

- the detailed internal package layout
- the exact implementation of validators, interpreters, or executors
- export or visualisation APIs
- migration strategy from older code
- packaging or release policy

Those concerns should be addressed in separate ADRs where required.

## Validation and Review Criteria

Future code and review work should check proposed changes against the following questions:

1. Does this change preserve the model of one dataset and two primary workflows?
2. Does this change clarify or weaken the PhosR-aligned workflow story?
3. Does this add a real public concept, or just expose an implementation detail?
4. Does this make invalid usage harder and valid usage clearer?
5. Does this keep the public API small and explainable?

If the answer to these questions is weak or negative, the change should be reconsidered.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR*. GitHub repository. https://github.com/PYangLab/PhosR
