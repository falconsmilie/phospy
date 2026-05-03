# ADR: Result Model Design for PhosPy

## Document Control

- **ADR ID:** ADR-005
- **Title:** Result Model Design for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines how public and internal result models should be designed in PhosPy. The package is being developed as a maintainable Python port of PhosR. To support that goal, result objects must remain honest, small, and aligned with the real scientific stages of the workflows.

The decision is that result objects should be primarily typed data containers. Public workflow results should expose nested stage results directly rather than duplicating them through mirrored convenience accessors. Export logic, plotting logic, compatibility façades, and broad ownership-management APIs should remain outside the core result contract unless a very strong product reason emerges.

## Status

Accepted.

This ADR defines the result-model rules that support the public API, workflow architecture, dataset boundary, and reference-resolution decisions established by earlier ADRs.

## Context and Problem Statement

Result models are a common place for architectural drift. Once a workflow produces data, there is often a temptation to make the result object "helpful" by adding:

- duplicated top-level aliases for nested stage outputs
- export helpers
- plotting adapters
- compatibility layers
- ownership and mutability variants for every field
- multiple alternate views over the same underlying data

These additions usually begin as convenience, but they quickly make result objects larger, harder to document, and harder to maintain. They also obscure the true scientific structure of the workflow by training users to depend on shortcuts instead of understanding what the workflow actually produced.

PhosPy is intended to expose one analysis-ready dataset model and two primary workflows. Its result models should reinforce that workflow story rather than becoming secondary service layers.

## Decision Drivers

The decision is driven by the following considerations:

1. **Honest workflow modelling.** Results should reflect the real stages of analysis, not conceal them behind redundant aliases.
2. **Maintainability.** Smaller result models are easier to evolve, test, and explain.
3. **Documentation quality.** The public result structure should be easy to describe without enumerating many near-duplicate access paths.
4. **Separation of concerns.** Data containers should remain separate from export, plotting, and orchestration behaviour.
5. **PhosR alignment.** Workflow outputs should map clearly onto the scientific stages users expect.
6. **Freedom from backwards-compatibility burden.** Historical result shapes and convenience mirrors do not need to be preserved.

## Proposed Decision

Public result models in PhosPy will be designed as typed data containers that expose the real workflow structure directly.

Top-level workflow result objects may contain nested stage-result objects, but they should not mirror those nested fields through duplicated convenience accessors unless a shortcut is clearly essential to the product experience.

Export logic, plotting logic, compatibility adapters, and large ownership-management APIs should remain outside the core result contract.

## Result Model Principles

### Principle 1: Results Are Data Containers First

A result object should primarily answer the question:

- what did this workflow or stage produce?

It should not try to become a general-purpose service object.

### Principle 2: Nested Structure Should Be Visible

If a workflow has meaningful internal stages, the result should expose those stages honestly.

For example, a kinase workflow result should expose distinct scoring, prediction, and activity results rather than flattening all of them into one large surface.

### Principle 3: Duplication Is Discouraged

If a value already exists in a nested stage result, the top-level workflow result should not expose another convenience property for the same value unless the shortcut is demonstrably central to the user experience.

### Principle 4: Results Should Not Own Unrelated Behaviour

Result objects should not become the place where the package accumulates:

- CSV export
- plotting conversion
- compatibility bridges
- rich formatting logic
- resource discovery
- implicit recalculation

Those behaviours belong in separate services or adapters.

### Principle 5: Mutability Must Be Explicit and Minimal

Result models should not grow broad families of near-duplicate accessors purely to manage copy semantics. The design should prefer clear and limited access patterns over large `to_*` / `to_owned_*` / `to_mutable_*_unsafe` surfaces.

## Public Workflow Results

### Kinase Workflow Result

The proposed public `KinaseWorkflowResult` should contain:

- `dataset`
- `references`
- `scoring_result`
- `prediction_result`
- optional `activity_result`

This is the intended top-level shape.

The nested structure is deliberate. Users should navigate stage outputs through the stage results themselves.

Examples of intended usage:

- `result.scoring_result.profile_scores`
- `result.prediction_result.pred_mat`
- `result.activity_result.activity_scores` (`weighted_activity` compatibility alias)

The top-level result should not mirror these through repeated aliases.

### Signalome Workflow Result

The proposed public `SignalomeWorkflowResult` should contain:

- `dataset`
- `kinase_result`
- `module_assignments`
- `signalome_modules`
- `kinase_network`
- optional `expanded_signalome`

This keeps the signalome workflow tied to the coherent upstream analysis lineage instead of reconstructing that lineage through separate convenience state.

## Stage Result Models

Stage result models are encouraged where they represent real scientific boundaries.

Examples include:

- `KinaseScoringResult`
- `KinasePredictionResult`
- `KinaseActivityResult`

These stage results should be typed and focused.

### `KinaseScoringResult`

This result should hold scoring-stage outputs such as:

- `profile_scores`
- optional `motif_scores`
- optional `rank_weighted_fusion_scores`
- optional `weights`

### `KinasePredictionResult`

This result should hold prediction-stage outputs such as:

- `pred_mat`
- `substrate_list`

### `KinaseActivityResult`

This result should hold activity-stage outputs such as:

- primary `activity_scores` matrix (`weighted_activity` compatibility alias)
- `thresholded_substrate_mean_activity`
- `thresholded_substrate_counts`
- `target_counts`
- `target_table`

These stage models should remain narrow and should not be wrapped repeatedly by additional façade types without a strong reason.

## Convenience Accessor Policy

Convenience accessors are not banned, but they are tightly constrained.

A convenience accessor is justified only if all of the following are true:

1. it represents a central product concept rather than an incidental implementation detail
2. it materially improves usability
3. it does not create ambiguity about the true source of the data
4. it does not multiply maintenance burden across similar fields

Top-level convenience properties for prediction output should be avoided. Prediction output should be accessed through `prediction_result` directly.

By contrast, stage-field mirrors such as top-level aliases for:

- `profile_scores`
- `rank_weighted_fusion_scores`
- `weights`
- `substrate_list`

should be avoided, because they simply duplicate data already available through nested stage results.

## Export and Visualisation Policy

Export and visualisation concerns are outside the core result contract.

Result models should not directly own:

- CSV export methods
- plotting methods
- map-conversion helpers
- network-rendering helpers
- format-conversion services

If those capabilities are needed, they should be implemented through separate services or adapters that consume result objects.

This keeps result models smaller and prevents scientific outputs from becoming mixed with presentation concerns.

## Ownership and Copy-Semantics Policy

Result models should not expose a wide public matrix of ownership accessors by default.

The preferred direction is:

- typed result fields are exposed directly
- any specialised copying or mutable escape hatches remain tightly controlled
- ownership semantics should be documented clearly where necessary, but should not dominate the public result surface

The project should avoid result APIs that grow into large families of:

- safe copy accessors
- owned-state accessors
- mutable unsafe accessors

unless there is a compelling and proven need.

## Internal Result Models

Internal execution stages may still use additional DTOs where needed, such as:

- interpreted execution result DTOs
- intermediate scoring/prediction outputs
- provider-specific internal data structures

However, these should not leak into the public contract unless intentionally promoted.

The public result story must remain smaller than the internal implementation.

## Consequences

### Positive Consequences

- Result objects stay smaller and easier to understand.
- The workflow structure remains visible and honest.
- Documentation becomes easier because there are fewer duplicate access paths.
- Export and plotting concerns stay decoupled from scientific data modelling.
- Maintenance burden is reduced because there are fewer mirrored fields and fewer façade layers.

### Negative Consequences

- Some users may need to navigate nested result objects rather than relying on top-level shortcuts.
- Convenience APIs will be more constrained than in a wrapper-heavy design.
- Separate adapter or exporter components may need to be introduced where presentation or serialisation support is required.

### Neutral Consequences

- Internal pipelines may still use richer temporary result DTOs where helpful.
- High-value shortcuts may still exist in carefully chosen cases.

## Rejected Alternatives

### Alternative 1: Flatten All Stage Outputs Onto the Top-Level Workflow Result

This option was rejected because it obscures the workflow structure, creates duplicate access paths, and encourages result-surface bloat.

### Alternative 2: Let Result Objects Accumulate Export and Plotting Behaviour

This option was rejected because it mixes scientific data modelling with presentation and transport concerns.

### Alternative 3: Preserve Historical Result Aliases and Compatibility Wrappers

This option was rejected because the project is not constrained by a backwards-compatibility goal and should not preserve architectural drift by default.

### Alternative 4: Introduce Broad Ownership and Mutability APIs on All Result Types

This option was rejected because it inflates the public surface and shifts result-model design toward defensive utility patterns rather than honest workflow modelling.

## Resolved Decisions

The following decisions are now resolved for this ADR.

1. `pred_mat_result` should not be retained as a top-level convenience property on `KinaseWorkflowResult`. Prediction output should be accessed through `prediction_result` directly.
2. `SignalomeWorkflowResult` should expose `kinase_result` directly.
3. Result-level provenance attachment is not a current concern and should remain outside result models for now.
4. The `dataset` field should be included on all public workflow result objects for now.

## Implementation Guidance

A likely healthy split is:

- public workflow result as a narrow typed container
- nested stage-result models for real stage boundaries
- separate exporter or adapter services when presentation concerns arise
- internal DTOs where execution stages need richer internal structure

Reviewers should reject changes that add mirrored aliases, export helpers, plotting helpers, or ownership-management sprawl to public result models without a strong and explicit justification.

## Scope Boundaries

This ADR defines result-model design only.

It does not define:

- the exact dataset contract
- reference resolution strategy
- transformer design
- exporter or visualisation service APIs
- migration strategy from current code

Those concerns should be addressed separately.

## Validation and Review Criteria

Future code and review work should check proposed changes against the following questions:

1. Does this result model expose the real workflow structure honestly?
2. Does this change add a real result concept, or just a duplicate access path?
3. Does this keep export and plotting behaviour out of core result models?
4. Does this keep convenience accessors rare and justified?
5. Does this make the public result contract clearer or more confusing?

If the answers are weak or negative, the design should be reconsidered.

## Relationship to Earlier ADRs

This ADR complements the earlier architecture decisions.

- ADR-001 defines the intended public API contract.
- ADR-002 defines the internal workflow architecture.
- ADR-003 defines the dataset and preprocessing boundary.
- ADR-004 defines the reference resolution strategy and `ReferenceBundle` contract.
- ADR-005 defines how workflow and stage results should be modelled.

Together, these ADRs establish:

- one public dataset model
- two public workflows
- one consistent internal workflow pattern
- one explicit reference-resolution path
- one disciplined approach to result-model design

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR*. GitHub repository. https://github.com/PYangLab/PhosR
