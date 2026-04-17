# ADR: Reference Resolution Strategy and `ReferenceBundle` Contract for PhosPy

## Document Control

- **ADR ID:** ADR-004
- **Title:** Reference Resolution Strategy and `ReferenceBundle` Contract for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Authors:** OpenAI ChatGPT with project direction from the PhosPy maintainer
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines how reference resources are represented, resolved, and validated in PhosPy. The package is being developed as a maintainable Python port of PhosR. To support that goal, the kinase workflow must operate against a predictable reference contract rather than a mix of implicit defaults, loosely structured resources, and workflow-local resolution logic.

The decision is to define a single public `ReferenceBundle` contract and a simple reference resolution strategy based on either an explicit bundle or a reference preset. Public workflows should not assemble references ad hoc. Instead, they should delegate reference resolution to a dedicated reference provider path during interpretation, with shared validation handled in the validation domain.

## Status

Accepted.

This ADR defines the reference boundary that supports the public API, dataset boundary, and internal workflow architecture established by earlier ADRs.

## Context and Problem Statement

The kinase workflow depends on reference resources such as kinase-substrate relationships and sequence-oriented assets used in scoring and prediction. These resources are central to the scientific behaviour of the package, but they are also a source of ambiguity if their resolution is not explicitly defined.

Without a clear reference strategy, the workflow boundary becomes hard to reason about. Questions such as the following remain underspecified:

- What does `ReferencePreset.AUTO` actually do?
- How should dataset organism influence reference choice?
- When should bundled references be used?
- When should user-supplied references be accepted?
- What is the minimum valid content of `ReferenceBundle`?
- Where should reference compatibility be validated?

If these decisions remain implicit, the interpreter stage becomes vague, error messages become inconsistent, and workflow behaviour becomes harder to document and maintain.

## Decision Drivers

The decision is driven by the following considerations:

1. **Predictability.** Users should be able to understand how PhosPy chooses reference resources.
2. **PhosR alignment.** Reference handling should support the expected kinase-analysis workflow without making users manage unnecessary internal detail.
3. **Maintainability.** Reference resolution logic should not be duplicated across workflows or hidden inside unrelated components.
4. **Validation quality.** Reference compatibility and completeness should be checked consistently.
5. **Extensibility.** The design should allow bundled references and user-supplied references without fragmenting the workflow contract.
6. **Clarity of responsibility.** Reference resolution belongs in interpretation; shared reference validity rules belong in the validation domain.

## Proposed Decision

PhosPy will support two public reference input forms for the kinase workflow:

1. a `ReferencePreset`
2. an explicit `ReferenceBundle`

`ReferencePreset.AUTO` remains a valid public option, but it must resolve through explicit, documented rules rather than implicit workflow-local behaviour.

Reference resolution will occur during the interpreter stage, through a dedicated reference provider path.

Reference validation will remain private and will live in the validation domain. Workflow-level validators will compose shared reference validation with workflow-specific compatibility checks.

## Public Reference Inputs

### `ReferencePreset`

`ReferencePreset` is a public enum-like input that allows users to request built-in references without supplying a full bundle manually.

The proposed public values are:

- `AUTO`
- `HUMAN`
- `MOUSE`
- `RAT`

These values define organism lanes for reference resolution. A given release may
bundle only a subset of those lanes (currently rat only).

Additional built-in organisms are future work.

### `ReferenceBundle`

`ReferenceBundle` is the only public structured reference model.

It represents a fully materialised reference package suitable for kinase scoring and prediction.

A workflow that receives a `ReferenceBundle` should not need to perform further discovery or resource assembly.

## Proposed `ReferenceBundle` Contract

The proposed public `ReferenceBundle` contains:

- `organism`
- `kinase_substrate_map`
- `site_sequences`

These are the only required public fields in the initial contract. Additional internal reference structures may exist, but they should be derived behind the contract rather than required as separate first-class public inputs.

### Field meanings

#### `organism`

A required organism label describing the biological scope of the bundle.

#### `kinase_substrate_map`

A non-empty structured table representing the kinase-substrate relationships required for scoring and related downstream logic.

#### `site_sequences`

A non-empty structured sequence resource required for the motif-aware scoring path.

## Bundle Invariants

The proposed `ReferenceBundle` invariants are:

- `kinase_substrate_map` must be non-empty
- `site_sequences` must be non-empty
- `site_sequences.index` must be unique
- `organism` is required and must be a supported organism enum value

Additional internal reference structures may exist behind this contract, but they should not be required as separate public inputs.

## Resolution Strategy

Reference resolution should follow a simple ordered strategy.

### Rule 1: Explicit `ReferenceBundle` wins

If the user supplies a `ReferenceBundle`, that bundle is used directly, subject to validation.

The workflow must not silently override it with bundled defaults.

### Rule 2: Explicit preset resolves to bundled references when available

If the user supplies a specific preset such as `HUMAN`, `MOUSE`, or `RAT`, the
interpreter resolves that preset to a bundled `ReferenceBundle` only when that
organism has packaged bundled support in the current release. Otherwise it fails
explicitly.

### Rule 3: `AUTO` resolves from dataset organism

If the user supplies `ReferencePreset.AUTO`, resolution should use the dataset organism.

The proposed rule is:

- if `dataset.organism` is present and has bundled support in the current
  release, resolve to the matching bundled reference bundle
- otherwise fail with a validation error rather than guessing

`AUTO` should not infer organism from fragile heuristics or partially trusted metadata.

### Rule 4: Unsupported organism fails explicitly

If organism is missing, unsupported, or incompatible with the requested preset, PhosPy should fail with a clear validation error.

## Organism Compatibility Rules

The following compatibility rules are proposed.

### Explicit bundle with organism

If a user-supplied `ReferenceBundle` organism and the dataset organism are both present, the two must be checked for compatibility.

A mismatch must fail clearly.

This is a scientific safety rule, not just a software preference. A dataset-organism and reference-organism mismatch can imply wrong sequence context, wrong substrate mappings, fragile gene-symbol assumptions, or unspoken cross-species translation assumptions. PhosPy should not silently proceed in those cases.

If cross-species support is needed in the future, it should be handled through an explicit adaptation or reference-building step rather than a permissive mismatch override inside the workflow.

### Explicit preset

If a user requests a specific preset and the dataset organism is present but does not match, the request must fail clearly.

PhosPy should not silently prioritise one source of truth over the other.

### `AUTO`

If `AUTO` is used, the dataset organism becomes the source of truth for bundled reference selection.

## Responsibility Boundaries

### Dataset responsibility

The dataset may declare `organism`, but it does not resolve references.

### Workflow request responsibility

The public request declares either a preset or an explicit bundle, but it does not perform resolution.

### Validator responsibility

The validation domain is responsible for reusable checks such as:

- whether a bundle is structurally valid
- whether organism values are supported
- whether dataset and reference organism values are compatible

Workflow-level validators are responsible for composing those checks with workflow-specific requirements.

### Interpreter responsibility

The interpreter is responsible for turning the request reference input into a concrete `ReferenceBundle`.

This is the stage that should call the provider path and resolve `AUTO`.

### Executor responsibility

The executor consumes a fully resolved `ReferenceBundle` and should not perform reference discovery or fallback resolution.

## Reference Provider Direction

The proposed internal direction is to introduce a dedicated reference provider seam.

A likely shape is:

- a private reference provider interface
- one bundled-reference provider implementation for the initial built-in organisms
- optional future provider implementations for additional organisms or external APIs
- organism-oriented contracts at the provider boundary so support can grow without pushing ad hoc organism logic into workflows

The provider should expose a uniform `run(...)` method in line with the project naming convention for stage and component methods.

The reference provider is not intended to be a broad public plugin system at this stage. It is an internal seam that supports interpretation cleanly.

## Error Strategy

Reference-related failures should be explicit and actionable.

Examples include:

- missing dataset organism when `AUTO` is requested
- unsupported organism value
- empty or malformed `ReferenceBundle`
- dataset and bundle organism mismatch
- dataset and preset mismatch

Errors should explain what was requested, what was found, and why resolution could not proceed.

## Consequences

### Positive consequences

- Reference resolution becomes predictable and easy to document.
- The kinase workflow remains simple at the public boundary.
- `AUTO` gains explicit meaning instead of behaving like a fuzzy convenience.
- Custom references remain possible without expanding the public workflow surface.
- Reference logic is kept out of the executor.

### Negative consequences

- Users must supply organism information when relying on `AUTO`.
- Unsupported organism scenarios fail earlier and more explicitly.
- Reference resolution logic must be built and maintained centrally rather than ad hoc.

### Neutral consequences

- Internal reference representations may still be richer than the public `ReferenceBundle` contract.
- Additional bundled reference sets may be introduced later without changing the workflow contract, as long as they fit the same resolution rules.

## Rejected Alternatives

### Alternative 1: Resolve references directly inside the workflow executor

This option was rejected because it mixes interpretation with execution and makes the executor responsible for decisions that should be settled earlier.

### Alternative 2: Allow raw reference fragments as first-class workflow inputs

This option was rejected because it weakens the public contract and shifts assembly burden onto users and workflow logic.

### Alternative 3: Make `AUTO` infer organism from heuristics

This option was rejected because implicit inference would be brittle, hard to explain, and likely to create silent scientific mistakes.

### Alternative 4: Expose the reference provider as a broad public plugin system immediately

This option was rejected because it would introduce extension-oriented complexity before the core workflow contract is fully stabilised.

## Resolved Decisions

The following decisions are now resolved for this ADR.

1. `ReferenceBundle.organism` is required for all public bundles.
2. The required public `ReferenceBundle` fields remain:
   - `organism`
   - `kinase_substrate_map`
   - `site_sequences`
3. Additional public built-in organisms beyond human, mouse, and rat are future work.
4. Future organism expansion should come through provider interfaces and organism-oriented contracts rather than ad hoc workflow conditionals.
5. Dataset-organism and reference-organism mismatch must always fail in the public workflow path.
6. Public signalome analysis should not accept direct reference input. It should rely on the resolved references already carried through the kinase workflow result.

## Implementation Guidance

A likely healthy split is:

- public request declares preset or explicit bundle
- workflow validator composes shared reference validation
- interpreter resolves to a concrete `ReferenceBundle`
- executor consumes only fully resolved references
- bundled resource loading stays behind the reference provider path

Reviewers should reject changes that push reference assembly, organism guessing, bundled-resource fallback logic, or mismatch-tolerance policies into the workflow executor.

## Scope Boundaries

This ADR defines the reference resolution strategy and public `ReferenceBundle` contract only.

It does not define:

- the full internal layout of the references package
- the full algorithmic details of scoring and prediction
- result-model design
- export or visualisation behaviour
- migration strategy from current code

Those concerns should be addressed separately.

## Validation and Review Criteria

Future code and review work should check proposed changes against the following questions:

1. Does this preserve a single clear path from reference input to concrete `ReferenceBundle`?
2. Does this keep `AUTO` explicit rather than heuristic?
3. Does this keep reference validation in the validation domain?
4. Does this keep reference discovery out of the executor?
5. Does this make bundled and custom references easier to reason about rather than harder?

If the answers are weak or negative, the design should be reconsidered.

## Relationship to Earlier ADRs

This ADR complements the earlier architecture decisions.

- ADR-001 defines the intended public API contract.
- ADR-002 defines the internal workflow architecture.
- ADR-003 defines the dataset and preprocessing boundary.
- ADR-004 defines how references are selected, resolved, and validated for workflow use.

Together, these ADRs establish:

- one public dataset model
- two public workflows
- one consistent internal workflow pattern
- one explicit path from reference request to resolved bundle

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR*. GitHub repository. https://github.com/PYangLab/PhosR
