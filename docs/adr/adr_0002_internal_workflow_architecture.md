# ADR: Internal Workflow Architecture for PhosPy

## Document Control

- **ADR ID:** ADR-0002
- **Title:** Internal Workflow Architecture for PhosPy
- **Status:** Accepted
- **Date:** 2026-04-16
- **Decision Type:** Architecture Decision Record

## Abstract

This Architecture Decision Record defines the intended internal workflow architecture for PhosPy. The package is being developed as a maintainable Python port of PhosR. To support that goal, internal workflow orchestration must remain simple, explicit, and easy to reason about.

The decision is to implement each public workflow through three internal stages:

- validator
- interpreter
- executor

These stages create a predictable path from user request to result and provide disciplined boundaries for validation, default resolution, domain preparation, and execution. This ADR also defines rules for where interfaces should exist, how DTOs should be used, and how result objects should be designed.

## Status

Accepted.

This ADR establishes the target internal workflow architecture for the ongoing architecture reset and rewrite effort.

Update note (2026-07-17, package DAG enforcement): Workflow and builder
composition may inject private validators into science-owned interpreters,
preprocessors, and resolvers through consumer-owned protocols. Science
collaborators must not import concrete `phospy.validation` implementations;
public API adapters and workflow orchestration own concrete validator wiring.

Update note (2026-07-30, preprocessing ownership decomposition): Dataset
preprocessing planning, state trace, and stage-result DTOs are owned by
`phospy.science.datasets.preprocessing.plan`,
`phospy.science.datasets.preprocessing.trace`, and
`phospy.science.datasets.preprocessing.results`. The legacy
`phospy.science.datasets.preprocessing.models` route is a compatibility facade
only.

Update note (2026-08-05, preprocessing stage extension contract):
Preprocessing stage composition separates static registration metadata,
plan-dependent contract interpretation, and executable stage behavior. This is
an internal extension seam, not a supported third-party public API. The
deprecated `PreprocessingStageMetadata` and `stage_metadata_registry` aliases
remain temporary compatibility routes and keep their documented PhosPy 1.8.0
removal target.

Update note (2026-08-06, dataset preprocessor lifecycle):
`DatasetBuildExecutor` composes one internal `DatasetPreprocessorContract`.
Every collaborator on that seam exposes `validate_quantitative_contracts(...)`
and `run(...)` with the same declared initial quantitative context. The builder
does not discover optional hooks, inspect collaborator signatures, or traverse
the preprocessing stage registry as a fallback; interpretation of preprocessing
stage quantitative contracts remains owned by `PreprocessingPipeline`.

Update note (2026-08-08, contracts facade dependency rule): Public request
DTOs remain passive transport objects. `phospy.contracts` may expose exact
science-owned domain values through stable facade routes. The executable rule
allows `phospy.science.configs.*` policy/configuration modules by prefix; other
science modules must declare the private module marker
`__phospy_contracts_facade_role__` with an approved role such as
`science_owned_public_model`, `science_owned_public_enum`,
`science_owned_public_constant`, `science_owned_public_table_contract`, or a
narrow public helper role. Contextual validation, execution, construction,
internal views, workflow implementation modules, reference builders/reference
validation, and private validation implementations remain outside the contracts
package even if a forbidden module declares that marker.

## Context and Problem Statement

PhosPy has accumulated orchestration patterns that are more complex than necessary for the product it is intended to be. The current direction has shown signs of wrapper-heavy execution paths, repeated validation across multiple layers, duplicated accessors, loose helper composition, and abstractions that are more "smart" than useful.

This complexity is not justified by the product goal. PhosPy is intended to be a maintainable Python port of PhosR. That requires internal architecture that preserves scientific correctness and conceptual clarity without turning workflow execution into a framework.

The package therefore needs a strict internal workflow model that is:

- easy to understand
- easy to test
- easy to extend at genuine seams
- resistant to wrapper proliferation and accidental complexity

## Decision Drivers

The decision is driven by the following considerations:

1. **Clarity of execution flow.** The path from request to result must be obvious.
2. **Separation of concerns.** Validation, interpretation, and execution are different responsibilities and should not be collapsed into ambiguous service objects.
3. **Maintainability.** Internal architecture should reduce drift, duplicated plumbing, and pass-through wrappers.
4. **Testability.** Each workflow stage should be testable independently.
5. **Extensibility at real seams.** Interfaces should support meaningful substitution, not exist for decoration.
6. **Discipline.** Internal structure should not expand simply because a helper or wrapper seems convenient.

## Decision

Each public workflow in PhosPy will be implemented through three internal stages:

- validator
- interpreter
- executor

The public workflow class is a thin coordinator that passes a request through these stages and returns a result.

### Workflow Pattern

The intended shape is:

1. accept a single public request DTO
2. validate the request
3. interpret the request after validation has succeeded into explicit executable inputs
4. execute domain logic
5. return a single public result DTO

The public workflow class should remain deliberately boring.

Illustrative shape:

```python
class KinaseWorkflow:
    def __init__(self, validator, interpreter, executor) -> None:
        self._validator = validator
        self._interpreter = interpreter
        self._executor = executor

    def run(self, request: KinaseWorkflowRequest) -> KinaseWorkflowResult:
        validated = self._validator.run(request)
        interpreted = self._interpreter.run(validated)
        return self._executor.run(interpreted)
```

The same pattern applies to `DifferentialAnalysisWorkflow` and
`SignalomeWorkflow`.

Configuration follows the same boundary. Public request/config DTOs record user
intent and remain stable transport objects. When a workflow needs contextual
interpretation, the interpreter converts those DTOs into distinctly named
resolved execution models before numerical science code is called. Executors and
science kernels consume the resolved model or primitive science-owned policy
values, not raw public workflow request DTOs.

## Method Naming Convention

All workflow-stage components must expose a `run(...)` method.

This includes:

- validators
- interpreters
- executors
- public workflows

The responsibility of a component is expressed by its class name, not by a verb-specific method name. This keeps stage call sites consistent and easy to remember.

Examples:

- `self._validator.run(request)`
- `self._interpreter.run(validated)`
- `self._executor.run(interpreted)`

Verb-specific stage methods such as `validate(...)`, `interpret(...)`, or `execute(...)` are discouraged for these pipeline components.

## Stage Definitions

### Validator

The validator is responsible for deciding whether a request is structurally and scientifically valid.

It may:

- validate types
- validate required inputs
- validate configuration ranges and thresholds
- validate dataset shape and metadata requirements
- validate compatibility between request components
- normalise simple input forms where necessary for validation only

It may not:

- run scoring, prediction, activity, or signalome algorithms
- perform broad default resolution that belongs to interpretation
- assemble final public result objects
- act as a general-purpose workflow service

Its normal output should be the original request unchanged, or otherwise simple successful completion where the surrounding pipeline does not require a returned value. The validator remains validation-focused rather than reshape-focused.

### Interpreter

The interpreter is responsible for turning a valid request into explicit executable domain inputs.

It may:

- resolve presets such as `auto`
- resolve references into concrete bundles
- make default policies explicit
- prepare execution-ready request structures
- derive domain inputs from validated user input

It may not:

- repeat validation already completed by the validator without a specific reason
- run scientific algorithms
- own workflow result assembly concerns

Its output should be an interpreted DTO that contains everything the executor needs, without requiring additional public-contract reasoning.

### Executor

The executor is responsible for running the actual domain logic on fully interpreted inputs.

It may:

- call preprocessing collaborators where still needed internally
- call scoring services
- call prediction services
- call activity analysis services
- call signalome analysis services
- assemble the public result model

It may not:

- perform broad public API validation
- interpret loose user-facing options that should already be resolved
- compensate for weak request modelling through additional ad hoc policy logic

The executor should be the stage where domain services are coordinated, not the stage where the public request is still being understood.

### Post-Resolution Validator Pattern

Some workflows have scientific eligibility checks that cannot be evaluated from
the raw public request. Those checks depend on references, policies, or
configuration that the interpreter must first resolve into executable inputs.
This is the **post-resolution validator** pattern.

Responsibility is split as follows:

- initial workflow validators validate request, dataset, and configuration
  invariants available before interpretation
- interpreters may resolve references, configuration, policies, and identifiers
  into execution-ready inputs
- resolved validators validate scientific eligibility that depends on those
  resolved inputs
- executors must not perform validation except defensive internal checks for
  impossible or programmer-error states

The kinase workflow is the motivating case. `KinaseWorkflowValidator` can check
request-level configuration policy, dataset identity requirements, localisation
policy, and reference-input compatibility before interpretation. It cannot know
whether a selected or caller-supplied reference bundle will still have enough
usable substrate overlap after display IDs are projected to dataset `site_key`
rows, site-sequence support is merged, and execution config is resolved.

`KinaseWorkflowInterpreter` therefore resolves the references, projects
reference substrate IDs onto dataset row identity, resolves sequence support and
execution config, and then delegates the resolved scientific eligibility checks
to `ResolvedKinaseEligibilityValidator` before the executor receives the
execution request. That resolved validator owns checks such as reference
coverage, eligible kinase counts after projection, sequence-supported scoring
sites, scored-site support, attrition policy, and kinase-library resource
usability.

A post-resolution validator is not a general validation hook. It must have a
clear resolved-input DTO, it must run after interpretation and before execution,
and it must not repeat request-boundary validation or run scoring, prediction,
activity, or signalome algorithms. The public workflow still coordinates a
validator, interpreter, and executor; the post-resolution validator is a narrow
validation seam used only when interpretation creates facts required for
scientific eligibility.

### Protocol Adapter Rule

Interpreters, executors, and preprocessing stages define the protocols they
consume when they need validation or external loading collaborators. Concrete
validators remain private to `phospy.validation` and are wired by API/workflow
adapter code. Local imports, dynamic imports, or service-locator style lookup
must not be used to hide package edges.

### Preprocessing Stage Extension Contract

Dataset preprocessing stages are composed through one internal contract path:
`PreprocessingStageContract` entries in
`phospy.science.datasets.preprocessing.stage_registry`. Each contract provides:

- static registration metadata: stage key, display label, provenance stage,
  table inputs/outputs, diagnostics metadata, backend, provenance inclusion, and
  a stage factory
- plan-dependent interpretation: operation name, serialized parameter payload,
  quantitative-operation contract, determinism declaration, inclusion predicate,
  random-seed resolver, and plan-specific validation
- executable behavior: a stage object implementing the explicit
  `PreprocessingStage` protocol

The executable stage protocol has two lifecycle hooks:

- `validate_before_quantitative_contract(state) -> None`
- `run(state) -> PreprocessingStageResult`

The first hook is explicit even when a stage has no work to do. Stages that do
not need data-dependent validation before quantitative contract enforcement must
implement it as a no-op. The pipeline invokes this hook directly; it must not
discover the hook through `getattr(...)`, `hasattr(...)`, or another
attribute-presence capability probe.

Registry construction validates only static structure: non-empty keys and
labels, duplicate stage keys, declared table-key shape, contract declaration
shape, callable contract fields, diagnostics metadata shape, and factory-backed
construction. It must not resolve operation names, parameters, inclusion
predicates, or quantitative contracts using `PreprocessingPlan.default()` or an
unrelated `PreprocessingPlan()`.

Plan-dependent validation runs against the actual `PreprocessingPlan` supplied
to the pipeline or builder. The pipeline interprets each scheduled stage against
that plan, validates declared stage-order table dependencies beyond initial
dataset inputs and optional preprocessing state sidecars, checks the
quantitative transition before numerical execution, validates emitted evidence
after execution, and applies the quantitative state transition. Quantitative
transition enforcement remains pipeline-owned; individual stage executors own
only their operation and any explicit pre-contract data validation they require.

The dataset builder calls preprocessing through one fixed internal lifecycle:

1. `validate_quantitative_contracts(...)` validates the actual plan's
   quantitative transitions from the declared initial scale kind and
   quantitative meaning.
2. `run(...)` receives the phospho table, site metadata, sample metadata, total
   table, actual plan, corrected preprocessing output including `None`, initial
   quantitative scale kind including `None`, and initial quantitative meaning
   including `None`.
3. The preprocessor emits typed preprocessed tables, trace records, stage
   evidence, and transformation events for downstream transformation-state and
   provenance assembly.

This is an internal maintenance seam. It is not a supported public extension API.

## Interface Rules

Interfaces are allowed and encouraged only at real extension seams.

### Interfaces That Are Appropriate

Good candidates for injected interfaces include:

- workflow validator
- workflow interpreter
- workflow executor
- reference provider
- local table/path reader
- reference source reader
- nested workflow runner used by orchestration
- dataset builder or preprocessing builder
- kinase activity analyser
- output publisher

These are real seams where alternative implementations are plausible and useful.

### Orchestration-Owned Concrete Adapters

Scientific packages must express external collaborators as protocols or narrow
constructor-injected contracts. They must not instantiate concrete I/O readers or
other workflows directly.

The public API and workflow orchestration layers are responsible for wiring
default adapters. Examples include local dataset table readers, local reference
source readers, Kinase Library resource readers, and the SPS/RUV-style batch
correction workflow runner used by dataset preprocessing. This keeps scientific
models and stages testable without importing `phospy.io` or `phospy.workflows`.

### Interfaces That Are Discouraged

Interfaces should not be introduced for:

- DTOs
- simple models
- value objects
- tiny helper functions
- wrappers that only delegate once
- helpers with one stable implementation and no realistic substitute

The existence of a class does not automatically justify an interface.

### Runtime Validation of Collaborators

Runtime collaborator validation must not be implemented through loose `hasattr(...)` checks on method names.

Where runtime validation is needed, use one of the following:

- explicit protocol types
- well-defined abstract base classes
- constructor-time assumptions backed by tests

Attribute-presence probing is too weak and too easy to drift.

## DTO Rules

DTOs are part of the architecture and must be used consistently.

### Public DTOs

Public DTOs define the stable product contract. These include:

- request DTOs
- result DTOs
- major public configuration DTOs

Public DTOs should be small, explicit, and easy to validate.

Request DTOs are passive. Their constructors store caller intent and may not
call workflow validators, dataset validators, construction services, executors,
or contextual scientific checks. When a request field needs a science-owned
enum, model, or policy object, the contracts facade imports and re-exports the
owned object identity instead of copying the definition into contracts.

The contracts-to-science import boundary is rule based, not maintained as an
exact exception list. `phospy.science.configs.*` is the only approved prefix
because configuration policies are direct science-owned values consumed by
workflow and numerical code. Other science owner modules opt in by declaring
`__phospy_contracts_facade_role__` with a constrained public role. The marker
means only that `phospy.contracts` may re-export the original object identity;
it does not transfer behavioural ownership to contracts and it cannot legalize
private modules, builders, construction services, executors, internal views,
validation implementations, reference builders/reference validation, or
workflow implementation modules.

### Internal DTOs

Internal DTOs define workflow stage boundaries. These include:

- interpreted request DTOs
- execution output DTOs if needed
- other stage-boundary DTOs where they clarify a real handoff

Internal DTOs are encouraged where they make stage boundaries clearer and reduce reliance on loosely coupled scalars. Validation itself should not manufacture dedicated DTOs as its normal job.

### DTO Usage Rules

The following rules apply:

1. Public workflows accept exactly one public request DTO.
2. Public workflows return exactly one public result DTO.
3. Workflow stages should pass typed DTOs, not long lists of unrelated scalar values.
4. Configuration should move through the workflow as config objects, not be repeatedly exploded into scalar parameters.
5. DTOs should model real stage boundaries, not act as generic bags for convenience.

## Result Design Rules

Result objects should be primarily data containers.

### What Result Objects Should Do

Result objects should:

- hold the domain outputs of a workflow or workflow stage
- expose nested stage results directly where relevant
- provide a very small number of high-value convenience properties only when clearly justified

### What Result Objects Should Not Do

Result objects should not:

- act as orchestration services
- duplicate nested state through mirrored accessors
- become export engines
- become plotting adapters
- accumulate compatibility aliases by default
- hide expensive copy operations behind innocent-looking properties

### Key Rule

Top-level workflow result objects must not mirror nested stage outputs through duplicated convenience accessors unless the shortcut is clearly essential to the product experience.

Users should be able to understand the real nested structure of results.

## Dataset Design Rules

The analysis-ready dataset model is a public data boundary, not an orchestration component.

The dataset model should:

- represent validated analysis-ready state
- expose clear data fields and simple derived properties
- avoid owning workflow logic

Creation, shaping, and validation concerns should be handled by separate builders or validators where helpful.

A healthy split is:

- dataset model
- dataset builder or preprocessing service
- dataset validator

## Anti-Patterns and Prohibited Directions

The following patterns are considered architectural anti-patterns for PhosPy unless there is a strong and explicit justification:

- redundant internal workflow wrappers
- composition graphs that exist mainly to delegate
- repeated validation across multiple workflow layers without clear stage ownership
- duplicated result accessor surfaces
- `hasattr(...)`-based collaborator validation
- magic export indirection for public modules
- long scalar-heavy workflow signatures
- result models that attempt to solve mutability, export, plotting, aliasing, and convenience access simultaneously
- abstractions introduced only to appear clean rather than to solve a real seam

## Consequences

### Positive Consequences

- Every workflow follows the same internal logic and is easier to understand.
- Validation concerns become explicit and testable.
- Default resolution and interpretation logic are no longer mixed with execution.
- Domain services can be coordinated cleanly without wrapper proliferation.
- Extension remains possible at real seams without turning the package into a framework.
- Result models stay smaller, clearer, and more honest.

### Negative Consequences

- Some current internal classes and helper layers will need to be removed or rewritten.
- Existing code that combines multiple concerns in one class will need to be split.
- The architecture imposes discipline that may make quick convenience shortcuts less acceptable.
- Some current abstractions will be judged unnecessary and deleted.

### Neutral Consequences

- Domain services may still be rich internally, but they no longer define the public story.
- Internal helper modules may continue to exist, but they should not shape workflow contracts.

## Rejected Alternatives

### Alternative 1: Single Large Workflow Service per Domain

This option was rejected because it encourages classes that mix validation, interpretation, execution, and result shaping in one place. That shape tends to become hard to test and hard to reason about over time.

### Alternative 2: Highly Abstract Composition Graph

This option was rejected because it introduces infrastructure-oriented complexity that is not justified by the size or product goal of the package.

### Alternative 3: Thin Wrappers Around Many Small Helpers with No Strict Stage Boundaries

This option was rejected because it tends to create pass-through layers, duplicated plumbing, and unclear ownership of decisions.

### Alternative 4: Expose Internal Service Boundaries as Public API

This option was rejected because the product contract should be workflow-oriented and PhosR-aligned, not shaped by internal implementation pieces.

## Implementation Guidance

### Workflow Classes

Each public workflow class should remain a simple coordinator. If the public workflow class grows beyond simple stage orchestration, its responsibilities should be reviewed.

Workflow stages should be invoked through a uniform `run(...)` method rather than verb-specific method names.

### Validation

Validation should happen at the workflow boundary. After validation, downstream stages should not repeatedly re-check the same public conditions unless there is a narrowly scoped reason.

### Interpretation

Interpretation is a distinct stage and must not be hidden inside validation or execution. This separation is important because a request can be valid without yet being fully resolved into concrete execution inputs.

### Execution

Executors should coordinate domain services and assemble results. They should not compensate for weak earlier stages.

### Oversized Stage Decomposition

When an interpreter or executor grows large enough that design construction,
scientific eligibility, model fitting, diagnostics, result shaping, and
provenance are mixed in one class, the stage must be decomposed into
domain-specific collaborators instead of being split into passive files.

For differential workflows, responsibility limits are:

- design assembly owns execution design matrices, contrast vectors, covariate
  column metadata, block metadata, sample order, formula, and design
  descriptions
- eligibility owns pre-fit feature eligibility, imputation-based withholding,
  post-fit numerical eligibility, and filtering to the tested feature family
- fitting owns calls into the scientific differential computation kernel only
- result assembly owns expansion back to the public feature index, identity
  metadata attachment, diagnostics, and public result construction
- provenance assembly owns row-attrition metrics, row-attrition reports, and
  workflow provenance payload updates

For enrichment workflows, responsibility limits are:

- validation and interpretation own identifier semantics, explicit background
  universe preparation, and method configuration
- set execution preparation owns min/max set-size filtering against the
  background-intersected set size
- the ORA kernel remains a pure scientific kernel with no pandas result-table,
  workflow, provenance, or output-serialization imports
- result assembly owns public enrichment records, summaries, diagnostics,
  caveats, and public result construction
- provenance assembly owns run provenance, row-attrition reports, and table
  fingerprints

Workflow executors coordinate these collaborators. They should not directly
construct every diagnostic payload, provenance object, and public result table.
This decomposition is not a license to introduce a generic service framework:
collaborators must remain concrete, domain-specific, and acyclic.

### Config Movement

Configuration should flow through the system as typed objects. Repeatedly unpacking config into long scalar argument lists is discouraged and should be treated as a maintainability smell.

Typed config flow does not mean public transport DTOs should cross every layer.
Public configs may be stored on validated/interpreted workflow records for
requested-value provenance, while resolved execution models carry effective
defaults, references, seeds, eligible sets, and contextual choices for execution.
Those resolved models must use different names from public DTOs and must be
owned by the workflow or science domain that interprets them.

## Scope Boundaries

This ADR defines internal workflow architecture only.

It does not define:

- the detailed public API contract beyond its interaction with workflow stages
- the exact module layout of every internal package
- the detailed structure of every domain service
- migration strategy for current internal code
- packaging or release management

Those concerns should be handled in separate ADRs or implementation plans.

## Validation and Review Criteria

Future code and review work should check future changes against the following questions:

1. Does this class clearly belong to validation, interpretation, execution, or data modelling?
2. Does this change make the request-to-result path clearer?
3. Does this interface represent a real extension seam?
4. Does this DTO define a real stage boundary?
5. Does this result model expose domain outputs honestly rather than through duplication?
6. Does this change reduce or increase wrapper and plumbing complexity?

If the answers are weak, the design should be reconsidered.

## Relationship to ADR-0001

This ADR complements ADR-0001, which defines the public API contract for PhosPy.

ADR-0001 defines what PhosPy exposes.
ADR-0002 defines how the internal workflow architecture should deliver that public contract.

The public API remains centred on one dataset model and three primary
downstream workflows (differential, kinase, signalome). The internal
architecture defined here is the implementation discipline that supports that
public shape.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14), i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub. https://github.com/PYangLab/PhosR
