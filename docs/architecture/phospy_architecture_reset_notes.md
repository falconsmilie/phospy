# PhosPy Architecture Reset Notes

## Purpose

This document captures the intended direction for PhosPy as a maintainable Python port of PhosR. It records the product goal, target public API, architectural rules, and the simplification principles that should guide future implementation work.

## Rewrite Boundary

The codebase now enforces a hard cutover boundary:

- `src/phospy/` is the only valid home for new implementation work.
- `src/phospy_legacy/` is retained only as migration reference material.
- The legacy tree is not a supported public API and must not be extended.

## Product Goal

PhosPy is intended to be a Python port of PhosR.

That means the software should prioritise:

- scientific correctness
- conceptual fidelity to PhosR workflows
- usability for phosphoproteomics analysis
- maintainability of the Python implementation

It should not prioritise historical compatibility with previous internal Python designs.

## Core Product Shape

PhosPy should expose a very small public API:

- one analysis-ready dataset model
- one kinase workflow
- one signalome workflow

The package should feel like a workflow-oriented scientific application, not a sprawling utility framework.

## Primary Architectural Principle

Each workflow should follow the same simple path:

1. accept user input
2. validate user input
3. interpret user input into executable domain inputs
4. execute domain logic
5. return a result

Any class or abstraction that does not make that path clearer should be questioned.

## Public API Contract

### Public types

The intended top-level public API is:

- `AnalysisReadyPhosphoDataset`
- `Organism`
- `ReferencePreset`
- `ReferenceBundle`
- `KinaseScoringConfig`
- `KinasePredictionConfig`
- `KinaseActivityConfig`
- `SimpleKinaseWorkflowRequest`
- `SimpleKinaseWorkflowResult`
- `SimpleKinaseWorkflow`
- `SignalomeConfig`
- `SignalomeWorkflowRequest`
- `SignalomeWorkflowResult`
- `SignalomeWorkflow`
- `PhosPyValidationError`

No additional public workflows, result facades, or compatibility wrappers should be added without a strong domain reason.

### Dataset

`AnalysisReadyPhosphoDataset` is the only public dataset model.

It represents validated, analysis-ready phosphoproteomics input data.

It should contain:

- `phospho`
- `site_metadata`
- optional `sample_metadata`
- optional `total`
- optional `organism`
- `intensity_scale`

Required `site_metadata` columns should be standardised rather than passed around as workflow arguments:

- `gene_symbol`
- `site`
- `site_sequence`

### Kinase workflow

`SimpleKinaseWorkflow` is the public entry point for PhosR-style kinase analysis.

It should accept exactly one request object:

- `SimpleKinaseWorkflowRequest`

It should return exactly one result object:

- `SimpleKinaseWorkflowResult`

The request should contain:

- dataset
- reference selection
- scoring config
- prediction config
- activity config

The result should contain:

- dataset
- resolved references
- scoring result
- prediction result
- optional activity result

### Signalome workflow

`SignalomeWorkflow` is the public entry point for PhosR-style signalome analysis.

It should accept exactly one request object:

- `SignalomeWorkflowRequest`

It should return exactly one result object:

- `SignalomeWorkflowResult`

The request should depend on the kinase workflow result so the public pipeline stays coherent and hard to misuse.

The request should contain:

- kinase workflow result
- signalome config

The result should contain:

- dataset
- kinase workflow result
- module assignments
- signalome modules
- kinase network
- optional expanded signalome

## Internal Workflow Design

Each workflow should be implemented through three internal stages:

- Validator
- Interpreter
- Executor

### Validator

The validator is responsible for deciding whether a request is structurally and scientifically valid.

It should:

- validate types
- validate required inputs
- validate ranges and thresholds
- validate dataset and reference compatibility

It should not:

- run domain algorithms
- resolve defaults beyond what is required for validation
- assemble public results

### Interpreter

The interpreter is responsible for turning valid user input into explicit executable domain inputs.

It should:

- resolve presets such as `auto`
- resolve reference bundles
- make defaults explicit
- transform request input into execution-ready models

It should not:

- repeat validation already completed
- run scientific analysis
- mirror public result concerns

### Executor

The executor is responsible for running the actual domain logic on fully interpreted inputs.

It should:

- invoke scoring
- invoke prediction
- invoke activity analysis
- invoke signalome analysis
- return the final public result model

It should not:

- perform public API validation
- perform policy interpretation that belongs earlier

## Interface Guidance

Interfaces should only be used for real extension seams.

Good candidates for injected interfaces:

- workflow validator
- workflow interpreter
- workflow executor
- reference provider
- dataset builder or preprocessing builder
- kinase activity analyser
- output publisher

Interfaces should not be created for every helper or model.

Avoid interfaces for:

- simple DTOs
- simple value translators
- wrappers that only delegate once
- helpers with one stable implementation and no realistic alternative

## DTO Guidance

DTOs are encouraged, but should be used in a disciplined way.

### Public DTOs

Public DTOs define the supported public contract:

- workflow request models
- workflow result models
- major configuration models

### Internal DTOs

Internal DTOs should define stage boundaries:

- validated request DTOs
- interpreted request DTOs
- execution result DTOs if needed

Avoid passing long lists of scalar arguments through orchestration layers.

## Result Design Rules

Result objects should be mostly data containers.

They should:

- expose real domain outputs
- expose nested stage results directly
- include only a very small number of high-value convenience properties if truly justified

They should not:

- duplicate nested state through alias properties
- act as service layers
- own export pipelines
- own plotting adapters
- implement broad mutable ownership APIs unless clearly necessary

### Important rule

Top-level result objects should not mirror nested fields through convenience accessors unless the shortcut is clearly essential.

Users should navigate the nested structure honestly.

## Dataset Design Rules

The dataset model should be the single public input boundary for analysis-ready phosphoproteomics data.

It should not contain workflow logic.

Creation and validation logic should live in separate builders or validators where appropriate.

Good split:

- `AnalysisReadyPhosphoDataset` as model
- dataset builder service
- dataset validator service

## Anti-Goals

The following should be treated as anti-patterns unless a strong case is made:

- redundant internal workflow wrappers
- composition graphs that exist mainly to delegate
- broad compatibility layers
- duplicated result accessors
- `hasattr(...)`-based collaborator validation
- magic export indirection
- hidden expensive copies behind innocent-looking properties
- long scalar-heavy workflow signatures

## Decision Rules for New Code

When adding or reviewing a class, ask:

- does this help preserve PhosR scientific behaviour?
- does this make the path from request to result clearer?
- does this represent a real domain seam?
- does this reduce maintenance cost?

If the answer is no, the class should probably not exist.

## Recommended Package Shape

A target structure could look like this:

```text
phospy/
  api/
    datasets.py
    workflows.py
    requests.py
    results.py

  workflows/
    kinase/
      contracts.py
      validator.py
      interpreter.py
      executor.py
    signalome/
      contracts.py
      validator.py
      interpreter.py
      executor.py

  datasets/
    model.py
    builder.py
    validator.py

  preprocessing/
    services.py

  references/
    provider.py
    models.py
    resources.py

  prediction/
    scoring_service.py
    prediction_service.py
    models.py

  activities/
    activity_service.py
    models.py

  signalomes/
    signalome_service.py
    models.py

  io/
    publishing.py
```

This keeps domain separation while removing unnecessary orchestration layers.

## Rewrite Guidance

Backward compatibility is not the priority.

The current codebase should be treated as:

- a source of working scientific logic
- a source of validation rules
- a source of tests and examples

It should not be treated as a structure that must be preserved.

The rewrite should preserve:

- correct scientific behaviour
- PhosR-aligned workflow concepts
- the new public API contract

It should not preserve old wrappers, aliases, or accidental internal structure.

## Immediate Next Steps

1. Finalise the public API contracts for the dataset and the two workflows.
2. Define validator, interpreter, and executor interfaces for each workflow.
3. Port scientific domain logic into the new seams.
4. Remove redundant internal workflow layers and duplicate result accessors.
5. Keep results and dataset models deliberately boring.

## Short Architecture Statement

PhosPy is a Python port of PhosR built around one analysis-ready dataset model and two primary workflow entry points: kinase analysis and signalome analysis. Each workflow accepts a single typed request, validates it, interprets it into executable domain inputs, executes domain services, and returns a typed result. Internal implementation uses injected interfaces only at real extension seams. Result objects are data containers, not orchestration or export services.
