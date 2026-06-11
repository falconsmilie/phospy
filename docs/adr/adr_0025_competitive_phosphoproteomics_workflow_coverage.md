# ADR: Competitive Phosphoproteomics Workflow Coverage Roadmap

## Document Control

- **ADR ID:** ADR-0025
- **Title:** Competitive Phosphoproteomics Workflow Coverage Roadmap
- **Status:** Accepted
- **Date:** 2026-06-11
- **Decision Type:** Architecture and Scientific Roadmap

Complements ADR-0013, ADR-0015, ADR-0019, ADR-0022, and ADR-0024.

## Status

Accepted.

## Context

PhosPy is compared against established phosphoproteomics tools and reference
surfaces such as PhosR, MSstatsPTM, and Kinase Library. Those comparisons are
useful for planning, but they can create implementation drift when the roadmap
lives only in reviews, tickets, or external commentary.

The project needs one internal source of truth that records intended direction
without turning future capabilities into present-tense feature claims.

## Decision

PhosPy will maintain competitive workflow coverage as a phased roadmap governed
by this ADR and surfaced through `docs/scientific-coverage.md`.

This ADR is not an implementation claim. A roadmap item becomes supported only
when code, public contracts, documentation, and tests exist, and when the
scientific coverage matrix is updated to the correct claim category.

## Current State

The supported scientific workflow interface is the Python API.

Current executable public lanes are:

- `AnalysisReadyDatasetBuilder`
- `DifferentialAnalysisWorkflow`
- `KinaseWorkflow`
- optional kinase activity tables within the kinase workflow
- optional `SignalomeWorkflow` after kinase prediction

Current differential support is scoped to two-condition unpaired simple
contrasts with explicit `ExperimentalDesign` and `Contrast` objects,
empirical-Bayes `standard` or `robust` modes, optional trend moderation, and
Benjamini-Hochberg adjustment. Batch-aware, block, paired, and repeated-measure
modeling are not executable in the current public workflow lane.

Current kinase support provides profile/motif scoring, rank-weighted fusion,
deterministic/adaptive prediction, and two explicit activity methods:
`simplified_weighted_substrate_activity_v1` and
`ksea_zscore_activity_v1`. These outputs are support summaries or
substrate-set enrichment statistics, not calibrated causal inference.

Current bundled runtime reference data is rat-only. Human and mouse analyses
can be run only when the caller supplies an explicit `ReferenceBundle`.

Current ingestion supports analysis-ready tables and generic table I/O
contracts used by Python workflows. PhosPy does not currently provide broad
semantic importers for vendor, search-engine, or upstream statistical tool
outputs.

Current core PhosPy does not provide a first-class visualization workflow/API.

Current PhosPy does not support command-line scientific workflow execution.
ADR-0022 remains authoritative: Python API is the supported scientific
workflow interface unless a future ADR changes that boundary.

PhosPy does not claim full PhosR, MSstatsPTM, or Kinase Library equivalence.

## Desired Direction

The desired direction is a competitive but scoped Python phosphoproteomics
workflow stack:

- richer reference handling with explicit provenance, compatibility checks,
  and user-supplied external reference bundles
- improved kinase inference and activity methods with method-specific
  validation and clear output meaning
- semantic importers for common phosphoproteomics table outputs, without
  bypassing dataset validation or site identity contracts
- richer differential designs, including batch-aware and repeated-measure
  designs, only after explicit design/result contracts and parity or validation
  evidence exist
- enrichment workflows that are separate from kinase scoring unless the method
  is explicitly a kinase activity or substrate-set activity method
- visualization adapters that consume validated result objects without becoming
  the source of scientific truth
- possible CLI workflow support only as a thin, validated wrapper over the
  Python API after ADR-0022's reintroduction criteria are satisfied

## Non-Goals

This roadmap does not make PhosPy a clone or full replacement for PhosR,
MSstatsPTM, Kinase Library, MaxQuant, FragPipe, Spectronaut, DIA-NN, or any
other upstream processing or analysis package.

The roadmap does not permit:

- describing planned capabilities as already supported
- broad global parity claims
- hidden sample-name inference for scientific design
- silent reference remapping or protein identity guessing
- bundling reference data without redistribution permission
- adding workflow logic to visualization or CLI layers that diverges from the
  Python API
- treating statistical association, enrichment, or scoring output as causal
  biological proof

## Implementation Phases

Phases are directional, not calendar commitments.

### Phase 0: Guardrails

- Keep `docs/scientific-coverage.md` as the user-facing current support matrix.
- Keep roadmap items separate from current executable support.
- Require tickets for roadmap work to name the target claim category and
  affected workflow contracts.

### Phase 1: Reference and Import Foundations

- Harden `ReferenceBundle` provenance, compatibility diagnostics, and
  external bundle validation.
- Add semantic importer contracts only when they emit typed tables or requests
  that still pass existing builder/workflow validation.
- Keep raw/vendor/search-engine format interpretation out of core workflows
  unless a dedicated importer contract owns it.

### Phase 2: Differential Design Depth

- Extend experimental-design contracts before adding batch, block, paired, or
  repeated-measure execution.
- Preserve explicit contrast definitions and provenance.
- Add validation and parity/evidence tests before public support claims.

### Phase 3: Kinase, Activity, and Enrichment Depth

- Add kinase inference/activity methods one method at a time with stable policy
  IDs, output-scale definitions, and tests.
- Separate broad pathway/gene-set enrichment from kinase scoring unless the
  method is explicitly a kinase activity method.
- Keep prediction scores, activity scores, and enrichment statistics labeled by
  their actual statistical meaning.

### Phase 4: Visualization and CLI

- Add visualization as result-consuming adapters, not as hidden analysis
  engines.
- Reintroduce CLI workflow execution only after a separate decision confirms
  Python API validation parity, complete configuration coverage, provenance,
  documentation, and tests.

## Reference-Data Redistribution Rule

PhosPy may bundle reference data only when redistribution is explicitly allowed
by the source license or written permission, and when source provenance is
documented.

If a useful reference source has absent, ambiguous, restrictive, or
non-redistributable terms, PhosPy must not bundle that data. Acceptable
alternatives are:

- document the required schema
- provide validators or adapters for user-supplied local files
- provide scripts that transform user-provided inputs locally
- include only synthetic, minimal, or otherwise redistributable fixtures

Derived reference tables inherit the redistribution limits of their sources
unless the source terms explicitly allow redistribution of derived data.

The rat-only bundled runtime reference status remains the current support
boundary until a future implementation changes the bundled data and passes
license, provenance, documentation, and test review.

## Scientific Claim Categories

Public scientific claims must use the categories maintained in
`docs/scientific-coverage.md`:

- `parity-gated`
- `validated PhosPy implementation`
- `experimental`
- `open gap`
- `deliberate scope difference`
- `not planned`

A roadmap item is not a claim category. A future direction remains an
`open gap` or `deliberate scope difference` until executable support exists. An
`experimental` claim requires executable behavior and explicit caveats. A
`parity-gated` claim requires maintained parity evidence.

## Architecture Responsibility Boundaries

Documentation boundaries:

- `docs/scientific-coverage.md` owns current user-facing support status.
- `docs/workflow_contracts.md` owns executable workflow contracts and known
  limitations.
- ADRs own rationale, direction, and governance.
- README and API guide pages may link to the roadmap but must not expand
  current support claims.

Code ownership boundaries:

- Dataset construction and semantic importers must feed the
  `AnalysisReadyDatasetBuilder` contract rather than bypassing validation.
- Differential extensions belong under differential design, workflow, result,
  and provenance modules; sample-name inference remains out of scope.
- Reference resolution belongs under reference models/resources/resolution and
  must not be embedded ad hoc in workflow code.
- Kinase scoring, prediction, activity, and enrichment methods must own stable
  scientific policy records in their domain modules.
- Visualization layers may read result models but must not mutate scientific
  outputs or implement alternate workflow semantics.
- CLI layers, if reintroduced, must delegate to Python API request/workflow
  objects and preserve the same validation, provenance, and failure behavior.

## Consequences

Future work can be planned against competitive workflow areas without
overclaiming current support.

Review, documentation, and release checks must reject roadmap language that
describes planned capabilities as available behavior.

Tickets that add scientific scope should update the coverage matrix, workflow
contracts, and tests in the same change that adds executable support.
