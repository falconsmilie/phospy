# ADR 0003: High-level adapter workflow for kinase inference

- Status: Proposed
- Date: 2026-04-10

## Context

PhosPy already has the core pieces needed for kinase inference:

- validated phospho and total dataset handling
- preprocessing and protein correction
- site-matrix construction
- prediction and activity workflows

What it does not have is one high-level public entry point for the common user path.

Today, users who start with study-shaped phospho data still need to translate those
inputs into workflow-shaped inputs before they can call the full kinase workflow.
That usually means doing some or all of the following outside the package:

- normalising phosphosite identifiers
- resolving one row per phosphosite
- cleaning or centring site sequences
- aligning phospho and total inputs
- building `site_sequences`
- providing `substrate_map`
- providing `motif_sequences`

This is too much work for the common path.

The problem is not that the current low-level workflow is wrong. The problem is that it
is currently the only real public lane for end-to-end kinase inference. That makes
PhosPy harder to use than a custom notebook for users who already have biologically
meaningful inputs and just want to run the package.

## Decision

PhosPy will add a high-level adapter workflow for kinase inference.

This workflow will accept user-shaped biological inputs and build the lower-level
workflow inputs required by the existing prediction and activity layers.

The adapter layer will introduce two supporting concepts:

- `AnalysisReadyPhosphoDataset`
- `ReferenceProvider`

### 1. `AnalysisReadyPhosphoDataset`

PhosPy will introduce a normalised intermediate dataset object for analysis-ready
phosphosite data.

This object should represent:

- one row per phosphosite
- stable phosphosite identifiers
- clean site-centred sequences
- aligned phospho values
- optional protein-corrected values
- row metadata
- preprocessing provenance

This becomes the boundary between raw input handling and kinase inference.

### 2. `ReferenceProvider`

PhosPy will introduce a reference-provider abstraction responsible for resolving kinase
prior inputs for the workflow.

At minimum, a provider must be able to supply:

- `substrate_map`
- `motif_sequences`
- reference metadata and provenance

This keeps reference resolution explicit and testable instead of pushing that burden
onto every caller.

### 3. `SimpleKinaseWorkflow`

PhosPy will add a new high-level workflow for the common kinase analysis path.

The intended shape is:

```python
result = SimpleKinaseWorkflow().run(
    phospho="phospho.tsv",
    total="total.tsv",
    species="mouse",
    reference="auto",
)
```

This workflow should:

1. load and normalise inputs
2. build an `AnalysisReadyPhosphoDataset`
3. resolve kinase reference inputs through a `ReferenceProvider`
4. call the existing workflow and activity layers
5. return a simple result object for downstream use

## Consequences

Positive consequences:

- the common path becomes shorter and easier to explain
- users no longer need to manually build workflow-shaped inputs in routine cases
- preprocessing and reference resolution become reusable package responsibilities
- PhosPy becomes easier to use as a replacement for notebook-style or PhosR-style workflows
- existing low-level workflow APIs can remain available for advanced use

Trade-offs:

- the public API surface grows
- PhosPy takes on explicit responsibility for a reference-resolution layer
- the adapter must stay narrow enough to avoid becoming a god interface
- some messy biological inputs will still require explicit user decisions

## Scope boundaries

This decision does not remove the existing low-level workflow APIs.

Advanced users should still be able to call the current workflow layers directly when
they need explicit control over:

- `substrate_map`
- `motif_sequences`
- `site_sequences`
- workflow settings
- intermediate results

This decision also does not require PhosPy to solve every raw import problem
immediately. The adapter should cover the common phosphosite-to-workflow path first.
Support for additional import styles can be added later through explicit adapters.

## Follow-on work

- define `AnalysisReadyPhosphoDataset`
- define `ReferenceProvider` and a first reference bundle contract
- add `SimpleKinaseWorkflow`
- document the simple lane and the advanced lane separately
- add migration examples from notebook-style phosphoproteomics workflows