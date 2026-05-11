# ADR: DataFrame and Series Ownership at Public Boundaries

## Document Control

- **ADR ID:** ADR-0016
- **Title:** DataFrame and Series Ownership at Public Boundaries
- **Status:** Accepted
- **Date:** 2026-05-02
- **Decision Type:** Architecture Decision Record

## Abstract

This ADR records the pandas ownership contract at PhosPy public boundaries.
The contract is already implemented and tested, but it must be explicit at the
architecture level because it governs API behaviour, provenance validity,
reproducibility guarantees, and expected copying costs.

## Status

Accepted.

This ADR is aligned with ADR-0003 (dataset boundary), ADR-0005 (result model
design), ADR-0006 (dataset state contract), and ADR-0014 (test policy).

## Context and Problem Statement

PhosPy workflows and datasets carry mutable pandas objects internally.
Without an explicit architecture rule, contributors can unintentionally weaken
public defensive-copy behaviour, introduce implicit borrow paths across public
boundaries, or create provenance drift between owned internal state and public
mutations.

This decision must cover both `DataFrame` and `Series`. A DataFrame-only policy
is incomplete because activity outputs and related result fields include
publicly exposed Series values.

The architecture must also state that pandas objects are not deeply immutable.
PhosPy controls ownership and export semantics; it does not claim deep
immutability of pandas internals.

## Decision

PhosPy adopts and enforces the following ownership rules:

1. Internal models and workflow objects own mutable pandas state.
2. Public pandas accessors return defensive snapshots for both DataFrames and
   Series.
3. Borrowed (non-copying) access is private/internal only and is permitted only
   inside controlled implementation boundaries.
4. Provenance fingerprints describe owned internal state at result creation
   time.
5. Mutating pandas objects returned by public accessors must not mutate internal
   dataset/result state.
6. High-volume persistence is an explicit export/publishing concern, not an
   implicit public borrow path.
7. Defensive-copy cost at public boundaries is accepted by design.

## Consequences

### Positive

- Public API semantics are stable and explicit across datasets and results.
- Provenance fingerprint validity is protected from post-export user mutation.
- Reproducibility assumptions are clearer for users and contributors.
- Internal mutation and public snapshot semantics remain decoupled.

### Negative

- Public accessor calls allocate copies and increase memory/CPU cost.
- Contributors must maintain stricter discipline around internal borrow helpers.
- Some high-throughput user paths need explicit export/publisher use for
  efficiency.

### Neutral

- Internal owned-state transfer between trusted components remains allowed where
  contractually controlled.

## Trade-offs

- **Safety and reproducibility over boundary performance:** public copies add
  overhead, but prevent external mutation from corrupting internal state.
- **Explicit persistence APIs over implicit zero-copy access:** exporting and
  publishing are clearer, auditable boundaries for high-volume output.
- **Strict boundary policy over convenience:** removing public `copy=...`
  toggles avoids drift and keeps one predictable semantics.

## Alternatives Considered

### Alternative 1: Optional public `copy=False` escape hatches

Rejected. This weakens contract clarity, complicates provenance expectations,
and invites accidental state aliasing across boundary calls.

### Alternative 2: Zero-copy public access by default

Rejected. This breaks isolation guarantees and allows user mutations to alter
internal state and invalidate provenance assumptions.

### Alternative 3: Claim deep immutability for public pandas objects

Rejected. pandas objects are mutable. PhosPy can enforce snapshot boundaries but
cannot honestly claim deep immutability.

## Affected Modules

- `src/phospy/_frame_ownership.py`
- `src/phospy/datasets/models.py`
- `src/phospy/api/results.py`
- `src/phospy/activities/models.py`
- `src/phospy/prediction/models.py`
- `src/phospy/signalomes/models.py`
- `src/phospy/provenance/`
- `src/phospy/io/`

## Testing and Governance Expectations

- `tests/unit/test_frame_ownership_policy.py` is the governance anchor for this
  ADR and must continue to enforce:
  - dataset DataFrame defensive-copy behaviour;
  - result DataFrame defensive-copy behaviour;
  - activity result Series defensive-copy behaviour.
- New or changed public pandas accessors must include boundary-mutation tests
  showing internal state isolation.
- Public accessors must not reintroduce public copy-semantics toggles.
- Provenance-sensitive tests must continue to validate that public exports do
  not alter owned-state fingerprints.
- Export/publisher documentation and implementations must remain the preferred
  persistence boundary for high-volume output.
