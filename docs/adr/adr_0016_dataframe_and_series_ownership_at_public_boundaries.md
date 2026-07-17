# ADR: DataFrame and Series Ownership at Public Boundaries

## Document Control

- **ADR ID:** ADR-0016
- **Title:** DataFrame and Series Ownership at Public Boundaries
- **Status:** Accepted
- **Date:** 2026-05-02
- **Decision Type:** Architecture Decision Record

## Abstract

This ADR defines pandas ownership semantics at PhosPy public boundaries.
The contract is already implemented and tested, but it is architecture-level
governance because it protects API stability, provenance validity, and replay
expectations.

## Status

Accepted.

This ADR is aligned with ADR-0003 (dataset boundary), ADR-0005 (result model
design), ADR-0006 (dataset state contract), and ADR-0014 (test policy).

Update note (2026-07-16, pandas global option isolation): PhosPy frame ownership
helpers must not set, restore, or otherwise mutate process-global pandas options.
Borrowing semantics are local to PhosPy-owned objects.

## Context and Problem Statement

PhosPy datasets and workflow results carry mutable pandas objects internally.
Without explicit boundary rules, public accessor changes can accidentally:

- leak writable references to internal state
- undermine provenance fingerprint meaning
- make replay behavior depend on caller mutation order

The policy must cover both `DataFrame` and `Series`. A DataFrame-only rule is
incomplete because some public result fields expose `Series` values.

PhosPy also cannot promise deep immutability of pandas internals. It can govern
ownership and export semantics at package boundaries.

## Decision Drivers

1. Keep public API semantics predictable and easy to reason about.
2. Preserve provenance and replay meaning after public data export.
3. Prevent boundary aliasing between caller-owned and package-owned pandas state.
4. Keep high-volume output workflows explicit rather than implicit mutable-frame paths.
5. Accept boundary copy cost in exchange for contract clarity.

## Decision

PhosPy adopts and enforces the following ownership rules:

1. Internal models and workflow objects own mutable pandas state.
2. Public pandas accessors return defensive, caller-writable snapshots for both
   DataFrames and Series.
3. Mutation-isolated internal frame access is package-private, read-only by
   contract, and limited to controlled internal boundaries.
4. Provenance fingerprints describe owned internal state at result creation
   time.
5. Mutating pandas objects returned by public accessors must not mutate internal
   dataset/result state.
6. High-volume persistence is an explicit export or publishing concern, not an
   implicit public borrow path.
7. Defensive-copy cost at public boundaries is accepted by design.
8. PhosPy must not change host-application pandas global options, including
   `mode.copy_on_write`.

## Boundary Rules

### Public Accessors

Public accessors must expose mutation-isolated snapshots that callers may mutate
locally without changing owned dataset/result state. Public APIs must not offer
user-facing `copy=False` toggles for owned boundary data.

### Internal Frame Snapshot Paths

Package-private `_borrow_*` helpers are allowed for trusted internal
collaboration only. They return mutation-isolated internal snapshots, are
read-only by contract, are not part of the public contract, and must stay out of
public API routes. Callers must not rely on successful mutation of a borrowed
object; writes may raise or detach locally depending on the pandas runtime, but
must not mutate the owner.

Borrowed snapshots are implemented without process-global pandas mutation:

- pandas runtimes with native copy-on-write isolation may use shallow pandas
  copies.
- NumPy-backed pandas 2.x frames use shallow pandas copies whose borrowed blocks
  are local read-only views.
- Unsupported pandas internals, including extension arrays that cannot be made
  read-only through the local helper, fall back to deep copies.

Implementation note (2026-06-14): workflow dataset access is mediated by the
dataset-owned `DatasetInternalView`. Workflows may depend on that defensive
internal view for the specific frame snapshots they require, but must not call
dataset `_borrow_*` methods directly. Workflow access to prediction and scoring
result frames follows the same domain-owned internal view pattern.

### Provenance

Provenance fields describe owned state at creation time. Caller mutation of
exported snapshots must not alter those recorded fingerprints.

Derived quantitative datasets created inside workflows follow the same rule.
Their fingerprints describe the owned derived tables at derived-object creation
time, not the parent dataset's tables. The dataset domain recomputes those
fingerprints from the actual derived phospho matrix, site metadata, optional
sample metadata, total matrix, comparisons, and imputation observation mask
before accepting lineage. Parent dataset fingerprints may appear only as
explicit lineage inputs. Source run provenance and preprocessing reports must
not be reused as the derived object's own provenance.

## Consequences

### Positive Consequences

- Public API boundary behavior is explicit and consistent.
- Provenance fingerprint validity is protected from post-export mutation.
- Replay assumptions are clearer for users and contributors.
- Internal mutation semantics stay decoupled from caller mutation semantics.
- Host applications keep control of pandas process-global configuration.

### Negative Consequences

- Public accessor calls allocate copies and add memory/CPU overhead.
- Contributors must maintain stricter discipline around internal frame helpers.
- Some high-throughput paths need explicit publisher/export APIs for
  performance.
- Extension-array-backed borrowed snapshots may allocate deep copies on pandas
  versions without native local copy-on-write guarantees.

### Neutral Consequences

- Internal owned-state transfer between trusted components remains allowed where
  contractually controlled.

## Rejected Alternatives

### Alternative 1: Optional Public `copy=False` Escape Hatches

Rejected. This weakens contract clarity and increases aliasing risk across
public boundaries.

### Alternative 2: Zero-Copy Public Access by Default

Rejected. This allows caller mutation to alter internal state and provenance
assumptions.

### Alternative 3: Claim Deep Immutability for Public Pandas Objects

Rejected. pandas objects are mutable. PhosPy can enforce snapshot boundaries
but cannot claim deep immutability.

## Affected Modules

- `src/phospy/frames/ownership.py`
- `src/phospy/science/datasets/models.py`
- `src/phospy/api/results.py`
- `src/phospy/science/activities/models.py`
- `src/phospy/science/prediction/models.py`
- `src/phospy/science/signalomes/models.py`
- `src/phospy/provenance/`
- `src/phospy/io/`

## Testing and Governance Expectations

- `tests/unit/test_frame_ownership_policy.py` is the governance anchor for this
  ADR and must continue to enforce:
  - dataset DataFrame defensive-copy behavior
  - result DataFrame defensive-copy behavior
  - activity-result Series defensive-copy behavior
  - pandas global-option preservation
  - borrowed-view mutation isolation without process-global option changes
- New or changed public pandas accessors must include boundary-mutation tests
  that demonstrate internal state isolation.
- Public accessors must not reintroduce public copy-semantics toggles.
- Provenance-sensitive tests must continue to validate that public exports do
  not alter owned-state fingerprints.
- Derived-data tests must validate both sides of lineage: parent fingerprints
  remain recorded as inputs, and derived fingerprints match the actual derived
  matrices, sample mapping, and optional masks/tables.
- Export/publisher implementations remain the preferred high-volume persistence
  boundary.
- `benchmarks/measure_dataframe_ownership_copy_policy.py` records representative
  shallow/deep copy counts for the internal borrow policy.

## Scope Boundaries

This ADR defines pandas ownership behavior at public boundaries.
It does not define:

- scientific scoring rules
- workflow stage architecture
- broader performance policy

Those concerns are governed by other ADRs.

## Validation and Review Criteria

Future changes must satisfy all of the following:

1. Do public accessors return mutation-isolated snapshots?
2. Are package-private internal frame paths still private/internal only?
3. Do provenance-sensitive paths avoid aliasing with exported objects?
4. Are new accessor behaviors covered by explicit boundary-mutation tests?
5. Are high-throughput persistence paths kept in explicit publisher/export APIs?

## Relationship to Earlier ADRs

This ADR complements the earlier architecture decisions.

- ADR-0003 defines the dataset boundary where these semantics apply.
- ADR-0005 defines result-model shape that consumes these semantics.
- ADR-0006 defines dataset state contracts that rely on stable owned state.
- ADR-0014 defines the testing policy that enforces these rules.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi,
R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub.
https://github.com/PYangLab/PhosR
