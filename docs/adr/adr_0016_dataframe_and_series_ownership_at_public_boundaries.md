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

Update note (2026-08-02, dataset-owned immutable snapshots): internal dataset
read paths use owner-detached immutable snapshots owned by the constructed
dataset. `DatasetInternalView` instances return shallow read-only pandas
wrappers over those private snapshots, so repeated workflow runs against an
unchanged dataset do not recreate full matrix deep copies for internal reads.
This does not change public export semantics.

Update note (2026-08-07, public equality and hashing): public pandas-bearing
containers must not rely on dataclass-generated equality. Stable dataset,
result, table-wrapper, request, and provenance-bearing containers use explicit
identity equality and identity hashing unless they define a named scientific
content-comparison method. Those named methods compare owned pandas leaves with
`Index.equals`, `DataFrame.equals`, `Series.equals`, typed scalar comparison, or
stable fingerprints selected by the owning domain model.

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
3. Mutation-isolated internal frame access is package-private, immutable/read-only by
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
9. Dataclass-generated equality is forbidden for public containers that own or
   may carry pandas/NumPy payloads. Table fields must not be hidden with
   `compare=False` while retaining misleading partial value equality, and
   `unsafe_hash=True` is forbidden. Scientific content equality belongs to the
   owning domain model through a named method.

## Boundary Rules

### Public Accessors

Public accessors must expose mutation-isolated snapshots that callers may mutate
locally without changing owned dataset/result state. Public APIs must not offer
user-facing `copy=False` toggles for owned boundary data.

### Internal Frame Snapshot Paths

Package-private `_borrow_*` helpers are allowed for trusted internal
collaboration only. They return owner-detached internal snapshots, are read-only
by contract, are not part of the public contract, and must stay out of public
API routes. Callers must not rely on successful mutation of a borrowed object;
writes to NumPy-backed internal snapshots are expected to raise because PhosPy
marks the detached snapshot blocks non-writeable. Extension-array-backed
snapshots that cannot be made non-writeable through NumPy flags must remain
owner-detached and may fall back to copying when exposed from a cached internal
snapshot.

Borrowed snapshots are implemented without process-global pandas mutation and
without relying on undocumented pandas copy-on-write behavior:

- `frames.ownership` owns the copy/freeze helpers.
- Internal immutable snapshots first make one owner-detached deep pandas copy.
- Supported mutable object cells are frozen once for internal snapshots
  (`list` to `tuple`, `dict` to read-only mapping, `set` to `frozenset`, and
  nested `numpy.ndarray` values to non-writeable arrays), so workflow metadata
  reads do not recursively copy those cells on every access.
- NumPy-backed pandas blocks are marked non-writeable on the detached snapshot.
- Validated dataset construction creates one private
  `DatasetInternalFrameStore` owned by the dataset domain. The store lazily
  constructs at most one immutable snapshot per dataset frame.
- `DatasetInternalView` does not cache frames or snapshots. It returns shallow
  pandas wrappers over the dataset-owned read-only snapshot blocks on each
  access. Mutating wrapper metadata such as adding columns is local to that
  wrapper and does not mutate the cached snapshot or dataset owner.
- Unsupported pandas internals, including extension arrays that cannot be made
  read-only through NumPy flags, stay owner-detached and fall back to detached
  copies when exposed from a cached internal snapshot.

Implementation note (2026-06-14): workflow dataset access is mediated by the
dataset-owned `DatasetInternalView`. Workflows may depend on that defensive
internal view for the specific frame snapshots they require, but must not call
dataset `_borrow_*` methods directly. Workflow access to prediction and scoring
result frames follows the same domain-owned internal view pattern.

Implementation note (2026-08-02): the differential workflow threads one
`DatasetInternalView` from request validation into interpretation when the
validated dataset is unchanged. If technical-replicate aggregation produces a
derived dataset, interpretation creates a new view for that independent derived
workflow state. The views are run-scoped, but the immutable frame snapshots are
dataset-scoped so repeated runs can reuse them without exposing mutable frames.
Representative differential workflow tests bound full phospho-matrix
`DataFrame.copy(deep=True)` calls and separately assert that repeated runs build
the dataset phospho/site-metadata internal snapshots once.

Implementation note (2026-08-04): ordinary Signalome workflow execution threads
one validator-owned `DatasetInternalView` from the private validated request
through interpretation, protein-group resolution, executor table construction,
and result identity validation. Signalome result reconstruction from persisted
bundles is intentionally separate: bundle loading builds an isolated validation
view so standalone reconstruction remains safe without depending on a workflow
instance cache.

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
- Dataset-owned internal snapshots retain one owner-detached copy per accessed
  frame for the lifetime of the dataset. This is intentional and must remain a
  private dataset-domain store, not a workflow-instance cache, public borrow
  path, global cache, or object-ID-keyed cache.

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
- Workflow copy-count instrumentation must cover at least one representative
  differential run and assert the full phospho matrix is not repeatedly
  deep-copied by validator/interpreter handoffs.
- API contract tests must statically audit public pandas/NumPy-bearing
  dataclasses for implicit equality and `unsafe_hash=True`, and must exercise
  same-instance comparison, independent equivalent objects, scientifically
  different tables, different provenance over equal tables, bundle
  reconstruction/round-trip comparison, and hash behavior.
- Release-scale memory profiling remains an explicit local benchmark concern.
  The 50,000 x 48 builder+differential workload is measured by
  `make benchmark-release-scale` / `benchmarks/measure_release_scale_builder_differential.py`
  rather than default CI.

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
6. Are dataset-owned immutable snapshots private, read-only, non-global, and
   exposed to workflows only through `DatasetInternalView` wrappers?
7. Do public pandas/NumPy-bearing containers avoid implicit dataclass equality,
   pandas Boolean coercion, partial `compare=False` equality, and unsafe hashes?

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
