# ADR 0001: Data ownership and mutability policy

- Status: Accepted
- Date: 2026-04-04

## Context

PhosPy mixes trusted request and configuration models with pandas-backed processing objects.
That is reasonable, but the project has not had one explicit policy for how ownership,
copying, and mutation should work across those different layers.

The clearest example is `PhosphoDataset`. The class owns validated pandas `DataFrame`
instances and exposes them through explicit `*_live` and `*_copy` accessors. That
means callers can mutate the dataset's owned state through those accessors. Without an
explicit policy, different modules can drift into conflicting assumptions about whether a
pandas-backed object is immutable, mutable, copied, or shared.

This confusion also affects:

- where boundary validation should end
- when external inputs should be copied
- when internal services may mutate owned tabular state
- whether outward-facing APIs should return shared frames or detached copies

## Decision

PhosPy adopts the following data ownership and mutability policy.

### 1. Configuration objects and small value records should be immutable

Objects that represent execution options or small metadata without owning pandas
tables should be immutable value objects where practical.

Validated request bundles that carry pandas tables are different: they are
trusted boundary bundles by convention, but they are not truly immutable value
objects just because they wrap mutable pandas state.

### 2. DataFrame-carrying workspace objects are mutable

Objects that own pandas `DataFrame` instances used for active processing are mutable
workspace objects, not immutable snapshots.

Examples include:

- `PhosphoDataset`
- preprocessing workspaces
- other internal processing objects that own working tabular state

### 3. Ownership transfers at construction and other boundary points

When caller-supplied pandas tables enter trusted application state, PhosPy should make the
ownership transfer explicit at the construction or validation boundary.

For workspace constructors such as `PhosphoDataset`, that means the workspace takes
ownership by isolating itself from caller-managed frames once at construction time.
After that boundary, internal code should treat the resulting tables as owned application
state.

### 4. Copy once when taking ownership

PhosPy should copy external caller-owned tables when ownership transfer requires isolation.
It should not repeatedly copy the same tables inside already trusted internal flows.

### 5. Snapshots must be explicit

When an API returns a detached copy for caller-owned mutation or inspection, that must be an
explicit operation. Snapshots should never be implied accidentally by vague naming or by a
class presenting itself as immutable when it is not.

## Consequences

Positive consequences:

- the mutability contract becomes honest
- internal processing code can work on owned tables without copy theatre
- future performance work can reduce duplicate copying without weakening boundary rules
- future API changes can distinguish explicit live-access paths from explicit copies more clearly

Trade-offs:

- public docs must be careful not to imply that pandas-backed objects are immutable by
  default
- APIs that expose owned frames must say so plainly
- tests should pin aliasing and ownership-transfer behaviour directly

## Current application to `PhosphoDataset`

`PhosphoDataset` should be understood as an owner of mutable tabular state.

- `dataset.total_df_live` returns the owned validated total table
- `dataset.phospho_df_live` returns the owned validated phospho table
- constructing a dataset isolates it from later caller mutation of the original input frames
- `dataset.total_df_copy` and `dataset.phospho_df_copy` return detached deep copies
- `dataset.copy_inputs()` returns detached deep copies for caller-owned mutation

In other words, `PhosphoDataset` is a mutable workspace owner, not an immutable snapshot.
Follow-on refactors should keep making that contract clearer in the API itself, especially
around how callers choose between explicit live-access paths and explicit copies.

## Follow-on work

Future tickets should align naming and access patterns with this policy by:

- auditing other DataFrame-carrying types for the same problem
- standardising ownership-transfer and copy behaviour across public boundaries
