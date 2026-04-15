# Preprocessing Owned-Path Model

This document is the contributor contract for preprocessing ownership and copy behaviour.

## Execution model

Preprocessing has two explicit entry lanes:

| Lane | Intended caller | Boundary method | Copy expectation |
| --- | --- | --- | --- |
| Public safe lane | External callers that may pass caller-owned frames | `CoreProcessor.process()` / `process_phospho_only()` | One detached copy at entry, then owned processing |
| Owned fast lane | Trusted internal workflows that already own mutable frames | `CoreProcessor.process_owned()` / `process_phospho_only_owned()` | No extra full-frame defensive boundary copy |

The main workflow wiring is:

1. `DatasetLoader` materializes validated dataset frames once.
2. `DatasetPreprocessing.from_owned(...)` marks those frames as already owned.
3. `DatasetPreprocessing.run(...)` routes to `CoreProcessor.process_owned(...)`.
4. Services and site-matrix builders stay on `*_owned` paths (`prepare_owned`, `correct_owned`, `build_owned`).

## Copying rules

Keep these rules stable across refactors:

1. Copy once when ownership crosses a public boundary.
2. Do not re-copy full DataFrames inside trusted owned paths.
3. Keep safe and owned methods explicit instead of inferring ownership from call shape.
4. Prefer `detached_frame_copy(...)` at public boundaries so pandas 3 Copy-on-Write can use shallow detaches safely.

## Contributor checklist

When touching preprocessing code:

1. If the call site is internal and already owns mutable frames, route to `*_owned` methods.
2. If adding a new public entrypoint, enforce one explicit boundary copy and document it.
3. Avoid adding `DataFrame.copy(...)` in service internals unless the copy is required for correctness.
4. Keep ownership signalling explicit (`from_owned`, owned-path method names, owned-frame markers).

## Benchmark and regression guard

Use the preprocessing benchmark harness to measure and enforce this model:

```bash
python benchmarks/measure_preprocessing_copy_churn.py --check
```

The harness benchmarks:

- public safe path (`CoreProcessor.process`)
- owned fast path (`CoreProcessor.process_owned`)
- large-matrix owned preprocessing

`--check` fails if the copy-churn guard expectations are violated (for example, if the public path no longer shows the expected boundary copy overhead or copy budgets are exceeded).

Useful local variants:

```bash
python benchmarks/measure_preprocessing_copy_churn.py
python benchmarks/measure_preprocessing_copy_churn.py --repeats 1 --stdout-only --check
```
