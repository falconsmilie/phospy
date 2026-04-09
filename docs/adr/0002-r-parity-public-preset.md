# ADR 0002: Public support decision for the `r_parity` preset

- Status: Accepted
- Date: 2026-04-09

## Context

PhosPy now exposes two public prediction presets through `svm_mode`:

- `default`
- `r_parity`

The codebase and parity tests already treat `r_parity` as intentional, but until now
there has not been one explicit decision record stating why the preset exists, what it
optimises for, and what trade-offs are accepted relative to `default`.

That ambiguity creates three problems:

- users cannot tell whether `r_parity` is a supported public choice or a temporary
  implementation detail
- contributors cannot tell what level of parity improvement is required to justify a
  second public preset
- future changes to prediction, sampling, and final-score aggregation are harder to
  review because the public support expectation is implicit instead of explicit

PhosPy already has fixture-backed parity coverage and a reproducible benchmark harness.
The benchmark command is:

```bash
python benchmarks/compare_prediction_modes.py
```

That harness compares the two public presets on the selected parity fixture families and
produces local review artifacts under `benchmarks/reports/latest/`.

## Decision

PhosPy retains `svm_mode="r_parity"` as a supported public preset.

### 1. Role of `default`

`default` remains the recommended stable native prediction preset.

It is the mode users should choose when they want the normal PhosPy path without
optimising specifically for the closest supported alignment to the protected R-backed
prediction seams.

### 2. Role of `r_parity`

`r_parity` remains the supported parity-oriented preset.

It exists for users and maintainers who need the closest supported alignment to the
protected learner, sampling, replay, and final-scoring contracts covered by the parity
fixtures. It is not a claim of full package equivalence to `PhosR`.

### 3. Benchmark evidence reviewed

The decision is grounded in the benchmark harness and the existing parity test suite.
The reviewed evidence source is the generated report from:

```bash
python benchmarks/compare_prediction_modes.py
```

The benchmark compares `default` and `r_parity` on:

- `tests/fixtures/r_reference_l6` for prediction ranking agreement
- `tests/fixtures/r_reference_l6` for replayed sampling-trace fidelity
- `tests/fixtures/public_workflow_reference` for the public `PredMatWorkflow` and
  `SignalomeWorkflow` demo outputs

The protected metric classes include:

- mean Spearman rank agreement
- top-N overlap
- exact or tolerance-based replay fidelity checks
- exact public workflow benchmark equality
- wall-clock runtime

### 4. Runtime trade-off

A separate public preset is only justified if its parity benefit is real enough to offset
its additional review and maintenance cost.

PhosPy therefore accepts that `r_parity` may be slower than `default`, provided that:

- the slower path materially improves or preserves the protected parity evidence
- the public workflow benchmarks still pass
- the extra preset remains narrowly documented and fixture-backed

### 5. Threshold for keeping a separate public preset

`r_parity` remains public only while all of the following stay true:

- it continues to meet the protected parity thresholds defined for the parity benchmark
  and parity tests
- it remains at least as strong as `default` on the parity-sensitive evidence that the
  preset is intended to protect
- it does not break the documented public workflow benchmarks
- its intent remains narrow, explicit, and maintainable

If future benchmark evidence no longer shows a meaningful parity benefit, the project
should demote or remove `r_parity` rather than keeping a second preset by inertia.

### 6. Support expectation going forward

`r_parity` is a supported secondary preset, not the recommended default.

Future changes to prediction policy, sampling, or final-score aggregation must preserve
its documented role and must be reviewed against the parity benchmark evidence. Public
communication should continue to describe:

- `default` as the recommended stable native mode
- `r_parity` as the supported parity-oriented mode

## Consequences

Positive consequences:

- the public meaning of both prediction presets is now explicit
- future benchmark and parity reviews have a named decision target
- contributors have a clear rule for when a second public preset is justified

Trade-offs:

- the project accepts the maintenance cost of a second public preset
- public docs must continue to explain that `r_parity` is narrow and evidence-driven
- future refactors must preserve benchmarkability, not just test pass/fail behaviour

## Follow-on work

- remove `r_parity` if future evidence no longer justifies public support
