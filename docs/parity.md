# Parity to PhosR

PhosPy is inspired by `PhosR`, but the parity claim is intentionally narrow.

## What Parity Means Here

Parity in this repository means:

- fixture-backed
- seam-level
- selective

It does not mean:

- full package equivalence with `PhosR`
- every `PhosR` feature is implemented
- every native Python path should match `PhosR` numerically

## Parity Contract Matrix

This page is the single contract map for parity-sensitive fixture families.
It answers the review question: “what does this fixture protect?”

Trace directories such as `prediction_trace/` belong to the fixture family listed below. They are not separate fixture families.

| Fixture family | Dataset origin | Protected workflow or seam | Protected metric or threshold class | Relevant mode or modes | Surface |
| --- | --- | --- | --- | --- | --- |
| `tests/fixtures/r_reference` | Small R-generated reference fixtures | Core preprocessing outputs, site-matrix construction, and downstream wrapper flow | Exact table equality or numeric equality within tight frame-comparison tolerances in parity tests | `default` and `r_parity` where applicable | Internal seam |
| `tests/fixtures/r_reference_l6` | Main R-generated L6 reference dataset | Downstream kinase-analysis outputs, native profile and combined-score seams, candidate-substrate selection, prediction ranking agreement, and replay against committed R sampling traces | Exact or tolerance-based table equality, candidate membership equality, ranking-agreement checks, top-N overlap checks, and replay-path agreement checks | Primarily `r_parity`, with explicit `default` versus `r_parity` comparison where tests protect mode intent | Mixed: internal seam and parity-sensitive workflow behaviour |
| `tests/fixtures/fragile_support_reference` | Curated L6-derived support-screening reference set | Boundary conditions around substrate support, profile construction, combined-score recomputation, and candidate-substrate inclusion cut-offs | Exact candidate-state expectations plus exact or tolerance-based seam recomputation checks | `default` and `r_parity` where applicable | Internal seam |
| `tests/fixtures/r_reference_l6_seam_stress` | Filtered L6-derived seam-stress subset plus filtered R trace tables | Combined-score seam, candidate-substrate boundary selection, and replay of smaller seam-stress sampling paths | Exact or tolerance-based table equality, exact candidate membership equality, and replay-path agreement checks | Primarily `r_parity` | Internal seam |
| `tests/fixtures/synthetic_adaptive_sampling_edge` | Fully synthetic deterministic edge-case data | Adaptive-sampling decision seams such as stable ordering, tiny pools, per-iteration overrides, and deterministic replay decisions | Exact trace-table equality and exact replay-decision checks | `r_parity` | Internal seam |
| `tests/fixtures/public_workflow_reference` | Small committed benchmark outputs generated from the public demos under `examples/` | Public `PredMatWorkflow` and `SignalomeWorkflow` example paths | Exact benchmark-table equality and fixed end-to-end workflow assertions | `default` and `r_parity` | Public workflow |

## What Is Covered

The current parity layer covers selected seams, including:

- core preprocessing outputs
- downstream kinase-analysis outputs
- selected native workflow seams
- selected prediction trace and replay checks
- end-to-end benchmark fixtures for the documented predMat and signalome workflow demos

For fixture generation and trace rebuild commands, see [`fixtures.md`](fixtures.md).

## `KinaseWorkflow` and `svm_mode`

`KinaseWorkflow` is a supported public API, but it is still a native Python workflow.

Use:

- `svm_mode="default"` for the recommended stable native path
- `svm_mode="r_parity"` for the closest supported parity-oriented learner, sampling, and final-scoring preset

Using `svm_mode="r_parity"` does not make the full workflow equivalent to `PhosR`, but it is the supported secondary preset used when parity-sensitive prediction checks need the closest supported reference-oriented path.

## Release Thresholds by Mode

The release bar is intentionally different for `default` and `r_parity`.
That is deliberate, because the two presets have different jobs.

### `default`

`default` is the recommended stable native mode.
A release is acceptable for `default` when it continues to satisfy all of the following:

- the documented public workflow benchmarks still match their committed fixtures in
  `tests/fixtures/public_workflow_reference`
- the L6 ranking-parity floor in `tests/test_parity-with_metrics.py` still holds:
  - mean Spearman rank agreement `>= 0.96`
  - mean top-20 overlap `>= 0.85`
  - mean top-30 overlap `>= 0.88`
  - kinases with top-10 overlap of at least 70%: `>= 20`
- non-parity tests and parity tests both continue to pass

`default` is **not** required to match the replay-trace parity bar used to justify
`r_parity` as a separate preset.

### `r_parity`

`r_parity` is the supported parity-oriented preset.
A release is acceptable for `r_parity` only when it continues to satisfy all of the following:

- the documented public workflow benchmarks still match their committed fixtures in
  `tests/fixtures/public_workflow_reference`
- on the L6 ranking benchmark, it remains at least as strong as `default` for:
  - mean Spearman rank agreement
  - mean top-10 overlap
  - mean top-20 overlap
  - mean top-30 overlap
- on the same ranking benchmark, mean top-10 overlap remains `>= 0.82`
- on the replayed L6 sampling-trace benchmark, it continues to meet the protected replay floor:
  - initial negative rows: exact match
  - iteration sample rows: exact match
  - iteration decision class-1 Pearson correlation `>= 0.999999`
  - iteration decision mean absolute difference `<= 1e-12`
  - iteration probability class-1 Pearson correlation `>= 0.998`
  - iteration probability mean absolute difference `<= 0.015`
  - final top-site matches: exact match

These are the release thresholds already enforced in parity tests. This page makes them
explicit so review does not depend on reading assertions in isolation.

### Threshold enforcement map

| Surface | Mode | Fixture family | Required outcome | Enforced by |
| --- | --- | --- | --- | --- |
| Documented `PredMatWorkflow` benchmark | `default`, `r_parity` | `tests/fixtures/public_workflow_reference` | Exact committed benchmark match | `tests/test_end_to_end_parity.py` |
| Documented `SignalomeWorkflow` benchmark | `default`, `r_parity` | `tests/fixtures/public_workflow_reference` | Exact committed benchmark match | `tests/test_end_to_end_parity.py` |
| L6 ranking parity floor | `default` | `tests/fixtures/r_reference_l6` | Meets the explicit ranking floor listed above | `tests/test_parity-with_metrics.py` |
| L6 ranking comparison | `r_parity` | `tests/fixtures/r_reference_l6` | Matches or exceeds `default` on protected ranking metrics and keeps top-10 overlap floor | `tests/test_parity-with_metrics.py` |
| L6 replayed trace fidelity | `r_parity` | `tests/fixtures/r_reference_l6` | Meets the explicit replay floor listed above | `tests/test_parity-with_metrics.py` |

### Parity-sensitive release review

When a change touches prediction policy, sampling, scoring, fixture generation, or the
public workflow examples, release review must check all of the following:

1. `pytest -m "not parity"` passes.
2. `pytest -m parity` passes.
3. `pytest tests/test_readme_smoke.py tests/test_end_to_end_parity.py` passes.
4. If prediction behaviour changed intentionally, regenerate the affected fixtures and explain
   the contract change in the pull request.
5. If prediction policy or sampling changed, regenerate and review the mode-comparison
   benchmark report before release:

   ```bash
   python benchmarks/compare_prediction_modes.py --repeats 1
   ```

6. Confirm the public docs still describe `default` as the recommended native mode and
   `r_parity` as the supported parity-oriented mode.

## Run the Parity Suite

```bash
pytest -m parity
```

Useful variants:

```bash
pytest -m parity -rs
pytest -m parity -vv
pytest -m parity -k l6
```

Make targets:

```bash
make test-parity
make test-seams
```

## Public Support Decision for `r_parity`

PhosPy retains `svm_mode="r_parity"` as a supported public preset.

Use:

- `svm_mode="default"` for the recommended stable native path
- `svm_mode="r_parity"` for the supported parity-oriented learner, sampling, and final-scoring preset

This decision is explicit and tracked in [ADR 0002](adr/0002-r-parity-public-preset.md).
The decision is evidence-driven: `r_parity` remains public because it is the mode used
when parity-sensitive prediction checks need the closest supported reference-oriented path.
It is not a claim of full package equivalence to `PhosR`.

## Benchmark the Public Prediction Modes

Run the reproducible mode-comparison harness from the repository root:

```bash
python benchmarks/compare_prediction_modes.py
```

By default this writes two review artifacts under `benchmarks/reports/latest/`:

- `compare_prediction_modes.json`
- `compare_prediction_modes.md`

The harness uses the selected parity fixture families rather than ad hoc inputs:

- `tests/fixtures/r_reference_l6` for ranking parity and replayed sampling-trace fidelity
- `tests/fixtures/public_workflow_reference` for the documented `PredMatWorkflow` and `SignalomeWorkflow` demo outputs

It records the same metric classes protected in parity tests, plus wall-clock runtime for both `default` and `r_parity`.

Useful variants:

```bash
python benchmarks/compare_prediction_modes.py --repeats 1
python benchmarks/compare_prediction_modes.py --stdout-only
```

## Optional Debug Output

Some parity tests can print extra summaries when `PHOSPY_SHOW_PARITY=1` is enabled.

Available flags:

- `PHOSPY_SHOW_PARITY`
- `PHOSPY_SHOW_PROFILE_CONSTRUCTION`
- `PHOSPY_SHOW_PREDICTION_MODE_COMPARISON`
- `PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON`
