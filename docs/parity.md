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
