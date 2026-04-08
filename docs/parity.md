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

## What Is Covered

The current parity layer covers selected seams, including:

- core preprocessing outputs
- downstream kinase-analysis outputs
- selected native workflow seams
- selected prediction trace and replay checks
- end-to-end benchmark fixtures for the documented predMat and signalome workflow demos

For fixture locations, see [`fixtures.md`](fixtures.md).

## `KinaseWorkflow` and `svm_mode`

`KinaseWorkflow` is a supported public API, but it is still a native Python workflow.

Use:

- `svm_mode="default"` for the recommended stable native path
- `svm_mode="r_parity"` for the closest supported parity-oriented learner, sampling, and final-scoring preset

Using `svm_mode="r_parity"` does not make the full workflow equivalent to `PhosR`, but it is the preset used when parity-sensitive prediction checks need the closest supported reference-oriented path.

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

## Optional Debug Output

Some parity tests can print extra summaries when `PHOSPY_SHOW_PARITY=1` is enabled.

Available flags:

- `PHOSPY_SHOW_PARITY`
- `PHOSPY_SHOW_PROFILE_CONSTRUCTION`
- `PHOSPY_SHOW_PREDICTION_MODE_COMPARISON`
- `PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON`
