# Parity to the R `PhosR` Package

This page explains what PhosPy means by parity to the R `PhosR` package.

## What the Claim Means

In this repository, a parity claim is narrow.

A seam is called parity-backed only when all three of these are true:

- a committed reference fixture exists for that seam
- an automated test compares the Python output with that fixture
- the claim stays limited to that seam

In practice, parity here is:

- seam-level
- fixture-backed
- intentionally narrow

It does **not** mean:

- the whole package is numerically identical to `PhosR`
- every `PhosR` workflow branch is implemented
- every native Python workflow path should match `PhosR`

## What is Covered

The current parity layer covers selected seams backed by committed fixtures, including:

- core preprocessing outputs
- downstream kinase-analysis outputs
- selected native workflow seam checks
- selected prediction trace and replay checks

For fixture locations, see [`docs/fixtures.md`](fixtures.md).

## `KinaseWorkflow` and `svm_mode`

`KinaseWorkflow` is a supported public API, but it is still a native Python workflow.

Use:

- `svm_mode="default"` for the normal native path
- `svm_mode="r_parity"` for the narrower learner-seam comparison used in parity checks

Example:

```python
from phospy import KinaseWorkflow

native = KinaseWorkflow(svm_mode="default")
comparison = KinaseWorkflow(svm_mode="r_parity")
```

Using `svm_mode="r_parity"` does not make the full workflow equivalent to `PhosR`.

## Run the Parity Tests

From the repository root:

```bash
pytest -m parity
```

Useful variants:

```bash
pytest -m parity -rs
pytest -m parity -vv
pytest -m parity --maxfail=1
pytest -m parity -k l6
```

Make targets:

```bash
make test-parity
make test-seams
```

## Optional Metrics Output

Some parity tests can print extra comparison summaries.

Available environment flags:

- `PHOSPY_SHOW_PARITY`
- `PHOSPY_SHOW_PROFILE_CONSTRUCTION`
- `PHOSPY_SHOW_PREDICTION_MODE_COMPARISON`
- `PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON`

Example:

```bash
PHOSPY_SHOW_PARITY=1 pytest -m parity -s
```

The extra flags only take effect when `PHOSPY_SHOW_PARITY` is enabled.
