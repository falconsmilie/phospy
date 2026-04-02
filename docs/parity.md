# Parity to the R `PhosR` Package

This page explains the repository's parity claim, what is covered today, and how to run the parity checks.

If you only want the short version, read [`docs/validation-and-parity.md`](validation-and-parity.md) first.

## What Parity Means in this Repository

In PhosPy, parity means parity to the R `PhosR` package for a specific seam backed by committed fixtures.

A parity claim is made only when all three of these are true:

- a committed fixture exists for the seam
- an automated test checks the Python output against that fixture
- the claim stays limited to that seam

Parity here is:

- seam-level
- fixture-backed
- intentionally narrow

It does **not** mean:

- the whole package is numerically identical to `PhosR`
- every `PhosR` branch, option, or edge case is implemented
- every native Python workflow path should match `PhosR`

## What is Covered Today

The current parity layer covers:

- core preprocessing outputs backed by committed R-generated fixtures
- downstream kinase-analysis outputs backed by committed R-generated fixtures
- selected native workflow seams backed by committed L6 reference tables
- selected prediction trace and replay checks backed by committed reference traces

For the fixture and trace layout, see [`docs/fixtures.md`](fixtures.md).

## `KinaseWorkflow` and `svm_mode`

`KinaseWorkflow` is a supported public API, but it is still a native Python workflow.

Use:

- `svm_mode="default"` for the normal native path
- `svm_mode="r_parity"` when you want the narrower learner-seam comparison used in the parity fixtures

Example:

```python
from phospy import KinaseWorkflow

native = KinaseWorkflow(svm_mode="default")
comparison = KinaseWorkflow(svm_mode="r_parity")
```

Using `svm_mode="r_parity"` does not widen the package claim. It only narrows one comparison seam.

## Run the Parity Suite

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

- `make test-parity` runs the parity suite with the standard parity summary output enabled
- `make test-seams` runs the seam-focused parity files only

## Optional Parity Metrics

Some parity tests can print extra comparison summaries.

Available environment flags:

- `PHOSPY_SHOW_PARITY`: master switch for parity metrics output
- `PHOSPY_SHOW_PROFILE_CONSTRUCTION`: adds the optional profile-construction summary
- `PHOSPY_SHOW_PREDICTION_MODE_COMPARISON`: compares default mode with `r_parity`
- `PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON`: adds replayed prediction comparison metrics

Example:

```bash
PHOSPY_SHOW_PARITY=1 pytest -m parity -s
```

The more specific flags do nothing unless `PHOSPY_SHOW_PARITY` is also enabled.

## How to Read the Metrics

These summaries help you interpret the seam-level parity claim. They do not widen it.

- profile-construction and score-matrix metrics show how closely numeric tables match the committed reference outputs
- prediction metrics focus on rank agreement and overlap rather than strict cell-by-cell equality
- replayed trace metrics check how closely Python follows the committed reference traces for selected kinases

## Related Docs

- [`docs/validation-and-parity.md`](validation-and-parity.md) for the short guide
- [`docs/api.md`](api.md) for the public API reference
- [`docs/fixtures.md`](fixtures.md) for fixture locations and trace directories
- [`CONTRIBUTING.md`](../CONTRIBUTING.md) for contributor workflow
