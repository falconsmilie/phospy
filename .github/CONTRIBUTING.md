# Contributing

Thanks for your interest in PhosPy.

The project is intentionally narrow. The best contributions are usually small, well-tested, and easy to review.

## Setup

```bash
pip install -e ".[test,dev]"
pre-commit install
```

## Usual Local Checks

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

## Style

PhosPy uses Ruff for linting and formatting.

```bash
ruff check --fix
ruff format
```

Repository rules live in:

- `.pre-commit-config.yaml`
- `pyproject.toml`

## Tests

Use:

```bash
pytest -m "not parity"
pytest -m parity
```

- the non-parity suite covers the normal Python tests
- the parity suite covers fixture-backed seam checks

For parity scope and release thresholds, see [`docs/parity.md`](../docs/parity.md).

## Release Gate

Before a release or substantial PR, run:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
pytest tests/test_readme_smoke.py tests/test_end_to_end_parity.py
```

For parity-sensitive changes to prediction policy, sampling, scoring, or public workflow fixtures, also review the lightweight mode-comparison benchmark:

```bash
python benchmarks/compare_prediction_modes.py --repeats 1
```

`make test-parity` prints the standard parity summary output.

## CI

GitHub Actions currently runs:

- `pre-commit`
- non-parity tests
- parity tests against committed fixtures

If you intentionally change fixture-producing behaviour, regenerate the affected fixtures and explain the change clearly in the pull request.


## Package Layout

PhosPy is organised by domain capability first. Use the package that owns the
behaviour rather than adding new root-level modules.

Start here:

- [`docs/architecture/package-layout.md`](../docs/architecture/package-layout.md)
- [`docs/adr/0004-reorganise-by-domain.md`](../docs/adr/0004-reorganise-by-domain.md)

A few rules matter most:

- keep `phospy.api` thin and orchestration-only
- put preprocessing in `phospy.preprocessing`
- put prediction and `predMat` execution in `phospy.prediction`
- put bundled biological reference handling in `phospy.references`
- keep `phospy.internal` narrow and internal-only

## Good Starting Areas

The public roadmap lives in [`docs/roadmap.md`](../docs/roadmap.md). Good contributions usually improve the supported surface without over-claiming: clearer diagnostics, better validation, tighter docs, and carefully scoped PhosR-inspired additions.


## Scientific policy defaults

PhosPy now treats key scientific preprocessing heuristics as explicit policies.
Do not bury these decisions in helper functions or silent defaults.

Current policy objects are:

- `SiteMatrixPolicy(duplicate_site_strategy="max_mean_signal")`
  - owns duplicate phosphosite collapse during site-matrix creation
  - `"max_mean_signal"` is the current PhosPy behaviour
  - that default is not claimed as exact PhosR parity unless a parity test proves it
- `KinaseProfilePolicy(missing_value_strategy="propagate_any_missing")`
  - owns missing-value handling when aggregating kinase substrate profiles
  - this keeps the current strict PhosPy behaviour and is the closest supported PhosR-style profile seam
- `SignalomeModuleSelectionPolicy(strategy="correlation_thresholds")`
  - owns automatic signalome module-count selection through the explicit `module_selection_strategy` field
  - this is an explicit PhosPy heuristic with diagnostics, not a hidden claim of exact PhosR parity

Keep these policies visible at the public boundary that owns them:

- duplicate-site policy stays with preprocessing / site-matrix construction
- missing-value aggregation policy stays with kinase profile aggregation
- signalome module-selection policy stays with signalome clustering
