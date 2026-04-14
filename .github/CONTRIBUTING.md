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

## Release Gate

Before a release or substantial PR, run:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
pytest tests/test_readme_smoke.py tests/test_end_to_end_parity.py
```

For parity-sensitive changes to prediction policy, sampling, scoring, or public workflow fixtures, also review:

```bash
python benchmarks/compare_prediction_modes.py --repeats 1
```

`make test-parity` prints the standard parity summary output.

## CI

GitHub Actions runs:

- `pre-commit`
- non-parity tests
- parity tests against committed fixtures

If you intentionally change fixture-producing behaviour, regenerate the affected fixtures and explain the change clearly in the pull request.

## Package Layout

PhosPy is organised by domain capability first. Put new code in the package that owns the behaviour instead of adding new root-level modules.

Start here:

- [`docs/architecture/package-layout.md`](../docs/architecture/package-layout.md)
- [`docs/adr/0004-reorganise-by-domain.md`](../docs/adr/0004-reorganise-by-domain.md)

Rules that matter most:

- keep `phospy.api` thin and orchestration-only
- put dataset ownership and dataset-shaped builders in `phospy.datasets`
- put preprocessing in `phospy.preprocessing`
- put prediction and `predMat` execution in `phospy.prediction`
- put kinase activity analysis in `phospy.activities`
- put bundled biological reference handling in `phospy.references`
- keep `phospy.io` for shared file and publishing concerns, not scientific logic
- keep `phospy.validation` focused on validation entry points and request models
- keep `phospy.internal` narrow and internal-only

## Good Starting Areas

The public roadmap lives in [`docs/roadmap.md`](../docs/roadmap.md). Good contributions usually improve the supported surface without over-claiming: clearer diagnostics, tighter validation, better docs, and carefully scoped PhosR-inspired additions.

## Scientific Policy Defaults

PhosPy keeps key scientific heuristics explicit. Do not bury these decisions in helper functions or silent defaults.

Current policy objects are:

- `SiteMatrixPolicy(duplicate_site_strategy="max_mean_signal", missing_data_policy="drop_any_missing")`
  - owns duplicate phosphosite collapse and missing-data filtering during site-matrix creation
- `KinaseProfilePolicy(missing_value_strategy="propagate_any_missing")`
  - owns missing-value handling during kinase substrate profile aggregation
- `SignalomeModuleSelectionPolicy(strategy="correlation_thresholds")`
  - owns automatic signalome module-count selection through the explicit `module_selection_strategy` field

Keep these policies visible at the public boundary that owns them:

- duplicate-site policy stays with preprocessing and site-matrix construction
- missing-value policy stays with kinase profile aggregation
- signalome module-selection policy stays with signalome clustering
