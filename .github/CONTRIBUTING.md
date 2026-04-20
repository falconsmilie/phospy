# Contributing

Thanks for your interest in PhosPy.

The project is intentionally narrow. The best contributions are usually small, well-tested, and easy to review.

## Setup

```bash
pip install -e ".[test,dev]"
pre-commit install
```

## Usual local checks

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest tests/parity -m parity -s
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

## Release gate

Before a release or substantial pull request, run:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
pytest tests/integration/test_public_examples_smoke.py
pytest tests/parity/test_public_predmat_parity.py -m parity -s
```

For performance-sensitive changes in kinase scoring lanes, also review:

```bash
python benchmarks/measure_kinase_scoring_runtime_alignment.py
```

For signalome/prediction science hot-path changes, also review:

```bash
python benchmarks/measure_signalome_prediction_hot_paths.py
```

For benchmark tooling edits, also run:

```bash
pytest tests/unit/test_benchmark_scripts_smoke.py
```

`make test-parity` runs `pytest tests/parity -m parity -s` and prints the
rewrite-owned parity summary chatter by default (no `PHOSPY_SHOW_*` flags).

Active parity authority is the rewrite suite in `tests/parity/` with fixture
inputs under:

- `tests/fixtures/rewrite_parity/**`
- `tests/fixtures/public_workflow_reference/**`

`tests_legacy/` remains archival/provenance only.

## CI

GitHub Actions runs:

- `pre-commit`
- non-parity tests
- parity tests against committed fixtures
- build validation for distributable packages

If you intentionally change fixture-producing behaviour, regenerate the affected fixtures and explain the contract change clearly in the pull request.

## Package layout

PhosPy is organised by domain capability first. Put new code in the package that owns the behaviour instead of adding new root-level modules.

Start here:

- [`docs/architecture/rewrite_cutover_boundary.md`](../docs/architecture/rewrite_cutover_boundary.md)
- [`docs/adr/adr_0010_internal_package_and_module_layout.md`](../docs/adr/adr_0010_internal_package_and_module_layout.md)
- [`docs/adr/adr_0002_internal_workflow_architecture.md`](../docs/adr/adr_0002_internal_workflow_architecture.md)

Rules that matter most:

- keep `phospy.api` thin and orchestration-only
- put dataset ownership and dataset-shaped builders in `phospy.datasets`
- put transformation-state establishment and transformer logic in `phospy.transformations`
- put prediction and `predMat` execution in `phospy.prediction`
- put kinase activity analysis in `phospy.activities`
- put signalome science domain logic in `phospy.signalomes`
- keep workflow validator/interpreter/executor staging in `phospy.workflows`
- put bundled biological reference handling in `phospy.references`
- keep `phospy.io` for shared file and publishing concerns, not scientific logic
- keep `phospy.validation` focused on internal validation domains
- keep exception taxonomy and user-handleable error classes in `phospy.errors`
- keep `phospy.data` for packaged runtime resource payloads only

## Good starting areas

The public roadmap lives in [`docs/roadmap.md`](../docs/roadmap.md). Good contributions usually improve the supported surface without over-claiming: clearer diagnostics, tighter validation, better docs, and carefully scoped PhosR-inspired additions.

## Scientific policy defaults

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
