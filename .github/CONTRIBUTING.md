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

`make test-parity` runs `pytest tests/parity -m parity -s` and prints parity
summary output by default (no `PHOSPY_SHOW_*` flags).

Core kinase L6 ranking release gates are enforced in
`tests/parity/test_l6_prediction_parity.py` with explicit hard thresholds
(`Spearman >= 0.96`, `top-20 overlap >= 0.85`, `top-30 overlap >= 0.88`,
`good_top10_count >= 20` on both prediction-matrix and top-k export surfaces).
Informational parity-reporting tests are marked `parity_diagnostic` and are not
threshold-bearing release blockers.

Active parity authority is the suite in `tests/parity/` with fixture
inputs under:

- `tests/fixtures/rewrite_parity/**`
- `tests/fixtures/public_workflow_reference/**`

## CI

GitHub Actions runs:

- `pre-commit`
- non-parity tests
- parity tests against committed fixtures
- build validation for distributable packages

If you intentionally change fixture-producing behaviour, regenerate the affected fixtures and explain the contract change clearly in the pull request.

For the active public signalome workflow fixture family
(`tests/fixtures/public_workflow_reference/signalome_rewrite_l6_*`), use:

```bash
make fixtures-public-workflow-reference
```

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

## Public config surface (current)

PhosPy keeps scientific decisions explicit through public request/config models.
Do not introduce new policy wrapper objects or hidden helper-level defaults.

Use the current public contract from `phospy.api`:

- `DatasetBuildRequest(..., preprocessing_config=DatasetPreprocessingConfig(...))`
- `KinaseWorkflowRequest(..., scoring_config=KinaseScoringConfig(...), prediction_config=KinasePredictionConfig(...), activity_config=KinaseActivityConfig(...) | None)`
- `SignalomeWorkflowRequest(..., config=SignalomeConfig(...))`

Builder preprocessing policy is grouped under `DatasetPreprocessingConfig`:

- `missing_data=DatasetMissingDataConfig(policy="forbid" | "impute_row_median", min_observed_values=...)`
- `total_protein_correction=DatasetTotalProteinCorrectionConfig(policy="none" | "ratio_to_total")`
- `site_matrix=DatasetSiteMatrixConfig(policy="as_input" | "build_from_metadata", duplicate_site_strategy=..., missing_data_policy="drop_any_missing")`
- `comparisons=DatasetComparisonBuildingConfig(policy="none" | "sample_metadata_pairs", sample_group_column=..., pairs=...)`

Workflow-stage scientific controls stay at their owning stage:

- kinase scoring controls stay in `KinaseScoringConfig` (`min_substrates`, `include_diagnostic_scoring_tables`, `profile_missing_value_strategy`)
- kinase prediction controls stay in `KinasePredictionConfig` (`mode`, `adaptive_policy`, `top_k`, `ensemble_size`, `n_iterations`, `random_state`)
- signalome controls stay in `SignalomeConfig` (`network_policy`, `assignment_policy`, `score_preconditioning_policy`, `module_count`, module-selection thresholds)

When contributing examples or new APIs, verify names and fields against
[`docs/api.md`](../docs/api.md) and keep `phospy.api` as the authoritative
contract namespace.
