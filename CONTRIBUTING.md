# Contributing

Thank you for your interest in PhosPy.

PhosPy is intentionally narrow. The most helpful contributions keep that scope clear: small changes, good tests,
and straightforward review.

## Get Set Up

Install the package with test and development tooling:

```bash
pip install -e ".[test,dev]"
```

Install the Git hooks:

```bash
pre-commit install
```

Run the usual local checks:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

## Code Style

PhosPy uses Ruff for linting and formatting.

- `ruff check --fix` handles linting, import sorting, and safe auto-fixes
- `ruff format` handles formatting
- `pre-commit` runs the repository checks before commit

The current repository rules live in:

- `.pre-commit-config.yaml`
- `pyproject.toml`

## Tests

PhosPy has two main test layers:

- **non-parity tests** for package behaviour that should stay stable without R
- **parity tests** for fixture-backed comparisons against committed R/PhosR reference outputs

Run them from the repository root:

```bash
pytest -m "not parity"
pytest -m parity
```

Use `pytest -m parity` when you want only the fixture-backed seam checks.

For the short explanation of what parity means, see
[`docs/validation-and-parity.md`](docs/validation-and-parity.md).

For the full parity guide, including metrics flags, cross-platform examples, and sample output, see
[`docs/parity.md`](docs/parity.md).

## Release Gate

The practical release gate is:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

The non-parity suite includes the documented example smoke workflow. The parity suite covers the committed R-backed
seams. The `make test-parity` shortcut enables `PHOSPY_SHOW_PARITY=1` and
`PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON=1`, so it prints the core parity summaries plus the replayed
mode-comparison block rather than every optional metrics block.

## CI Expectations

GitHub Actions currently runs three checks:

- `pre-commit`
- non-parity tests
- parity tests against the committed fixture snapshots

If you intentionally change fixture-producing behaviour, regenerate the affected fixtures, commit them, and explain the
change clearly in the pull request.

## Good First Directions

The public roadmap lives in [`docs/roadmap.md`](docs/roadmap.md). The most useful contributions are the ones that
extend the current supported surface without over-claiming: native-workflow CLI coverage, better seam-level validation,
clearer diagnostics, and carefully scoped PhosR-inspired additions.
