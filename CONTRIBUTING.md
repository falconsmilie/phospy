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

Run the usual test commands from the repository root:

```bash
pytest -m "not parity"
pytest -m parity
```

- `pytest -m "not parity"` runs the regular Python test suite
- `pytest -m parity` runs the fixture-backed parity checks

For parity details, see [`docs/parity.md`](docs/parity.md).

## Release Gate

The practical release gate is:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

The `make test-parity` shortcut enables the standard parity summary output.

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
