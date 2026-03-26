# Contributing

Thank you for your interest in contributing to PhosPy.

PhosPy 1.0.0 is intentionally narrow: it supports the documented preprocessing, `predMat`-driven kinase-analysis, and native `KinaseWorkflow` paths. The clearest way to help is to keep changes small, well-tested, and easy to review.

## Development Setup

Install the package with test and development tooling:

```bash
pip install -e ".[test,dev]"
```

Install the Git hooks:

```bash
pre-commit install
```

Run the local quality checks on demand:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

## Code Conventions

This repository uses Ruff for both linting and formatting.

- `ruff check --fix` handles linting, import sorting, and safe auto-fixes
- `ruff format` handles code formatting
- `pre-commit` enforces both before commit

The current repository policy lives in:

- `.pre-commit-config.yaml`
- `pyproject.toml`

## Test Layers

### Unit Tests

The regular non-parity run covers package logic that should stay stable regardless of the R reference fixtures.

### Parity Tests

The parity layer compares Python outputs against CSV fixtures generated from real R/PhosR runs.

Run it with:

```bash
pytest -m parity
```

See [`docs/parity.md`](docs/parity.md) for the fixture model, scope boundaries, and current parity claims.

## Release Gate

The practical 1.0.0 release gate is:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

The non-parity suite includes the documented example smoke workflow. The parity suite covers the committed R-backed seams.

## CI Expectations

GitHub Actions currently runs three quality gates:

- Ruff via `pre-commit`
- unit tests
- parity tests against the committed fixture snapshots

If you intentionally change fixture-producing behaviour, regenerate the affected R fixtures, commit them, and explain
the change clearly in the pull request.