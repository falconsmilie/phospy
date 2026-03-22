# Contributing

## Development setup

Install the package with test and development tooling:

```bash
pip install -e ".[test,dev]"
```

Install the git hooks:

```bash
pre-commit install
```

Run the local quality checks on demand:

```bash
pre-commit run --all-files
pytest
pytest -m parity
```

## Code conventions

This repository uses Ruff for both linting and formatting.

- `ruff check --fix` handles linting, import sorting, and safe auto-fixes.
- `ruff format` handles code formatting.
- `pre-commit` enforces both before commit.

The current repository policy lives in:

- `.pre-commit-config.yaml`
- `pyproject.toml`

## Test layers

### Unit tests

The regular `pytest` run covers package logic that should stay stable regardless of the R reference fixtures.

### Parity tests

The parity layer compares Python outputs against CSV fixtures generated from real R/PhosR runs.

Run it with:

```bash
pytest -m parity
```

See [`docs/parity.md`](docs/parity.md) for the fixture model, scope boundaries, and current claims.

## CI expectations

GitHub Actions runs three quality gates:

- Ruff via `pre-commit`
- unit tests
- parity tests against the committed fixture snapshots

If you change fixture-producing behaviour intentionally, regenerate the affected R fixtures, commit them, and explain the change in the pull request.
