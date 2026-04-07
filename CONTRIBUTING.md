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

```bash
pytest -m "not parity"
pytest -m parity
```

- the non-parity suite covers the normal Python tests
- the parity suite covers fixture-backed seam checks

For parity scope, see [`docs/parity.md`](docs/parity.md).

## Release Gate

Before a release or substantial PR, run:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

`make test-parity` prints the standard parity summary output.

## CI

GitHub Actions currently runs:

- `pre-commit`
- non-parity tests
- parity tests against committed fixtures

If you intentionally change fixture-producing behaviour, regenerate the affected fixtures and explain the change clearly in the pull request.

## Good Starting Areas

The public roadmap lives in [`docs/roadmap.md`](docs/roadmap.md). Good contributions usually improve the supported surface without over-claiming: clearer diagnostics, better validation, tighter docs, and carefully scoped PhosR-inspired additions.
