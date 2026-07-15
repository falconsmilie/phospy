# Contributing to PhosPy

Thanks for helping improve PhosPy.

## Before You Open a Pull Request

Please make sure your change is:

- accurate for the current code
- clear for scientists and first-time users
- small enough to review without detective work
- honest about the supported public lane

## Documentation Expectations

When editing docs or examples:

- keep the beginner path small and runnable
- avoid repeating the same rule in many places
- keep top-level user guides in `docs/`, workflow API pages in `docs/api/`,
  architecture decisions in `docs/adr/`, and testing-audit material in
  `docs/testing/`
- use `phospy.api` for request, config, result, enum, reference, and error imports
- keep beginner docs focused on the rat bundled-reference lane unless the page is explicitly advanced

## Local Setup

```bash
pip install -e ".[dev]"
pip install -e ".[dev,parquet]"  # optional parquet support
```

## Tests to Run

For most changes:

```bash
pytest -m "not parity"
```

For public docs and examples:

```bash
pytest tests/unit/test_public_contract_import_routes.py
pytest tests/unit/test_public_examples_contract.py
pytest tests/integration/test_public_examples_smoke.py
```

Run parity tests when scientific logic or fixture-backed behaviour changes:

```bash
pytest tests/parity -m parity -s
```

Run full release validation when changing scientific/parity/provenance/performance
behavior or before preparing a release. The authoritative release-gate command
is `make test-release-gate`; default `pytest` is not sufficient for release, and
release tests, reproducibility/golden checks, parity tests, and performance
contracts are release-blocking:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet]"
make test-release-gate
```

## Style

Use Ruff for linting and formatting:

```bash
ruff check --fix
ruff format
```

For a deeper architecture rationale, review
[ADR 0002: Internal Workflow Architecture for PhosPy](adr/adr_0002_internal_workflow_architecture.md).

## Type Checking

Run Pyright locally with the same scope used in CI:

```bash
pyright
```

The checked scope is explicitly configured in `pyproject.toml` under
`[tool.pyright].include` and covers:

- `src/phospy/api`
- `src/phospy/errors`
- `src/phospy/frames`
- `src/phospy/io`
- `src/phospy/policies`
- `src/phospy/provenance`
- `src/phospy/release`
- `src/phospy/science`
- `src/phospy/tables`
- `src/phospy/validation`
- `src/phospy/workflows`

Expectations for typing changes:

- `src/phospy/science/datasets/models.py` is listed under
  `[tool.pyright].strict` and is already strict-checked
- prefer precise public/model/protocol boundary types over `Any`
- keep suppressions narrow and local
- avoid broad ignores and config-wide suppression

Suppression policy:

> Avoid suppressions by default. Use them only when Pyright cannot model correct runtime behaviour. Every suppression must be narrow, error-code-specific, commented, and justified by tests where practical.
