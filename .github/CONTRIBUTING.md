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
- keep `docs/` flat except for `docs/adr/`
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
pytest tests/integration/test_cli_smoke.py
```

Run parity tests when scientific logic or fixture-backed behaviour changes:

```bash
pytest tests/parity -m parity -s
```

Run full release validation when changing scientific/parity/provenance/performance
behavior or before preparing a release. The authoritative release-gate command
is `make test-release-gate`; default `pytest` is not sufficient for release, and
release tests, reproducibility/golden checks, parity tests, and performance
contracts are release-blocking. CI runs this same command on Python 3.10, 3.11,
and 3.12:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet]"
make test-release-gate
```

Before publishing distributions, use `make build`. It validates reference-bundle
files from the actual Git index, builds the wheel and sdist, and checks both
archives against the committed reference manifests. Release-gate duration and
JUnit reports are written under `build/reports/`.

## Style

Use Ruff for linting and formatting:

```bash
ruff check --fix
ruff format
```
