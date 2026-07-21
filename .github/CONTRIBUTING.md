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
pytest tests/parity -m "parity and not parity_diagnostic" -s
```

Run full release checks when changing scientific/parity/provenance/performance
behavior or before preparing a release. The maintainer release command is
`make release-check`; default `pytest` is not sufficient for publishing, and
blocking parity tests, performance contracts, release/golden/reproducibility
tests, checked-in reference validation, metadata checks, packaged-reference
checks, and the wheel smoke test are release-blocking. Default pytest `testpaths`
omit `tests/release`, `tests/golden`, and `tests/performance`; the
`test-release-gates` Make target selects release/golden tests explicitly. This
provides normal CI/build confidence, not formal exact-source/exact-artifact
attestation:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet]"
make release-check
```

Before publishing distributions, use `make build`. It starts from an empty
`dist/`, builds one wheel and one sdist, runs metadata checks, and validates the
packaged reference manifests and declared file hashes in both archives.

## Style

Use Ruff for linting and formatting:

```bash
ruff check --fix
ruff format
```
