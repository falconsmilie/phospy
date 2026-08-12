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

Run full release checks before preparing a release, and after changes to
scientific/parity/provenance/performance, distribution, reference-bundle, or
public-contract behavior. The maintainer release command is
`make release-check`; default `pytest` is useful for normal contributor
confidence but is not sufficient for publishing:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet,docs]"
make release-check
```

See [Maintenance](../docs/maintenance.md) for the detailed process. Final
release verification requires a Git-backed checkout for staged-byte validation
and verifies the freshly built wheel and sdist; a source-tree test run alone
does not prove the built distributions are valid.

For local scale profiling, use `make benchmark-release-scale`. That optional
50,000 x 48 benchmark is machine-dependent and informational; it is not part of
`make release-check` or CI.

Before publishing distributions, use `make release-check`; do not substitute a
source-tree test pass for the aggregate release command.

## Style

Use Ruff for linting and formatting:

```bash
ruff check --fix
ruff format
```
