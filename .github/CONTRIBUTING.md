# Contributing to PhosPy

Thanks for helping improve PhosPy.

## Before you open a pull request

Please make sure your change is:

- accurate for the current code
- clear for scientists and first-time users
- small enough to review without detective work
- honest about the supported public lane

## Documentation expectations

When editing docs or examples:

- prefer one clear beginner path over several overlapping overview pages
- avoid repeating the same rule in many places
- keep examples runnable
- use `phospy.api` for request, config, result, enum, reference, and error imports
- keep beginner docs focused on the rat bundled-reference lane unless the page is explicitly advanced

## Local setup

```bash
pip install -e ".[dev]"
pip install -e ".[dev,parquet]"  # optional parquet support
```

## Tests to run

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

Run parity tests when your change touches scientific logic or fixture-backed behaviour.

## Pull request tips

A good pull request usually includes:

- a short problem statement
- the public effect of the change
- any contract or example updates needed to keep docs aligned

## Style

Use Ruff for linting and formatting:

```bash
ruff check --fix
ruff format
```
