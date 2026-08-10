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
pytest tests/parity -m "parity and not parity_diagnostic" -s
```

Run full release checks when changing scientific/parity/provenance/performance
behavior or before preparing a release. The maintainer release command is
`make release-check`; default `pytest` is not sufficient for publishing, and
blocking parity tests, performance contracts, release/golden/reproducibility
tests, public-consumer contract tests, strict documentation build, checked-in
reference validation, metadata checks, packaged-reference checks, and installed
wheel/sdist verification are release-blocking. Default pytest `testpaths` omit
`tests/contract`, `tests/release`, `tests/golden`, and `tests/performance`; the
`test-contract` and `test-release-gates` Make targets select those suites
explicitly. This
provides normal CI/build confidence, not formal exact-source/exact-artifact
attestation:

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet,docs]"
make release-check
```

For local scale profiling, use `make benchmark-release-scale`. That optional
50,000 x 48 benchmark is machine-dependent and informational; it is not part of
`make release-check` or CI.

Before publishing distributions, use the documented build command:

```bash
make build
```

It starts from an empty `dist/`, builds one wheel and one sdist, runs metadata
checks, validates the packaged reference manifests and declared file hashes in
both archives, and installs/executes both artifacts outside the checkout. To
inspect only the archives, run `make build`; to rebuild and run the
installed-artifact verifier, run `make verify-installed-distributions`.

## Style

Use Ruff for linting and formatting:

```bash
ruff check --fix
ruff format
```

For a deeper architecture rationale, review
[ADR 0002: Internal Workflow Architecture for PhosPy](adr/adr_0002_internal_workflow_architecture.md).

## Type Checking

Run type checking locally with the same scope and strict-coverage policy used in CI:

```bash
python tools/testing/pyright_strict_coverage.py --check
python scripts/run_pyright.py
```

The checked scope is explicitly configured in `pyproject.toml` under
`[tool.pyright].include` and covers:

- `src/phospy/api`
- `src/phospy/contracts`
- `src/phospy/errors`
- `src/phospy/frames`
- `src/phospy/io`
- `src/phospy/policies`
- `src/phospy/provenance`
- `src/phospy/science`
- `src/phospy/tables`
- `src/phospy/validation`
- `src/phospy/workflows`

Expectations for typing changes:

- `src/phospy/science/datasets/models.py` is listed under
  `[tool.pyright].strict` and is already strict-checked
- declared strict paths must exist and remain inside `[tool.pyright].include`;
  `python tools/testing/pyright_strict_coverage.py --check` enforces this
- prefer precise public/model/protocol boundary types over `Any`
- keep suppressions narrow and local
- avoid broad ignores and config-wide suppression

Suppression policy:

> Avoid suppressions by default. Use them only when Pyright cannot model correct runtime behaviour. In strict files, every Pyright suppression must use this inline format:
>
> `# pyright: ignore[reportRuleName] - concrete technical rationale`
>
> The rule list must contain one or more explicit Pyright `report...` diagnostic names, and the rationale must explain the concrete typing limitation and why the runtime operation is safe. Blanket `# pyright: ignore`, placeholder reasons such as `TODO`, file-wide diagnostic downgrades such as `# pyright: reportUnknownMemberType=false`, file-wide strictness downgrades, and `[tool.pyright].ignore` entries that intersect effective strict files are rejected by `python tools/testing/pyright_strict_coverage.py --check`.
