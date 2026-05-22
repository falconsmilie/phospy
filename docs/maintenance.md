# Maintenance

This page describes the maintainer material.

## Development Setup

```bash
pip install -e ".[dev]"
pip install -e ".[dev,parquet]"  # optional parquet support
```

For CI-aligned dependency resolution:

```bash
pip install -c constraints/ci.txt -e ".[dev,test]"
```

For full release-gate validation (includes reproducibility/golden and performance):

```bash
pip install -c constraints/ci.txt -e ".[dev,test,parquet]"
```

If `make test-release-gate` fails with import errors for optional engines, install
the optional extras above and rerun.

## Common Checks

```bash
ruff check .
ruff format --check .
pyright
pytest -m "not parity"
```

Run parity tests when scientific logic, fixture data, reference handling, or
scoring behaviour changes:

```bash
pytest tests/parity -m parity -s
```

Run performance checks when preprocessing, scoring, prediction, or signalome hot
paths change:

```bash
pytest tests/performance -m performance
```

Run full scientific release validation before tagging a release:

```bash
make test-release-gate
```

The publish pipeline (`.github/workflows/publish.yml`) runs this same release
gate before building and uploading tagged distributions.

## Type Checking

Pyright is the configured type checker. The checked scope is listed in
`pyproject.toml` under `[tool.pyright]` and includes:

- `src/phospy/api`
- `src/phospy/science/datasets`
- `src/phospy/io` (including publishing/export paths)
- `src/phospy/science/prediction`
- `src/phospy/provenance`
- `src/phospy/science/references`
- `src/phospy/science/signalomes/clustering`
- `src/phospy/tables`
- `src/phospy/validation`
- `src/phospy/workflows`

Strict checking is enabled for selected stable scientific/core modules listed
under `[tool.pyright].strict` (for example provenance/reference model modules),
and can be expanded incrementally.

Avoid suppressions by default. Use them only when Pyright cannot model correct
runtime behaviour. Every suppression must be narrow, error-code-specific,
commented, and justified by tests where practical.

## Fixture Policy

Active fixture roots:

- `tests/fixtures/`
- `tests/support/`
- `scripts/active/`

Regeneration scripts should be deterministic and should say which fixture family
they update. Generated benchmark reports belong in `benchmarks/reports/`, which
is ignored by git.

## Documentation Policy

Docs should stay flat, beginner-friendly, and tested against the code. Prefer one
clear beginner path over several overlapping overview pages. Keep examples small
and runnable.

Docs subdirectories are intentional:

- `docs/api/` for workflow API references
- `docs/adr/` for architecture decision records
- `docs/testing/` for testing-audit and consolidation material

## Frame Ownership Policy

PhosPy treats DataFrames as owned mutable state internally.

Input DataFrames are copied when accepted into validated dataset/table objects.
Workflow internals may pass owned DataFrames without repeated defensive copies.
Public result/table access should either return a safe copy or clearly mark the
returned object as an internal-only borrowed reference.

Provenance fingerprints describe the owned internal state at creation time.

Exposure categories:

- `owned_internal`: DataFrames stored in dataset/result/table dataclass fields.
- `safe_public_copy`: `to_dataframe(...)`, `to_pandas(...)`, and
  `*_dataframe(...)` helpers (always defensive snapshots).
- `borrowed_internal_view`: private/internal helpers only (`_borrow_dataframe`,
  `_borrow_optional_dataframe`).
- `export_snapshot`: persisted outputs and provenance fingerprints.

## Release Notes

- Current release notes: [PhosPy Release Notes](release-notes.md)
- Changelog: [`CHANGELOG.md`](https://github.com/falconsmilie/phospy/blob/main/CHANGELOG.md)
- Citation metadata: [`CITATION.cff`](https://github.com/falconsmilie/phospy/blob/main/CITATION.cff)

## ADRs

Architecture and governance decisions live in [ADR Index](adr/index.md). ADRs are
advanced maintainer documents; day-to-day users should start with the
[Quickstart](quickstart.md).

