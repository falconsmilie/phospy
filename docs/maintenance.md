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

## Type Checking

Pyright is the configured type checker. The checked scope is listed in
`pyproject.toml` under `[tool.pyright]` and includes:

- `src/phospy/api`
- `src/phospy/datasets`
- `src/phospy/io` (including publishing/export paths)
- `src/phospy/prediction`
- `src/phospy/provenance`
- `src/phospy/references`
- `src/phospy/signalomes/clustering`
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

The only docs subdirectory should be `docs/adr/` for decision records.

## Frame Ownership Policy

PhosPy treats DataFrames as owned mutable state internally.

Input DataFrames are copied when accepted into validated dataset/table objects.
Workflow internals may pass owned DataFrames without repeated defensive copies.
Public result/table access should either return a safe copy or clearly mark the
returned object as borrowed and unsafe to mutate.

Provenance fingerprints describe the owned internal state at creation time.

Exposure categories:

- `owned_internal`: DataFrames stored in dataset/result/table dataclass fields.
- `safe_public_copy`: `to_dataframe(...)`, `to_pandas(...)`, and
  `*_dataframe(...)` helpers with `copy=True` (default).
- `borrowed_public_view`: same helpers with `copy=False`; borrowed and unsafe to
  mutate unless intentional owner mutation is desired.
- `export_snapshot`: persisted outputs and provenance fingerprints.

## Release Notes

- Current release notes: [PhosPy 1.5.0](release-notes-1.5.0.md)
- Changelog: [`CHANGELOG.md`](../CHANGELOG.md)
- Citation metadata: [`CITATION.cff`](../CITATION.cff)

## ADRs

Architecture and governance decisions live in [ADR Index](adr/index.md). ADRs are
advanced maintainer documents; day-to-day users should start with the
[Quickstart](quickstart.md).
