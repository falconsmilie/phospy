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
`pyproject.toml` under `[tool.pyright]`. Some modules are marked strict as the
first stronger gate. Avoid broad ignores; prefer narrower types at public and
workflow boundaries.

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

Public boundaries should own their DataFrames unless a clearly internal transfer
has already established ownership. Avoid exposing mutable internals through
public result objects. Internal DTOs may pass owned frames between adjacent stages
when the ownership transfer is obvious and tested.

## Release Notes

- Current release notes: [PhosPy 1.5.0](release-notes-1.5.0.md)
- Changelog: [`CHANGELOG.md`](../CHANGELOG.md)
- Citation metadata: [`CITATION.cff`](../CITATION.cff)

## ADRs

Architecture and governance decisions live in [ADR Index](adr/index.md). ADRs are
advanced maintainer documents; day-to-day users should start with the
[Quickstart](quickstart.md).
