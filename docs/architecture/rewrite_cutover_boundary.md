# Rewrite Cutover Boundary

This repository now uses a hard package boundary for the rewrite.

> Audience: contributors and maintainers working on code placement and migration boundaries.
> For end-user workflow usage, see [Getting started](../getting-started/index.md).

## Package Layout

```text
src/phospy/                    # new architecture and supported public package
legacy_archive/phospy_legacy/  # old implementation, reference-only during migration
```

## Contributor Rules

- All new architecture and implementation work must land under `src/phospy/`.
- `phospy_legacy` is internal migration reference material only.
- Do not extend old structures in `legacy_archive/phospy_legacy/`.
- No compatibility bridge layer from `phospy` into `phospy_legacy` is maintained by default.
- Temporary functional incompleteness in `src/phospy/` is expected at this stage.

## Where Next

- Maintainer navigation: [Contributor and maintainer docs](../contributor/index.md)
- Architecture context: [Architecture index](index.md)
- Decision rationale: [ADR Index](../adr/index.md)
