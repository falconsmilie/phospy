# Fixtures

This page documents fixture roots that are part of the current supported test
and parity contract.

## Active Fixture Roots

- `tests/fixtures/rewrite_parity/`
- `tests/fixtures/public_workflow_reference/`

These paths are the default source for active unit/parity/integration lanes.

## Regeneration Paths

Supported fixture generators:

- `scripts/active/generate_r_l6_fixtures.R`
- `scripts/active/generate_signalome_public_workflow_reference.py`

Default bootstrap:

```bash
make fixtures-all
```

## Policy

- Fixture data needed for active regression gates stays in-tree.
- Historical-only fixture artifacts are not retained in live archive trees.
- Historical context belongs in commit history and focused provenance notes
  inside active fixture directories.

## Where Next

- Parity governance and confidence tiers: [Parity to PhosR](parity.md)
- Maintainer navigation hub: [Contributor and maintainer docs](contributor/index.md)
