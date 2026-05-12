# Parity to PhosR

PhosPy parity is intentionally narrow and fixture-backed. Passing a parity test
for one lane does not mean the whole PhosR package is implemented.

Scope ownership split:

- [Scientific Coverage](scientific-coverage.md) owns user-facing coverage status,
  intended parity scope, and interpretation limits.
- [Parity](parity.md) owns fixture-level comparison evidence, test locations,
  and comparison mechanics.

## What Parity Means Here

A parity claim must say:

- which input fixture was used
- which output table or metric was compared
- which tolerance or acceptance rule was used
- which PhosPy workflow or stage produced the output

The strongest label is `PARITY_GATED_ACTIVE_SCIENCE`: behaviour protected by
active parity-focused tests in `tests/parity/`.

## Active Parity Areas

Current active parity coverage includes:

- differential phosphorylation (`tests/parity/test_differential_analysis_parity.py`)
- kinase scoring and prediction surfaces (L6/public fixture lanes)
- selected preprocessing and activity-stage behaviours with explicit fixtures
- signalome workflow and clustering backend fixture lanes

Run the parity suite with:

```bash
pytest tests/parity -m parity -s
```

Release decisions should run the full release gate (`make test-release-gate`).
Parity failures in that gate are release-blocking.

Some diagnostic parity tests are informational. Release decisions should use the
threshold-bearing gates and the documented fixture expectations, not visual
inspection alone.

## Fixture Locations

| Purpose | Location |
| --- | --- |
| Parity tests | `tests/parity/` |
| Shared parity helpers | `tests/support/` |
| Public workflow reference fixtures | `tests/fixtures/public_workflow_reference/` |
| Regeneration scripts | `scripts/active/` |

## Open Gaps

Open gaps should be described as open gaps, not as partial equivalence. Common
examples include broader organism-specific bundled references, additional PhosR
workflow surfaces, and any method not protected by fixture-backed comparison.
