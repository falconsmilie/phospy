# Parity to PhosR

PhosPy parity is intentionally narrow and fixture-backed. Passing a parity test
for one lane does not mean the whole PhosR package is implemented.

## What Parity Means Here

A parity claim must say:

- which input fixture was used
- which output table or metric was compared
- which tolerance or acceptance rule was used
- which PhosPy workflow or stage produced the output

The strongest label is `PARITY_GATED_ACTIVE_SCIENCE`: behaviour protected by
active parity-focused tests in `tests/parity/`.

## Active Parity Areas

Current active parity coverage focuses on kinase scoring and prediction surfaces,
including L6-style fixture-backed ranking checks and selected preprocessing or
activity-stage behaviours where explicit fixtures exist.

Run the parity suite with:

```bash
pytest tests/parity -m parity -s
```

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
