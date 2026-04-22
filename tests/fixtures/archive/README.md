# Fixture Archive

This tree contains committed fixture outputs that are kept for historical
provenance and forensic parity investigation only.

These paths are not part of the active rewrite fixture lane and are excluded
from routine maintainer bootstrap targets such as `make fixtures-all`.

## Current Archived Families

- `adaptive_sampling_edge_trace_debug/`
  - source: legacy synthetic adaptive-sampling seam-debug exports promoted into
    rewrite fixtures during early parity investigation
  - moved from `tests/fixtures/rewrite_parity/adaptive_sampling_edge/` on
    2026-04-22
- `adaptive_sampling_replay_trace_debug/`
  - source: donor replay trace debug exports from the legacy L6 prediction
    trace family
  - moved from `tests/fixtures/rewrite_parity/adaptive_sampling_replay/` on
    2026-04-22

## Maintainer Note

Keep active parity gates and active fixture regeneration centered on:

- `tests/fixtures/rewrite_parity/`
- `tests/fixtures/public_workflow_reference/`

Archive material should only be touched for explicit provenance/audit work.
