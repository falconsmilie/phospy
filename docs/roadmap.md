# Roadmap

This page is directional and non-committal; it is not a release promise.

## Current Supported Product Contract

Supported today in `src/phospy/`:

- `AnalysisReadyPhosphoDataset` as the workflow dataset boundary
- one public dataset builder story with two supported input routes
  (`DataFrame` and file path)
- `KinaseWorkflow.run(KinaseWorkflowRequest(...))`
- `SignalomeWorkflow.run(SignalomeWorkflowRequest(...))`
- nested stage result models (no flattened convenience surface)
- rat bundled-reference runtime lane with `ReferencePreset.AUTO`/`RAT`
- external bundle persistence in `phospy.io`

Historical code under `legacy_archive/phospy_legacy/` is migration reference only.

## Supported Science (Current Lane)

- Kinase scoring publishes:
  authoritative `profile_scores` and `combined_scores`, with optional
  diagnostic `motif_scores` and `weights`
- Prediction uses the legacy-matching downstream score lane:
  `combined_scores` when present, otherwise `profile_scores`
- Adaptive ensemble prediction mode is supported (`mode="adaptive_ensemble"`)
  in the standard package install (no extra dependency step required)
- Kinase activity stage is supported and optional
- Signalome workflow consumes the same downstream score lane and returns module
  assignments, module table, kinase network, and populated
  `expanded_signalome`
- Signalome supports `assignment_policy="weighted_top"` and network-policy
  variants (`positive_only`, `absolute_threshold`, `signed`)

## Science-Parity Ticket Status (2026-04-21)

Audit reference:
[Legacy Science Gap Audit](architecture/legacy_science_gap_audit.md)

- `SCI-GAP-01`: profile missing-value strategy (`strict` + `median_skipna`) -
  completed
- `SCI-GAP-12`: core kinase downstream scoring/prediction parity restoration -
  completed (L6 ranking/candidate/replay gates are enforced in rewrite-owned
  parity lanes)
- `SCI-GAP-05`: adaptive ensemble prediction lane - completed
- `SCI-GAP-06`: signalome clustering + module-count diagnostics - completed
- `SCI-GAP-08`: weighted-top assignment policy + fractional module shares -
  completed
- `SCI-GAP-09`: signalome network policy expansion (`positive_only`,
  `absolute_threshold`, `signed`) - completed
- `SCI-GAP-10`: `expanded_signalome` output population - completed
- `SCI-GAP-11`: activity/KSEA parity lock - completed (regression lock remains
  active)
- `SCI-GAP-07`: deterministic tie-metadata hardening for current lexicographic
  assignment policy - completed

Tracked `SCI-GAP-*` tickets above are closed in the supported rewrite lane as
of 2026-04-20, but this does not imply full legacy-science parity.
See the full legacy-science inventory in:
`docs/architecture/legacy_science_gap_audit.md` and `docs/parity.md`.

## Remaining Roadmap (Real Next Steps)

The remaining roadmap is governance/doc synchronization and parity-maintenance
hygiene for already-landed legacy-science areas.

1. Keep contract docs aligned with landed science status across `docs/roadmap.md`,
   `docs/api.md`, and `docs/parity.md` so deferred wording does not reappear for
   completed lanes.
2. Keep governance truth sources synchronized in the same change window:
   `docs/architecture/legacy_science_gap_audit.md`, parity donor inventory, and
   release notes/changelog entries.
3. Maintain parity-lock hygiene for landed science lanes by updating fixture
   provenance and regression references when parity fixtures are promoted or
   refreshed.

## Not a Near-Term Goal

- broad parity claims across all legacy lanes
- public expansion of unsupported legacy-facing APIs
