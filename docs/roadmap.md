# Roadmap

This page is directional and non-committal; it is not a release promise.

> Audience: advanced users and maintainers tracking project direction and governance follow-ons.
> If you need onboarding, begin with [Getting started](getting-started/index.md).

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

## Supported Science (Current Lane)

- Kinase scoring publishes:
  authoritative `profile_scores` and `combined_scores`, with optional
  diagnostic `motif_scores` and `weights`
- Prediction uses the established downstream score lane:
  `combined_scores` when present, otherwise `profile_scores`
- Adaptive ensemble prediction mode is supported (`mode="adaptive_ensemble"`)
  in the standard package install (no extra dependency step required)
- Kinase activity stage is supported and optional
- Signalome workflow consumes the same downstream score lane and returns module
  assignments, module table, kinase network, and populated
  `expanded_signalome`
- Signalome supports `assignment_policy="weighted_top"` and network-policy
  variants (`positive_only`, `absolute_threshold`, `signed`)

## Science-Parity Ticket Status (2026-04-22)

Audit reference:
[Legacy Science Gap Audit](architecture/legacy_science_gap_audit.md)

- `SCI-GAP-01`: profile missing-value strategy (`strict` + `median_skipna`) -
  completed
- `SCI-GAP-12`: core kinase downstream scoring/prediction parity restoration -
  completed with repaired like-for-like ranking comparison surfaces and active
  ranking closure gates in parity CI
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
of 2026-04-20. Ticket closure labels alone still do not imply full
legacy-science parity outside the audited inventory and explicit coverage tiers.
See the full legacy-science inventory in:
`docs/architecture/legacy_science_gap_audit.md` and `docs/parity.md`.

## Remaining Roadmap (Real Next Steps)

The remaining roadmap focuses on sustaining parity governance quality rather
than reopening repaired ranking-surface work.

1. Keep governance truth sources synchronized in the same change window:
   `docs/architecture/legacy_science_gap_audit.md`, `docs/parity.md`, and
   release notes/changelog entries.
2. Keep ranking surface contracts explicit in parity metrics helpers so
   closure-grade ranking assertions remain source-consistent and
   policy-consistent.
3. Revisit ranking thresholds only through evidence-backed fixture updates and
   explicit release-governance review.

## Not a Near-Term Goal

- broad parity claims across all legacy lanes
- public expansion of unsupported legacy-facing APIs

## Where Next

- Scientific confidence tiers and gates: [Parity to PhosR](parity.md)
- Detailed evidence inventory: [Legacy science gap audit](architecture/legacy_science_gap_audit.md)
- Current public contract: [API Guide](api.md)
