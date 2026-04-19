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
  `profile_scores`, `motif_scores`, `combined_scores`, `weights`
- Prediction uses the legacy-matching downstream score lane:
  `combined_scores` when present, otherwise `profile_scores`
- Kinase activity stage is supported and optional
- Signalome workflow consumes the same downstream score lane and returns module
  assignments, module table, and kinase network

## Deferred / Experimental / Not Yet Ported

- `SignalomeWorkflowResult.expanded_signalome` population
- Adaptive sampling ensemble prediction lane
- Additional legacy science components tracked in gap-audit follow-ons

## Science Gap Follow-ons

Audit reference:
[Legacy Science Gap Audit](architecture/legacy_science_gap_audit.md)

- `SCI-GAP-01`: profile missing-value strategy lane (`median_skipna` optional)
- `SCI-GAP-05`: adaptive sampling ensemble prediction lane
- `SCI-GAP-06`: signalome clustering and module-count diagnostics
- `SCI-GAP-07`: deterministic tie-metadata hardening
- `SCI-GAP-08`: weighted-top assignment policy and fractional module shares
- `SCI-GAP-09`: signalome network policy expansion
- `SCI-GAP-10`: `expanded_signalome` output population
- `SCI-GAP-11`: activity/KSEA parity fixture lock while prediction lane evolves

## Not a Near-Term Goal

- broad parity claims across all legacy lanes
- public expansion of unsupported legacy-facing APIs
