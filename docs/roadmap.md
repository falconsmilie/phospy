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
  authoritative `profile_scores` and `rank_weighted_fusion_scores`, with optional
  diagnostic `motif_scores` and `weights`
- Prediction uses the established downstream score lane:
  `rank_weighted_fusion_scores` when present, otherwise `profile_scores`
- Adaptive ensemble prediction mode is supported (`mode="adaptive_ensemble"`)
  in the standard package install (no extra dependency step required)
- Kinase activity stage is supported and optional
- Signalome workflow consumes the same downstream score lane and returns module
  assignments, module table, kinase network, and populated
  `expanded_signalome`
- Signalome supports `assignment_policy="weighted_top"` and network-policy
  variants (`positive_only`, `absolute_threshold`, `signed`)

## Remaining Roadmap (Real Next Steps)

The remaining roadmap focuses on sustaining parity governance quality rather
than reopening repaired ranking-surface work.

1. Keep governance truth sources synchronized in the same change window:
   `docs/architecture/science_gap_audit.md`, `docs/parity.md`, and
   release notes/changelog entries.
2. Keep ranking surface contracts explicit in parity metrics helpers so
   closure-grade ranking assertions remain source-consistent and
   policy-consistent.
3. Revisit ranking thresholds only through evidence-backed fixture updates and
   explicit release-governance review.

## Not a Near-Term Goal

- broad parity claims across all legacy lanes
- public expansion of unsupported legacy-facing APIs
