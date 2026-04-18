# Roadmap

This page is direction, not a release promise.

## Current Rewrite Contract

- Dataset builder: supported
- Kinase workflow: supported
- Signalome workflow: first real vertical slice implemented

Everything above applies to `src/phospy/` only.
`legacy_archive/phospy_legacy/` remains migration reference material.

## Landed and Supported

- `AnalysisReadyPhosphoDataset` as the public dataset boundary
- `DatasetBuildRequest` and `AnalysisReadyDatasetBuilder.run(request)` for dataset construction
- `KinaseWorkflow.run(request)` with nested stage outputs:
  `scoring_result`, `prediction_result`, `activity_result`
- Rewrite CLI support for dataset build and kinase workflow from files
- `SignalomeWorkflow.run(request)` with module assignments, signalome modules, and kinase network outputs
- `ReferencePreset`/`ReferenceBundle` routing with rat bundled references
- Explicit cutover policy: non-rat presets are enum-level lanes only, not bundled
  runtime support; non-rat execution requires caller-supplied `ReferenceBundle`
- External save/load output-bundle services for `KinaseWorkflowResult` with
  versioned manifest (`phospy.io`)
- Rewrite-only examples and tests for the supported route

## Not Yet Implemented in Signalome

- `SignalomeWorkflowResult.expanded_signalome` population in the rewrite path

## Likely Next Steps

- Expand signalome scientific depth beyond the first vertical slice
- Add broader fixture-backed regression coverage for the kinase lane
- Expand bundled references beyond rat only after provenance and validation are in place
- Improve workflow diagnostics and error-action guidance

## Science Gap Audit Follow-ons (P1)

- Audit reference: [Legacy Science Gap Audit](architecture/legacy_science_gap_audit.md)
- `SCI-GAP-01`: profile missing-value strategy lane (`median_skipna` optional; default unchanged)
- `SCI-GAP-02`: motif scoring stage for rewrite kinase workflow (landed in scoring outputs)
- `SCI-GAP-03`: profile+motif weighted combination stage and outputs (landed)
- `SCI-GAP-04`: candidate selection (`top`, threshold, inclusion) before prediction (ported seam helpers)
- `SCI-GAP-05`: adaptive sampling ensemble prediction lane (deferred pending rewrite dependency/reproducibility contract)
- `SCI-GAP-06`: signalome clustering and module-count selection diagnostics
- `SCI-GAP-07`: deterministic tie metadata parity hardening for assignments
- `SCI-GAP-08`: weighted-top assignment policy and fractional module shares
- `SCI-GAP-09`: signalome network policy expansion (`positive_only` and `absolute_threshold`)
- `SCI-GAP-10`: `expanded_signalome` population in rewrite result model
- `SCI-GAP-11`: activity/KSEA parity fixture lock while prediction lane evolves

## Not a Near-term Goal

- Broad parity claims across all legacy lanes
- Public expansion of unsupported legacy-facing APIs
