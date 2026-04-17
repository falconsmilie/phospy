# Roadmap

This page is direction, not a release promise.

## Current Rewrite Contract

- Dataset builder: supported
- Simple kinase workflow: supported
- Signalome workflow: first real vertical slice implemented

Everything above applies to `src/phospy/` only.
`src/phospy_legacy/` remains migration reference material.

## Landed and Supported

- `AnalysisReadyPhosphoDataset` as the public dataset boundary
- `DatasetBuildRequest` and `AnalysisReadyDatasetBuilder.run(request)` for dataset construction
- `SimpleKinaseWorkflow.run(request)` with nested stage outputs:
  `scoring_result`, `prediction_result`, `activity_result`
- Rewrite CLI support for dataset build and simple kinase workflow from files
- `SignalomeWorkflow.run(request)` with module assignments, signalome modules, and kinase network outputs
- `ReferencePreset`/`ReferenceBundle` routing with rat bundled references
- Rewrite-only examples and tests for the supported route

## Not Yet Implemented in Signalome

- `SignalomeWorkflowResult.expanded_signalome` population in the rewrite path

## Likely Next Steps

- Expand signalome scientific depth beyond the first vertical slice
- Add broader fixture-backed regression coverage for the simple kinase lane
- Expand bundled references beyond rat only after provenance and validation are in place
- Improve workflow diagnostics and error-action guidance

## Not a Near-term Goal

- Broad parity claims across all legacy lanes
- Public expansion of unsupported legacy-facing APIs
