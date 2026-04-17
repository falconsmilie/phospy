# Roadmap

This page is direction, not a release promise.

## Current Rewrite Contract

- Dataset builder: supported
- Simple kinase workflow: supported
- Signalome workflow: placeholder until implemented

Everything above applies to `src/phospy/` only.
`src/phospy_legacy/` remains migration reference material.

## Landed and Supported

- `AnalysisReadyPhosphoDataset` as the public dataset boundary
- `DatasetBuildRequest` and `AnalysisReadyDatasetBuilder.run(request)` for dataset construction
- `SimpleKinaseWorkflow.run(request)` with nested stage outputs:
  `scoring_result`, `prediction_result`, `activity_result`
- `ReferencePreset`/`ReferenceBundle` routing with rat bundled references
- Rewrite-only examples and tests for the supported route

## Placeholder (Not Yet Implemented)

- `SignalomeWorkflow.run(request)` public shell and request/result model
- Scientific signalome computation and non-empty signalome outputs

## Likely Next Steps

- Implement real signalome workflow internals behind the existing request/result shell
- Add broader fixture-backed regression coverage for the simple kinase lane
- Expand bundled references beyond rat only after provenance and validation are in place
- Improve workflow diagnostics and error-action guidance

## Not a Near-term Goal

- Broad parity claims across all legacy lanes
- Public expansion of unsupported legacy-facing APIs
