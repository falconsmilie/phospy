# Roadmap

This page is direction, not a release promise.

## Landed Recently

These now exist in the supported public surface:

- a first-class `PredMatWorkflow`
- a first-class `SignalomeWorkflow`
- stable `PredMatResult`, `SignalomeResult`, `SignalomeMapData`, and `SignalomeNetworkData` contracts
- public examples for `predMat` and signalome workflows
- explicit `svm_mode` guidance for `default` and `r_parity`

## Most Likely Next Steps

### Native Workflow CLI Support

The Python API supports `KinaseWorkflow`, but the CLI still covers the file-based preprocessing path plus optional `predMat` analysis.

### Better Diagnostics

The next useful improvements are clearer validation messages, easier debugging around overlap failures, and tighter guidance when prediction thresholds are too strict.

### Broader Native Workflow Validation

The repository already includes L6 references, prediction traces, and public workflow benchmarks. Extending that coverage would raise confidence without widening the parity claim.

### Smoother Workflow I/O

The package would benefit from easier movement between validated in-memory results, reproducible configuration, and saved output bundles.

## Likely PhosR-Inspired Ports

### Site- and Gene-Centric Downstream Analysis

Likely areas include:

- phosphosite-to-gene collapsing
- over-representation pathway analysis
- rank-based pathway enrichment

### Wider Preprocessing Helpers

The current public surface focuses on the main path from total and phospho inputs to corrected phosphosite matrices. Smaller targeted helpers remain a natural next layer.

### Richer Downstream Kinase Summaries

The downstream layer already covers weighted activity, KSEA-style summaries, and target counts. Reporting and export can continue to grow from there.

## Less Likely Near-Term Work

These are still unlikely near-term directions:

- package-wide parity claims
- trying to mirror all of `PhosR`
- broad public API expansion without fixture-backed evidence
- deeper SPS- and RUV-style normalisation work before the narrower public workflows are fully hardened
