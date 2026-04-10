# Roadmap

This page is direction, not a release promise.

## Landed Recently

These now exist in the supported public surface:

- a first-class `PredMatWorkflow`
- a first-class `SignalomeWorkflow`
- stable `PredMatResult`, `SignalomeResult`, `SignalomeMapData`, and `SignalomeNetworkData` contracts
- `AnalysisReadyPhosphoDataset` as the preprocessing-to-inference boundary
- `ReferenceBundle`, `ReferenceProvider`, and a first `BundledReferenceProvider` lane for `rat` / `l6_native`
- `SimpleKinaseWorkflow` as the supported common end-to-end kinase inference lane
- separate documentation and examples for the simple lane and the advanced native lane
- explicit `svm_mode` guidance for `default` and `r_parity`
- public examples for predMat, simple workflow, native workflow, and signalome workflows

## Most Likely Next Steps

### Fixture-Backed Hardening of the Simple Lane

The simple workflow lane now exists as a supported public path. The next useful step is deeper fixture-backed regression coverage so the common path stays stable and reviewable as the package grows.

### Native Workflow CLI Support

The Python API supports `KinaseWorkflow`, but the CLI still covers the file-based preprocessing path plus optional `predMat` analysis.

### Better Diagnostics

The next useful improvements are clearer validation messages, easier debugging around overlap failures, and tighter guidance when prediction thresholds are too strict.

### Broader Native Workflow Validation

The repository already includes L6 references, prediction traces, and public workflow benchmarks. Extending that coverage would raise confidence without widening the parity claim.

### Broader Bundled Reference Coverage

The first bundled provider lane is intentionally narrow. Expanding beyond `rat` / `l6_native` needs the same level of explicit provenance, validation, and fixture-backed evidence before it should be presented as supported.

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
