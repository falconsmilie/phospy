# Roadmap

This page is direction, not a release promise.

## Landed Recently

These now exist in the supported public surface:

- `PhosphoDataset` as the supported preprocessing entry point
- a first-class `SimpleKinaseWorkflow`
- a first-class `SignalomeWorkflow`
- stable `PredMatResult`, `SignalomeResult`, `SignalomeMapData`, and `SignalomeNetworkData` contracts
- `AnalysisReadyPhosphoDataset` as the preprocessing-to-inference boundary
- `ReferenceBundle`, `ReferenceProvider`, and a first `BundledReferenceProvider` lane for `rat` / `l6_native`
- separate documentation and examples for preprocessing, kinase scoring, and signalome lanes
- explicit `svm_mode` guidance for `default` and `r_parity`
- public examples for simple workflow and signalome workflows

## Most Likely Next Steps

### Fixture-Backed Hardening of the Simple Lane

The simple workflow lane now exists as a supported public path. The next useful step is deeper fixture-backed regression coverage so the common path stays stable and reviewable as the package grows.

### CLI Coverage

The CLI covers file-based preprocessing plus optional kinase activity analysis from an existing `predMat`.

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
