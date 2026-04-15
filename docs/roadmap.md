# Roadmap

This page is direction, not a release promise.

## Already Landed

These now exist in the supported public surface:

- `PhosphoDataset` as the supported preprocessing entry point
- `SimpleKinaseWorkflow` as the main end-to-end workflow
- `SignalomeWorkflow` as the supported signalome workflow
- stable `PredMatResult`, `SignalomeResult`, `SignalomeMapData`, and `SignalomeNetworkData` contracts
- `AnalysisReadyPhosphoDataset` as the preprocessing-to-inference boundary
- `ReferenceBundle`, `ReferenceProvider`, and a first `BundledReferenceProvider` lane for `rat` / `l6_native`
- separate documentation and examples for preprocessing, kinase scoring, and signalome analysis
- explicit `svm_mode` guidance for `default` and `r_parity`

## Most Likely Next Steps

### Better Diagnostics

Clearer validation messages, better overlap debugging, and more direct guidance when prediction thresholds are too strict.

### Hardening the Simple Workflow Lane

The main end-to-end path exists. The next useful step is deeper fixture-backed regression coverage so it stays stable as the package grows.

### CLI Coverage

The CLI already covers file-based preprocessing and optional kinase activity analysis from an existing `predMat`. A natural next step is stronger test and example coverage around that path.

### Broader Native Workflow Validation

The repository already includes L6 references, prediction traces, and public workflow benchmarks. Extending that coverage would improve confidence without widening the parity claim.

### Broader Bundled Reference Coverage

The first bundled provider lane is intentionally narrow. Expanding beyond `rat` / `l6_native` needs the same level of provenance, validation, and fixture-backed evidence before it should be presented as supported.

### Smoother Workflow I/O

The package would benefit from easier movement between validated in-memory results, reproducible configuration, and saved output bundles.

## Likely PhosR-inspired Additions

- phosphosite-to-gene collapsing
- over-representation pathway analysis
- rank-based pathway enrichment
- smaller targeted preprocessing helpers
- richer kinase summary and export layers

## Less Likely Near-term Work

- package-wide parity claims
- trying to mirror all of `PhosR`
- broad public API expansion without fixture-backed evidence
- deeper SPS- and RUV-style normalisation work before the narrower public workflows are fully hardened
