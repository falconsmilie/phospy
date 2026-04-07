# Roadmap

This page is direction, not a release promise.

## Most Likely Next Steps

### First-Class `predMat` Workflow

User feedback is strongest around creating `predMat` in a clear, supported way. The next priority is to make that a
first-class workflow rather than a lower-level assembly of scoring and prediction pieces.

This likely includes:

- one obvious public entry point for `predMat` generation
- a canonical result contract for `predMat`
- validation tied to the actual scoring and prediction subset requirements
- a short end-to-end example showing the recommended path

### Native `KinaseWorkflow` CLI Support

The Python API already supports the native workflow. The CLI does not yet expose it.

### Broader Native Workflow Validation

The repository already contains L6 references, prediction traces, and curated seam fixtures. Extending that coverage
would increase confidence without widening the parity claim.

### Better Diagnostic Tooling

The current trace-export scripts are useful, but they could be easier to run and compare.

## Likely PhosR-Inspired Ports

### Signalome Construction and Outputs

After `predMat` becomes a first-class workflow, the most likely next PhosR-inspired feature is signalome support.

This likely includes:

- a dedicated public workflow for signalome construction
- a stable result model for signalome outputs
- signalome map-ready output generation
- kinase-network outputs derived from signalome results
- concise user-facing examples for the recommended workflow

### Site- and Gene-Centric Downstream Analysis

Once `predMat` and signalome workflows are in place, the next likely layer is broader interpretation support, such as:

- phosphosite-to-gene collapsing
- over-representation pathway analysis
- rank-based pathway enrichment

### Wider Preprocessing Helpers

The current release line focuses on the core path from total and phospho inputs to corrected phosphosite matrices. More
targeted preprocessing helpers remain a natural future direction, but they are not the current headline priority.

### Richer Downstream Kinase Summaries

The downstream layer already covers weighted activity, KSEA-style summaries, and target counts. Reporting and export can
continue to grow from there.

### Smoother Workflow Inputs and Outputs

The native workflow would benefit from easier movement between raw tables, validated requests, reproducible
configuration, and saved output bundles. This remains an ongoing direction, especially around prediction and export
contracts.

## Less Likely Near-Term Work

Large package-wide parity claims, a rush to mirror all of PhosR, or a wide public API expansion without fixture-backed
evidence are all unlikely near-term directions.

Deep SPS- and RUV-style normalisation work is also less likely near-term. It is valuable, but it is broader, harder to
validate well, and lower priority than making `predMat` generation and downstream signalome analysis feel native in
Python.