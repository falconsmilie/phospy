# Roadmap

This roadmap is a guide to direction, not a release promise. It highlights the next steps that fit the current
codebase most naturally.

## Most Likely Next Steps

### Native `KinaseWorkflow` CLI Support

`KinaseWorkflow` is already part of the supported Python API, but the CLI currently stops at the core preprocessing and
`predMat` path. A natural next step is to expose the native workflow from the command line with explicit workflow
inputs and prediction outputs.

### Broader Native Workflow Validation

The repository already includes L6 references, prediction traces, and a curated fragile-support dataset. Extending that
evidence to more kinases, more fixture shapes, and more seam-specific assertions would strengthen confidence while
keeping the parity claim narrow and honest.

### Better Diagnostic Tooling

The current trace-export and comparison scripts are already useful when Python and R diverge. Another likely step is to
make them easier to run, compare, and interpret.

## Likely PhosR-Inspired Ports

### Wider Preprocessing Helpers

Version 1.0.0 focuses on the core path from total and phospho tables to corrected phosphosite matrices. A sensible next
wave would be carefully chosen preprocessing helpers that fit the current dataset model.

### Richer Downstream Kinase Summaries

The current downstream analysis is centred on `predMat`, weighted activity, KSEA-style summaries, and target counts. A
natural extension would be a broader reporting layer for interpretation and export.

### Smoother Workflow Inputs and Outputs

The native workflow already covers profile construction, motif scoring, score combination, candidate selection, and
adaptive SVM prediction. A practical next step is to make it easier to move between raw tables, validated request
objects, reproducible configuration, and saved output bundles.

## Less Likely Near-Term Work

PhosPy is still intentionally narrow. Large package-wide parity claims, a rushed attempt to mirror all of PhosR, or a
wide public API expansion without fixture-backed evidence are all less likely than incremental, well-tested growth.

## Reading This Roadmap

Treat the items above as direction rather than commitment. The most credible next work is the work that keeps PhosPy
small, test-backed, and clear about which seams are validated.
