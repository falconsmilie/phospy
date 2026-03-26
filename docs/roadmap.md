# Roadmap

This roadmap explains where PhosPy is most likely to go after 1.0.0. It is a guide to direction, not a release
promise. The aim is to show the most natural next steps from the current codebase without pretending that every idea is
scheduled.

## First Expansion Areas

The first expansion areas are the ones that sit closest to the supported 1.0.0 API, examples, fixtures, and tests.

### Native `KinaseWorkflow` CLI Support

`KinaseWorkflow` is already part of the supported Python API, but the CLI currently stops at the core preprocessing and
`predMat` path. A natural next step is to expose the native workflow from the command line with explicit inputs for
substrate maps, site sequences, motif sequences, and prediction outputs.

### Broader Native Workflow Parity Evidence

The repository already includes L6 references, prediction traces, and a curated fragile-support dataset. Extending that
evidence to more kinases, more fixture shapes, and more seam-specific assertions would strengthen confidence while
keeping the parity claim narrow and honest.

### Better Diagnostic and Replay Tooling

The current trace-export and comparison scripts are already useful when Python and R diverge. Another likely step is to
make that tooling easier to run, easier to compare, and easier to interpret during native-workflow debugging.

## Likely PhosR-Inspired Ports

These are the areas that make sense once the current release surface has settled and the validation story stays clear.

### Wider Preprocessing Helpers

Version 1.0.0 focuses on the core path from total and phospho tables to corrected phosphosite matrices. A sensible next
wave would be carefully selected ports from the broader PhosR preprocessing stack, especially extra filtering,
transformation, and preparation helpers that fit the current dataset model.

### Richer Downstream Kinase Summaries

The current downstream analysis is centred on `predMat`, weighted activity, KSEA-style summaries, and target counts. A
natural extension would be to broaden that reporting layer with more kinase-centric summaries and more convenient output
surfaces for interpretation and export.

### Smoother Workflow Inputs and Outputs

The native workflow already covers profile construction, motif scoring, score combination, candidate selection, and
adaptive SVM prediction. A practical next step is to make it easier to move between raw tables, validated request
objects, reproducible configuration, and saved output bundles.

## What Is Less Likely in the Near Term

PhosPy is still intentionally narrow. Large package-wide parity claims, a rushed attempt to mirror all of PhosR, or a
wide public API expansion without fixture-backed evidence are all less likely than incremental, well-tested growth.

## How to Read This Roadmap

The items above are best read as direction rather than commitment. In practice, the most credible next work is the work
that keeps PhosPy small, test-backed, and explicit about which seams are validated.
