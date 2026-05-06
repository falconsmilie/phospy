# PhosPy Release Notes

## Release Overview

PhosPy keeps the supported public shape clear: build an
`AnalysisReadyPhosphoDataset`, run `KinaseWorkflow`, and optionally run
`SignalomeWorkflow` when protein identifiers are available.

This release strengthens preprocessing provenance, signalome clustering
configuration, activity-method reporting, and performance guardrail
communication.

## Added

- KSEA-style kinase activity inference through `KinaseActivityConfig(method="ksea_zscore")`.
- Grouped public configuration models for dataset preprocessing, kinase scoring and prediction, and signalome workflows.
- Explicit signalome clustering engine selection with `"scipy_hierarchical"` as the production default and `"exact_python"` for reference/debug checks.
- Public output and bundle metadata for activity method identity, quantitative meaning, preprocessing state, and signalome scale guards.

## Changed

- User-facing documentation now describes the current API and CLI option names.
- Signalome candidate scoring is documented separately from tree construction so `candidate_scoring_policy="sampled"` is not mistaken for a general large-input bypass.
- Public output guidance distinguishes simple CLI/publisher directories from reloadable workflow bundles.

## Fixed

- Documentation links now point to the current release notes page.
- The CLI guide now documents the current `--clustering-engine` option.
- Citation metadata now matches the package version in `pyproject.toml`.

## Scientific Scope

Bundled runtime references in this release are rat-only. Human and mouse remain
valid enum values, but they require a caller-supplied `ReferenceBundle` for
workflow execution.

KSEA-style activity output is supported as a PhosPy activity method, but it is
not reported as equivalent to PhosR kinase activity inference.

Next: [Quickstart](quickstart.md) or [API Guide](api.md).
