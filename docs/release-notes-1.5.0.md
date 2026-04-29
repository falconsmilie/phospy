# PhosPy 1.5.0 Release Notes

Released: 2026-04-22

## Release Overview

PhosPy `1.5.0` clarifies the supported public product shape: build an
analysis-ready dataset, run kinase scoring/prediction, and optionally run
signalome analysis.

## Supported Public Shape

- top-level `phospy` exposes only the main workflow entrypoints
- `phospy.api` exposes requests, configs, results, enums, references, and public
  exceptions
- all public workflow classes use `run(request)`
- the CLI supports file-based `dataset-build`, `kinase`, and `signalome` commands

## Scientific Scope

- Rat bundled runtime references are supported in this release.
- `ReferencePreset.AUTO` resolves from `dataset.organism` and works in the rat
  bundled-reference lane.
- Human and mouse work requires an explicit `ReferenceBundle`.
- Signalome requires explicit `site_metadata.protein_id` for every interpreted
  site.
- Activity output is documented as thresholded substrate-mean activity and
  weighted activity, not full KSEA enrichment.

## Documentation Improvements

- Beginner docs were consolidated into one quickstart.
- API and CLI docs now describe the current public parameters directly.
- Output docs now separate simple publisher output from reloadable bundles.
- Performance docs now call out signalome guard behaviour.

Next: [Quickstart](quickstart.md) or [API Guide](api.md).
