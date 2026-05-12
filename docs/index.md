# PhosPy Documentation

## Welcome

PhosPy helps you turn phosphosite intensity tables into kinase scoring,
kinase prediction, and optional signalome analysis in Python.

## Start Here

1. [Quickstart](quickstart.md): copy/paste a complete first workflow.
2. [API Guide](api/guide.md): see the public Python classes and config options.
3. [CLI Guide](cli.md): run the file-based lane from a terminal.
4. [Validation Guide](validation.md): fix common input and boundary errors.

## What PhosPy Does

PhosPy supports clear public workflow lanes:

1. build an `AnalysisReadyPhosphoDataset`
2. run `DifferentialAnalysisWorkflow` with explicit design and contrasts
3. run `KinaseWorkflow`
4. optionally run `SignalomeWorkflow`

The package does not provide HTTP endpoints. Use the Python API for DataFrame
work and explicit references. Use the CLI for simple file-based runs.

## Page Map

| Need | Page                                                                  |
| --- |-----------------------------------------------------------------------|
| First successful run | [Quickstart](quickstart.md)                                           |
| Public Python classes and parameters | [API Guide](api/guide.md)                                             |
| Terminal commands | [CLI Guide](cli.md)                                                   |
| Input rules and common errors | [Validation Guide](validation.md)                                     |
| Workflow expectations, assumptions, and result interpretation | [Workflow Contracts](workflow_contracts.md)                           |
| Written output files and reloadable bundles | [Output Bundles](output_bundles.md)                                   |
| Runtime limits and larger datasets | [Performance Contracts](performance.md)                               |
| Scientific coverage and PhosR comparison | [Scientific Coverage](scientific-coverage.md) and [Parity](parity.md) |
| Development, CI, fixtures, and release notes | [Maintenance](maintenance.md)                                         |

A good first outcome is modest: build a two-site rat dataset, run kinase with
`ReferencePreset.AUTO`, and add signalome once `protein_id` is present.
