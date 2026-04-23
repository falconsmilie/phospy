# PhosPy Documentation

Welcome. The docs are intentionally split into two lanes:

- **new users and scientists**: start with the quickstart and troubleshooting
- **advanced users and maintainers**: use the API, validation, and architecture pages when needed

## Start here

1. [Quickstart: first workflow](getting-started/quickstart-first-workflow.md)
2. [Troubleshooting: first run](getting-started/troubleshooting-first-run.md)
3. [CLI Guide](cli.md) or [API Guide](api.md)

## What PhosPy does

PhosPy supports one clear workflow:

1. build an analysis-ready dataset
2. run kinase scoring and prediction
3. optionally run signalome analysis

For most first runs, use the rat bundled-reference lane:

- `organism=Organism.RAT`
- `references=ReferencePreset.AUTO`

## When to use which page

| Need | Best page |
| --- | --- |
| My first successful run | [Quickstart](getting-started/quickstart-first-workflow.md) |
| My run failed | [Troubleshooting](getting-started/troubleshooting-first-run.md) |
| File-based command line usage | [CLI Guide](cli.md) |
| Python classes, requests, configs, and results | [API Guide](api.md) |
| Exact validation rules | [Validation Guide](validation.md) |
| Saved output layout | [Output Bundles](output_bundles.md) |
| Release framing | [Release Notes 1.5.0](release_notes/1.5.0.md) |

## A good first outcome

After the quickstart, you should be able to:

- prepare a valid `phospho` table and matching `site_metadata`
- build `AnalysisReadyPhosphoDataset`
- run kinase in the supported bundled rat lane
- understand when signalome needs `protein_id`
