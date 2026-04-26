# CLI Guide

The `phospy` CLI is the file-first public lane.

It supports three commands:

- `dataset-build`
- `kinase`
- `signalome`

Use the CLI when you want reproducible file-based runs. Use the Python API when
you need DataFrames, explicit `ReferenceBundle` injection, or full config
control.

## Before you run a command

Check these first:

- your files are `.csv`, `.tsv`, `.txt`, or `.parquet`
- `site_metadata.index` matches `phospho.index`
- `--organism rat` is set for bundled first runs
- `protein_id` is present if you plan to run `phospy signalome`

## Command summary

| Command | What it does |
| --- | --- |
| `phospy dataset-build` | Builds and writes an analysis-ready dataset |
| `phospy kinase` | Builds a dataset, runs kinase workflow, writes outputs |
| `phospy signalome` | Builds a dataset, runs kinase, then signalome, and writes outputs |

## Common arguments

Dataset input arguments:

- `--phospho`
- `--site-metadata`
- `--sample-metadata` (optional)
- `--total` (optional)
- `--organism {human,mouse,rat}`

Output arguments:

- `--outdir`
- `--output-format {csv,tsv,parquet}`

Supported read and write formats:

- `.csv`
- `.tsv`
- `.txt` (tab-separated)
- `.parquet` with `pip install "phospy[parquet]"`

## `dataset-build`

```bash
phospy dataset-build   --phospho ./input/phospho.csv   --site-metadata ./input/site_metadata.csv   --organism rat   --outdir ./out
```

This writes a `dataset/` directory and a short summary of written paths.

## `kinase`

```bash
phospy kinase   --phospho ./input/phospho.csv   --site-metadata ./input/site_metadata.csv   --organism rat   --reference auto   --outdir ./out
```

Additional kinase options:

- `--reference {auto,human,mouse,rat}`
- `--scoring-min-substrates`
- `--prediction-top-k`
- `--prediction-deterministic-max-selected-kinases`
- `--prediction-adaptive-ensemble-runs`
- `--prediction-mode {deterministic_ranking,adaptive_ensemble}`
- `--prediction-adaptive-policy {stable,r_parity}`
- `--prediction-n-iterations`
- `--prediction-random-state`
- `--skip-activity`
- `--activity-threshold`
- `--activity-min-substrates`
- `--activity-top-n-substrates`

Notes:

- bundled runtime references are rat-only in this release
- `--reference auto` is the recommended beginner lane when `--organism rat` is set
- `--skip-activity` disables the activity stage completely
- `--prediction-ensemble-size` is accepted as a legacy alias and maps to both mode-specific prediction-size options
- `--prediction-mode adaptive_ensemble` works in the normal install; no extra dependency lane is required

## `signalome`

```bash
phospy signalome   --phospho ./input/phospho.csv   --site-metadata ./input/site_metadata.csv   --organism rat   --reference auto   --outdir ./out
```

Additional signalome options:

- `--substrate-support-cutoff`
- `--network-correlation-threshold`
- `--network-policy {positive_only,absolute_threshold,signed}`
- `--assignment-policy {cutoff_binary,weighted_top}`
- `--score-preconditioning-policy {allow_and_report,error_on_drop}`

Signalome requires explicit, non-empty `protein_id` values in `site_metadata`.

## Output layout

`dataset-build` writes:

```text
<outdir>/dataset/
  phospho.<fmt>
  site_metadata.<fmt>
  sample_metadata.<fmt>   # optional
  total.<fmt>             # optional
  manifest.json
```

`kinase` also writes:

```text
<outdir>/kinase/
  scoring/
  prediction/
  activity/               # optional
  references/
  manifest.json
```

`signalome` also writes:

```text
<outdir>/signalome/
  module_assignments.<fmt>
  signalome_modules.<fmt>
  kinase_network_nodes.<fmt>   # optional
  kinase_network_edges.<fmt>
  kinase_network_candidate_correlations.<fmt>   # optional traceability sidecar
  expanded_signalome.<fmt>     # optional by contract
  manifest.json
```

`<fmt>` is `csv`, `tsv`, or `parquet`.

## When the CLI is not enough

Use the Python API when you need:

- DataFrame inputs
- `DatasetPreprocessingConfig`
- explicit `ReferenceBundle` injection
- advanced scoring or signalome config beyond the CLI surface

## Where next

- [Quickstart](getting-started/quickstart-first-workflow.md)
- [Troubleshooting](getting-started/troubleshooting-first-run.md)
- [API Guide](api.md)
- [Output Bundles](output_bundles.md)
