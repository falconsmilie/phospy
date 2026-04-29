# CLI Guide

The `phospy` CLI is the file-based lane for supported public workflows. Use the
Python API when you need DataFrame inputs, custom `ReferenceBundle` objects, or
advanced preprocessing configuration.

## Commands

```bash
phospy dataset-build --phospho phospho.csv --site-metadata site_metadata.csv --outdir out
phospy kinase --phospho phospho.csv --site-metadata site_metadata.csv --organism rat --reference auto --outdir out
phospy signalome --phospho phospho.csv --site-metadata site_metadata.csv --organism rat --reference auto --outdir out
```

Supported input formats are `.csv`, `.tsv`, `.txt` as tab-separated text, and
`.parquet`. Supported output formats are `csv`, `tsv`, and `parquet`.

## Shared Dataset Arguments

| Option | Required | Meaning |
| --- | --- | --- |
| `--phospho` | yes | phosphosite intensity table |
| `--site-metadata` | yes | site metadata table |
| `--sample-metadata` | no | optional sample metadata table |
| `--total` | no | optional total-protein table |
| `--organism` | no | `human`, `mouse`, or `rat` |
| `--outdir` | no | output root; default `phospy-output` |
| `--output-format` | no | `csv`, `tsv`, or `parquet`; default `csv` |

CSV/TSV/TXT inputs are read with the first column as the row index.

## `dataset-build`

Builds and writes an analysis-ready dataset.

```bash
phospy dataset-build \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --outdir ./out
```

Written files include `dataset/phospho.*`, `dataset/site_metadata.*`, optional
metadata tables, and `dataset/manifest.json`.

## `kinase`

Builds the dataset, resolves references, then runs kinase scoring and prediction.

```bash
phospy kinase \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --reference auto \
  --skip-activity \
  --outdir ./out
```

Kinase-specific options:

| Option | Default | Meaning |
| --- | --- | --- |
| `--reference` | `auto` | `auto`, `human`, `mouse`, or `rat` |
| `--scoring-min-substrates` | `2` | minimum quantified substrates per kinase |
| `--prediction-top-k` | `30` | top predicted substrate sites per kinase |
| `--prediction-deterministic-max-selected-kinases` | `10` | retained kinases in deterministic mode |
| `--prediction-adaptive-ensemble-runs` | `10` | ensemble runs in adaptive mode |
| `--prediction-mode` | `deterministic_ranking` | `deterministic_ranking` or `adaptive_ensemble` |
| `--prediction-adaptive-policy` | `stable` | `stable` or `r_parity` |
| `--prediction-n-iterations` | `5` | adaptive sampling iterations |
| `--prediction-random-state` | none | optional random state |
| `--skip-activity` | false | disable activity outputs |
| `--activity-threshold` | `0.6` | prediction-score threshold for activity |
| `--activity-min-substrates` | `3` | minimum selected substrates per kinase |
| `--activity-top-n-substrates` | `20` | top predicted substrates for weighted activity |

`ReferencePreset.AUTO` works with bundled rat references when the dataset
organism is rat. Human and mouse reference presets are enum values, but bundled
runtime data is not shipped for them in `1.5.0`; use the Python API with an
explicit `ReferenceBundle` for those organisms.

## `signalome`

Builds the dataset, runs kinase, then runs signalome. It requires non-empty
`site_metadata.protein_id` values for all interpreted sites.

```bash
phospy signalome \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --reference auto \
  --skip-activity \
  --outdir ./out
```

Signalome-specific options:

| Option | Default | Meaning |
| --- | --- | --- |
| `--substrate-support-cutoff` | `0.5` | prediction support cutoff |
| `--network-correlation-threshold` | `0.5` | edge threshold used by the network policy |
| `--network-policy` | `signed` | `positive_only`, `absolute_threshold`, or `signed` |
| `--assignment-policy` | `cutoff_binary` | `cutoff_binary` or `weighted_top` |
| `--score-preconditioning-policy` | `allow_and_report` | allow/report or reject all-missing score rows |
| `--tree-engine` | `exact` | exact tree construction; only supported value |
| `--candidate-scoring-policy` | `full` | `full` or `sampled` |
| `--max-exact-tree-sites` | `2000` | hard guard for exact tree construction |
| `--max-full-candidate-scoring-sites` | `2000` | hard guard for full candidate scoring |

The CLI does not expose every Python config field. Use the API for explicit
`module_count`, comparison building, total-protein correction, custom references,
or full signalome config control.

## Output Layout

A kinase run writes:

```text
out/
  dataset/
  kinase/
    scoring/
    prediction/
    activity/        # only when activity is enabled
    references/
    manifest.json
```

A signalome run also writes:

```text
out/
  signalome/
    module_assignments.*
    signalome_modules.*
    kinase_network_edges.*
    kinase_network_nodes.*                  # when present
    kinase_network_candidate_correlations.* # when present
    expanded_signalome.*                    # when present
    manifest.json
```
