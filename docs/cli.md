# CLI Guide

The supported rewrite CLI lane is intentionally narrow:

1. Build an analysis-ready dataset from files.
2. Run the simple kinase workflow from files.

## Commands

### Dataset Build

```bash
phospy dataset-build \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --outdir ./out
```

### Simple Kinase Workflow

```bash
phospy simple-kinase \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --reference auto \
  --prediction-top-k 6 \
  --prediction-ensemble-size 8 \
  --outdir ./out
```

Bundled preset policy in the rewrite cutover:

- bundled runtime coverage is rat-only
- `--reference auto` and `--reference rat` are the supported bundled paths
- `--reference human` / `--reference mouse` fail with
  `UnsupportedOrganismError` until those lanes have provenance and fixture-backed
  validation

## Input Formats

Supported file formats for read and write:

- `.csv`
- `.tsv` / `.txt`
- `.parquet`

## Output Layout

`dataset-build` writes:

```text
<outdir>/
  dataset/
    phospho.<fmt>
    site_metadata.<fmt>
    sample_metadata.<fmt>   # optional
    total.<fmt>             # optional
    manifest.json
```

`simple-kinase` writes:

```text
<outdir>/
  dataset/
    ...
  simple_kinase/
    scoring/
      profile_scores.<fmt>
      motif_scores.<fmt>      # optional
      combined_scores.<fmt>   # optional
      weights.<fmt>           # optional
    prediction/
      pred_mat.<fmt>
      substrate_list.<fmt>    # optional
    activity/
      activity_scores.<fmt>   # optional
    references/
      kinase_substrate_map.<fmt>
      site_sequences.<fmt>
    manifest.json
```

`<fmt>` is selected with `--output-format` (`csv`, `tsv`, or `parquet`).
