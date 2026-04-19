# CLI Guide

The supported public CLI lane is intentionally narrow:

1. Build an analysis-ready dataset from files.
2. Run the kinase workflow from files.
3. Run the signalome workflow from files via dataset -> kinase -> signalome.

## Commands

### Dataset Build

```bash
phospy dataset-build \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --outdir ./out
```

### Kinase Workflow

```bash
phospy kinase \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --reference auto \
  --scoring-min-substrates 2 \
  --prediction-top-k 6 \
  --prediction-ensemble-size 8 \
  --activity-threshold 0.6 \
  --activity-min-substrates 3 \
  --activity-top-n-substrates 20 \
  --outdir ./out
```

`--prediction-mode adaptive_ensemble` is supported in the normal CLI lane and
works after standard package installation with no extra dependency setup.

### Signalome Workflow

```bash
phospy signalome \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --reference auto \
  --scoring-min-substrates 2 \
  --prediction-top-k 6 \
  --prediction-ensemble-size 12 \
  --substrate-support-cutoff 0.5 \
  --network-correlation-threshold 0.5 \
  --network-policy signed \
  --outdir ./out
```

Signalome threshold knobs:

- `--substrate-support-cutoff`: prediction-score cutoff used to select
  kinase-supported substrates.
- `--network-correlation-threshold`: correlation threshold used by
  `--network-policy` to keep network edges.
- `--network-policy`:
  - `positive_only`: positive correlations only
  - `absolute_threshold`: absolute-threshold inclusion with unsigned edge values
  - `signed`: absolute-threshold inclusion with signed edge values

Kinase scoring support floor:

- `--scoring-min-substrates` defaults to `2`.
- single-substrate kinase profiles are rejected (`min_substrates` must be `>= 2`).

Bundled preset policy in the current public lane:

- bundled runtime coverage is rat-only
- `--reference auto` and `--reference rat` are the supported bundled paths
- `--reference human` / `--reference mouse` fail with
  `UnsupportedOrganismError` until those lanes have provenance and fixture-backed
  validation
- transformation state is established by PhosPy's supported transformer path
  during dataset build; no user-declared transformation flag is part of this lane

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

`kinase` writes:

```text
<outdir>/
  dataset/
    ...
  kinase/
    scoring/
      profile_scores.<fmt>
      motif_scores.<fmt>      # optional by contract
      combined_scores.<fmt>   # optional by contract
      weights.<fmt>           # optional by contract
    prediction/
      pred_mat.<fmt>
      substrate_list.<fmt>    # optional
    activity/
      weighted_activity.<fmt> # optional
      ksea_scores.<fmt>       # optional
      ksea_counts.<fmt>       # optional
      target_counts.<fmt>     # optional
      target_table.<fmt>      # optional
    references/
      kinase_substrate_map.<fmt>
      site_sequences.<fmt>
    manifest.json
```

`signalome` writes:

```text
<outdir>/
  dataset/
    ...
  kinase/
    ...
  signalome/
    module_assignments.<fmt>
    signalome_modules.<fmt>
    kinase_network_nodes.<fmt>
    kinase_network_edges.<fmt>
    expanded_signalome.<fmt>   # optional
    manifest.json
```

`<fmt>` is selected with `--output-format` (`csv`, `tsv`, or `parquet`).

Notes:

- In the default supported kinase lane, scoring emits
  `profile_scores` and `combined_scores`.
- Diagnostic scoring tables (`motif_scores`, `weights`) are emitted only when
  `include_diagnostic_scoring_tables=True` is set in the scoring config.
- Optional entries above are contract-level optional fields that may be absent in
  custom/future lanes or when an upstream optional stage is disabled.
