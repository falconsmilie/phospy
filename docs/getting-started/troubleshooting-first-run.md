# Troubleshooting: first run and supported-lane failures

Use this page when your first PhosPy run fails. It is written for the supported
beginner lane first, not for every advanced edge case.

## Fast sanity check

Before reading anything deeper, confirm these basics:

- `phospho.index` contains site IDs like `GENE;SITE;`
- `site_metadata.index` exactly matches `phospho.index`
- `site_metadata` includes `gene_symbol` and `site`
- you set `organism=Organism.RAT` or `--organism rat` for bundled first runs
- you added `protein_id` only if you are running signalome

## Common failures

| What you saw | Usually means | What to do |
| --- | --- | --- |
| `ModuleNotFoundError: No module named 'phospy'` | Wrong environment | Activate the environment where you installed PhosPy |
| `phospy: command not found` | CLI not on your path | Install the package in the active environment and reopen the shell |
| parquet read/write failure | Optional dependency missing | Install `phospy[parquet]` |
| `ReferencePreset.AUTO requires dataset.organism` | Organism was not set | Set `organism=Organism.RAT` or `--organism rat` |
| human/mouse bundled-reference failure | Bundled references are rat-only | Use an explicit `ReferenceBundle` in Python |
| `seam=kinase.activity.valid_candidates` with tiny toy inputs | Default activity filters are stricter than a 2-site demo matrix | Set `activity_config=None` for first-run toy examples, or lower `activity_config.min_substrates`/`activity_config.threshold` |
| signalome error mentioning `protein_id` | Protein identity is required | Add a non-empty `protein_id` column |
| rows disappeared during site-matrix building | Some rows could not be kept in that preprocessing lane | Check sequence support and chosen preprocessing policy |
| overlap/support boundary error | Dataset and references do not overlap enough | Read the seam details and adjust inputs or references |

## Import or CLI command not available

Examples:

```text
ModuleNotFoundError: No module named 'phospy'
phospy: command not found
```

Fix:

- reinstall in the active environment
- confirm `python -m pip show phospy`
- if needed, run the CLI as `python -m phospy.cli --help`

## Input and file-loading failures

Supported file formats are:

- `.csv`
- `.tsv`
- `.txt` (tab-separated)
- `.parquet` with the optional parquet extra installed

Also check:

- `site_metadata.index` matches `phospho.index`
- `site_metadata` has `gene_symbol` and `site`
- values in `phospho` are numeric

## AUTO reference resolution fails

If you see an error about `ReferencePreset.AUTO` or missing organism, the usual
cause is simple: `dataset.organism` was never set.

Beginner fix:

- Python: `organism=Organism.RAT`
- CLI: `--organism rat`

## Bundled reference organism failure

Bundled runtime references are rat-only in this release.

That means:

- `ReferencePreset.AUTO` is the easiest bundled first-run lane when the dataset organism is rat
- human and mouse work need an explicit `ReferenceBundle`
- enum presence does not mean bundled data ships for every organism

## Signalome fails on `protein_id`

Signalome requires explicit protein identity.

Fix:

- add `site_metadata.protein_id`
- make sure it is non-empty for every interpreted site

Important: a site ID such as `TSC2;S939;` gives site identity, not protein
identity. It is not a fallback for `protein_id`.

## Site-matrix row drops

If you use site-matrix building from metadata, row count can become smaller than
the original metadata table.

This is not always a bug. Common reasons are:

- missing sequence support for some rows
- duplicate-site handling collapsing rows
- the public lane keeping only rows that can end in a strict dataset

Check the input row count against `dataset.phospho.shape[0]` and review the
preprocessing policy you selected.

If site-matrix duplicate handling ran, inspect
`dataset.preprocessing_report.duplicate_site_resolution` and
`dataset.preprocessing_report.metadata_conflicts` to see which source rows were
retained, dropped, or aggregated and whether duplicate rows disagreed on key
metadata fields.

If comparison building ran (`comparisons.policy="sample_metadata_pairs"`),
inspect `dataset.preprocessing_report.comparison_group_stats` and
`dataset.preprocessing_report.comparison_pair_stats` to see the replicate-level
summaries and pairwise values behind each `dataset.comparisons` effect-size
column.

## Kinase or signalome boundary errors

Some workflow failures include seam names, counts, and a `next_action` hint.
Those messages are trying to be useful, not dramatic.

Typical causes:

- too little overlap between dataset sites and reference sites
- too little support after scoring thresholds
- signalome inputs not suitable for downstream network/module stages

Read the error details first. They usually tell you which boundary failed.

## Still stuck?

Use the stricter pages only after this one:

- [Validation Guide](../validation.md)
- [API Guide](../api.md)
- [CLI Guide](../cli.md)
