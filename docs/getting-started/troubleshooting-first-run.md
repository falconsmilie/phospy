# Troubleshooting: first-run and supported-lane failures

Start here when your first PhosPy run fails, or when a supported dataset -> kinase -> signalome run stops at a contract boundary.

This page is organised by symptom:

- what you saw
- what it usually means
- what to check next
- where to go for deeper contract detail

Use it before diving into the full [Validation Guide](../validation.md).

## Fast Sanity Check

Before chasing a deeper bug, confirm these basics:

- you installed `phospy` into the same Python environment you are using
- your phospho row IDs look like `GENE;SITE;` such as `TSC2;S939;`
- `site_metadata.index` lines up exactly with `phospho.index`
- you are using `organism=Organism.RAT` or `--organism rat` for the bundled first-run lane
- you added `protein_id` only if you are running signalome

## Confirm You Are in the Supported Lane

The quickest supported first run in 1.5.0 is:

1. install `phospy` (or `phospy[parquet]` if you need parquet)
2. build a dataset with `organism=Organism.RAT` or `--organism rat`
3. run kinase with `ReferencePreset.AUTO` or `--reference auto`
4. run signalome only if `site_metadata.protein_id` is present and non-empty

Also remember:

- bundled runtime references are rat-only in this release
- human and mouse runs need an explicit `ReferenceBundle` in Python
- the CLI is intentionally file-first and does not expose the full Python API surface

## Jump to the Symptom You Saw

| What you saw | Go to |
| --- | --- |
| Import fails or the `phospy` command is missing | [Import or CLI command is not available](#import-or-cli-command-is-not-available) |
| File/path/format/parquet error | [Input and file-loading failures](#input-and-file-loading-failures) |
| `ReferencePreset.AUTO requires dataset.organism` | [AUTO reference resolution fails](#auto-reference-resolution-fails) |
| Human or mouse bundled-reference failure | [Bundled reference organism failure](#bundled-reference-organism-failure) |
| Signalome complains about `protein_id` | [Signalome fails on protein_id](#signalome-fails-on-protein_id) |
| Dataset row count drops during site-matrix preprocessing | [Site-matrix row drops](#site-matrix-row-drops) |
| Kinase workflow fails with overlap/support seam details | [Kinase overlap or support boundary failure](#kinase-overlap-or-support-boundary-failure) |
| Signalome fails later with support/network seam details | [Signalome support or network boundary failure](#signalome-support-or-network-boundary-failure) |

## Import or CLI Command Is Not Available

### What you saw

Examples:

```text
ModuleNotFoundError: No module named 'phospy'
phospy: command not found
```

### What it usually means

You are not running in the environment where PhosPy was installed, or you followed local-clone install instructions for a situation that only needed a normal package install.

### What to check next

- For normal use, install with `pip install phospy`.
- If you need parquet input/output, install with `pip install "phospy[parquet]"`.
- For a local clone, use `pip install -e ".[dev]"` or `pip install -e ".[dev,parquet]"`.
- Make sure the Python interpreter that runs your script is the same one where you installed the package.
- Make sure the shell that runs `phospy ...` can see the same environment.

### Where to go deeper

- [Quickstart: first workflow](quickstart-first-workflow.md)
- [CLI Guide](../cli.md)

## Input and File-Loading Failures

### What you saw

Common public error text includes:

```text
input file does not exist: ...
unsupported table file format for '...'. supported formats: csv (.csv), tsv (.tsv), txt as tab-separated tsv (.txt), parquet (.parquet)
parquet input requires optional parquet dependencies (for example pyarrow)
failed to parse table input '...': ...
```

### What it usually means

The CLI or builder received a missing path, an unsupported suffix, a file that does not parse as a table, or a parquet file without the optional parquet dependency installed.

### What to check next

- Confirm the file path is correct.
- Use only `.csv`, `.tsv`, `.txt`, or `.parquet`.
- Treat `.txt` as tab-separated input, not arbitrary plain text.
- Install the parquet extra before using `.parquet`:
  `pip install "phospy[parquet]"`.
- If you are using the Python API, remember the builder accepts either a `DataFrame` or a supported file path at the request boundary.

### Where to go deeper

- [CLI Guide](../cli.md#input-formats)
- [Validation Guide: Builder flexibility vs dataset strictness](../validation.md#builder-flexibility-vs-dataset-strictness)

## AUTO Reference Resolution Fails

### What you saw

```text
ReferencePreset.AUTO requires dataset.organism
```

### What it usually means

You asked PhosPy to resolve bundled references from the dataset organism, but the dataset does not carry an organism.

### What to check next

- In Python, set `organism` when building the dataset:
  `DatasetBuildRequest(..., organism=Organism.RAT)`.
- In the CLI, pass `--organism rat` before using `--reference auto`.
- Rebuild the dataset first, then rerun kinase.

### Where to go deeper

- [Quickstart: first workflow](quickstart-first-workflow.md)
- [Validation Guide: Reference validation](../validation.md#reference-validation)

## Bundled Reference Organism Failure

### What you saw

Common public error text includes:

```text
no bundled references are available for organism 'human' in the current release; supported bundled organisms: rat; non-bundled organism lanes require a caller-supplied ReferenceBundle
```

You may also see:

```text
dataset.organism and requested reference preset must match
references.organism must match dataset.organism when both are present
```

### What it usually means

You are trying to use a bundled lane that is not shipped in 1.5.0, or your dataset organism and reference selection disagree.

### What to check next

- Use rat for the bundled first-run lane.
- For human or mouse work, use the Python API and provide an explicit `ReferenceBundle`.
- Keep `dataset.organism`, requested preset, and explicit bundle organism aligned.
- Do not expect `ReferencePreset.HUMAN` or `ReferencePreset.MOUSE` to load bundled runtime references in this release.

### Where to go deeper

- [Quickstart: first workflow](quickstart-first-workflow.md#2-know-required-data-and-reference-scope)
- [Validation Guide: Reference validation](../validation.md#reference-validation)
- [API Guide](../api.md)

## Signalome Fails on protein_id

### What you saw

Common public error text includes:

```text
site_metadata is missing required columns: protein_id
site_metadata.protein_id must contain non-empty string values
```

The public error text also makes the contract explicit:

```text
Supported signalome execution requires explicit site_metadata.protein_id; gene-symbol site-ID prefixes are not a protein-identity substitute.
```

### What it usually means

Signalome needs explicit protein identity for the interpreted dataset sites. A gene-symbol-prefixed site ID such as `"MAPK14;T180;"` is not enough.

### What to check next

- Add a `protein_id` column to `site_metadata` before dataset build.
- Make sure every value is a non-empty string.
- Do not rely on `gene_symbol` or the site-ID prefix as a substitute for protein identity.
- Rebuild the dataset, rerun kinase, then rerun signalome.

### Where to go deeper

- [Quickstart: first workflow](quickstart-first-workflow.md#4-optional-run-signalome)
- [Validation Guide: Workflow validation](../validation.md#workflow-validation)
- [Workflow guides](../workflow-guides/index.md)

## Site-Matrix Row Drops

### What you saw

This is often not an exception. Instead, the built dataset contains fewer rows than the source input when you use:

```text
site_matrix.policy='build_from_metadata'
```

You may also find diagnostics on the built phospho matrix:

```python
row_drop_stats = dataset.phospho.attrs.get("site_matrix_row_drop_stats")
print(row_drop_stats)
```

### What it usually means

The supported site-matrix build lane keeps only rows that can be constructed cleanly from metadata. Rows without usable `site_sequence` support are excluded from that path. In the strict public lane, the supported site-matrix missing-data policy is `drop_any_missing`.

### What to check next

- Inspect `dataset.phospho.attrs["site_matrix_row_drop_stats"]`.
- Check whether `site_metadata.site_sequence` is missing, blank, or unresolved for dropped rows.
- Check whether incomplete phospho values are being dropped by `missing_data_policy='drop_any_missing'`.
- Treat the row loss as a preprocessing-policy effect, not a silent bug.

### Where to go deeper

- [Validation Guide: Builder flexibility vs dataset strictness](../validation.md#builder-flexibility-vs-dataset-strictness)
- [Validation Guide: Builder preprocessing policy rules](../validation.md#builder-preprocessing-policy-rules)
- [API Guide](../api.md)

## Kinase Overlap or Support Boundary Failure

### What you saw

Typical public boundary errors include seam details such as:

```text
kinase workflow boundary validation failed at seam=kinase.interpreter.reference_coverage
kinase workflow boundary validation failed at seam=kinase.interpreter.eligible_kinases
kinase workflow boundary validation failed at seam=kinase.interpreter.sequence_support
```

The details usually include counts such as `overlap_sites=0`, `eligible_kinases=0`, or `sequence_supported_sites=0`, plus a `next_action=...` hint.

### What it usually means

One of the core support checks failed:

- the dataset phosphosite IDs do not overlap the resolved reference substrate sites
- overlap exists, but not enough kinases meet `scoring_config.min_substrates`
- resolved `references.site_sequences` does not support any dataset sites for kinase scoring

### What to check next

- Verify dataset site IDs are canonical and look like `"<gene_symbol>;<site>;"`.
- Confirm dataset organism and reference selection are compatible.
- If overlap is shallow, use references with better coverage for the dataset.
- If the failure is about eligible kinases, remember the scientific floor remains `min_substrates >= 2`.
- Read the `seam=...`, count fields, and `next_action=...` text closely. Those details are part of the public recovery path.

### Where to go deeper

- [Validation Guide: Workflow validation](../validation.md#workflow-validation)
- [API Guide](../api.md)

## Signalome Support or Network Boundary Failure

### What you saw

Typical public boundary errors include seam details such as:

```text
seam=signalome.executor.kinase_support
seam=signalome.executor.network
seam=signalome.interpreter.site_alignment
seam=signalome.interpreter.kinase_overlap
```

These errors also include concrete counts and a `next_action=...` hint.

### What it usually means

The interpreted kinase result is not usable for the requested signalome stage. Common causes are:

- no kinases pass `substrate_support_cutoff`
- prediction and downstream score matrices do not align on sites or kinase columns
- downstream score signal is too weak or not variable enough for network construction

### What to check next

- Confirm the kinase result is from the same dataset you intend to analyse.
- Lower `substrate_support_cutoff` only if that still makes scientific sense for your run.
- Inspect `kinase_result.prediction_result.pred_mat` and the downstream score matrix being used.
- Read the boundary counts and `next_action=...` text rather than treating the failure as a generic crash.

### Where to go deeper

- [Validation Guide: Workflow validation](../validation.md#workflow-validation)
- [API Guide](../api.md)

## Still Stuck?

Capture the exact error text, the request or CLI command you used, and which stage failed:

- dataset build
- reference resolution
- kinase
- signalome

Then compare the failure with the deeper contract pages:

- [Validation Guide](../validation.md)
- [API Guide](../api.md)
- [CLI Guide](../cli.md)
