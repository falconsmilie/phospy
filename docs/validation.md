# Validation Guide

PhosPy validates early and loudly so scientific assumptions are not hidden. Most
errors are fixable once you know which boundary rejected the input.

## Dataset Input Rules

`phospho` must be a non-empty numeric pandas `DataFrame` or supported file path.
Rows are phosphosites and columns are samples. The index must use standard site
IDs such as `MAPK14;Y182;`. Missing values are rejected by default.

`site_metadata` must be a non-empty table aligned to `phospho.index`. It must
include non-empty `gene_symbol` and `site` columns. `site_sequence` is optional
but needed for motif-based scoring. `protein_id` is optional for kinase but
required for signalome.

`sample_metadata`, when provided, must align to the phospho sample columns.

`total`, when provided, must be numeric, missing-value-free, and aligned to the
phospho sample columns.

## Site Metadata Conventions

Accepted column aliases are narrow:

| Accepted alias | Normalised column |
| --- | --- |
| `gene_name` | `gene_symbol` |
| `centralized_sequence` | `site_sequence` |

If `gene_symbol` or `site` is missing, the builder can derive them from an index
like `TSC2;S939;`. It does not derive `protein_id`.

## Analysis-Ready Dataset Boundary

A built `AnalysisReadyPhosphoDataset` must have:

- numeric, non-empty, missing-value-free `phospho`
- unique sample columns
- unique site IDs
- `site_metadata.index` exactly matching `phospho.index`
- required non-empty `gene_symbol` and `site`
- `sample_metadata.index` exactly matching `phospho.columns` when provided
- `total.columns` exactly matching `phospho.columns` when provided
- an `Organism` enum value or `None`
- explicit intensity-scale and processing-state metadata

## Preprocessing Rules

Defaults are intentionally strict: no transform, no normalisation, no imputation,
no total-protein correction, and no comparison construction.

Common cross-field checks:

- `subtract_log_total` requires `total` input data.
- `subtract_log_total` requires `intensity_transform.policy="log2"`.
- `sample_metadata_pairs` requires `sample_metadata`.
- site-matrix construction may drop incomplete rows because the public output
  dataset must be complete.
- duplicate-site resolution is recorded in `dataset.preprocessing_report` when it runs.

## Reference Validation

`ReferenceBundle` requires:

- `organism` as an `Organism` enum value
- `kinase_substrate_map` with non-empty `kinase` and `substrate_site`
- `site_sequences` indexed by site ID with non-empty `site_sequence`
- no duplicate `(kinase, substrate_site)` pairs

`ReferencePreset.AUTO` uses `dataset.organism`. In `1.5.0`, bundled runtime
references are rat-only.

## Workflow Validation

### Kinase Workflow

`KinaseWorkflowRequest.dataset` must be an `AnalysisReadyPhosphoDataset`.
References must be compatible with the dataset organism when organism information
is present.

`KinaseScoringConfig.min_substrates` must be at least `2`. The activity stage can
be disabled with `activity_config=None`, which is useful for tiny examples.

### Signalome Workflow

`SignalomeWorkflowRequest.kinase_result` must be a `KinaseWorkflowResult`.
Signalome also requires explicit `protein_id` values for every interpreted site.
Gene-symbol prefixes in site IDs are not treated as protein identity.

Signalome scale guards protect expensive clustering work:

- `max_exact_tree_sites` limits exact tree construction.
- `max_full_candidate_scoring_sites` limits full candidate correlation scoring.
- `candidate_scoring_policy="sampled"` can reduce candidate-scoring cost but still
  needs exact tree construction.

After a successful run, `result.provenance.workflow_parameters["scale_guard"]`
shows exact tree-generation details and candidate-scoring details separately.

## Quick Fix Table

| Error shape | What to check first |
| --- | --- |
| unsupported file format | Use `.csv`, `.tsv`, `.txt`, or `.parquet`; install parquet support for `.parquet`. |
| missing `gene_symbol` or `site` | Add those columns or use site IDs formatted as `GENE;SITE;`. |
| signalome protein identity error | Add non-empty `protein_id` for every site. |
| reference resolution error | Use rat with `AUTO`, or pass an explicit `ReferenceBundle`. |
| total-protein correction error | Provide `total`, set `intensity_transform.policy="log2"`, and configure identity mapping. |
| activity error on a tiny example | Disable activity or provide enough supported substrates. |
| signalome scale error | Reduce sites, use `sampled` candidate scoring where appropriate, or raise guards deliberately. |

## Error Families

Common public exception families are:

- `PhosPyInputError`: file, table, or request input problem
- `PhosPyValidationError`: validated object does not satisfy its contract
- `PhosPyReferenceError`: reference resolution or compatibility problem
- `PhosPyWorkflowError`: workflow boundary or execution problem
