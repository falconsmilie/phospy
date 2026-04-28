# Validation Guide

This guide explains the main public validation rules.

Use [Troubleshooting](getting-started/troubleshooting-first-run.md) first if you
just want the fastest fix.

## What validation is protecting

PhosPy validation is there to keep the public workflow honest:

- builder inputs must be readable and aligned
- dataset objects must be structurally valid
- references must match the dataset
- workflow requests must contain the right types and usable values

## Builder input rules

`DatasetBuildRequest` supports two public input routes:

- pandas `DataFrame` values
- file paths for supported table formats

Required fields:

- `phospho`
- `site_metadata`

Optional fields:

- `sample_metadata`
- `total`
- `organism`
- `preprocessing_config`

Main checks:

- `phospho` and `site_metadata` must be DataFrames or supported file paths
- `phospho` must be numeric
- `site_metadata.index` must align to `phospho.index`
- required site metadata information must be available through columns or supported derivation

## Site metadata conventions

Supported aliases are intentionally narrow:

- `gene_symbol`: `gene_symbol`, `gene_name`
- `site`: `site`
- `site_sequence`: `site_sequence`, `centralized_sequence`
- `protein_id`: `protein_id`

Unsupported historical aliases are rejected, including:

- `gene`
- `residue`
- `phosphosite`
- `site_position`
- `sequence`
- `protein`

If `gene_symbol` and/or `site` are missing, PhosPy can derive them only from
index values exactly matching `"<gene_symbol>;<site>;"`.

## Final dataset boundary

`AnalysisReadyPhosphoDataset` is the strict workflow-facing boundary.

`AnalysisReadyPhosphoDataset` itself is strict, missing-value-free, and intended
for workflow execution rather than loose exploratory ingestion.

Main expectations:

- `phospho` and `site_metadata` are DataFrames
- site identity is coherent between row IDs and metadata
- required metadata values are non-empty
- intensity scale state is established and coherent
- processing state is explicit and coherent with intensity scale state
- the supported builder lane hands workflows a missing-value-free dataset

## Internal table schemas

Public APIs still expose pandas `DataFrame` tables, but internal workflow
boundaries now use lightweight schema wrappers in `phospy.tables`.

These wrappers are intentionally small:

- they validate one concrete scientific table shape
- they preserve frame ownership rules
- they expose the validated table through `.frame`
- they raise the same domain error type used at that boundary

Current internal wrappers cover core dataset/reference/prediction/activity and
signalome sidecar tables, including:

- `PhosphoIntensityMatrix`, `SiteMetadataTable`, `SampleMetadataTable`, `TotalProteinMatrix`
- `KinaseSubstrateReference`, `SiteSequenceReference`
- `KinaseScoreMatrix`, `KinasePredictionMatrix`
- `ActivityMatrix`, `ActivityCountSeries`, `ActivityTargetTable`
- `SignalomeSiteContext`, `SignalomeProteinSiteContext`

When adding a new internal scientific `DataFrame` contract, prefer adding one
small wrapper in `phospy.tables` instead of passing raw frames across executor
seams.

## Preprocessing rules

`DatasetPreprocessingConfig` groups six policy areas:

- `intensity_transform`
- `normalisation`
- `missing_data`
- `total_protein_correction`
- `site_matrix`
- `comparisons`

Public preprocessing config dataclasses validate local policy/value constraints
at construction time. Unsupported config-local state is rejected immediately.
`DatasetPreprocessingConfigValidator` still validates the request boundary type
and nested config object types; cross-object request requirements remain in
`DatasetBuildRequestValidator`.

Key public rules:

- `intensity_transform.policy="identity"` is the default
- `intensity_transform.policy="log2"` requires non-negative `pseudocount` and values where `value + pseudocount > 0`
- `normalisation.policy="none"` is the default
- `normalisation.policy="median_center"` and `normalisation.policy="quantile"` are explicit opt-in methods
- `missing_data.policy="forbid"` is the strict default
- `missing_data.policy="impute_row_median"` requires `min_observed_values`
- `total_protein_correction.policy="subtract_log_total"` requires:
  - `intensity_transform.policy="log2"`
  - a `total` table aligned to phospho sample columns and site-to-protein mapping
- total-protein correction formula in the supported lane:
  `log2(phospho + pseudocount) - log2(total + pseudocount)`
- `site_matrix.policy="build_from_metadata"` may reduce row count when rows cannot be supported in that lane
- `site_matrix.minimum_observed_values` is internal-only state and must remain `None` in the public lane
- the public builder lane still ends in a missing-value-free `AnalysisReadyPhosphoDataset`
- `comparisons.policy="sample_metadata_pairs"` requires matching `sample_metadata` and a usable sample-group column
- `dataset.preprocessing_report.row_counts` reports stage-level row counts, and `dataset.preprocessing_report.operations` reports stage-level policy/parameter summaries
- `dataset.preprocessing_report.row_audit` reports row-level preprocessing actions (`dropped`, `imputed`, `retained`, `collapsed`, `aggregated`) with reasons and parameter snapshots
- when comparison building runs, `dataset.comparisons` stays the compact workflow matrix and comparison provenance is exposed in `dataset.preprocessing_report.comparison_group_stats` and `dataset.preprocessing_report.comparison_pair_stats`

Not yet supported in this public lane: `knn` imputation, `min_prob` imputation,
and `combat` batch correction.

## Reference validation

Reference rules are simple but strict:

- `ReferencePreset.AUTO` requires `dataset.organism`
- explicit preset and dataset organism must agree when both are set
- explicit `ReferenceBundle.organism` and dataset organism must agree when both are set
- bundled runtime references are rat-only in this release
- `ReferencePreset.HUMAN` and `ReferencePreset.MOUSE` are valid enum values, but they are not bundled runtime lanes here

`ReferenceBundle` itself must contain non-empty, internally consistent tables.

## Workflow validation

### Kinase workflow

`KinaseWorkflowRequest` validation checks:

- `dataset` is `AnalysisReadyPhosphoDataset`
- `references` is `ReferencePreset` or `ReferenceBundle`
- config objects are the right public dataclass types
- scoring support floor is respected (`min_substrates >= 2`)

Motif scoring sequence-context validation is part of the kinase scoring-stage
scientific contract:

- default expected motif window length is `15` residues (centre index `7`)
- default scoring expects centred windows and does not midpoint-crop longer input
  sequences
- longer centred sequences are accepted only in explicit centred-sequence
  extraction mode
- centre residue is validated against site identity when available
- supported residue characters are the 20 canonical amino acids
- query/target windows that are missing, short, non-centred,
  site-residue-mismatched, non-phospho-centre (`S/T/Y`), or unsupported are
  excluded from motif scoring
- reference/library windows are validated with the same rules before motif
  profile construction; invalid reference windows are excluded from motif model
  construction
- explicit motif-sequence libraries accept both:
  - bare sequence entries, and
  - structured entries with metadata equivalent to
    `reference_id`, `site_id`, `kinase`, and `sequence`
- bare explicit entries can validate motif-window quality, alphabet support, and
  phospho-compatible centre residue (`S/T/Y`), but cannot always assert intended
  phosphosite identity
- structured explicit entries are preferred because supplied `site_id` can be
  validated and centre-residue matching can be checked deterministically
- site-residue mismatch checks require site-like metadata (for example `site_id`
  or a reference identifier that encodes the site)
- query-side sequence validation diagnostics are exposed on
  `result.scoring_result.motif_sequence_validation`
- reference-library validation diagnostics are exposed on
  `result.scoring_result.motif_library_validation`

`KinaseScoringConfig`, `KinasePredictionConfig`, and `KinaseActivityConfig`
validate local policy/range rules at object construction.

### Signalome workflow

`SignalomeWorkflowRequest` validation checks:

- `kinase_result` is `KinaseWorkflowResult`
- `config` is `SignalomeConfig`
- upstream matrices are usable for signalome execution
- `kinase_result.dataset.site_metadata.protein_id` exists and is non-empty for all interpreted sites

`SignalomeConfig` validates local policy/range rules at object construction.
This includes:

- `cluster_tree_backend` policy values
- `candidate_scoring_backend` policy values
- `max_exact_cluster_tree_sites >= 1`
- `max_full_correlation_sites >= 1`

Signalome execution enforces runtime scale guards at the expensive call sites:

- exact cluster-tree construction is guarded by `max_exact_cluster_tree_sites`
- full candidate-correlation scoring is guarded by `max_full_correlation_sites`

If either guard is exceeded, execution fails with `SignalomeScaleError`.

Signalome is intentionally strict about protein identity. A site ID such as
`TSC2;S939;` is not a substitute for `protein_id`.

## Boundary errors during workflow execution

Some failures happen after request validation, during workflow interpretation or
execution. These often raise `WorkflowBoundaryError` with:

- a seam name
- counts or other details
- a `next_action` hint

Typical examples are overlap failures, low-support failures, or signalome
network/module preconditions.

## Quick fix table

| Problem | Common fix |
| --- | --- |
| input format rejected | pass a DataFrame or supported file path |
| dataset organism missing for `AUTO` | set `organism=Organism.RAT` for bundled first runs |
| bundled human/mouse preset fails | use an explicit `ReferenceBundle` |
| signalome fails on `protein_id` | add a non-empty `protein_id` column |
| signalome scale guard fails with `SignalomeScaleError` | reduce interpreted sites, use `candidate_scoring_backend="sampled"` only for candidate module-count scoring cost (it does not bypass exact cluster-tree construction), and/or deliberately raise `max_exact_cluster_tree_sites` / `max_full_correlation_sites` |
| rows dropped in site-matrix building | review sequence support and preprocessing policy |

## Validation ownership summary

| Invariant | Owner |
| --- | --- |
| builder input source checks | `DatasetBuildRequestValidator` |
| preprocessing config policy | `DatasetPreprocessingConfigValidator` |
| analysis-ready dataset structure/content | `AnalysisReadyDatasetValidator` |
| intensity-scale-state coherence | `IntensityScaleStateValidator` |
| reference compatibility | `ReferenceCompatibilityValidator` |
| reference bundle structure/content | `ReferenceBundleValidator` |
| kinase workflow request/config validity | `KinaseWorkflowValidator` |
| signalome workflow request/config validity | `SignalomeWorkflowValidator` |
| runtime seam diagnostics | workflow interpreters/executors via `WorkflowBoundaryError` |

## Where next

- [Troubleshooting](getting-started/troubleshooting-first-run.md)
- [API Guide](api.md)
- [CLI Guide](cli.md)
