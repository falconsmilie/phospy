# Dataset Builders

PhosPy's dataset-builder entrypoint is `AnalysisReadyDatasetBuilder`.
Detailed API usage lives in
[Dataset Build Workflow](dataset-build-workflow.md).

Identity boundary summary:

- Analysis-ready row identity is `site_key`.
- `display_id` (for example `GENE;SITE;`) is a human-readable label and may
  repeat after `site_key` becomes row identity.
- `site_key` must be unique.
- `AnalysisReadyPhosphoDataset.phospho.index` is `site_key`.
- `AnalysisReadyPhosphoDataset.site_metadata.index` is `site_key`.
- `AnalysisReadyPhosphoDataset.site_metadata["display_id"]` is required.
- `AnalysisReadyPhosphoDataset.site_metadata["site_key"]` must exactly match
  `site_metadata.index`.
- Direct analysis-ready construction requires auditable protein context metadata:
  `organism`, `protein_namespace`, `protein_identifier`, `gene_symbol`, `site`,
  and `site_sequence`.
- Direct analysis-ready construction requires encoded `site_key` indexes and
  does not silently fall back to display-site identity.
- Builder input may accept legacy display-indexed shape only when enough
  protein context exists to derive `site_key`.
- Workflows operate on `site_key`; site-level outputs that materialize row
  identity include both `site_key` and `display_id`.
- Kinase references may remain display-ID keyed only at the reference boundary;
  the kinase workflow projects them through an explicit `display_id` ->
  `site_key` mapping before scoring.
- Duplicate rows that resolve to the same `site_key` are a scientific ambiguity
  and fail by default during site-matrix preprocessing.
- Non-error duplicate-site policies (`max_mean_signal`, `first`,
  `aggregate_mean`, `aggregate_median`) are deliberate scientific choices.
- When a non-error duplicate-site policy is used, inspect
  `dataset.preprocessing_report.duplicate_site_resolution` and
  `metadata_conflicts`.
- Duplicate `display_id` values remain valid when the corresponding `site_key`
  values differ.
- Duplicate rows are not automatically renamed.
- Peptide-evidence modelling scope is unchanged in this identity migration.

See
[ADR-0024: Protein-Scoped Phosphosite Row Identity](../adr/adr_0024_protein_scoped_phosphosite_row_identity.md)
for the formal architecture decision.
