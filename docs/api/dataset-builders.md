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
- Protein context is required to construct `site_key`.
- Direct analysis-ready datasets must not silently fall back to display-site
  identity.
- Builder input may accept legacy display-indexed shape only when enough
  protein context exists to derive `site_key`.
- Duplicate rows are not automatically aggregated beyond explicit duplicate
  policies.
- Duplicate rows are not automatically renamed.
- Peptide-evidence modelling scope is unchanged in this identity migration.

See
[ADR-0024: Protein-Scoped Phosphosite Row Identity](../adr/adr_0024_protein_scoped_phosphosite_row_identity.md)
for the formal architecture decision.
