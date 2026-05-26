# Dataset Builders

PhosPy's dataset-builder entrypoint is `AnalysisReadyDatasetBuilder`.
Detailed API usage lives in
[Dataset Build Workflow](dataset-build-workflow.md).

Identity boundary summary:

- PhosPy currently supports one analysis-ready row per normalised phosphosite
  display identifier such as `GENE;SITE;`.
- Duplicate display-site rows are rejected during dataset construction.
- Duplicate rows are not automatically aggregated, collapsed, or renamed.
- Protein/isoform/source fields can be stored in metadata, but they do not
  currently define row identity.
- Resolve duplicate display-site rows upstream before constructing
  `AnalysisReadyPhosphoDataset`.
