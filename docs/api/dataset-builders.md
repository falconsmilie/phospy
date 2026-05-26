# Dataset Builders

PhosPy's dataset-builder entrypoint is `AnalysisReadyDatasetBuilder`.
Detailed API usage lives in
[Dataset Build Workflow](dataset-build-workflow.md).

Identity boundary summary:

- PhosPy currently supports one analysis-ready row per normalised phosphosite
  display identifier such as `GENE;SITE;`.
- Duplicate display-site rows are rejected during dataset construction.
- Duplicate display-site rows are rejected even when duplicate rows carry
  identical protein metadata.
- Duplicate rows are not automatically aggregated.
- Duplicate rows are not automatically renamed.
- Protein and isoform context can be retained as metadata, but protein and
  isoform context does not currently define row identity.
- Source and mapping context can be retained as metadata, but source context
  does not currently define row identity.
- Peptide-evidence-scoped row identity is not currently supported.
- Resolve duplicated display-site rows upstream before constructing
  `AnalysisReadyPhosphoDataset`.

See
[ADR-0023: Supported Phosphosite Display-Site Identity Scope](../adr/adr_0023_supported_phosphosite_display_site_identity_scope.md)
for the formal architecture decision.
