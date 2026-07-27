# Large Differential Limma Trend Fixture Provenance

Generated with R version 4.5.2 (2025-10-31 ucrt)
limma version: 3.66.0
Seed: 20260724
Generation timestamp (UTC): 2026-07-24T00:00:00Z
Command: `Rscript scripts/active/generate_large_differential_limma_trend_fixture.R --outdir tests/fixtures/rewrite_parity/differential_limma_trend_large --seed 20260724 --timestamp 2026-07-24T00:00:00Z --n_features 1600`
Source policy: deterministic synthetic fixture generated locally without network access; limma is the external authority for exported parity quantities.
Classification: external parity for limma result columns; simulation diagnostics are fixture sanity metadata.
Design: ~0 + condition with groups A/B and unbalanced 5/7 replicates
Contrast: B_vs_A = B - A
Rows: 1600 phosphosites/features; columns: 12 samples
Mean-variance trend: expected residual variance increases smoothly with mean intensity plus deterministic sinusoidal structure.
Shifted features: 152 total; positive=100; negative=52
Generator SHA-256: 17eb99447b2b7cc33e8b76df3e2481031f795c96a7628ae68af824310e35a0e2

Output files are listed with SHA-256 digests in `MANIFEST.json`.
