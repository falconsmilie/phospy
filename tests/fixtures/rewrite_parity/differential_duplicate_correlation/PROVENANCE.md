# Differential Duplicate-Correlation Limma Fixture Provenance

Generated with R version 4.5.2 (2025-10-31 ucrt)
Bioconductor version: 3.22
limma version: 3.66.0
Seed: 20260818
Generation timestamp (UTC): 2026-08-18T00:00:00Z
Command: `Rscript scripts/active/generate_differential_duplicate_correlation_limma_fixtures.R --outdir tests/fixtures/rewrite_parity/differential_duplicate_correlation --seed 20260818 --timestamp 2026-08-18T00:00:00Z --allow-unpinned-environment false`
Generator SHA-256: d12f2c97ba889530665588118b0229a0502c1266437f4bc361a0a217a8b47d3e
Byte policy: utf-8 LF with final newline
Serialization policy: CSV uses comma separators, a header row, row.names=FALSE, UTF-8, LF line endings, a final newline, options(digits=17, scipen=999), and the literal NA token for missing numeric values. R NaN and +/-Inf numeric outputs are serialized as NA; feature-correlation failure details are retained in companion status and missing-kind columns. JSON manifests use stable key ordering and UTF-8 LF bytes.
Source policy: deterministic synthetic fixture generated locally without network access; limma duplicateCorrelation, lmFit, contrasts.fit, and eBayes outputs are the external scientific authority for expected numerical columns.
Classification: external parity for R/limma duplicate-correlation intermediate and final outputs.
Expected outputs come only from the pinned R/limma run. PhosPy is not imported or executed by this generator.
Redistribution metadata: synthetic deterministic inputs and black-box limma numeric outputs are repository test fixtures; limma source code is not redistributed.

Scientific citations:
- Smyth GK (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology 3(1), Article 3.
- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.

Fixtures:
- Fixture A - complete two-condition pairs (`fixture_a_complete_pairs`)
- Fixture B - more than two observations per block (`fixture_b_three_observation_blocks`)
- Fixture C - incomplete and unequal blocks (`fixture_c_incomplete_unequal_blocks`)
- Fixture D - feature-level estimator failures (`fixture_d_feature_level_failures`)

Each fixture uses a fixed-effects design matrix without block dummy variables and supplies block IDs only to limma duplicateCorrelation/lmFit.
Output files and SHA-256 digests are listed in `MANIFEST.json`.
