# Differential DEqMS Depth Fixture Provenance

Generated with R version 4.5.2 (2025-10-31 ucrt)
Bioconductor version: 3.22
limma version: 3.66.0
DEqMS version: 1.28.0
Seed: 20260909
Generation timestamp (UTC): 2026-09-09T00:00:00Z
Command: `Rscript tests/fixtures/rewrite_parity/differential_deqms_depth/generate_fixture.R --outdir tests/fixtures/rewrite_parity/differential_deqms_depth --seed 20260909 --timestamp 2026-09-09T00:00:00Z --allow-unpinned-environment false`
Generator SHA-256: 691393efe7854c320c76893775392eab3bbd921c5cc417371bac92494eba279b
Byte policy: utf-8 LF with final newline
Serialization policy: CSV uses comma separators, a header row, row.names=FALSE, UTF-8, LF line endings, a final newline, options(digits=17, scipen=999), and the literal NA token for missing numeric values. R NaN and +/-Inf numeric outputs are serialized as NA. JSON manifests use stable key ordering and UTF-8 LF bytes.
Source policy: deterministic synthetic fixture generated locally without network access; DEqMS::spectraCounteBayes and limma outputs are the external scientific authority for expected numerical columns.
Classification: external parity for DEqMS count-aware empirical-Bayes outputs; feature metadata are fixture sanity metadata.
Expected outputs come only from the pinned R/DEqMS and R/limma run. PhosPy is not imported or executed by this generator.
Redistribution metadata: synthetic deterministic inputs and black-box DEqMS/limma numeric outputs are repository test fixtures; DEqMS and limma source code is not redistributed.
Design: ~0 + condition with groups A/B and unbalanced 5/7 replicates
Contrast: B_vs_A = B - A
Rows: 144 phosphosites/features; columns: 12 samples
Quantification depth: feature-level synthetic PSM count supplied as fit$count before DEqMS::spectraCounteBayes.
DEqMS model: spectraCounteBayes(fit.method='loess') uses log2(count) and loess span 0.75 to estimate count-dependent prior variance.
Depth groups: lower_depth=49; middle_depth=48; higher_depth=47
Comparable depth pair: QD021;S21; has depth 2; QD117;S117; has depth 32; their intensity rows are byte-identical before model fitting.
Terminology conclusion: PhosPy's validated feature should be described as quantification-depth-aware empirical Bayes moderation inspired by DEqMS, not exact DEqMS-compatible numerical equivalence, because PhosPy preserves its existing deterministic trend smoother while DEqMS::spectraCounteBayes uses R loess span 0.75.

Scientific citations:
- Zhu Y, Orre LM, Zhou Tran Y, Mermelekas G, Johansson HJ, Malyutina A, Anders S, Lehtio J (2020). DEqMS: A Method for Accurate Variance Estimation in Differential Protein Expression Analysis. Molecular & Cellular Proteomics 19(6), 1047-1057.
- Smyth GK (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology 3(1), Article 3.
- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.

Output files and SHA-256 digests are listed in `MANIFEST.json`.
