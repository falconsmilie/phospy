suppressPackageStartupMessages(library(limma))

samples <- c("A_1", "A_2", "B_1", "B_2")
group <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
design <- model.matrix(~0 + group)
colnames(design) <- c("A", "B")
rownames(design) <- samples

mat <- matrix(
  c(
    0.10, 0.20, 1.10, 1.20,  # SITE_01 ordinary positive effect
    0.20, 0.10, 0.90, 1.00,  # SITE_02 ordinary positive effect
    0.30, 0.20, 0.35, 0.15,  # SITE_03 near-null effect
   -0.10, 0.00, 1.20, -0.80, # SITE_04 unequal variance
    2.00, 2.00, 2.00, 2.00,  # SITE_05 zero variance
    1.10, 1.00, 0.10, 0.20,  # SITE_06 negative effect for B_vs_A
    0.00, 0.05, 0.10, 0.00,  # SITE_07 tiny effect
   -0.40, -0.30, 0.40, 0.60  # SITE_08 larger positive effect
  ),
  byrow = TRUE,
  nrow = 8,
  ncol = 4
)
rownames(mat) <- sprintf("SITE_%02d", seq_len(nrow(mat)))
colnames(mat) <- samples

contrast_mat <- makeContrasts(B_vs_A = B - A, A_vs_B = A - B, levels = design)

fit <- lmFit(mat, design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

tt_ba <- topTable(
  fit2,
  coef = "B_vs_A",
  number = Inf,
  sort.by = "none",
  adjust.method = "BH"
)
tt_ab <- topTable(
  fit2,
  coef = "A_vs_B",
  number = Inf,
  sort.by = "none",
  adjust.method = "BH"
)

write.csv(
  data.frame(site_id = rownames(mat), mat, check.names = FALSE),
  "tests/fixtures/rewrite_parity/differential_limma_envelope/matrix.csv",
  row.names = FALSE
)
write.csv(
  data.frame(sample = rownames(design), design, check.names = FALSE),
  "tests/fixtures/rewrite_parity/differential_limma_envelope/design.csv",
  row.names = FALSE
)
write.csv(
  data.frame(
    coefficient = rownames(contrast_mat),
    contrast_mat,
    check.names = FALSE
  ),
  "tests/fixtures/rewrite_parity/differential_limma_envelope/contrasts.csv",
  row.names = FALSE
)
write.csv(
  data.frame(
    site_id = rownames(tt_ba),
    tt_ba[, c("logFC", "t", "P.Value", "adj.P.Val"), drop = FALSE],
    check.names = FALSE
  ),
  "tests/fixtures/rewrite_parity/differential_limma_envelope/limma_B_vs_A.csv",
  row.names = FALSE
)
write.csv(
  data.frame(
    site_id = rownames(tt_ab),
    tt_ab[, c("logFC", "t", "P.Value", "adj.P.Val"), drop = FALSE],
    check.names = FALSE
  ),
  "tests/fixtures/rewrite_parity/differential_limma_envelope/limma_A_vs_B.csv",
  row.names = FALSE
)

provenance <- c(
  "# Differential Limma Envelope Fixture Provenance",
  "",
  paste0("Generated with ", R.version$version.string),
  paste0("limma version: ", as.character(packageVersion("limma"))),
  "Seed: none (deterministic matrix literals)",
  "Design: ~0 + group with groups A/B and 2 replicates each",
  "Contrasts (column order): B_vs_A then A_vs_B",
  "Rows: 8 sites; columns: 4 samples",
  "Contains edge rows: zero variance (SITE_05), unequal variance (SITE_04)"
)
writeLines(
  provenance,
  con = "tests/fixtures/rewrite_parity/differential_limma_envelope/PROVENANCE.md"
)
