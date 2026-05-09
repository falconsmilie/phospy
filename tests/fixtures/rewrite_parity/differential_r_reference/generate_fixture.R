suppressPackageStartupMessages(library(limma))

set.seed(42)

samples <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2")
group <- factor(c("A", "A", "B", "B", "C", "C"), levels = c("A", "B", "C"))
design <- model.matrix(~0 + group)
colnames(design) <- c("A", "B", "C")
rownames(design) <- samples

mat <- matrix(
  rnorm(12 * length(samples), mean = 0, sd = 0.2),
  nrow = 12,
  ncol = length(samples)
)
mat[1:4, group == "B"] <- mat[1:4, group == "B"] + 1.0
mat[5:8, group == "C"] <- mat[5:8, group == "C"] + 1.5
rownames(mat) <- sprintf("SITE_%02d", seq_len(nrow(mat)))
colnames(mat) <- samples

contrast_mat <- makeContrasts(B_vs_A = B - A, C_vs_A = C - A, levels = design)

fit <- lmFit(mat, design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

tt_b <- topTable(
  fit2,
  coef = "B_vs_A",
  number = Inf,
  sort.by = "none",
  adjust.method = "BH"
)
tt_c <- topTable(
  fit2,
  coef = "C_vs_A",
  number = Inf,
  sort.by = "none",
  adjust.method = "BH"
)

write.csv(
  data.frame(site_id = rownames(mat), mat, check.names = FALSE),
  "tests/fixtures/rewrite_parity/differential_r_reference/matrix.csv",
  row.names = FALSE
)
write.csv(
  data.frame(sample = rownames(design), design, check.names = FALSE),
  "tests/fixtures/rewrite_parity/differential_r_reference/design.csv",
  row.names = FALSE
)
write.csv(
  data.frame(
    coefficient = rownames(contrast_mat),
    contrast_mat,
    check.names = FALSE
  ),
  "tests/fixtures/rewrite_parity/differential_r_reference/contrasts.csv",
  row.names = FALSE
)
write.csv(
  data.frame(
    site_id = rownames(tt_b),
    tt_b[, c("logFC", "t", "P.Value", "adj.P.Val"), drop = FALSE],
    check.names = FALSE
  ),
  "tests/fixtures/rewrite_parity/differential_r_reference/limma_B_vs_A.csv",
  row.names = FALSE
)
write.csv(
  data.frame(
    site_id = rownames(tt_c),
    tt_c[, c("logFC", "t", "P.Value", "adj.P.Val"), drop = FALSE],
    check.names = FALSE
  ),
  "tests/fixtures/rewrite_parity/differential_r_reference/limma_C_vs_A.csv",
  row.names = FALSE
)

provenance <- c(
  "# Differential Parity Fixture Provenance",
  "",
  paste0("Generated with ", R.version$version.string),
  paste0("limma version: ", as.character(packageVersion("limma"))),
  "Seed: 42",
  "Design: ~0 + group with groups A/B/C and 2 replicates each",
  "Contrasts: B_vs_A and C_vs_A",
  "Rows: 12 sites; columns: 6 samples"
)
writeLines(
  provenance,
  con = "tests/fixtures/rewrite_parity/differential_r_reference/PROVENANCE.md"
)
