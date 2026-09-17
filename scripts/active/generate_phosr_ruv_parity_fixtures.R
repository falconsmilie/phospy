#!/usr/bin/env Rscript

# Generate synthetic external-reference evidence for SPS discovery and RUV-III.
#
# This script is a maintainer-only regeneration tool.  Normal Python tests read
# the committed outputs and neither start R nor access the network.

options(stringsAsFactors = FALSE, digits = 17, scipen = 999)

EXPECTED <- list(
  r = "R version 4.5.2 (2025-10-31 ucrt)",
  bioconductor = "3.22",
  BiocManager = "1.30.27",
  PhosR = "1.13.1",
  PhosR_remote_sha = "1be74902b775833c64f5833e70538eaf843cf6a5",
  ruv = "0.9.7.1",
  SummarizedExperiment = "1.40.0",
  S4Vectors = "0.48.0",
  BiocGenerics = "0.56.0",
  Matrix = "1.7.4",
  jsonlite = "2.0.0",
  digest = "0.6.39"
)

parse_args <- function(args) {
  result <- list(
    outdir = "tests/fixtures/rewrite_parity/phosr_ruv",
    timestamp = "2026-09-17T00:00:00Z",
    allow_unpinned_environment = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (i == length(args)) stop("missing value for argument: ", key)
    value <- args[[i + 1L]]
    if (key == "--outdir") result$outdir <- value
    else if (key == "--timestamp") result$timestamp <- value
    else if (key == "--allow-unpinned-environment") {
      result$allow_unpinned_environment <- identical(tolower(value), "true")
    } else stop("unknown argument: ", key)
    i <- i + 2L
  }
  result
}

required_pkgs <- c(
  "BiocManager", "PhosR", "ruv", "SummarizedExperiment", "S4Vectors",
  "BiocGenerics", "Matrix", "jsonlite", "digest"
)
for (package in required_pkgs) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("required pinned reference package is unavailable: ", package)
  }
}

observed <- list(
  r = R.version.string,
  bioconductor = as.character(BiocManager::version()),
  BiocManager = as.character(packageVersion("BiocManager")),
  PhosR = as.character(packageVersion("PhosR")),
  PhosR_remote_sha = packageDescription("PhosR")$RemoteSha,
  ruv = as.character(packageVersion("ruv")),
  SummarizedExperiment = as.character(packageVersion("SummarizedExperiment")),
  S4Vectors = as.character(packageVersion("S4Vectors")),
  BiocGenerics = as.character(packageVersion("BiocGenerics")),
  Matrix = as.character(packageVersion("Matrix")),
  jsonlite = as.character(packageVersion("jsonlite")),
  digest = as.character(packageVersion("digest"))
)
if (!isTRUE(all.equal(EXPECTED, observed)) &&
    !parse_args(commandArgs(trailingOnly = TRUE))$allow_unpinned_environment) {
  stop(
    "reference environment does not match the pinned specification\nexpected: ",
    paste(names(EXPECTED), unlist(EXPECTED), sep = "=", collapse = ", "),
    "\nobserved: ",
    paste(names(observed), unlist(observed), sep = "=", collapse = ", ")
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
outdir <- normalizePath(args$outdir, winslash = "/", mustWork = FALSE)
if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "sps"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "ruv_iii"), recursive = TRUE, showWarnings = FALSE)

write_lf <- function(path, text) {
  payload <- paste0(paste(text, collapse = "\n"), "\n")
  connection <- file(path, open = "wb")
  on.exit(close(connection))
  writeBin(charToRaw(enc2utf8(payload)), connection)
}

write_csv_lf <- function(frame, path) {
  connection <- textConnection("lines", "w", local = TRUE)
  write.table(
    frame,
    connection,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE,
    na = "NA",
    qmethod = "double"
  )
  close(connection)
  write_lf(path, lines)
}

write_json_lf <- function(value, path) {
  text <- jsonlite::toJSON(
    value,
    auto_unbox = TRUE,
    pretty = TRUE,
    digits = 17,
    null = "null",
    na = "string"
  )
  write_lf(path, text)
}

site_key <- function(protein, residue, position) {
  paste0(
    "phospy:v1|organism=rat|protein_namespace=protein_id|protein_identifier=",
    protein,
    "|residue=", residue,
    "|position=", position
  )
}

matrix_frame <- function(matrix, keys) {
  data.frame(site_key = keys, matrix, check.names = FALSE)
}

# SPS: source absolute log2 abundance is explicitly centred by subtracting the
# per-site mean of the two baseline samples.  The resulting condition-relative
# log2 matrix is supplied unchanged to both getSPS and PhosPy.
sps_count <- 210L
sps_index <- seq_len(sps_count)
sps_proteins <- sprintf("SPS%04d", sps_index)
sps_genes <- sps_proteins
sps_residue <- rep(c("S", "T", "Y"), length.out = sps_count)
sps_position <- 100L + sps_index
sps_keys <- mapply(site_key, sps_proteins, sps_residue, sps_position)
sps_ids <- paste0(sps_genes, ";", sps_residue, sps_position)
sps_reference_ids <- c("reference_a", "reference_b", "reference_c")
sps_samples <- c("baseline_1", "baseline_2", "stimulated_1", "stimulated_2")
sps_conditions <- c("baseline", "baseline", "stimulated", "stimulated")
sps_relative <- list()
sps_absolute <- list()
sps_phospho <- list()
sps_conds <- list()

for (dataset_position in seq_along(sps_reference_ids)) {
  dataset_id <- sps_reference_ids[[dataset_position]]
  base <- 18 + dataset_position * 0.4 + sps_index / 500
  change <- 0.025 + (sps_index - 1L) * 0.004
  change[sps_index > 180L] <- 2 + (sps_index[sps_index > 180L] - 180L) * 0.08
  change <- change + (dataset_position - 1L) * ((sps_index %% 7L) - 3L) * 0.0007
  # A deliberate tie group exercises reference midranks.  The external final
  # ordering within that group is treated as unspecified by the Python tests.
  change[c(20L, 21L)] <- 0.125 + dataset_position * 0.0003
  baseline_noise <- 0.006 + (sps_index %% 5L) * 0.0002
  treatment_noise <- 0.008 + (sps_index %% 3L) * 0.0003
  # Make the deliberate tie byte-identical after centring and serialization;
  # every other site's score is separated from its neighbours.
  base[21L] <- base[20L]
  baseline_noise[21L] <- baseline_noise[20L]
  treatment_noise[21L] <- treatment_noise[20L]
  absolute <- cbind(
    baseline_1 = base - baseline_noise,
    baseline_2 = base + baseline_noise,
    stimulated_1 = base + change - treatment_noise,
    stimulated_2 = base + change + treatment_noise
  )
  baseline_mean <- rowMeans(absolute[, sps_conditions == "baseline", drop = FALSE])
  relative <- sweep(absolute, 1L, baseline_mean, FUN = "-")
  row_order <- switch(
    dataset_position,
    sps_index,
    rev(sps_index),
    c(sps_index[71:210], sps_index[1:70])
  )
  absolute <- absolute[row_order, , drop = FALSE]
  relative <- relative[row_order, , drop = FALSE]
  genes <- sps_genes[row_order]
  residues <- sps_residue[row_order]
  positions <- sps_position[row_order]
  keys <- sps_keys[row_order]
  sps_absolute[[dataset_id]] <- absolute
  sps_relative[[dataset_id]] <- relative
  sps_phospho[[dataset_id]] <- PhosR::PhosphoExperiment(
    assays = list(Quantification = relative),
    UniprotID = sps_proteins[row_order],
    GeneSymbol = genes,
    Site = positions,
    Residue = residues,
    Sequence = rep("AAAAAAASAAAAAAA", sps_count),
    Localisation = rep(1, sps_count)
  )
  sps_conds[[dataset_id]] <- sps_conditions
  write_csv_lf(
    matrix_frame(absolute, keys),
    file.path(outdir, "sps", paste0(dataset_id, "_source_absolute_log2.csv"))
  )
  write_csv_lf(
    matrix_frame(relative, keys),
    file.path(outdir, "sps", paste0(dataset_id, "_condition_relative_log2.csv"))
  )
}

sps_top_n <- 25L
sps_selected_raw <- PhosR::getSPS(
  phosData = unname(sps_phospho),
  assays = "Quantification",
  conds = unname(sps_conds),
  num = sps_top_n
)
sps_selected_ids <- sub(";$", "", sps_selected_raw)
sps_key_by_id <- setNames(sps_keys, sps_ids)
sps_selected_keys <- unname(sps_key_by_id[sps_selected_ids])

# Capture the pinned getSPS intermediates explicitly and assert that their final
# list agrees with the package call.  This mirrors getSPS 1.13.1 solely for
# diagnostic evidence; getSPS itself remains the external selection authority.
sps_scores <- list()
for (dataset_id in sps_reference_ids) {
  matrix <- sps_relative[[dataset_id]]
  ids_in_order <- paste0(
    PhosR::GeneSymbol(sps_phospho[[dataset_id]]), ";",
    PhosR::Residue(sps_phospho[[dataset_id]]),
    PhosR::Site(sps_phospho[[dataset_id]])
  )
  condition_means <- sapply(
    unique(sps_conditions),
    function(condition) rowMeans(matrix[, sps_conditions == condition, drop = FALSE])
  )
  score <- apply(condition_means, 1L, function(values) max(abs(values)))
  names(score) <- ids_in_order
  sps_scores[[dataset_id]] <- score
}
common_ids <- sort(Reduce(intersect, lapply(sps_scores, names)))
score_matrix <- sapply(sps_reference_ids, function(id) sps_scores[[id]][common_ids])
colnames(score_matrix) <- sps_reference_ids
rank_quantiles <- apply(-abs(score_matrix), 2L, rank, ties.method = "average")
rank_quantiles <- (rank_quantiles - 0.5) / nrow(score_matrix)
fisher_statistic <- -2 * rowSums(log(rank_quantiles))
consensus_score <- pchisq(
  fisher_statistic,
  df = (length(sps_reference_ids) - 1L) * 2L,
  lower.tail = FALSE
)
package_order <- names(sort(consensus_score, decreasing = TRUE))
stopifnot(identical(package_order[seq_len(sps_top_n)], sps_selected_ids))
sps_trace <- data.frame(
  phosr_site_id = common_ids,
  site_key = unname(sps_key_by_id[common_ids]),
  reference_a_stability = score_matrix[, "reference_a"],
  reference_a_rank_quantile = rank_quantiles[, "reference_a"],
  reference_b_stability = score_matrix[, "reference_b"],
  reference_b_rank_quantile = rank_quantiles[, "reference_b"],
  reference_c_stability = score_matrix[, "reference_c"],
  reference_c_rank_quantile = rank_quantiles[, "reference_c"],
  fisher_statistic = fisher_statistic,
  consensus_score = consensus_score,
  package_rank = match(common_ids, package_order),
  selected_top_n = common_ids %in% sps_selected_ids,
  check.names = FALSE
)
write_csv_lf(sps_trace, file.path(outdir, "sps", "phosr_getsps_trace.csv"))
write_csv_lf(
  data.frame(
    rank = seq_len(sps_top_n),
    phosr_site_id = sps_selected_ids,
    site_key = sps_selected_keys
  ),
  file.path(outdir, "sps", "phosr_getsps_selected.csv")
)
write_csv_lf(
  data.frame(
    dataset_id = rep(sps_reference_ids, each = length(sps_samples)),
    sample_id = rep(sps_samples, times = length(sps_reference_ids)),
    condition = rep(sps_conditions, times = length(sps_reference_ids)),
    is_baseline = rep(sps_conditions == "baseline", times = length(sps_reference_ids))
  ),
  file.path(outdir, "sps", "sample_metadata.csv")
)

# RUV-III: controls are already selected.  Replicate-set means carry biological
# variation; two independent within-set factors carry unwanted variation.
ruv_samples <- paste0("sample_", seq_len(8L))
ruv_sets <- rep(paste0("replicate_", letters[1:4]), each = 2L)
ruv_w1 <- c(-2, 2, -1, 1, -3, 3, -1.5, 1.5)
ruv_w2 <- c(-1, 1, 2, -2, -0.5, 0.5, 3, -3)
ruv_biology <- rep(c(-4, -1, 2, 5), each = 2L)
ruv_values <- rbind(
  control_1 = 40 + 1.0 * ruv_w1 + 0.2 * ruv_w2,
  control_2 = 20 - 1.7 * ruv_w1 + 0.4 * ruv_w2,
  control_3 = 30 + 0.3 * ruv_w1 - 1.4 * ruv_w2,
  control_4 = 12 + 1.2 * ruv_w1 + 1.1 * ruv_w2,
  biological_1 = 50 + ruv_biology + 0.8 * ruv_w1 - 0.2 * ruv_w2,
  biological_2 = 10 - 1.5 * ruv_biology + 0.1 * ruv_w1 + 0.7 * ruv_w2,
  mixed = 25 + 0.5 * ruv_biology - 0.6 * ruv_w1 + 0.9 * ruv_w2,
  invariant = rep(7, length(ruv_samples))
)
colnames(ruv_values) <- ruv_samples
ruv_keys <- mapply(
  site_key,
  sprintf("RUV%04d", seq_len(nrow(ruv_values))),
  rep(c("S", "T", "Y"), length.out = nrow(ruv_values)),
  200L + seq_len(nrow(ruv_values))
)
ruv_controls <- seq_len(4L)
ruv_mapping <- model.matrix(~ factor(ruv_sets) - 1)
colnames(ruv_mapping) <- unique(ruv_sets)
rownames(ruv_mapping) <- ruv_samples
ruv_y <- t(ruv_values)

write_csv_lf(
  matrix_frame(ruv_values, ruv_keys),
  file.path(outdir, "ruv_iii", "input_matrix.csv")
)
write_csv_lf(
  data.frame(
    sample_id = ruv_samples,
    replicate_set = ruv_sets,
    planted_biological_effect = ruv_biology,
    planted_unwanted_factor_1 = ruv_w1,
    planted_unwanted_factor_2 = ruv_w2
  ),
  file.path(outdir, "ruv_iii", "sample_metadata.csv")
)
write_csv_lf(
  data.frame(site_key = ruv_keys[ruv_controls], control_status = "negative_control"),
  file.path(outdir, "ruv_iii", "negative_controls.csv")
)

ruv_diagnostics <- list()
for (k in c(0L, 1L, 2L)) {
  result <- ruv::RUVIII(
    Y = ruv_y,
    M = ruv_mapping,
    ctl = seq_len(ncol(ruv_y)) %in% ruv_controls,
    k = k,
    return.info = TRUE
  )
  corrected <- t(result$newY)
  write_csv_lf(
    matrix_frame(corrected, ruv_keys),
    file.path(outdir, "ruv_iii", paste0("ruv_expected_k", k, ".csv"))
  )
  residual <- (diag(nrow(ruv_y)) - ruv_mapping %*%
    solve(t(ruv_mapping) %*% ruv_mapping) %*% t(ruv_mapping)) %*% ruv_y
  singular_values <- svd(residual, nu = 0L, nv = 0L)$d
  ruv_diagnostics[[paste0("k", k)]] <- list(
    k = k,
    input_matrix_rank = qr(ruv_y)$rank,
    replicate_mapping_rank = qr(ruv_mapping)$rank,
    replicate_residual_rank = qr(residual)$rank,
    residual_singular_values = unname(singular_values[seq_len(k)]),
    corrected_matrix_file = paste0("ruv_iii/ruv_expected_k", k, ".csv")
  )
}

ruv_permutation <- c(7L, 2L, 5L, 4L, 1L, 8L, 3L, 6L)
permuted_result <- ruv::RUVIII(
  Y = ruv_y[ruv_permutation, , drop = FALSE],
  M = ruv_mapping[ruv_permutation, , drop = FALSE],
  ctl = seq_len(ncol(ruv_y)) %in% ruv_controls,
  k = 2L,
  return.info = TRUE
)
permuted_corrected <- t(permuted_result$newY)
colnames(permuted_corrected) <- ruv_samples[ruv_permutation]
write_csv_lf(
  matrix_frame(permuted_corrected, ruv_keys),
  file.path(outdir, "ruv_iii", "ruv_expected_k2_permuted.csv")
)
write_json_lf(
  list(
    cases = ruv_diagnostics,
    sample_permutation = ruv_samples[ruv_permutation],
    factor_comparison_policy = paste(
      "estimated factor coordinates are not compared directly because SVD factors",
      "are invariant to sign and, for tied subspaces, rotation"
    )
  ),
  file.path(outdir, "ruv_iii", "reference_diagnostics.json")
)

environment_payload <- c(
  EXPECTED,
  list(
    platform = R.version$platform,
    blas = if (nzchar(unname(extSoftVersion()[["BLAS"]]))) {
      unname(extSoftVersion()[["BLAS"]])
    } else {
      "R default BLAS"
    },
    lapack = La_version(),
    phosr_source = paste0(
      "https://github.com/SlavovLab/PhosR/commit/", EXPECTED$PhosR_remote_sha
    ),
    ruv_source = paste0(
      "https://cran.r-project.org/src/contrib/Archive/ruv/ruv_", EXPECTED$ruv,
      ".tar.gz"
    ),
    bioconductor_release = EXPECTED$bioconductor
  )
)
write_json_lf(environment_payload, file.path(outdir, "REFERENCE_ENVIRONMENT.json"))

metadata <- list(
  fixture_schema_version = "phosr-ruv-parity-v1",
  classification = "external_parity",
  generation_timestamp_utc = args$timestamp,
  scientific_citations = c(
    paste(
      "Kim HJ et al. (2021). PhosR enables processing and functional analysis",
      "of phosphoproteomic data. Cell Reports 34(8):108771.",
      "doi:10.1016/j.celrep.2021.108771"
    ),
    paste(
      "Gagnon-Bartsch J (2019). ruv: Detect and Remove Unwanted Variation",
      "using Negative Controls. R package version 0.9.7.1.",
      "doi:10.32614/CRAN.package.ruv"
    )
  ),
  quantitative_contract = list(
    sps_source_quantitative_meaning = "absolute processed phosphosite log2 abundance",
    sps_source_scale = "log2",
    selected_baseline_condition = "baseline",
    centering_operation = paste(
      "for each site and reference dataset, subtract the arithmetic mean of",
      "baseline_1 and baseline_2 from every sample on the log2 scale"
    ),
    centered_quantitative_meaning = "condition-relative log2 fold change",
    zero_semantics = "the per-site baseline-condition mean is exactly zero",
    matrix_passed_to_getSPS = "sps/*_condition_relative_log2.csv",
    matrix_passed_to_phospy = "the identical sps/*_condition_relative_log2.csv bytes",
    raw_absolute_input_is_not_valid_phospy_sps_input = TRUE
  ),
  sps_case = list(
    reference_dataset_count = length(sps_reference_ids),
    replicated_conditions = TRUE,
    deliberately_different_input_row_order = TRUE,
    shared_site_count = sps_count,
    requested_top_n = sps_top_n,
    selected_site_count = length(sps_selected_keys),
    minimum_datasets_per_site = 3L,
    all_sites_contribute_to_every_reference = TRUE,
    tie_policy = paste(
      "compare consensus scores and tie-group membership; PhosPy additionally",
      "guarantees site_key ascending order within exact ties"
    )
  ),
  contract_differences = list(
    partial_reference_contributions = paste(
      "PhosPy can rank a site with exactly minimum_datasets_per_site references.",
      "PhosR getSPS 1.13.1 forms one rectangular score table and propagates NA",
      "for sites absent from a reference, so external parity is asserted only",
      "for sites contributing to every reference. The PhosPy-only boundary is",
      "covered by native contract/science tests and is not weakened here."
    ),
    overlap_floor = paste(
      "PhosR getSPS rejects fewer than 200 sites shared by at least two datasets;",
      "PhosPy exposes minimum_shared_sites explicitly. This fixture uses 210."
    ),
    missing_values = paste(
      "The low-level PhosPy RUV-III kernel and this parity lane require complete",
      "finite input. ruv::RUVIII warns that missing values are unsupported. No",
      "missing-value parity is claimed."
    ),
    ambiguous_sps_input = paste(
      "PhosR leaves scale and baseline preparation to caller convention; PhosPy",
      "requires established condition-relative log2 evidence and rejects an",
      "unasserted absolute-abundance matrix."
    )
  ),
  comparison_policy = list(
    sps_stability_and_consensus_absolute_tolerance = 1e-12,
    sps_selection = "exact site set and exact ordering outside exact score ties",
    ruv_corrected_matrix_absolute_tolerance = 1e-10,
    ruv_corrected_matrix_relative_tolerance = 1e-12,
    diagnostics_absolute_tolerance = 1e-10,
    rationale = paste(
      "rank/selection is discrete; scalar SPS probability implementations agree",
      "to near machine precision; RUV-III uses SVD and linear solves whose last",
      "digits can vary across BLAS/LAPACK implementations"
    )
  ),
  redistribution = list(
    status = "approved_for_repository_test_fixture_redistribution",
    origin = paste(
      "small deterministic synthetic inputs created by this project; no PhosR",
      "or external biological reference dataset is included"
    ),
    package_policy = paste(
      "PhosR and ruv are invoked as installed black-box implementations; package",
      "source code is not copied into the fixture"
    )
  ),
  ruv_case = list(
    selected_negative_control_count = length(ruv_controls),
    replicate_set_count = length(unique(ruv_sets)),
    k_values = c(0L, 1L, 2L),
    unwanted_factor_count = 2L,
    biological_signal_present = TRUE,
    sample_order_permutation = TRUE
  )
)
write_json_lf(metadata, file.path(outdir, "METADATA.json"))

provenance_lines <- c(
  "# PhosR SPS and ruv RUV-III external reference evidence",
  "",
  paste("Generated:", args$timestamp),
  paste("R:", observed$r),
  paste("PhosR:", observed$PhosR, "commit", observed$PhosR_remote_sha),
  paste("ruv:", observed$ruv),
  "",
  "The inputs are deterministic synthetic values authored for PhosPy. No external",
  "biological dataset or bundled PhosR reference dataset is redistributed.",
  "",
  "SPS source files contain absolute processed log2 abundance. For every site, the",
  "mean of baseline_1 and baseline_2 was subtracted from all four samples. The",
  "resulting condition-relative log2 files, not the absolute files, were supplied",
  "unchanged to both PhosR getSPS and the PhosPy parity test.",
  "",
  "RUV-III consumes the checked-in, already-selected negative controls. It never",
  "reruns SPS. Corrected matrices are compared, while latent factors are not",
  "compared directly because SVD coordinates admit sign/rotation invariances.",
  "",
  "Regeneration is separate from CI and requires the exact environment recorded in",
  "REFERENCE_ENVIRONMENT.json. The generator refuses version drift by default.",
  "Normal test execution neither starts R, installs packages, nor accesses a network.",
  "PhosPy has no runtime dependency on R.",
  "",
  "Scientific citations:",
  "- Kim HJ et al. (2021), Cell Reports 34(8):108771, doi:10.1016/j.celrep.2021.108771.",
  "- Gagnon-Bartsch J (2019), ruv 0.9.7.1, doi:10.32614/CRAN.package.ruv."
)
write_lf(file.path(outdir, "PROVENANCE.md"), provenance_lines)

generator_path <- normalizePath(
  "scripts/active/generate_phosr_ruv_parity_fixtures.R",
  winslash = "/",
  mustWork = TRUE
)
all_files <- list.files(outdir, recursive = TRUE, full.names = TRUE)
all_files <- all_files[basename(all_files) != "MANIFEST.json"]
relative_files <- substring(all_files, nchar(outdir) + 2L)
relative_files <- gsub("\\\\", "/", relative_files)
order_index <- order(relative_files)
relative_files <- relative_files[order_index]
all_files <- all_files[order_index]
file_entries <- lapply(seq_along(all_files), function(i) {
  list(
    relative_path = relative_files[[i]],
    sha256 = digest::digest(file = all_files[[i]], algo = "sha256", serialize = FALSE)
  )
})
manifest <- list(
  manifest_schema_version = "fixture-manifest-v1",
  fixture_family = "phosr_sps_ruv_iii",
  classification = "external_parity",
  external_implementation = list(
    name = "PhosR getSPS and ruv RUVIII",
    r_version = observed$r,
    PhosR_version = observed$PhosR,
    PhosR_commit = observed$PhosR_remote_sha,
    ruv_version = observed$ruv
  ),
  generator = "scripts/active/generate_phosr_ruv_parity_fixtures.R",
  generator_sha256 = digest::digest(
    file = generator_path, algo = "sha256", serialize = FALSE
  ),
  command = paste(
    "Rscript scripts/active/generate_phosr_ruv_parity_fixtures.R",
    "--outdir tests/fixtures/rewrite_parity/phosr_ruv",
    paste("--timestamp", args$timestamp),
    "--allow-unpinned-environment false"
  ),
  generation_timestamp_utc = args$timestamp,
  byte_policy = "utf-8 LF with final newline",
  source_policy = paste(
    "deterministic synthetic project-authored inputs; external package outputs",
    "are serialized directly from the pinned R run"
  ),
  files = file_entries
)
write_json_lf(manifest, file.path(outdir, "MANIFEST.json"))
