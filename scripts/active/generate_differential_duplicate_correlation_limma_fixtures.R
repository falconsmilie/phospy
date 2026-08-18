#!/usr/bin/env Rscript

required_pkgs <- c("BiocManager", "limma")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, then rerun this script.",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(limma))

FIXTURE_FAMILY <- "differential_duplicate_correlation_limma"
DEFAULT_OUTDIR <- "tests/fixtures/rewrite_parity/differential_duplicate_correlation"
DEFAULT_SEED <- "20260818"
DEFAULT_TIMESTAMP <- "2026-08-18T00:00:00Z"
PINNED_R_VERSION <- "R version 4.5.2 (2025-10-31 ucrt)"
PINNED_BIOCONDUCTOR_VERSION <- "3.22"
PINNED_LIMMA_VERSION <- "3.66.0"
CANONICAL_TEXT_BYTE_POLICY <- "utf-8 LF with final newline"
SERIALIZATION_POLICY <- paste(
  "CSV uses comma separators, a header row, row.names=FALSE, UTF-8, LF line",
  "endings, a final newline, options(digits=17, scipen=999), and the literal",
  "NA token for missing numeric values. R NaN and +/-Inf numeric outputs are",
  "serialized as NA; feature-correlation failure details are retained in",
  "companion status and missing-kind columns. JSON manifests use stable key",
  "ordering and UTF-8 LF bytes."
)
SCIENTIFIC_CITATION <- c(
  "Smyth GK (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology 3(1), Article 3.",
  "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47."
)
REDISTRIBUTION_METADATA <- list(
  status = "approved_for_repository_test_fixture_redistribution",
  fixture_data_origin = "synthetic deterministic inputs generated locally; no third-party biological dataset is included",
  external_package_source_policy = "limma is invoked as an installed black-box scientific implementation; limma source code is not copied into this repository",
  scope = "exact generated CSV/JSON/Markdown fixture files listed in MANIFEST.json",
  attribution = "scientific citations recorded in scientific_citation and PROVENANCE.md"
)

read_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    outdir = DEFAULT_OUTDIR,
    seed = DEFAULT_SEED,
    timestamp = DEFAULT_TIMESTAMP,
    "manifest-outdir-label" = NA_character_,
    "allow-unpinned-environment" = "false"
  )

  if (length(args) == 0) {
    return(defaults)
  }
  if (length(args) %% 2 != 0) {
    stop("Arguments must be provided as --key value pairs.", call. = FALSE)
  }

  parsed <- defaults
  i <- 1
  while (i <= length(args)) {
    key <- gsub("^--", "", args[[i]])
    value <- args[[i + 1]]
    if (!(key %in% names(defaults))) {
      stop("Unknown argument: --", key, call. = FALSE)
    }
    parsed[[key]] <- value
    i <- i + 2
  }
  parsed
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(NA_character_)
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
}

canonical_text <- function(text) {
  normalized <- gsub("\r\n", "\n", text, fixed = TRUE)
  normalized <- gsub("\r", "\n", normalized, fixed = TRUE)
  if (!endsWith(normalized, "\n")) {
    normalized <- paste0(normalized, "\n")
  }
  enc2utf8(normalized)
}

write_canonical_text <- function(text, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(canonical_text(text)), con, useBytes = TRUE)
}

write_canonical_lines <- function(lines, path) {
  write_canonical_text(paste(lines, collapse = "\n"), path)
}

canonicalize_numeric_missing <- function(data) {
  for (column in names(data)) {
    if (is.numeric(data[[column]])) {
      values <- data[[column]]
      values[!is.finite(values)] <- NA_real_
      data[[column]] <- values
    }
  }
  data
}

write_canonical_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  utils::write.table(
    canonicalize_numeric_missing(data),
    con,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE,
    qmethod = "double",
    na = "NA",
    eol = "\n"
  )
}

json_escape <- function(value) {
  value <- enc2utf8(as.character(value))
  value <- gsub("\\\\", "\\\\\\\\", value)
  value <- gsub("\"", "\\\\\"", value)
  value <- gsub("\n", "\\\\n", value, fixed = TRUE)
  value <- gsub("\r", "\\\\r", value, fixed = TRUE)
  value <- gsub("\t", "\\\\t", value, fixed = TRUE)
  paste0("\"", value, "\"")
}

json_scalar <- function(value) {
  if (length(value) == 0 || is.null(value)) {
    return("null")
  }
  if (length(value) != 1) {
    stop("json_scalar received a non-scalar value", call. = FALSE)
  }
  if (is.na(value)) {
    return("null")
  }
  if (is.logical(value)) {
    return(ifelse(value, "true", "false"))
  }
  if (is.numeric(value)) {
    if (!is.finite(value)) {
      return("null")
    }
    return(format(value, digits = 17, scientific = FALSE, trim = TRUE))
  }
  json_escape(value)
}

json_value <- function(value, indent = 0) {
  indent_text <- paste(rep(" ", indent), collapse = "")
  child_indent <- paste(rep(" ", indent + 2), collapse = "")
  if (is.null(value)) {
    return("null")
  }
  if (is.atomic(value) && is.null(names(value))) {
    if (length(value) == 1) {
      return(json_scalar(value))
    }
    if (length(value) == 0) {
      return("[]")
    }
    return(paste0(
      "[",
      paste(vapply(value, json_scalar, character(1)), collapse = ", "),
      "]"
    ))
  }
  if (is.list(value) && is.null(names(value))) {
    if (length(value) == 0) {
      return("[]")
    }
    items <- vapply(
      value,
      function(item) paste0(child_indent, json_value(item, indent + 2)),
      character(1)
    )
    return(paste0("[\n", paste(items, collapse = ",\n"), "\n", indent_text, "]"))
  }
  if (is.list(value)) {
    keys <- names(value)
    if (is.null(keys) || any(!nzchar(keys))) {
      stop("named JSON objects require non-empty names", call. = FALSE)
    }
    items <- vapply(
      seq_along(value),
      function(i) {
        paste0(
          child_indent,
          json_escape(keys[[i]]),
          ": ",
          json_value(value[[i]], indent + 2)
        )
      },
      character(1)
    )
    return(paste0("{\n", paste(items, collapse = ",\n"), "\n", indent_text, "}"))
  }
  json_scalar(value)
}

write_canonical_json <- function(payload, path) {
  write_canonical_text(json_value(payload, 0), path)
}

sha256_text <- function(text) {
  path <- tempfile(pattern = "phospy-duplicate-correlation-hash-", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  write_canonical_text(text, path)
  unname(tools::sha256sum(path))
}

aggregate_role_hash <- function(file_table, roles) {
  selected <- file_table[file_table$role %in% roles, , drop = FALSE]
  selected <- selected[order(selected$relative_path), , drop = FALSE]
  text <- paste(
    paste(selected$relative_path, selected$sha256, sep = "\t"),
    collapse = "\n"
  )
  sha256_text(text)
}

matrix_csv <- function(mat) {
  data.frame(feature_id = rownames(mat), mat, check.names = FALSE)
}

matrix_with_row_id_csv <- function(row_id_name, mat) {
  data.frame(setNames(list(rownames(mat)), row_id_name), mat, check.names = FALSE)
}

feature_vector_csv <- function(feature_ids, ...) {
  data.frame(feature_id = feature_ids, ..., check.names = FALSE)
}

sample_metadata_frame <- function(sample_ids, block_ids, conditions, extra = NULL) {
  data <- data.frame(
    sample_order = seq_along(sample_ids),
    sample_id = sample_ids,
    block_id = block_ids,
    condition = as.character(conditions),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(extra)) {
    data <- cbind(data, extra, stringsAsFactors = FALSE)
  }
  rownames(data) <- sample_ids
  data
}

block_vector_frame <- function(sample_metadata) {
  block_sizes <- table(sample_metadata$block_id)
  data.frame(
    sample_order = sample_metadata$sample_order,
    sample_id = sample_metadata$sample_id,
    block_id = sample_metadata$block_id,
    block_size = as.integer(block_sizes[sample_metadata$block_id]),
    singleton_block = as.integer(block_sizes[sample_metadata$block_id]) == 1L,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

design_without_block_dummies <- function(condition, sample_ids, levels) {
  condition <- factor(condition, levels = levels)
  design <- model.matrix(~0 + condition)
  colnames(design) <- levels
  rownames(design) <- sample_ids
  design
}

feature_correlation_table <- function(fixture, mat, design, dc) {
  feature_ids <- rownames(mat)
  atanh_values <- as.numeric(dc$atanh.correlations)
  missing_kind <- ifelse(
    is.nan(atanh_values),
    "NaN",
    ifelse(is.na(atanh_values), "NA", "finite")
  )
  status <- ifelse(is.finite(atanh_values), "estimated", "missing")
  observed_count <- rowSums(!is.na(mat))
  finite_atanh <- ifelse(is.finite(atanh_values), atanh_values, NA_real_)
  correlation <- ifelse(is.finite(atanh_values), tanh(atanh_values), NA_real_)
  data.frame(
    feature_order = seq_along(feature_ids),
    feature_id = feature_ids,
    atanh_correlation = finite_atanh,
    correlation = correlation,
    status = status,
    atanh_correlation_missing_kind = missing_kind,
    observed_value_count = as.integer(observed_count),
    design_rank = qr(design)$rank,
    fixture = fixture,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

duplicate_correlation_summary <- function(fixture, mat, sample_metadata, dc) {
  block_sizes <- as.integer(table(sample_metadata$block_id))
  max_block_size <- max(block_sizes)
  lower_bound <- ifelse(max_block_size > 1, -1 / (max_block_size - 1), NA_real_)
  atanh_values <- as.numeric(dc$atanh.correlations)
  finite_values <- atanh_values[is.finite(atanh_values)]
  feature_correlations <- tanh(finite_values)
  data.frame(
    field = c(
      "fixture",
      "sample_count",
      "feature_count",
      "block_count",
      "repeated_block_count",
      "singleton_block_count",
      "min_block_size",
      "max_block_size",
      "compound_symmetry_lower_bound",
      "trim_fraction_each_tail",
      "estimated_feature_correlation_count",
      "missing_feature_correlation_count",
      "min_feature_correlation",
      "max_feature_correlation",
      "consensus_correlation",
      "consensus_atanh_correlation"
    ),
    value = c(
      fixture,
      ncol(mat),
      nrow(mat),
      length(block_sizes),
      sum(block_sizes > 1L),
      sum(block_sizes == 1L),
      min(block_sizes),
      max_block_size,
      lower_bound,
      0.15,
      length(finite_values),
      length(atanh_values) - length(finite_values),
      ifelse(length(feature_correlations), min(feature_correlations), NA_real_),
      ifelse(length(feature_correlations), max(feature_correlations), NA_real_),
      as.numeric(dc$consensus.correlation),
      atanh(as.numeric(dc$consensus.correlation))
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

fit_sigma_df_table <- function(fit) {
  sigma <- as.numeric(fit$sigma)
  data.frame(
    feature_id = rownames(fit$coefficients),
    sigma = sigma,
    residual_variance = sigma^2,
    residual_degrees_of_freedom = as.numeric(fit$df.residual),
    average_expression = as.numeric(fit$Amean),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

contrast_statistics_table <- function(eb_fit) {
  feature_ids <- rownames(eb_fit$coefficients)
  contrast_names <- colnames(eb_fit$coefficients)
  rows <- list()
  row_index <- 1L
  for (contrast in contrast_names) {
    p_values <- as.numeric(eb_fit$p.value[, contrast])
    adjusted <- p.adjust(p_values, method = "BH")
    standard_error <- as.numeric(
      eb_fit$stdev.unscaled[, contrast] * sqrt(eb_fit$s2.post)
    )
    rows[[row_index]] <- data.frame(
      feature_id = feature_ids,
      contrast = contrast,
      logFC = as.numeric(eb_fit$coefficients[, contrast]),
      AveExpr = as.numeric(eb_fit$Amean),
      SE = standard_error,
      t = as.numeric(eb_fit$t[, contrast]),
      P.Value = p_values,
      adj.P.Val = adjusted,
      B = as.numeric(eb_fit$lods[, contrast]),
      stdev.unscaled = as.numeric(eb_fit$stdev.unscaled[, contrast]),
      sigma = as.numeric(eb_fit$sigma),
      residual_degrees_of_freedom = as.numeric(eb_fit$df.residual),
      s2.prior = recycle_to_features(eb_fit$s2.prior, length(feature_ids)),
      df.prior = recycle_to_features(eb_fit$df.prior, length(feature_ids)),
      s2.post = as.numeric(eb_fit$s2.post),
      df.total = recycle_to_features(eb_fit$df.total, length(feature_ids)),
      var.prior = recycle_to_contrasts(eb_fit$var.prior, contrast_names)[[contrast]],
      proportion = as.numeric(eb_fit$proportion),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    row_index <- row_index + 1L
  }
  do.call(rbind, rows)
}

recycle_to_features <- function(value, n_features) {
  numeric_value <- as.numeric(value)
  if (length(numeric_value) == n_features) {
    return(numeric_value)
  }
  rep(numeric_value[[1]], n_features)
}

recycle_to_contrasts <- function(value, contrast_names) {
  numeric_value <- as.numeric(value)
  if (length(numeric_value) == length(contrast_names)) {
    return(setNames(as.list(numeric_value), contrast_names))
  }
  setNames(as.list(rep(numeric_value[[1]], length(contrast_names))), contrast_names)
}

ebayes_priors_table <- function(eb_fit) {
  feature_ids <- rownames(eb_fit$coefficients)
  sigma <- as.numeric(eb_fit$sigma)
  data.frame(
    feature_id = feature_ids,
    sigma = sigma,
    residual_variance = sigma^2,
    residual_degrees_of_freedom = as.numeric(eb_fit$df.residual),
    s2.prior = recycle_to_features(eb_fit$s2.prior, length(feature_ids)),
    df.prior = recycle_to_features(eb_fit$df.prior, length(feature_ids)),
    s2.post = as.numeric(eb_fit$s2.post),
    df.total = recycle_to_features(eb_fit$df.total, length(feature_ids)),
    average_expression = as.numeric(eb_fit$Amean),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

ebayes_global_table <- function(eb_fit) {
  data.frame(
    field = c("method", "trend", "robust", "proportion", "contrast_count"),
    value = c(
      "standard",
      "FALSE",
      "FALSE",
      format(as.numeric(eb_fit$proportion), digits = 17, scientific = FALSE),
      ncol(eb_fit$coefficients)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

make_cs_residual <- function(size, rho, sigma) {
  covariance <- matrix(rho, nrow = size, ncol = size)
  diag(covariance) <- 1
  decomposition <- eigen(covariance, symmetric = TRUE)
  root <- decomposition$vectors %*% diag(sqrt(pmax(decomposition$values, 0)))
  as.numeric(root %*% rnorm(size, mean = 0, sd = sigma))
}

build_fixture_a <- function() {
  blocks <- sprintf("pair_%02d", 1:6)
  sample_ids <- as.vector(rbind(paste0(blocks, "_A"), paste0(blocks, "_B")))
  block_ids <- rep(blocks, each = 2)
  condition <- rep(c("A", "B"), times = length(blocks))
  sample_metadata <- sample_metadata_frame(sample_ids, block_ids, condition)
  design <- design_without_block_dummies(condition, sample_ids, c("A", "B"))
  contrast_mat <- makeContrasts(B_vs_A = B - A, levels = design)

  feature_ids <- sprintf("FXA_%02d", 1:10)
  base <- c(9.1, 10.2, 8.8, 11.0, 7.5, 9.7, 10.6, 8.2, 9.9, 11.5)
  effect <- c(0.0, 0.0, 0.85, 1.10, -0.70, 0.0, 1.35, -1.05, 0.45, 0.0)
  block_offset <- c(-0.30, 0.12, 0.28, -0.18, 0.08, -0.06)
  mat <- matrix(NA_real_, nrow = length(feature_ids), ncol = length(sample_ids))
  for (i in seq_along(feature_ids)) {
    for (j in seq_along(sample_ids)) {
      block_index <- match(block_ids[[j]], blocks)
      deterministic_noise <- 0.025 * sin(i * 1.7 + j * 0.9) +
        0.015 * cos(i + block_index)
      mat[i, j] <- base[[i]] +
        ifelse(condition[[j]] == "B", effect[[i]], 0.0) +
        block_offset[[block_index]] +
        deterministic_noise
    }
  }
  rownames(mat) <- feature_ids
  colnames(mat) <- sample_ids
  feature_metadata <- data.frame(
    feature_order = seq_along(feature_ids),
    feature_id = feature_ids,
    controlled_case = ifelse(effect == 0.0, "null", "non_null"),
    true_B_minus_A = effect,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(
    id = "fixture_a_complete_pairs",
    title = "Fixture A - complete two-condition pairs",
    matrix = mat,
    sample_metadata = sample_metadata,
    design = design,
    contrasts = contrast_mat,
    feature_metadata = feature_metadata,
    contrasts_description = c("B_vs_A = B - A"),
    requirement_notes = c(
      "six complete paired blocks",
      "two observations per block",
      "complete quantitative matrix",
      "one simple treatment contrast",
      "contains null and non-null features"
    )
  )
}

build_fixture_b <- function(seed) {
  set.seed(seed + 101L)
  blocks <- sprintf("subject_%02d", 1:8)
  condition_levels <- c("T0", "T1", "T2")
  sample_ids <- as.vector(
    unlist(lapply(blocks, function(block) paste0(block, "_", condition_levels)))
  )
  block_ids <- rep(blocks, each = length(condition_levels))
  condition <- rep(condition_levels, times = length(blocks))
  sample_metadata <- sample_metadata_frame(
    sample_ids,
    block_ids,
    condition,
    extra = data.frame(
      time_point = condition,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
  design <- design_without_block_dummies(condition, sample_ids, condition_levels)
  contrast_mat <- makeContrasts(
    T1_vs_T0 = T1 - T0,
    T2_vs_T0 = T2 - T0,
    levels = design
  )

  feature_ids <- sprintf("FXB_%02d", 1:14)
  base <- seq(7.6, 11.5, length.out = length(feature_ids))
  effect_t1 <- c(0, 0.35, 0.0, 0.55, -0.30, 0.0, 0.80, -0.45, 0.15, 0.0, 0.65, -0.20, 0.0, 0.45)
  effect_t2 <- c(0, 0.70, 0.0, 1.10, -0.55, 0.0, 1.25, -0.85, 0.35, 0.0, 0.95, -0.40, 0.0, 0.75)
  sigma <- seq(0.10, 0.22, length.out = length(feature_ids))
  rho <- -0.47
  mat <- matrix(NA_real_, nrow = length(feature_ids), ncol = length(sample_ids))
  for (i in seq_along(feature_ids)) {
    for (block_index in seq_along(blocks)) {
      residual <- make_cs_residual(length(condition_levels), rho, sigma[[i]])
      for (condition_index in seq_along(condition_levels)) {
        sample_index <- (block_index - 1L) * length(condition_levels) + condition_index
        condition_name <- condition_levels[[condition_index]]
        effect <- ifelse(
          condition_name == "T1",
          effect_t1[[i]],
          ifelse(condition_name == "T2", effect_t2[[i]], 0.0)
        )
        mat[i, sample_index] <- base[[i]] + effect + residual[[condition_index]]
      }
    }
  }
  rownames(mat) <- feature_ids
  colnames(mat) <- sample_ids
  feature_metadata <- data.frame(
    feature_order = seq_along(feature_ids),
    feature_id = feature_ids,
    controlled_case = ifelse(effect_t1 == 0.0 & effect_t2 == 0.0, "null", "time_effect"),
    true_T1_minus_T0 = effect_t1,
    true_T2_minus_T0 = effect_t2,
    simulated_within_block_correlation = rho,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(
    id = "fixture_b_three_observation_blocks",
    title = "Fixture B - more than two observations per block",
    matrix = mat,
    sample_metadata = sample_metadata,
    design = design,
    contrasts = contrast_mat,
    feature_metadata = feature_metadata,
    contrasts_description = c("T1_vs_T0 = T1 - T0", "T2_vs_T0 = T2 - T0"),
    requirement_notes = c(
      "eight blocks with three observations per block",
      "three time points",
      "two contrasts",
      "simulated near the valid compound-symmetry lower bound for block size three"
    )
  )
}

build_fixture_c <- function() {
  sample_ids <- c(
    "ib01_A",
    "ib01_B",
    "ib02_A",
    "ib02_B",
    "ib02_C",
    "ib03_A",
    "ib04_B",
    "ib04_C",
    "ib05_A",
    "ib05_C1",
    "ib05_C2",
    "ib06_C"
  )
  block_ids <- c(
    "ib01",
    "ib01",
    "ib02",
    "ib02",
    "ib02",
    "ib03",
    "ib04",
    "ib04",
    "ib05",
    "ib05",
    "ib05",
    "ib06"
  )
  condition <- c("A", "B", "A", "B", "C", "A", "B", "C", "A", "C", "C", "C")
  condition_levels <- c("A", "B", "C")
  sample_metadata <- sample_metadata_frame(sample_ids, block_ids, condition)
  design <- design_without_block_dummies(condition, sample_ids, condition_levels)
  contrast_mat <- makeContrasts(B_vs_A = B - A, C_vs_A = C - A, levels = design)

  feature_ids <- sprintf("FXC_%02d", 1:9)
  base <- c(8.5, 9.0, 10.2, 7.8, 9.7, 8.9, 10.8, 7.4, 9.3)
  effect_b <- c(0.0, 0.45, -0.25, 0.0, 0.80, -0.55, 0.0, 0.30, -0.20)
  effect_c <- c(0.0, 0.70, -0.40, 0.55, 1.05, -0.75, 0.0, 0.50, -0.35)
  block_offset <- c(0.20, -0.15, 0.05, -0.28, 0.18, -0.08)
  names(block_offset) <- unique(block_ids)
  mat <- matrix(NA_real_, nrow = length(feature_ids), ncol = length(sample_ids))
  for (i in seq_along(feature_ids)) {
    for (j in seq_along(sample_ids)) {
      condition_name <- condition[[j]]
      condition_effect <- ifelse(
        condition_name == "B",
        effect_b[[i]],
        ifelse(condition_name == "C", effect_c[[i]], 0.0)
      )
      mat[i, j] <- base[[i]] +
        condition_effect +
        block_offset[[block_ids[[j]]]] +
        0.035 * sin(i * j / 2.0) +
        0.018 * cos(i + j)
    }
  }
  rownames(mat) <- feature_ids
  colnames(mat) <- sample_ids
  feature_metadata <- data.frame(
    feature_order = seq_along(feature_ids),
    feature_id = feature_ids,
    controlled_case = ifelse(effect_b == 0.0 & effect_c == 0.0, "null", "condition_effect"),
    true_B_minus_A = effect_b,
    true_C_minus_A = effect_c,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(
    id = "fixture_c_incomplete_unequal_blocks",
    title = "Fixture C - incomplete and unequal blocks",
    matrix = mat,
    sample_metadata = sample_metadata,
    design = design,
    contrasts = contrast_mat,
    feature_metadata = feature_metadata,
    contrasts_description = c("B_vs_A = B - A", "C_vs_A = C - A"),
    requirement_notes = c(
      "unequal repeated blocks",
      "singleton blocks retained",
      "incomplete condition coverage by block",
      "full-rank fixed-effects design without block columns"
    )
  )
}

build_fixture_d <- function() {
  blocks <- sprintf("failure_block_%02d", 1:5)
  sample_ids <- as.vector(rbind(paste0(blocks, "_A"), paste0(blocks, "_B")))
  block_ids <- rep(blocks, each = 2)
  condition <- rep(c("A", "B"), times = length(blocks))
  sample_metadata <- sample_metadata_frame(sample_ids, block_ids, condition)
  design <- design_without_block_dummies(condition, sample_ids, c("A", "B"))
  contrast_mat <- makeContrasts(B_vs_A = B - A, levels = design)

  mat <- rbind(
    D_valid_null = c(10.0, 10.1, 10.2, 10.0, 9.9, 10.0, 10.1, 10.2, 9.8, 9.9),
    D_valid_effect = c(8.0, 9.0, 8.2, 9.3, 7.9, 8.8, 8.1, 9.2, 8.0, 9.1),
    D_missing_still_estimable = c(6.0, 6.4, 6.1, NA, 5.9, 6.3, 6.2, 6.5, 6.1, 6.4),
    D_constant_all = rep(5.0, 10),
    D_near_constant = 4 + c(0, 1e-10, -1e-10, 2e-10, 0, -2e-10, 1e-10, 0, -1e-10, 2e-10),
    D_rank_loss_only_A = c(7.0, NA, 7.1, NA, 7.2, NA, 7.1, NA, 7.0, NA),
    D_insufficient_one = c(NA, NA, NA, NA, 3.0, NA, NA, NA, NA, NA),
    D_all_missing = rep(NA_real_, 10)
  )
  colnames(mat) <- sample_ids
  feature_ids <- rownames(mat)
  feature_metadata <- data.frame(
    feature_order = seq_along(feature_ids),
    feature_id = feature_ids,
    controlled_case = c(
      "valid_null",
      "valid_non_null",
      "missing_observations_but_estimable",
      "constant_feature",
      "near_constant_feature",
      "feature_specific_rank_loss",
      "insufficient_observed_values",
      "all_observations_missing"
    ),
    true_B_minus_A = c(0.0, 1.0, NA, 0.0, 0.0, NA, NA, NA),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(
    id = "fixture_d_feature_level_failures",
    title = "Fixture D - feature-level estimator failures",
    matrix = mat,
    sample_metadata = sample_metadata,
    design = design,
    contrasts = contrast_mat,
    feature_metadata = feature_metadata,
    contrasts_description = c("B_vs_A = B - A"),
    requirement_notes = c(
      "valid features coexist with invalid feature-level correlation estimates",
      "missing feature observations",
      "constant and near-constant features",
      "feature-specific rank loss",
      "insufficient observed values",
      "consensus remains available from valid feature estimates"
    )
  )
}

write_fixture <- function(fixture, outdir) {
  fixture_dir <- file.path(outdir, fixture$id)
  dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)

  mat <- fixture$matrix
  sample_metadata <- fixture$sample_metadata
  design <- fixture$design
  contrast_mat <- fixture$contrasts
  block_vector <- sample_metadata$block_id

  dc <- duplicateCorrelation(mat, design, block = block_vector)
  fit <- lmFit(mat, design, block = block_vector, correlation = dc$consensus.correlation)
  contrast_fit <- contrasts.fit(fit, contrast_mat)
  eb_fit <- eBayes(contrast_fit, trend = FALSE, robust = FALSE)

  files <- list(
    list(path = file.path(fixture$id, "matrix.csv"), role = "input"),
    list(path = file.path(fixture$id, "sample_metadata.csv"), role = "input"),
    list(path = file.path(fixture$id, "blocks.csv"), role = "input"),
    list(path = file.path(fixture$id, "design.csv"), role = "input"),
    list(path = file.path(fixture$id, "contrasts.csv"), role = "input"),
    list(path = file.path(fixture$id, "feature_metadata.csv"), role = "input"),
    list(path = file.path(fixture$id, "feature_correlations.csv"), role = "intermediate"),
    list(path = file.path(fixture$id, "duplicate_correlation_summary.csv"), role = "intermediate"),
    list(path = file.path(fixture$id, "fit_coefficients.csv"), role = "intermediate"),
    list(path = file.path(fixture$id, "fit_stdev_unscaled.csv"), role = "intermediate"),
    list(path = file.path(fixture$id, "fit_cov_coefficients.csv"), role = "intermediate"),
    list(path = file.path(fixture$id, "fit_sigma_df.csv"), role = "intermediate"),
    list(path = file.path(fixture$id, "contrast_fit_coefficients.csv"), role = "intermediate"),
    list(path = file.path(fixture$id, "contrast_fit_stdev_unscaled.csv"), role = "intermediate"),
    list(path = file.path(fixture$id, "contrast_fit_cov_coefficients.csv"), role = "intermediate"),
    list(path = file.path(fixture$id, "ebayes_priors.csv"), role = "output"),
    list(path = file.path(fixture$id, "ebayes_statistics.csv"), role = "output"),
    list(path = file.path(fixture$id, "ebayes_global.csv"), role = "output")
  )

  write_canonical_csv(matrix_csv(mat), file.path(fixture_dir, "matrix.csv"))
  write_canonical_csv(sample_metadata, file.path(fixture_dir, "sample_metadata.csv"))
  write_canonical_csv(block_vector_frame(sample_metadata), file.path(fixture_dir, "blocks.csv"))
  write_canonical_csv(
    matrix_with_row_id_csv("sample_id", design),
    file.path(fixture_dir, "design.csv")
  )
  write_canonical_csv(
    matrix_with_row_id_csv("coefficient", contrast_mat),
    file.path(fixture_dir, "contrasts.csv")
  )
  write_canonical_csv(fixture$feature_metadata, file.path(fixture_dir, "feature_metadata.csv"))
  write_canonical_csv(
    feature_correlation_table(fixture$id, mat, design, dc),
    file.path(fixture_dir, "feature_correlations.csv")
  )
  write_canonical_csv(
    duplicate_correlation_summary(fixture$id, mat, sample_metadata, dc),
    file.path(fixture_dir, "duplicate_correlation_summary.csv")
  )
  write_canonical_csv(
    matrix_with_row_id_csv("feature_id", fit$coefficients),
    file.path(fixture_dir, "fit_coefficients.csv")
  )
  write_canonical_csv(
    matrix_with_row_id_csv("feature_id", fit$stdev.unscaled),
    file.path(fixture_dir, "fit_stdev_unscaled.csv")
  )
  write_canonical_csv(
    matrix_with_row_id_csv("coefficient", fit$cov.coefficients),
    file.path(fixture_dir, "fit_cov_coefficients.csv")
  )
  write_canonical_csv(fit_sigma_df_table(fit), file.path(fixture_dir, "fit_sigma_df.csv"))
  write_canonical_csv(
    matrix_with_row_id_csv("feature_id", contrast_fit$coefficients),
    file.path(fixture_dir, "contrast_fit_coefficients.csv")
  )
  write_canonical_csv(
    matrix_with_row_id_csv("feature_id", contrast_fit$stdev.unscaled),
    file.path(fixture_dir, "contrast_fit_stdev_unscaled.csv")
  )
  write_canonical_csv(
    matrix_with_row_id_csv("contrast", contrast_fit$cov.coefficients),
    file.path(fixture_dir, "contrast_fit_cov_coefficients.csv")
  )
  write_canonical_csv(ebayes_priors_table(eb_fit), file.path(fixture_dir, "ebayes_priors.csv"))
  write_canonical_csv(contrast_statistics_table(eb_fit), file.path(fixture_dir, "ebayes_statistics.csv"))
  write_canonical_csv(ebayes_global_table(eb_fit), file.path(fixture_dir, "ebayes_global.csv"))

  block_sizes <- as.integer(table(sample_metadata$block_id))
  max_block_size <- max(block_sizes)
  lower_bound <- ifelse(max_block_size > 1, -1 / (max_block_size - 1), NA_real_)
  atanh_values <- as.numeric(dc$atanh.correlations)
  finite_count <- sum(is.finite(atanh_values))
  missing_count <- length(atanh_values) - finite_count
  list(
    fixture = fixture,
    files = files,
    manifest_entry = list(
      id = fixture$id,
      title = fixture$title,
      sample_ids = as.character(sample_metadata$sample_id),
      feature_ids = rownames(mat),
      block_ids = unique(as.character(sample_metadata$block_id)),
      condition_levels = colnames(design),
      design_columns = colnames(design),
      contrast_names = as.list(as.character(colnames(contrast_mat))),
      contrasts = as.list(as.character(fixture$contrasts_description)),
      block_count = length(block_sizes),
      repeated_block_count = sum(block_sizes > 1L),
      singleton_block_count = sum(block_sizes == 1L),
      min_block_size = min(block_sizes),
      max_block_size = max_block_size,
      compound_symmetry_lower_bound = lower_bound,
      consensus_correlation = as.numeric(dc$consensus.correlation),
      estimated_feature_correlation_count = finite_count,
      missing_feature_correlation_count = missing_count,
      requirement_notes = fixture$requirement_notes
    )
  )
}

validate_pinned_environment <- function(allow_unpinned) {
  r_version <- R.version$version.string
  bioconductor_version <- as.character(BiocManager::version())
  limma_version <- as.character(packageVersion("limma"))
  if (!allow_unpinned) {
    mismatches <- c()
    if (!identical(r_version, PINNED_R_VERSION)) {
      mismatches <- c(mismatches, paste0("R version: expected ", PINNED_R_VERSION, " actual ", r_version))
    }
    if (!identical(bioconductor_version, PINNED_BIOCONDUCTOR_VERSION)) {
      mismatches <- c(
        mismatches,
        paste0(
          "Bioconductor version: expected ",
          PINNED_BIOCONDUCTOR_VERSION,
          " actual ",
          bioconductor_version
        )
      )
    }
    if (!identical(limma_version, PINNED_LIMMA_VERSION)) {
      mismatches <- c(
        mismatches,
        paste0("limma version: expected ", PINNED_LIMMA_VERSION, " actual ", limma_version)
      )
    }
    if (length(mismatches) > 0) {
      stop(
        "Pinned fixture environment mismatch:\n",
        paste(mismatches, collapse = "\n"),
        "\nUse the pinned environment or pass --allow-unpinned-environment true for exploratory local runs that must not be committed.",
        call. = FALSE
      )
    }
  }
  list(
    r_version = r_version,
    bioconductor_version = bioconductor_version,
    limma_version = limma_version
  )
}

main <- function() {
  options(digits = 17, scipen = 999)
  args <- read_args()
  outdir <- args$outdir
  seed <- as.integer(args$seed)
  timestamp <- args$timestamp
  allow_unpinned <- tolower(args[["allow-unpinned-environment"]]) %in% c("true", "1", "yes")
  if (is.na(seed)) {
    stop("--seed must be an integer", call. = FALSE)
  }
  env <- validate_pinned_environment(allow_unpinned)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  fixtures <- list(
    build_fixture_a(),
    build_fixture_b(seed),
    build_fixture_c(),
    build_fixture_d()
  )
  written <- lapply(fixtures, write_fixture, outdir = outdir)

  manifest_outdir_label <- args[["manifest-outdir-label"]]
  if (is.na(manifest_outdir_label) || !nzchar(manifest_outdir_label)) {
    manifest_outdir_label <- outdir
  }
  command <- paste(
    "Rscript scripts/active/generate_differential_duplicate_correlation_limma_fixtures.R",
    "--outdir", manifest_outdir_label,
    "--seed", seed,
    "--timestamp", timestamp,
    "--allow-unpinned-environment", tolower(as.character(allow_unpinned))
  )
  script <- script_path()
  script_sha256 <- ifelse(is.na(script), NA_character_, unname(tools::sha256sum(script)))

  provenance <- c(
    "# Differential Duplicate-Correlation Limma Fixture Provenance",
    "",
    paste0("Generated with ", env$r_version),
    paste0("Bioconductor version: ", env$bioconductor_version),
    paste0("limma version: ", env$limma_version),
    paste0("Seed: ", seed),
    paste0("Generation timestamp (UTC): ", timestamp),
    paste0("Command: `", command, "`"),
    paste0("Generator SHA-256: ", script_sha256),
    paste0("Byte policy: ", CANONICAL_TEXT_BYTE_POLICY),
    paste0("Serialization policy: ", SERIALIZATION_POLICY),
    "Source policy: deterministic synthetic fixture generated locally without network access; limma duplicateCorrelation, lmFit, contrasts.fit, and eBayes outputs are the external scientific authority for expected numerical columns.",
    "Classification: external parity for R/limma duplicate-correlation intermediate and final outputs.",
    "Expected outputs come only from the pinned R/limma run. PhosPy is not imported or executed by this generator.",
    "Redistribution metadata: synthetic deterministic inputs and black-box limma numeric outputs are repository test fixtures; limma source code is not redistributed.",
    "",
    "Scientific citations:",
    paste0("- ", SCIENTIFIC_CITATION),
    "",
    "Fixtures:",
    unlist(lapply(written, function(item) paste0("- ", item$fixture$title, " (`", item$fixture$id, "`)"))),
    "",
    "Each fixture uses a fixed-effects design matrix without block dummy variables and supplies block IDs only to limma duplicateCorrelation/lmFit.",
    "Output files and SHA-256 digests are listed in `MANIFEST.json`."
  )
  write_canonical_lines(provenance, file.path(outdir, "PROVENANCE.md"))

  file_rows <- list()
  row_index <- 1L
  for (item in written) {
    for (file in item$files) {
      file_rows[[row_index]] <- data.frame(
        fixture = item$fixture$id,
        role = file$role,
        relative_path = file$path,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  file_rows[[row_index]] <- data.frame(
    fixture = "family",
    role = "governance",
    relative_path = "PROVENANCE.md",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  file_table <- do.call(rbind, file_rows)
  file_table$sha256 <- unname(tools::sha256sum(file.path(outdir, file_table$relative_path)))

  manifest_files <- lapply(
    seq_len(nrow(file_table)),
    function(i) {
      list(
        fixture = file_table$fixture[[i]],
        role = file_table$role[[i]],
        relative_path = file_table$relative_path[[i]],
        sha256 = file_table$sha256[[i]]
      )
    }
  )
  fixture_entries <- lapply(written, function(item) item$manifest_entry)
  manifest <- list(
    manifest_schema_version = "fixture-manifest-v1",
    fixture_family = FIXTURE_FAMILY,
    classification = "external_parity",
    external_implementation = list(
      name = "R limma",
      r_version = env$r_version,
      bioconductor_version = env$bioconductor_version,
      limma_version = env$limma_version
    ),
    pinned_environment = list(
      r_version = PINNED_R_VERSION,
      bioconductor_version = PINNED_BIOCONDUCTOR_VERSION,
      limma_version = PINNED_LIMMA_VERSION
    ),
    generator = "scripts/active/generate_differential_duplicate_correlation_limma_fixtures.R",
    generator_sha256 = script_sha256,
    command = command,
    seed = seed,
    generation_timestamp_utc = timestamp,
    scientific_citation = SCIENTIFIC_CITATION,
    redistribution_metadata = REDISTRIBUTION_METADATA,
    source_policy = "synthetic deterministic fixture generated locally; no network access; limma outputs are the external scientific authority for expected numerical columns",
    numeric_authority = "expected numeric values are serialized directly from the pinned R/limma run and are not generated by PhosPy or manually edited after generation",
    byte_policy = CANONICAL_TEXT_BYTE_POLICY,
    serialization_policy = SERIALIZATION_POLICY,
    model_policy = list(
      duplicate_correlation_trim_fraction_each_tail = 0.15,
      fixed_effects_design_excludes_block_dummy_variables = TRUE,
      block_correlation_supplied_to_lmFit = TRUE,
      empirical_bayes = "limma::eBayes(trend = FALSE, robust = FALSE)"
    ),
    hashes = list(
      input_files_sha256 = aggregate_role_hash(file_table, c("input")),
      intermediate_files_sha256 = aggregate_role_hash(file_table, c("intermediate")),
      output_files_sha256 = aggregate_role_hash(file_table, c("output")),
      governance_files_sha256 = aggregate_role_hash(file_table, c("governance"))
    ),
    fixtures = fixture_entries,
    files = manifest_files
  )
  write_canonical_json(manifest, file.path(outdir, "MANIFEST.json"))

  message("Done. Duplicate-correlation limma fixtures written to: ", outdir)
}

main()
