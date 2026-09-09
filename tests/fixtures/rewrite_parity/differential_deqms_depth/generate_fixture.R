#!/usr/bin/env Rscript

required_pkgs <- c("BiocManager", "limma", "DEqMS")
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
suppressPackageStartupMessages(library(DEqMS))

FIXTURE_FAMILY <- "differential_deqms_depth"
DEFAULT_OUTDIR <- "tests/fixtures/rewrite_parity/differential_deqms_depth"
DEFAULT_SEED <- "20260909"
DEFAULT_TIMESTAMP <- "2026-09-09T00:00:00Z"
PINNED_R_VERSION <- "R version 4.5.2 (2025-10-31 ucrt)"
PINNED_BIOCONDUCTOR_VERSION <- "3.22"
PINNED_LIMMA_VERSION <- "3.66.0"
PINNED_DEQMS_VERSION <- "1.28.0"
CANONICAL_TEXT_BYTE_POLICY <- "utf-8 LF with final newline"
SERIALIZATION_POLICY <- paste(
  "CSV uses comma separators, a header row, row.names=FALSE, UTF-8, LF line",
  "endings, a final newline, options(digits=17, scipen=999), and the literal",
  "NA token for missing numeric values. R NaN and +/-Inf numeric outputs are",
  "serialized as NA. JSON manifests use stable key ordering and UTF-8 LF bytes."
)
SCIENTIFIC_CITATION <- c(
  "Zhu Y, Orre LM, Zhou Tran Y, Mermelekas G, Johansson HJ, Malyutina A, Anders S, Lehtio J (2020). DEqMS: A Method for Accurate Variance Estimation in Differential Protein Expression Analysis. Molecular & Cellular Proteomics 19(6), 1047-1057.",
  "Smyth GK (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology 3(1), Article 3.",
  "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47."
)
REDISTRIBUTION_METADATA <- list(
  status = "approved_for_repository_test_fixture_redistribution",
  fixture_data_origin = "synthetic deterministic inputs generated locally; no third-party biological dataset is included",
  external_package_source_policy = "DEqMS and limma are invoked as installed black-box scientific implementations; package source code is not copied into this repository",
  scope = "exact generated CSV/JSON/Markdown/R fixture files listed in MANIFEST.json",
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
  path <- tempfile(pattern = "phospy-deqms-depth-hash-", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  write_canonical_text(text, path)
  unname(tools::sha256sum(path))
}

aggregate_role_hash <- function(file_table, roles) {
  selected <- file_table[file_table$role %in% roles, , drop = FALSE]
  sort_keys <- vapply(
    selected$relative_path,
    function(path) {
      paste(sprintf("%03d", as.integer(charToRaw(enc2utf8(path)))), collapse = ".")
    },
    character(1)
  )
  selected <- selected[order(sort_keys), , drop = FALSE]
  text <- paste(
    paste(selected$relative_path, selected$sha256, sep = "\t"),
    collapse = "\n"
  )
  sha256_text(text)
}

copy_generator_to_outdir <- function(script, outdir) {
  if (is.na(script)) {
    return(NA_character_)
  }
  target <- file.path(outdir, "generate_fixture.R")
  source_norm <- normalizePath(script, winslash = "/", mustWork = TRUE)
  target_norm <- normalizePath(target, winslash = "/", mustWork = FALSE)
  source_text <- readChar(script, file.info(script)$size, useBytes = TRUE)
  canonical_source_text <- canonical_text(source_text)
  if (
    !identical(source_norm, target_norm) ||
      !identical(source_text, canonical_source_text)
  ) {
    write_canonical_text(canonical_source_text, target)
  }
  unname(tools::sha256sum(target))
}

matrix_csv <- function(mat) {
  data.frame(site_id = rownames(mat), mat, check.names = FALSE)
}

matrix_with_row_id_csv <- function(row_id_name, mat) {
  data.frame(setNames(list(rownames(mat)), row_id_name), mat, check.names = FALSE)
}

build_depth_fixture <- function(seed) {
  set.seed(seed)

  a_samples <- sprintf("A_%02d", seq_len(5))
  b_samples <- sprintf("B_%02d", seq_len(7))
  samples <- c(a_samples, b_samples)
  condition <- factor(
    c(rep("A", length(a_samples)), rep("B", length(b_samples))),
    levels = c("A", "B")
  )
  design <- model.matrix(~0 + condition)
  colnames(design) <- c("A", "B")
  rownames(design) <- samples

  n_features <- 144L
  feature_index <- seq_len(n_features)
  site_ids <- sprintf("QD%03d;S%d;", feature_index, feature_index)
  low_depth <- c(1L, 2L, 3L, 4L)
  high_depth <- c(16L, 24L, 32L, 48L)
  middle_depth <- c(6L, 8L, 10L, 12L)
  quantification_depth <- rep(
    c(low_depth, middle_depth, high_depth),
    length.out = n_features
  )

  comparable_pair_low <- 21L
  comparable_pair_high <- 117L
  quantification_depth[[comparable_pair_low]] <- 2L
  quantification_depth[[comparable_pair_high]] <- 32L

  log2_depth <- log2(quantification_depth)
  base_intensity <- 9.4 +
    0.42 * sin(feature_index / 8.0) +
    0.18 * cos(feature_index / 19.0)
  residual_sd <- exp(0.18 - 0.23 * log2_depth + 0.06 * sin(feature_index / 11.0))

  true_shift <- rep(0.0, n_features)
  true_shift[feature_index %% 17L == 0L] <- 0.58
  true_shift[feature_index %% 31L == 0L] <- -0.46

  mat <- matrix(NA_real_, nrow = n_features, ncol = length(samples))
  for (i in seq_along(feature_index)) {
    for (j in seq_along(samples)) {
      condition_shift <- ifelse(condition[[j]] == "B", true_shift[[i]], 0.0)
      replicate_noise <- rnorm(1, mean = 0.0, sd = residual_sd[[i]])
      deterministic_offset <- 0.012 * sin(i * 1.3 + j * 0.7)
      mat[i, j] <- base_intensity[[i]] + condition_shift +
        deterministic_offset + replicate_noise
    }
  }

  mat[comparable_pair_high, ] <- mat[comparable_pair_low, ]
  base_intensity[[comparable_pair_high]] <- base_intensity[[comparable_pair_low]]
  residual_sd[[comparable_pair_high]] <- residual_sd[[comparable_pair_low]]
  true_shift[[comparable_pair_high]] <- true_shift[[comparable_pair_low]]

  rownames(mat) <- site_ids
  colnames(mat) <- samples

  contrast_mat <- makeContrasts(B_vs_A = B - A, levels = design)
  list(
    matrix = mat,
    design = design,
    contrasts = contrast_mat,
    quantification_depth = quantification_depth,
    log2_depth = log2(quantification_depth),
    base_intensity = base_intensity,
    residual_sd = residual_sd,
    true_shift = true_shift,
    comparable_pair_low = site_ids[[comparable_pair_low]],
    comparable_pair_high = site_ids[[comparable_pair_high]],
    depth_groups = ifelse(
      quantification_depth <= 4L,
      "lower_depth",
      ifelse(quantification_depth >= 16L, "higher_depth", "middle_depth")
    )
  )
}

deqms_reference_table <- function(fit) {
  coef_col <- "B_vs_A"
  coef_index <- 1L
  feature_ids <- rownames(fit$coefficients)
  sca_p <- as.numeric(fit$sca.p[, coef_index])
  data.frame(
    site_id = feature_ids,
    logFC = as.numeric(fit$coefficients[, coef_col]),
    AveExpr = as.numeric(fit$Amean),
    stdev.unscaled = as.numeric(fit$stdev.unscaled[, coef_col]),
    sigma = as.numeric(fit$sigma),
    residual_variance = as.numeric(fit$sigma)^2,
    quantification_depth = as.numeric(fit$count[feature_ids]),
    log2_quantification_depth = log2(as.numeric(fit$count[feature_ids])),
    limma.t = as.numeric(fit$t[, coef_col]),
    limma.P.Value = as.numeric(fit$p.value[, coef_col]),
    limma.adj.P.Val = p.adjust(as.numeric(fit$p.value[, coef_col]), method = "BH"),
    sca.priorvar = as.numeric(fit$sca.priorvar),
    sca.dfprior = rep(as.numeric(fit$sca.dfprior), length(feature_ids)),
    sca.postvar = as.numeric(fit$sca.postvar),
    sca.t = as.numeric(fit$sca.t[, coef_index]),
    sca.P.Value = sca_p,
    sca.adj.P.Val = p.adjust(sca_p, method = "BH"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

deqms_model_table <- function(fit) {
  feature_ids <- rownames(fit$coefficients)
  data.frame(
    site_id = feature_ids,
    quantification_depth = as.numeric(fit$count[feature_ids]),
    log2_quantification_depth = log2(as.numeric(fit$count[feature_ids])),
    log_residual_variance = log(as.numeric(fit$sigma)^2),
    deqms_fitted_log_prior_variance = log(as.numeric(fit$sca.priorvar)),
    loess_span = rep(0.75, length(feature_ids)),
    fit_method = rep(as.character(fit$fit.method), length(feature_ids)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

deqms_global_table <- function(fit) {
  data.frame(
    field = c(
      "method",
      "fit_method",
      "trend_covariate",
      "trend_covariate_transformation",
      "quantification_depth_kind",
      "sca.dfprior",
      "contrast_count",
      "feature_count"
    ),
    value = c(
      "spectraCounteBayes",
      as.character(fit$fit.method),
      "quantification_depth",
      "log2",
      "psm_count",
      format(as.numeric(fit$sca.dfprior), digits = 17, scientific = FALSE),
      ncol(fit$coefficients),
      nrow(fit$coefficients)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

feature_metadata_table <- function(fixture) {
  feature_ids <- rownames(fixture$matrix)
  comparable_pair_id <- rep("", length(feature_ids))
  comparable_pair_id[feature_ids %in% c(
    fixture$comparable_pair_low,
    fixture$comparable_pair_high
  )] <- "same_intensity_different_depth_pair"
  data.frame(
    feature_order = seq_along(feature_ids),
    site_id = feature_ids,
    quantification_depth = as.numeric(fixture$quantification_depth),
    quantification_depth_kind = "psm_count",
    depth_group = fixture$depth_groups,
    log2_quantification_depth = as.numeric(fixture$log2_depth),
    true_B_minus_A = fixture$true_shift,
    controlled_case = ifelse(
      comparable_pair_id != "",
      comparable_pair_id,
      ifelse(fixture$true_shift == 0.0, "null", "shifted")
    ),
    comparable_pair_id = comparable_pair_id,
    synthetic_base_intensity = fixture$base_intensity,
    synthetic_residual_sd = fixture$residual_sd,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

depth_table <- function(fixture) {
  data.frame(
    site_id = rownames(fixture$matrix),
    quantification_depth = as.numeric(fixture$quantification_depth),
    quantification_depth_kind = "psm_count",
    depth_group = fixture$depth_groups,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

validate_pinned_environment <- function(allow_unpinned) {
  r_version <- R.version$version.string
  bioconductor_version <- as.character(BiocManager::version())
  limma_version <- as.character(packageVersion("limma"))
  deqms_version <- as.character(packageVersion("DEqMS"))
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
    if (!identical(deqms_version, PINNED_DEQMS_VERSION)) {
      mismatches <- c(
        mismatches,
        paste0("DEqMS version: expected ", PINNED_DEQMS_VERSION, " actual ", deqms_version)
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
    limma_version = limma_version,
    deqms_version = deqms_version
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

  script <- script_path()
  script_sha256 <- copy_generator_to_outdir(script, outdir)
  fixture <- build_depth_fixture(seed)
  fit <- lmFit(fixture$matrix, fixture$design)
  contrast_fit <- contrasts.fit(fit, fixture$contrasts)
  limma_fit <- eBayes(contrast_fit, trend = FALSE, robust = FALSE)
  limma_fit$count <- fixture$quantification_depth
  deqms_fit <- spectraCounteBayes(limma_fit, fit.method = "loess", coef_col = 1)

  write_canonical_csv(matrix_csv(fixture$matrix), file.path(outdir, "matrix.csv"))
  write_canonical_csv(
    matrix_with_row_id_csv("sample", fixture$design),
    file.path(outdir, "design.csv")
  )
  write_canonical_csv(
    matrix_with_row_id_csv("coefficient", fixture$contrasts),
    file.path(outdir, "contrasts.csv")
  )
  write_canonical_csv(depth_table(fixture), file.path(outdir, "quantification_depth.csv"))
  write_canonical_csv(feature_metadata_table(fixture), file.path(outdir, "feature_metadata.csv"))
  write_canonical_csv(deqms_reference_table(deqms_fit), file.path(outdir, "deqms_B_vs_A.csv"))
  write_canonical_csv(deqms_model_table(deqms_fit), file.path(outdir, "deqms_model.csv"))
  write_canonical_csv(deqms_global_table(deqms_fit), file.path(outdir, "deqms_global.csv"))

  manifest_outdir_label <- args[["manifest-outdir-label"]]
  if (is.na(manifest_outdir_label) || !nzchar(manifest_outdir_label)) {
    manifest_outdir_label <- outdir
  }
  command <- paste(
    "Rscript tests/fixtures/rewrite_parity/differential_deqms_depth/generate_fixture.R",
    "--outdir", manifest_outdir_label,
    "--seed", seed,
    "--timestamp", timestamp,
    "--allow-unpinned-environment", tolower(as.character(allow_unpinned))
  )

  lower_count <- sum(fixture$depth_groups == "lower_depth")
  middle_count <- sum(fixture$depth_groups == "middle_depth")
  higher_count <- sum(fixture$depth_groups == "higher_depth")
  low_pair_row <- which(rownames(fixture$matrix) == fixture$comparable_pair_low)
  high_pair_row <- which(rownames(fixture$matrix) == fixture$comparable_pair_high)

  provenance <- c(
    "# Differential DEqMS Depth Fixture Provenance",
    "",
    paste0("Generated with ", env$r_version),
    paste0("Bioconductor version: ", env$bioconductor_version),
    paste0("limma version: ", env$limma_version),
    paste0("DEqMS version: ", env$deqms_version),
    paste0("Seed: ", seed),
    paste0("Generation timestamp (UTC): ", timestamp),
    paste0("Command: `", command, "`"),
    paste0("Generator SHA-256: ", script_sha256),
    paste0("Byte policy: ", CANONICAL_TEXT_BYTE_POLICY),
    paste0("Serialization policy: ", SERIALIZATION_POLICY),
    "Source policy: deterministic synthetic fixture generated locally without network access; DEqMS::spectraCounteBayes and limma outputs are the external scientific authority for expected numerical columns.",
    "Classification: external parity for DEqMS count-aware empirical-Bayes outputs; feature metadata are fixture sanity metadata.",
    "Expected outputs come only from the pinned R/DEqMS and R/limma run. PhosPy is not imported or executed by this generator.",
    "Redistribution metadata: synthetic deterministic inputs and black-box DEqMS/limma numeric outputs are repository test fixtures; DEqMS and limma source code is not redistributed.",
    "Design: ~0 + condition with groups A/B and unbalanced 5/7 replicates",
    "Contrast: B_vs_A = B - A",
    paste0("Rows: ", nrow(fixture$matrix), " phosphosites/features; columns: ", ncol(fixture$matrix), " samples"),
    "Quantification depth: feature-level synthetic PSM count supplied as fit$count before DEqMS::spectraCounteBayes.",
    "DEqMS model: spectraCounteBayes(fit.method='loess') uses log2(count) and loess span 0.75 to estimate count-dependent prior variance.",
    paste0("Depth groups: lower_depth=", lower_count, "; middle_depth=", middle_count, "; higher_depth=", higher_count),
    paste0("Comparable depth pair: ", fixture$comparable_pair_low, " has depth ", fixture$quantification_depth[[low_pair_row]], "; ", fixture$comparable_pair_high, " has depth ", fixture$quantification_depth[[high_pair_row]], "; their intensity rows are byte-identical before model fitting."),
    "Terminology conclusion: PhosPy's validated feature should be described as quantification-depth-aware empirical Bayes moderation inspired by DEqMS, not exact DEqMS-compatible numerical equivalence, because PhosPy preserves its existing deterministic trend smoother while DEqMS::spectraCounteBayes uses R loess span 0.75.",
    "",
    "Scientific citations:",
    paste0("- ", SCIENTIFIC_CITATION),
    "",
    "Output files and SHA-256 digests are listed in `MANIFEST.json`."
  )
  write_canonical_lines(provenance, file.path(outdir, "PROVENANCE.md"))

  file_rows <- data.frame(
    role = c(
      "generator",
      "input",
      "input",
      "input",
      "input",
      "metadata",
      "output",
      "output",
      "output",
      "governance"
    ),
    relative_path = c(
      "generate_fixture.R",
      "matrix.csv",
      "design.csv",
      "contrasts.csv",
      "quantification_depth.csv",
      "feature_metadata.csv",
      "deqms_B_vs_A.csv",
      "deqms_model.csv",
      "deqms_global.csv",
      "PROVENANCE.md"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  file_rows$sha256 <- unname(tools::sha256sum(file.path(outdir, file_rows$relative_path)))
  manifest_files <- lapply(
    seq_len(nrow(file_rows)),
    function(i) {
      list(
        role = file_rows$role[[i]],
        relative_path = file_rows$relative_path[[i]],
        sha256 = file_rows$sha256[[i]]
      )
    }
  )
  manifest <- list(
    manifest_schema_version = "fixture-manifest-v1",
    fixture_family = FIXTURE_FAMILY,
    classification = "external_parity",
    external_implementation = list(
      name = "R DEqMS",
      r_version = env$r_version,
      bioconductor_version = env$bioconductor_version,
      limma_version = env$limma_version,
      deqms_version = env$deqms_version
    ),
    pinned_environment = list(
      r_version = PINNED_R_VERSION,
      bioconductor_version = PINNED_BIOCONDUCTOR_VERSION,
      limma_version = PINNED_LIMMA_VERSION,
      deqms_version = PINNED_DEQMS_VERSION
    ),
    generator = "tests/fixtures/rewrite_parity/differential_deqms_depth/generate_fixture.R",
    generator_sha256 = script_sha256,
    command = command,
    seed = seed,
    generation_timestamp_utc = timestamp,
    scientific_citation = SCIENTIFIC_CITATION,
    redistribution_metadata = REDISTRIBUTION_METADATA,
    source_policy = "synthetic deterministic fixture generated locally; no network access; DEqMS and limma outputs are the external scientific authority for expected numerical columns",
    numeric_authority = "expected numeric values are serialized directly from the pinned R/DEqMS run and are not generated by PhosPy or manually edited after generation",
    byte_policy = CANONICAL_TEXT_BYTE_POLICY,
    serialization_policy = SERIALIZATION_POLICY,
    model_policy = list(
      upstream_lm_fit = "limma::lmFit",
      upstream_contrast_fit = "limma::contrasts.fit",
      upstream_pre_deqms_ebayes = "limma::eBayes(trend = FALSE, robust = FALSE)",
      deqms_function = "DEqMS::spectraCounteBayes",
      deqms_fit_method = "loess",
      deqms_loess_span = 0.75,
      trend_covariate = "log2(fit$count)",
      quantification_depth_kind = "psm_count"
    ),
    design = list(
      formula = "~0 + condition",
      contrast = "B_vs_A = B - A",
      condition_levels = c("A", "B"),
      sample_ids = colnames(fixture$matrix),
      sample_counts = list(A = 5L, B = 7L),
      n_features = nrow(fixture$matrix),
      lower_depth_feature_count = lower_count,
      middle_depth_feature_count = middle_count,
      higher_depth_feature_count = higher_count,
      comparable_depth_pair = list(
        low_depth_site_id = fixture$comparable_pair_low,
        high_depth_site_id = fixture$comparable_pair_high,
        low_depth = fixture$quantification_depth[[low_pair_row]],
        high_depth = fixture$quantification_depth[[high_pair_row]],
        intensity_rows_identical = TRUE
      )
    ),
    terminology_conclusion = "quantification-depth-aware empirical Bayes moderation inspired by DEqMS; no exact DEqMS-compatible numerical equivalence claim",
    hashes = list(
      input_files_sha256 = aggregate_role_hash(file_rows, c("input")),
      metadata_files_sha256 = aggregate_role_hash(file_rows, c("metadata")),
      output_files_sha256 = aggregate_role_hash(file_rows, c("output")),
      governance_files_sha256 = aggregate_role_hash(file_rows, c("governance", "generator"))
    ),
    files = manifest_files
  )
  write_canonical_json(manifest, file.path(outdir, "MANIFEST.json"))

  message("Done. DEqMS depth fixture written to: ", outdir)
}

main()
