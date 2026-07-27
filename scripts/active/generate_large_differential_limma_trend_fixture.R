#!/usr/bin/env Rscript

required_pkgs <- c("limma")
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

read_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    outdir = "tests/fixtures/rewrite_parity/differential_limma_trend_large",
    seed = "20260724",
    timestamp = "2026-07-24T00:00:00Z",
    n_features = "1600",
    "manifest-outdir-label" = NA_character_
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

CANONICAL_TEXT_BYTE_POLICY <- "utf-8 LF with final newline"

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

write_canonical_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  write.csv(data, con, row.names = FALSE, eol = "\n")
}

json_escape <- function(value) {
  value <- gsub("\\\\", "\\\\\\\\", value)
  value <- gsub("\"", "\\\\\"", value)
  value <- gsub("\n", "\\\\n", value)
  paste0("\"", value, "\"")
}

json_scalar <- function(value) {
  if (is.na(value)) {
    return("null")
  }
  if (is.logical(value)) {
    return(ifelse(value, "true", "false"))
  }
  if (is.numeric(value)) {
    return(format(value, scientific = FALSE, trim = TRUE))
  }
  json_escape(as.character(value))
}

json_array <- function(values) {
  paste0("[", paste(vapply(values, json_scalar, character(1)), collapse = ", "), "]")
}

write_manifest <- function(
  outdir,
  files,
  seed,
  timestamp,
  command,
  script_sha256,
  r_version,
  limma_version,
  n_features,
  samples
) {
  hashes <- tools::sha256sum(file.path(outdir, files))
  file_lines <- c()
  for (i in seq_along(files)) {
    file_lines <- c(
      file_lines,
      paste0(
        "    {\n",
        "      \"relative_path\": ", json_escape(files[[i]]), ",\n",
        "      \"sha256\": ", json_escape(unname(hashes[[i]])), "\n",
        "    }"
      )
    )
  }
  manifest <- c(
    "{",
    "  \"manifest_schema_version\": \"fixture-manifest-v1\",",
    "  \"fixture_family\": \"differential_limma_trend_large\",",
    "  \"classification\": \"external_parity\",",
    "  \"external_implementation\": {",
    "    \"name\": \"R limma\",",
    paste0("    \"r_version\": ", json_escape(r_version), ","),
    paste0("    \"limma_version\": ", json_escape(limma_version)),
    "  },",
    paste0("  \"generator\": ", json_escape("scripts/active/generate_large_differential_limma_trend_fixture.R"), ","),
    paste0("  \"generator_sha256\": ", json_escape(script_sha256), ","),
    paste0("  \"command\": ", json_escape(command), ","),
    paste0("  \"seed\": ", json_scalar(seed), ","),
    paste0("  \"generation_timestamp_utc\": ", json_escape(timestamp), ","),
    "  \"source_policy\": \"synthetic deterministic fixture generated locally; no network access; limma outputs are the external scientific authority for parity columns\",",
    paste0("  \"byte_policy\": ", json_escape(CANONICAL_TEXT_BYTE_POLICY), ","),
    "  \"design\": {",
    "    \"formula\": \"~0 + condition\",",
    "    \"contrast\": \"B_vs_A = B - A\",",
    "    \"condition_levels\": [\"A\", \"B\"],",
    paste0("    \"sample_ids\": ", json_array(samples), ","),
    "    \"sample_counts\": {\"A\": 5, \"B\": 7},",
    paste0("    \"n_features\": ", json_scalar(n_features), ","),
    "    \"mean_variance_trend\": \"feature residual variance increases smoothly with mean intensity plus sinusoidal structure\",",
    "    \"truth_sets\": {\"positive_shift_every\": 16, \"negative_shift_every\": 29, \"negative_shift_offset\": 8}",
    "  },",
    "  \"files\": [",
    paste(file_lines, collapse = ",\n"),
    "  ]",
    "}"
  )
  write_canonical_lines(manifest, file.path(outdir, "MANIFEST.json"))
}

recycle_to_features <- function(value, n_features) {
  if (length(value) == n_features) {
    return(as.numeric(value))
  }
  rep(as.numeric(value[[1]]), n_features)
}

main <- function() {
  options(digits = 17, scipen = 999)
  args <- read_args()
  outdir <- args$outdir
  seed <- as.integer(args$seed)
  n_features <- as.integer(args$n_features)
  timestamp <- args$timestamp
  if (is.na(seed)) {
    stop("--seed must be an integer", call. = FALSE)
  }
  if (is.na(n_features) || n_features < 1500) {
    stop("--n_features must be an integer >= 1500", call. = FALSE)
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)

  a_samples <- sprintf("A_%02d", seq_len(5))
  b_samples <- sprintf("B_%02d", seq_len(7))
  samples <- c(a_samples, b_samples)
  condition <- factor(c(rep("A", length(a_samples)), rep("B", length(b_samples))), levels = c("A", "B"))
  design <- model.matrix(~0 + condition)
  colnames(design) <- c("A", "B")
  rownames(design) <- samples

  feature <- seq_len(n_features)
  site_ids <- sprintf("LG%04d;S%d;", feature, feature)
  mean_intensity <- seq(7.25, 13.75, length.out = n_features) +
    0.15 * sin(feature / 23.0)
  log_sigma <- -1.35 +
    0.12 * (mean_intensity - mean(mean_intensity)) +
    0.22 * sin(mean_intensity * 1.7) +
    0.04 * cos(feature / 19.0)
  sigma <- exp(log_sigma)

  true_shift <- rep(0.0, n_features)
  positive <- (feature %% 16L) == 0L
  negative <- ((feature + 8L) %% 29L) == 0L
  negative <- negative & !positive
  true_shift[positive] <- 0.65
  true_shift[negative] <- -0.50

  mat <- matrix(NA_real_, nrow = n_features, ncol = length(samples))
  for (j in seq_along(samples)) {
    shift <- if (condition[[j]] == "B") true_shift else rep(0.0, n_features)
    mat[, j] <- rnorm(n_features, mean = mean_intensity + shift, sd = sigma)
  }
  rownames(mat) <- site_ids
  colnames(mat) <- samples

  contrast_mat <- makeContrasts(B_vs_A = B - A, levels = design)
  fit <- lmFit(mat, design)
  fit2 <- contrasts.fit(fit, contrast_mat)
  fit2 <- eBayes(fit2, trend = TRUE, robust = FALSE)

  coef_name <- "B_vs_A"
  p_values <- as.numeric(fit2$p.value[, coef_name])
  expected <- data.frame(
    site_id = rownames(mat),
    logFC = as.numeric(fit2$coefficients[, coef_name]),
    AveExpr = as.numeric(fit2$Amean),
    SE = as.numeric(fit2$stdev.unscaled[, coef_name] * sqrt(fit2$s2.post)),
    t = as.numeric(fit2$t[, coef_name]),
    P.Value = p_values,
    adj.P.Val = p.adjust(p_values, method = "BH"),
    s2.prior = recycle_to_features(fit2$s2.prior, n_features),
    df.prior = recycle_to_features(fit2$df.prior, n_features),
    s2.post = as.numeric(fit2$s2.post),
    df.total = recycle_to_features(fit2$df.total, n_features),
    sigma = as.numeric(fit2$sigma),
    stdev.unscaled = as.numeric(fit2$stdev.unscaled[, coef_name]),
    check.names = FALSE
  )

  diagnostics <- data.frame(
    site_id = rownames(mat),
    mean_intensity = mean_intensity,
    expected_sigma = sigma,
    true_shift = true_shift,
    is_shifted = true_shift != 0.0,
    shift_direction = ifelse(true_shift > 0.0, "positive", ifelse(true_shift < 0.0, "negative", "none")),
    check.names = FALSE
  )

  write_canonical_csv(
    data.frame(site_id = rownames(mat), mat, check.names = FALSE),
    file.path(outdir, "matrix.csv")
  )
  write_canonical_csv(
    data.frame(sample = rownames(design), design, check.names = FALSE),
    file.path(outdir, "design.csv")
  )
  write_canonical_csv(
    data.frame(coefficient = rownames(contrast_mat), contrast_mat, check.names = FALSE),
    file.path(outdir, "contrasts.csv")
  )
  write_canonical_csv(expected, file.path(outdir, "limma_B_vs_A.csv"))
  write_canonical_csv(diagnostics, file.path(outdir, "simulation_diagnostics.csv"))

  out_files <- c(
    "matrix.csv",
    "design.csv",
    "contrasts.csv",
    "limma_B_vs_A.csv",
    "simulation_diagnostics.csv",
    "PROVENANCE.md"
  )

  r_version <- R.version$version.string
  limma_version <- as.character(packageVersion("limma"))
  script <- script_path()
  script_sha256 <- ifelse(is.na(script), NA_character_, unname(tools::sha256sum(script)))
  manifest_outdir_label <- args[["manifest-outdir-label"]]
  if (is.na(manifest_outdir_label) || !nzchar(manifest_outdir_label)) {
    manifest_outdir_label <- outdir
  }
  command <- paste(
    "Rscript scripts/active/generate_large_differential_limma_trend_fixture.R",
    "--outdir", manifest_outdir_label,
    "--seed", seed,
    "--timestamp", timestamp,
    "--n_features", n_features
  )

  provenance <- c(
    "# Large Differential Limma Trend Fixture Provenance",
    "",
    paste0("Generated with ", r_version),
    paste0("limma version: ", limma_version),
    paste0("Seed: ", seed),
    paste0("Generation timestamp (UTC): ", timestamp),
    paste0("Command: `", command, "`"),
    "Source policy: deterministic synthetic fixture generated locally without network access; limma is the external authority for exported parity quantities.",
    "Classification: external parity for limma result columns; simulation diagnostics are fixture sanity metadata.",
    "Design: ~0 + condition with groups A/B and unbalanced 5/7 replicates",
    "Contrast: B_vs_A = B - A",
    paste0("Rows: ", n_features, " phosphosites/features; columns: ", length(samples), " samples"),
    "Mean-variance trend: expected residual variance increases smoothly with mean intensity plus deterministic sinusoidal structure.",
    paste0("Shifted features: ", sum(true_shift != 0.0), " total; positive=", sum(true_shift > 0.0), "; negative=", sum(true_shift < 0.0)),
    paste0("Generator SHA-256: ", script_sha256),
    "",
    "Output files are listed with SHA-256 digests in `MANIFEST.json`."
  )
  write_canonical_lines(provenance, file.path(outdir, "PROVENANCE.md"))

  write_manifest(
    outdir = outdir,
    files = out_files,
    seed = seed,
    timestamp = timestamp,
    command = command,
    script_sha256 = script_sha256,
    r_version = r_version,
    limma_version = limma_version,
    n_features = n_features,
    samples = samples
  )

  message("Done. Large limma trend fixture written to: ", outdir)
}

main()
