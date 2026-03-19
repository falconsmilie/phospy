#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(janitor)
})

required_pkgs <- c("PhosR")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, then rerun this script.",
    call. = FALSE
  )
}

cols_total <- paste0("group", 1:6)
cols_phospho <- paste0("p_group", 1:6)
cols_corrected <- paste0("phospho_corrected_", 1:6)

read_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    total = "examples/data/total.tsv",
    phospho = "examples/data/phospho.tsv",
    outdir = "tests/fixtures/r_reference"
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
    parsed[[key]] <- value
    i <- i + 2
  }
  parsed
}

replace_sentinel <- function(df, columns, sentinel) {
  df %>%
    mutate(across(all_of(columns), ~ ifelse(. == sentinel, NA_real_, as.numeric(.))))
}

collapse_duplicate_genes <- function(df, gene_col, value_cols) {
  df %>%
    group_by(.data[[gene_col]]) %>%
    mutate(.mean_signal = rowMeans(across(all_of(value_cols)), na.rm = TRUE)) %>%
    slice_max(order_by = .mean_signal, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-.mean_signal)
}

filter_min_observed <- function(df, columns, min_observed = 4L) {
  df %>%
    filter(rowSums(!is.na(select(., all_of(columns)))) >= min_observed)
}

build_weighted_kinase_activity <- function(pred_mat, phospho_matrix, top_n = 20L, min_substrates = 3L) {
  kinase_names <- colnames(pred_mat)
  sample_names <- colnames(phospho_matrix)
  kinase_mat <- matrix(NA_real_, nrow = length(kinase_names), ncol = length(sample_names))
  rownames(kinase_mat) <- kinase_names
  colnames(kinase_mat) <- sample_names

  for (kinase in kinase_names) {
    ordered_sites <- names(sort(pred_mat[, kinase], decreasing = TRUE))[seq_len(min(top_n, nrow(pred_mat)))]
    ordered_sites <- ordered_sites[ordered_sites %in% rownames(phospho_matrix)]

    if (length(ordered_sites) < min_substrates) {
      next
    }

    weights <- pred_mat[ordered_sites, kinase]
    kinase_mat[kinase, ] <- apply(
      phospho_matrix[ordered_sites, , drop = FALSE],
      2,
      function(x) weighted.mean(x, weights, na.rm = TRUE)
    )
  }

  kinase_mat[apply(is.na(kinase_mat), 1, all) == FALSE, , drop = FALSE]
}

build_ksea_scores <- function(pred_mat, phospho_matrix, threshold = 0.6, min_substrates = 3L) {
  kinase_names <- colnames(pred_mat)
  keep <- vapply(
    kinase_names,
    function(kinase) sum(pred_mat[, kinase] > threshold, na.rm = TRUE) >= min_substrates,
    logical(1)
  )
  kinase_names <- kinase_names[keep]

  if (length(kinase_names) == 0) {
    return(list(scores = matrix(numeric(0), nrow = 0, ncol = ncol(phospho_matrix)), counts = integer(0)))
  }

  ksea_scores <- sapply(kinase_names, function(kinase) {
    sites <- rownames(pred_mat)[pred_mat[, kinase] > threshold]
    sites <- intersect(sites, rownames(phospho_matrix))
    if (length(sites) < min_substrates) {
      return(rep(NA_real_, ncol(phospho_matrix)))
    }
    colMeans(phospho_matrix[sites, , drop = FALSE], na.rm = TRUE)
  })

  if (is.null(dim(ksea_scores))) {
    ksea_scores <- matrix(ksea_scores, nrow = 1)
    rownames(ksea_scores) <- colnames(pred_mat)[keep][1]
    colnames(ksea_scores) <- colnames(phospho_matrix)
  } else {
    ksea_scores <- t(ksea_scores)
    colnames(ksea_scores) <- colnames(phospho_matrix)
  }

  ksea_counts <- vapply(kinase_names, function(kinase) {
    sum(rownames(pred_mat)[pred_mat[, kinase] > threshold] %in% rownames(phospho_matrix))
  }, integer(1))

  list(scores = ksea_scores, counts = sort(ksea_counts, decreasing = TRUE))
}

write_session_info <- function(outdir) {
  sink(file.path(outdir, "sessionInfo.txt"))
  print(sessionInfo())
  sink()
}

load_phosphosite_mouse <- function() {
  env <- new.env(parent = emptyenv())
  tryCatch({
    data("PhosphoSite.mouse", package = "PhosR", envir = env)
    if (exists("PhosphoSite.mouse", envir = env, inherits = FALSE)) {
      return(get("PhosphoSite.mouse", envir = env, inherits = FALSE))
    }
    NULL
  }, error = function(e) NULL)
}

main <- function() {
  args <- read_args()
  outdir <- args$outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  total_path <- args$total
  phospho_path <- args$phospho

  message("Reading input tables...")
  df_total <- read_delim(total_path, delim = "\t", show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(
      genes = as.character(genes),
      across(all_of(cols_total), as.numeric)
    ) %>%
    replace_sentinel(cols_total, sentinel = 10)

  df_total_unique <- collapse_duplicate_genes(df_total, gene_col = "genes", value_cols = cols_total) %>%
    mutate(genes = toupper(genes))

  df_total_filtered <- filter_min_observed(df_total_unique, cols_total, min_observed = 4)

  df_phospho <- read_delim(
    phospho_path,
    delim = "\t",
    locale = locale(encoding = "UTF-16LE"),
    show_col_types = FALSE
  ) %>%
    clean_names() %>%
    mutate(
      gene_names = toupper(as.character(gene_names)),
      gene_p_site = as.character(gene_p_site),
      localization_prob = as.numeric(localization_prob),
      centralized_sequence = as.character(centralized_sequence),
      across(all_of(cols_phospho), as.numeric)
    ) %>%
    filter(!is.na(uid), !is.na(gene_names)) %>%
    filter(localization_prob >= 0.75) %>%
    replace_sentinel(cols_phospho, sentinel = 12)

  df_phospho_filtered <- filter_min_observed(df_phospho, cols_phospho, min_observed = 4)

  message("Building corrected phosphosite table...")
  df_phospho_corrected <- df_phospho_filtered %>%
    inner_join(
      df_total_filtered %>% select(genes, all_of(cols_total)),
      by = c("gene_names" = "genes")
    )

  for (i in seq_along(cols_total)) {
    df_phospho_corrected[[cols_corrected[[i]]]] <-
      df_phospho_corrected[[cols_phospho[[i]]]] - df_phospho_corrected[[cols_total[[i]]]]
  }

  message("Constructing PhosR input matrix...")
  phosr_input <- df_phospho_corrected %>%
    separate(gene_p_site, into = c("gene", "p_site"), sep = "_", remove = FALSE) %>%
    mutate(
      gene = as.character(gene),
      p_site = as.character(p_site),
      site_id = paste0(toupper(gene), ";", toupper(p_site), ";")
    ) %>%
    select(gene_names, gene, p_site, uid, centralized_sequence, site_id, all_of(cols_corrected)) %>%
    filter(!is.na(centralized_sequence)) %>%
    filter(if_all(all_of(cols_corrected), ~ !is.na(.))) %>%
    group_by(site_id) %>%
    mutate(.mean_signal = rowMeans(across(all_of(cols_corrected)), na.rm = TRUE)) %>%
    slice_max(order_by = .mean_signal, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-.mean_signal)

  mat_phospho_corrected <- phosr_input %>%
    select(site_id, all_of(cols_corrected)) %>%
    column_to_rownames("site_id") %>%
    as.matrix()

  seqs <- phosr_input$centralized_sequence
  names(seqs) <- phosr_input$site_id

  write.csv(df_total_unique, file.path(outdir, "df_total_unique.csv"), row.names = FALSE)
  write.csv(df_total_filtered, file.path(outdir, "df_total_filtered.csv"), row.names = FALSE)
  write.csv(df_phospho_filtered, file.path(outdir, "df_phospho_filtered.csv"), row.names = FALSE)
  write.csv(df_phospho_corrected, file.path(outdir, "df_phospho_corrected.csv"), row.names = FALSE)
  write.csv(phosr_input, file.path(outdir, "phosr_input.csv"), row.names = FALSE)
  write.csv(mat_phospho_corrected, file.path(outdir, "mat_phospho_corrected.csv"), row.names = TRUE)
  write.csv(
    data.frame(site_id = names(seqs), centralized_sequence = unname(seqs)),
    file.path(outdir, "site_sequences.csv"),
    row.names = FALSE
  )

  substrate_list <- load_phosphosite_mouse()
  if (is.null(substrate_list)) {
    stop(
      "Could not load PhosphoSite.mouse from PhosR. The package is installed, but the dataset was not found.",
      call. = FALSE
    )
  }

  message("Running PhosR kinase-substrate scoring...")
  set.seed(1)
  kss_mat <- PhosR::kinaseSubstrateScore(
    substrate.list = substrate_list,
    mat = mat_phospho_corrected,
    seqs = seqs,
    numMotif = 5,
    numSub = 1,
    species = "mouse",
    verbose = TRUE
  )

  pred_mat <- PhosR::kinaseSubstratePred(
    phosScoringMatrices = kss_mat,
    ensembleSize = 5,
    top = 10,
    cs = 0.6,
    inclusion = 5,
    iter = 3
  )

  kinase_activity <- build_weighted_kinase_activity(pred_mat, mat_phospho_corrected)
  ksea <- build_ksea_scores(pred_mat, mat_phospho_corrected)

  write.csv(pred_mat, file.path(outdir, "predMat.csv"), row.names = TRUE)
  write.csv(kinase_activity, file.path(outdir, "kinase_activity_matrix.csv"), row.names = TRUE)
  write.csv(ksea$scores, file.path(outdir, "ksea_scores.csv"), row.names = TRUE)
  write.csv(
    data.frame(kinase = names(ksea$counts), n_substrates = unname(ksea$counts)),
    file.path(outdir, "ksea_counts.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(
      kinase = colnames(pred_mat),
      n_targets = vapply(colnames(pred_mat), function(k) sum(pred_mat[, k] > 0.6, na.rm = TRUE), integer(1))
    ),
    file.path(outdir, "kinase_target_counts.csv"),
    row.names = FALSE
  )

  write_session_info(outdir)
  message("Done. R reference fixtures written to: ", outdir)
}

main()
