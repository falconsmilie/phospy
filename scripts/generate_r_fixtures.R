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

# Relaxed thresholds for the tiny synthetic fixture dataset so that
# PhosR scoring produces non-empty outputs more reliably.
phosr_num_motif <- 1L
phosr_num_sub <- 1L

comparisons <- list(
  c("group1", "group4"),
  c("group2", "group5"),
  c("group3", "group6"),
  c("group1", "group2"),
  c("group1", "group3"),
  c("group2", "group3"),
  c("group4", "group5"),
  c("group4", "group6"),
  c("group5", "group6")
)

read_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    total = "examples/data/total.tsv",
    phospho = "examples/data/phospho.tsv",
    outdir = "tests/fixtures/rewrite_parity/r_reference"
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

add_pairwise_comparisons <- function(df, comparisons, corrected_cols, output_prefix = "p_") {
  result <- df
  group_map <- stats::setNames(corrected_cols, paste0("group", seq_along(corrected_cols)))

  for (comparison in comparisons) {
    left <- comparison[[1]]
    right <- comparison[[2]]
    result[[paste0(output_prefix, left, "_", right)]] <- result[[group_map[[left]]]] - result[[group_map[[right]]]]
  }

  result
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
    empty_scores <- matrix(numeric(0), nrow = 0, ncol = ncol(phospho_matrix))
    colnames(empty_scores) <- colnames(phospho_matrix)
    return(list(scores = empty_scores, counts = integer(0)))
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
    rownames(ksea_scores) <- kinase_names[1]
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
    data("PhosphoSitePlus", package = "PhosR", envir = env)
    if (exists("PhosphoSite.mouse", envir = env, inherits = FALSE)) {
      return(get("PhosphoSite.mouse", envir = env, inherits = FALSE))
    }
    NULL
  }, error = function(e) NULL)
}

load_motif_mouse_list <- function() {
  env <- new.env(parent = emptyenv())
  tryCatch({
    data("KinaseMotifs", package = "PhosR", envir = env)
    if (exists("motif.mouse.list", envir = env, inherits = FALSE)) {
      return(get("motif.mouse.list", envir = env, inherits = FALSE))
    }
    NULL
  }, error = function(e) NULL)
}

score_phosphosites_motifs_compat <- function(mat, motif.mouse.list.filtered, seqs) {
  motifScoreMatrix <- matrix(
    NA_real_,
    nrow = nrow(mat),
    ncol = length(motif.mouse.list.filtered)
  )
  rownames(motifScoreMatrix) <- rownames(mat)
  colnames(motifScoreMatrix) <- names(motif.mouse.list.filtered)

  seqWin <- vapply(seqs, function(x) {
    mid <- (nchar(x) + 1) / 2
    substr(x, start = (mid - 7), stop = (mid + 7))
  }, character(1))

  message("Scoring phosphosites against kinase motifs:")
  for (i in seq_len(length(motif.mouse.list.filtered))) {
    motifScoreMatrix[, i] <- PhosR::frequencyScoring(
      seqWin,
      motif.mouse.list.filtered[[i]]
    )
    cat(paste0(i, ".\n"))
  }

  message("done.")
  PhosR::minmax(motifScoreMatrix)
}

score_phosphosite_profile_compat <- function(mat, ks.profile.list.filtered) {
  profileScoreMatrix <- (
    t(
      apply(
        mat,
        1,
        cor,
        t(do.call(rbind, ks.profile.list.filtered))
      )
    ) + 1
  ) / 2

  if (is.null(dim(profileScoreMatrix))) {
    profileScoreMatrix <- matrix(profileScoreMatrix, ncol = 1)
  }

  rownames(profileScoreMatrix) <- rownames(mat)
  colnames(profileScoreMatrix) <- names(ks.profile.list.filtered)
  profileScoreMatrix
}

kinase_substrate_score_compat <- function(
    substrate.list,
    mat,
    seqs,
    motif.mouse.list,
    numMotif = 5,
    numSub = 1
) {
  ks.profile.list <- PhosR::kinaseSubstrateProfile(substrate.list, mat)

  message("Number of kinases passed motif size filtering: ",
          sum(motif.mouse.list$NumInputSeq >= numMotif))
  message("Number of kinases passed profile size filtering: ",
          sum(ks.profile.list$NumSub >= numSub))

  motif.mouse.list.filtered <- motif.mouse.list[
    which(motif.mouse.list$NumInputSeq >= numMotif)
  ]
  ks.profile.list.filtered <- ks.profile.list[
    which(ks.profile.list$NumSub >= numSub)
  ]

  motifScoreMatrix <- score_phosphosites_motifs_compat(
    mat = mat,
    motif.mouse.list.filtered = motif.mouse.list.filtered,
    seqs = seqs
  )

  message("Scoring phosphosites against kinase-substrate profiles:")
  profileScoreMatrix <- score_phosphosite_profile_compat(
    mat,
    ks.profile.list.filtered
  )
  message("done.")

  o <- intersect(colnames(motifScoreMatrix), colnames(profileScoreMatrix))

  if (length(o) == 0) {
    warning(
      "No overlapping kinases between motif and profile score matrices. ",
      "Falling back to profile-only combined scores for fixture generation."
    )

    ksActivityMatrix <- do.call(rbind, ks.profile.list.filtered)
    rownames(ksActivityMatrix) <- names(ks.profile.list.filtered)
    ksActivityMatrix <- ksActivityMatrix[colnames(profileScoreMatrix), , drop = FALSE]

    return(list(
      motifScoreMatrix = motifScoreMatrix,
      profileScoreMatrix = profileScoreMatrix,
      combinedScoreMatrix = profileScoreMatrix,
      ksActivityMatrix = ksActivityMatrix,
      weights = setNames(rep(1, ncol(profileScoreMatrix)), colnames(profileScoreMatrix))
    ))
  }

  message("Generating combined scores for phosphosites")
  message("by motifs and phospho profiles:")

  combinedScoreMatrix <- matrix(
    NA_real_,
    nrow = nrow(motifScoreMatrix),
    ncol = length(o)
  )
  colnames(combinedScoreMatrix) <- o
  rownames(combinedScoreMatrix) <- rownames(motifScoreMatrix)

  w1 <- log(rank(motif.mouse.list$NumInputSeq[o]) + 1)
  w2 <- log(rank(ks.profile.list$NumSub[o]) + 1)
  w3 <- w1 + w2

  for (i in seq_len(length(o))) {
    combinedScoreMatrix[, i] <-
      (w1[i] / (w1[i] + w2[i]) * motifScoreMatrix[, o[i]]) +
      (w2[i] / (w1[i] + w2[i]) * profileScoreMatrix[, o[i]])
  }
  message("done.")

  ksActivityMatrix <- do.call(rbind, ks.profile.list.filtered)
  rownames(ksActivityMatrix) <- names(ks.profile.list.filtered)
  ksActivityMatrix <- ksActivityMatrix[o, , drop = FALSE]

  list(
    motifScoreMatrix = motifScoreMatrix,
    profileScoreMatrix = profileScoreMatrix,
    combinedScoreMatrix = combinedScoreMatrix,
    ksActivityMatrix = ksActivityMatrix,
    weights = w3
  )
}

can_run_kinase_pred_compat <- function(phosScoringMatrices, top = 10L, cs = 0.6, inclusion = 5L) {
  featureMat <- phosScoringMatrices$combinedScoreMatrix

  if (is.null(featureMat) || ncol(featureMat) == 0 || nrow(featureMat) == 0) {
    return(FALSE)
  }

  valid <- vapply(seq_len(ncol(featureMat)), function(i) {
    scores <- sort(featureMat[, i], decreasing = TRUE, na.last = NA)
    scores <- head(scores, min(top, length(scores)))
    sum(scores > cs, na.rm = TRUE) >= inclusion
  }, logical(1))

  any(valid)
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

  df_phospho_corrected <- add_pairwise_comparisons(
    df_phospho_corrected,
    comparisons = comparisons,
    corrected_cols = cols_corrected,
    output_prefix = "p_"
  )

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
      "Could not load PhosphoSite.mouse from the PhosphoSitePlus dataset in PhosR.",
      call. = FALSE
    )
  }

  motif.mouse.list <- load_motif_mouse_list()
  if (is.null(motif.mouse.list)) {
    stop(
      "Could not load motif.mouse.list from the KinaseMotifs dataset in PhosR.",
      call. = FALSE
    )
  }

  message("Running PhosR kinase-substrate scoring...")
  set.seed(1)
  kss_mat <- kinase_substrate_score_compat(
    substrate.list = substrate_list,
    mat = mat_phospho_corrected,
    seqs = seqs,
    motif.mouse.list = motif.mouse.list,
    numMotif = phosr_num_motif,
    numSub = phosr_num_sub
  )

  pred_top <- 10L
  pred_cs <- 0.6
  pred_inclusion <- 5L
  pred_iter <- 3L
  pred_ensemble <- 5L

  message("Predicting kinases for phosphosites:")

  if (!can_run_kinase_pred_compat(
    phosScoringMatrices = kss_mat,
    top = pred_top,
    cs = pred_cs,
    inclusion = pred_inclusion
  )) {
    warning(
      "Synthetic fixture dataset is too small for stable kinaseSubstratePred(); ",
      "using combinedScoreMatrix as predMat fallback."
    )
    pred_mat <- as.matrix(kss_mat$combinedScoreMatrix)
  } else {
    pred_mat <- tryCatch(
      PhosR::kinaseSubstratePred(
        phosScoringMatrices = kss_mat,
        ensembleSize = pred_ensemble,
        top = pred_top,
        cs = pred_cs,
        inclusion = pred_inclusion,
        iter = pred_iter
      ),
      error = function(e) {
        warning(
          "kinaseSubstratePred() failed on the synthetic fixture dataset: ",
          conditionMessage(e),
          ". Falling back to combinedScoreMatrix as predMat."
        )
        as.matrix(kss_mat$combinedScoreMatrix)
      }
    )
  }

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
