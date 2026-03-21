#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

required_pkgs <- c("PhosR", "SummarizedExperiment")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, then rerun this script.",
    call. = FALSE
  )
}

read_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    outdir = "tests/fixtures/r_reference_l6"
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

write_session_info <- function(outdir) {
  sink(file.path(outdir, "sessionInfo.txt"))
  print(sessionInfo())
  sink()
}

load_required_objects <- function() {
  env <- new.env(parent = emptyenv())
  data("phospho_L6_ratio_pe", package = "PhosR", envir = env)
  data("SPSs", package = "PhosR", envir = env)
  data("PhosphoSitePlus", package = "PhosR", envir = env)
  data("KinaseMotifs", package = "PhosR", envir = env)

  required <- c("phospho.L6.ratio.pe", "SPSs", "PhosphoSite.mouse", "motif.mouse.list")
  missing <- required[!vapply(required, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(missing) > 0) {
    stop("Missing required PhosR data objects: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  list(
    ppe = get("phospho.L6.ratio.pe", envir = env, inherits = FALSE),
    SPSs = get("SPSs", envir = env, inherits = FALSE),
    PhosphoSite.mouse = get("PhosphoSite.mouse", envir = env, inherits = FALSE),
    motif.mouse.list = get("motif.mouse.list", envir = env, inherits = FALSE)
  )
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
    stop("No overlapping kinases between motif and profile score matrices.", call. = FALSE)
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

main <- function() {
  args <- read_args()
  outdir <- args$outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message("Loading bundled PhosR L6 dataset...")
  objs <- load_required_objects()
  ppe <- objs$ppe
  SPSs <- objs$SPSs
  substrate_list <- objs$PhosphoSite.mouse
  motif.mouse.list <- objs$motif.mouse.list

  sites <- paste(
    sapply(PhosR::GeneSymbol(ppe), function(x) x),
    ";",
    sapply(PhosR::Residue(ppe), function(x) x),
    sapply(PhosR::Site(ppe), function(x) x),
    ";",
    sep = ""
  )
  grps <- gsub("_.+", "", colnames(ppe))
  design <- model.matrix(~ grps - 1)
  ctl <- which(sites %in% SPSs)

  message("Running PhosR normalization and filtering on L6 data...")
  ppe <- PhosR::RUVphospho(ppe, M = design, k = 3, ctl = ctl)
  phosphoL6 <- SummarizedExperiment::assay(ppe, "normalised")
  phosphoL6.mean <- PhosR::meanAbundance(phosphoL6, grps = grps)
  aov <- PhosR::matANOVA(mat = phosphoL6, grps = grps)
  idx <- (aov < 0.05) & (rowSums(phosphoL6.mean > 0.5) > 0)
  phosphoL6.reg <- phosphoL6[idx, , drop = FALSE]
  L6.phos.std <- PhosR::standardise(phosphoL6.reg)
  rownames(L6.phos.std) <- paste0(
    PhosR::GeneSymbol(ppe), ";",
    PhosR::Residue(ppe),
    PhosR::Site(ppe), ";"
  )[idx]
  L6.phos.seq <- PhosR::Sequence(ppe)[idx]

  message("Running PhosR kinase-substrate scoring on L6 data...")
  L6.matrices <- kinase_substrate_score_compat(
    substrate.list = substrate_list,
    mat = L6.phos.std,
    seqs = L6.phos.seq,
    motif.mouse.list = motif.mouse.list,
    numMotif = 5,
    numSub = 1
  )

  set.seed(1)
  L6.predMat <- PhosR::kinaseSubstratePred(L6.matrices, top = 30)

  kinase_activity <- build_weighted_kinase_activity(L6.predMat, L6.phos.std)
  ksea <- build_ksea_scores(L6.predMat, L6.phos.std)

  write.csv(L6.phos.std, file.path(outdir, "l6_phospho_matrix.csv"), row.names = TRUE)
  write.csv(
    data.frame(site_id = names(L6.phos.seq), centralized_sequence = unname(L6.phos.seq)),
    file.path(outdir, "l6_site_sequences.csv"),
    row.names = FALSE
  )
  write.csv(L6.predMat, file.path(outdir, "predMat.csv"), row.names = TRUE)
  write.csv(kinase_activity, file.path(outdir, "kinase_activity_matrix.csv"), row.names = TRUE)
  write.csv(ksea$scores, file.path(outdir, "ksea_scores.csv"), row.names = TRUE)
  write.csv(
    data.frame(kinase = names(ksea$counts), n_substrates = unname(ksea$counts)),
    file.path(outdir, "ksea_counts.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(
      kinase = colnames(L6.predMat),
      n_targets = vapply(colnames(L6.predMat), function(k) sum(L6.predMat[, k] > 0.6, na.rm = TRUE), integer(1))
    ),
    file.path(outdir, "kinase_target_counts.csv"),
    row.names = FALSE
  )

  write_session_info(outdir)
  message("Done. R L6 reference fixtures written to: ", outdir)
}

main()
