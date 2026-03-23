#!/usr/bin/env Rscript

required_pkgs <- c("PhosR", "SummarizedExperiment", "e1071")
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

multiAdaSampling_fn <- getFromNamespace("multiAdaSampling", "PhosR")

read_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    outdir = "tests/fixtures/r_reference_l6",
    trace_kinases = "PRKAA1,MAPK1,MAPK9,IRAK1,TBK1,LCK",
    trace_top_n = "10"
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

parse_csv_values <- function(value) {
  parts <- trimws(unlist(strsplit(value, ",", fixed = TRUE)))
  parts[nzchar(parts)]
}

write_session_info <- function(outdir) {
  sink(file.path(outdir, "sessionInfo.txt"))
  print(sessionInfo())
  sink()
}

write_trace_readme <- function(trace_dir, trace_kinases, trace_top_n) {
  lines <- c(
    "# R prediction trace fixtures",
    "",
    "These files are generated from the bundled PhosR L6 example path.",
    "They are intended for direct comparison with Python-side prediction debug traces.",
    "",
    paste0("Trace kinases: ", paste(trace_kinases, collapse = ", ")),
    paste0("Per-ensemble top-N export: ", trace_top_n),
    "",
    "Files:",
    "- trace_candidates.csv: ranked combined-score candidates for the traced kinases",
    "- trace_initial_negatives.csv: initial negative draw for each ensemble member",
    "- trace_iteration_probabilities.csv: per-iteration class probabilities on the base train set",
    "- trace_iteration_resampling_weights.csv: per-iteration class-specific resampling weights",
    "- trace_iteration_samples.csv: resampled site identities for each iteration and class",
    "- trace_final_ensemble_predictions.csv: final per-ensemble prediction probabilities for all sites",
    "- trace_final_ensemble_top.csv: final per-ensemble top-ranked sites"
  )
  writeLines(lines, file.path(trace_dir, "README.md"))
}



build_quantified_substrate_map <- function(
  substrate.list,
  observed_sites,
  min_substrates = 1L
) {
  substrate_map <- list()
  map_idx <- 0L

  for (kinase in names(substrate.list)) {
    quantified <- c()
    seen <- character(0)

    for (site in substrate.list[[kinase]]) {
      if ((site %in% observed_sites) && !(site %in% seen)) {
        quantified <- c(quantified, site)
        seen <- c(seen, site)
      }
    }

    if (length(quantified) >= min_substrates) {
      map_idx <- map_idx + 1L
      substrate_map[[map_idx]] <- quantified
      names(substrate_map)[map_idx] <- kinase
    }
  }

  substrate_map
}

flatten_grouped_mapping <- function(
  grouped_map,
  group_col = "kinase",
  value_col = "site_id"
) {
  rows <- list()
  row_idx <- 0L

  for (group_name in names(grouped_map)) {
    values <- grouped_map[[group_name]]
    if (length(values) == 0) {
      next
    }
    row_idx <- row_idx + 1L
    rows[[row_idx]] <- data.frame(
      kinase = rep(group_name, length(values)),
      site_id = values,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) {
    empty <- data.frame(
      kinase = character(0),
      site_id = character(0),
      stringsAsFactors = FALSE
    )
    names(empty) <- c(group_col, value_col)
    return(empty)
  }

  result <- do.call(rbind, rows)
  names(result) <- c(group_col, value_col)
  result
}

build_combined_weight_table <- function(
  overlap_kinases,
  motif_sizes,
  profile_sizes
) {
  motif_rank_weight <- log(rank(motif_sizes[overlap_kinases]) + 1)
  profile_rank_weight <- log(rank(profile_sizes[overlap_kinases]) + 1)
  total_weight <- motif_rank_weight + profile_rank_weight

  data.frame(
    kinase = overlap_kinases,
    motif_weight = as.numeric(motif_rank_weight / total_weight),
    profile_weight = as.numeric(profile_rank_weight / total_weight),
    motif_rank_weight = as.numeric(motif_rank_weight),
    profile_rank_weight = as.numeric(profile_rank_weight),
    stringsAsFactors = FALSE
  )
}

build_prediction_top_table <- function(pred_mat, top_n = 30L) {
  rows <- list()
  row_idx <- 0L

  for (kinase in colnames(pred_mat)) {
    ordered <- sort(pred_mat[, kinase], decreasing = TRUE)
    top_count <- min(top_n, length(ordered))
    if (top_count == 0) {
      next
    }
    top_sites <- names(ordered)[seq_len(top_count)]
    row_idx <- row_idx + 1L
    rows[[row_idx]] <- data.frame(
      kinase = kinase,
      site_id = top_sites,
      pred_score = as.numeric(ordered[seq_len(top_count)]),
      rank = seq_len(top_count),
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) {
    return(data.frame(
      kinase = character(0),
      site_id = character(0),
      pred_score = numeric(0),
      rank = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, rows)
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
    stop(
      "Missing required PhosR data objects: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
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

  message(
    "Number of kinases passed motif size filtering: ",
    sum(motif.mouse.list$NumInputSeq >= numMotif)
  )
  message(
    "Number of kinases passed profile size filtering: ",
    sum(ks.profile.list$NumSub >= numSub)
  )

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

build_weighted_kinase_activity <- function(
  pred_mat,
  phospho_matrix,
  top_n = 20L,
  min_substrates = 3L
) {
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

trace_substrate_list <- function(phosScoringMatrices, top, cs, inclusion, trace_kinases) {
  combined <- phosScoringMatrices$combinedScoreMatrix
  substrate.list <- list()
  kinaseSel <- c()
  candidate_rows <- list()
  candidate_count <- 0
  substrate_count <- 0

  for (i in seq_len(ncol(combined))) {
    kinase <- colnames(combined)[i]
    ordered <- sort(combined[, i], decreasing = TRUE)
    top_count <- min(top, length(ordered))
    top_sites <- ordered[seq_len(top_count)]
    selected_sites <- names(which(top_sites > cs))

    if (length(selected_sites) >= inclusion) {
      substrate_count <- substrate_count + 1
      substrate.list[[substrate_count]] <- selected_sites
      kinaseSel <- c(kinaseSel, kinase)
    }

    if (kinase %in% trace_kinases) {
      candidate_count <- candidate_count + 1
      candidate_rows[[candidate_count]] <- data.frame(
        kinase = kinase,
        rank = seq_along(ordered),
        site = names(ordered),
        combined_score = as.numeric(ordered),
        within_top = seq_along(ordered) <= top,
        above_threshold = as.numeric(ordered) > cs,
        selected_candidate = names(ordered) %in% selected_sites,
        stringsAsFactors = FALSE
      )
    }
  }

  names(substrate.list) <- kinaseSel
  candidates <- if (length(candidate_rows) > 0) {
    do.call(rbind, candidate_rows)
  } else {
    data.frame(
      kinase = character(0),
      rank = integer(0),
      site = character(0),
      combined_score = numeric(0),
      within_top = logical(0),
      above_threshold = logical(0),
      selected_candidate = logical(0),
      stringsAsFactors = FALSE
    )
  }

  list(substrate_list = substrate.list, candidates = candidates)
}

trace_multi_ada_sampling <- function(
  train.mat,
  test.mat,
  label,
  kernelType,
  iter = 5,
  kinase,
  ensemble_idx,
  trace_top_n = 10
) {
  X <- train.mat
  Y <- label

  model <- c()
  prob.mat <- c()
  iteration_prob_rows <- list()
  iteration_decision_rows <- list()
  iteration_weight_rows <- list()
  iteration_sample_rows <- list()
  probability_count <- 0
  decision_count <- 0
  weight_count <- 0
  sample_count <- 0

  for (i in seq_len(iter)) {
    tmp <- X
    rownames(tmp) <- NULL
    model <- e1071::svm(tmp, factor(Y), kernel = kernelType, probability = TRUE)
    train_pred <- predict(model, train.mat, decision.values = TRUE, probability = TRUE)
    prob.mat <- attr(train_pred, "probabilities")
    decision_values <- attr(train_pred, "decision.values")

    if (is.null(rownames(prob.mat))) {
      rownames(prob.mat) <- rownames(train.mat)
    }

    prob_col_1 <- rep(NA_real_, nrow(prob.mat))
    prob_col_2 <- rep(NA_real_, nrow(prob.mat))
    if ("1" %in% colnames(prob.mat)) {
      prob_col_1 <- prob.mat[, "1"]
    }
    if ("2" %in% colnames(prob.mat)) {
      prob_col_2 <- prob.mat[, "2"]
    }

    train_decision <- as.numeric(decision_values)
    if (length(train_decision) != nrow(prob.mat)) {
      stop("Unexpected decision.values length in training trace.", call. = FALSE)
    }
    if (stats::cor(train_decision, prob_col_1, use = "pairwise.complete.obs") < 0) {
      train_decision <- -train_decision
    }

    probability_count <- probability_count + 1
    iteration_prob_rows[[probability_count]] <- data.frame(
      kinase = kinase,
      ensemble = ensemble_idx,
      iteration = i,
      site = rownames(prob.mat),
      label = as.character(label[rownames(prob.mat)]),
      prob_class_1 = as.numeric(prob_col_1),
      prob_class_2 = as.numeric(prob_col_2),
      stringsAsFactors = FALSE
    )

    decision_count <- decision_count + 1
    iteration_decision_rows[[decision_count]] <- data.frame(
      kinase = kinase,
      ensemble = ensemble_idx,
      iteration = i,
      site = rownames(prob.mat),
      label = as.character(label[rownames(prob.mat)]),
      decision_value_class_1 = train_decision,
      stringsAsFactors = FALSE
    )

    X <- c()
    Y <- c()
    for (j in seq_len(ncol(prob.mat))) {
      class_name <- colnames(prob.mat)[j]
      voteClass <- prob.mat[label == class_name, , drop = FALSE]
      raw_weights <- as.numeric(voteClass[, j])
      normalized_weights <- raw_weights / sum(raw_weights)
      weight_count <- weight_count + 1
      iteration_weight_rows[[weight_count]] <- data.frame(
        kinase = kinase,
        ensemble = ensemble_idx,
        iteration = i,
        class_label = class_name,
        site = rownames(voteClass),
        raw_weight = raw_weights,
        normalized_weight = normalized_weights,
        stringsAsFactors = FALSE
      )
      idx <- sample(
        seq_len(nrow(voteClass)),
        size = nrow(voteClass),
        replace = TRUE,
        prob = raw_weights
      )
      sampled_sites <- rownames(voteClass)[idx]
      sample_count <- sample_count + 1
      iteration_sample_rows[[sample_count]] <- data.frame(
        kinase = kinase,
        ensemble = ensemble_idx,
        iteration = i,
        class_label = class_name,
        draw = seq_along(sampled_sites),
        site = sampled_sites,
        stringsAsFactors = FALSE
      )
      X <- rbind(X, train.mat[sampled_sites, , drop = FALSE])
      Y <- c(Y, label[sampled_sites])
    }
  }

  pred_obj <- predict(model, newdata = test.mat, probability = TRUE)
  pred <- attr(pred_obj, "probabilities")
  if (is.null(rownames(pred))) {
    rownames(pred) <- rownames(test.mat)
  }

  pred_col_1 <- rep(NA_real_, nrow(pred))
  pred_col_2 <- rep(NA_real_, nrow(pred))
  if ("1" %in% colnames(pred)) {
    pred_col_1 <- pred[, "1"]
  }
  if ("2" %in% colnames(pred)) {
    pred_col_2 <- pred[, "2"]
  }

  if (length(pred_decision) != nrow(pred)) {
    stop("Unexpected decision.values length in final prediction trace.", call. = FALSE)
  }
  if (stats::cor(pred_decision, pred_col_1, use = "pairwise.complete.obs") < 0) {
    pred_decision <- -pred_decision
  }

  ordered_sites <- rownames(pred)[order(pred_col_1, decreasing = TRUE)]
  top_count <- min(trace_top_n, length(ordered_sites))
  top_sites <- ordered_sites[seq_len(top_count)]

  list(
    pred = pred,
    iteration_probabilities = do.call(rbind, iteration_prob_rows),
    iteration_decision_values = do.call(rbind, iteration_decision_rows),
    iteration_resampling_weights = do.call(rbind, iteration_weight_rows),
    iteration_samples = do.call(rbind, iteration_sample_rows),
    final_predictions = data.frame(
      kinase = kinase,
      ensemble = ensemble_idx,
      site = rownames(pred),
      prob_class_1 = as.numeric(pred_col_1),
      prob_class_2 = as.numeric(pred_col_2),
      stringsAsFactors = FALSE
    ),
    final_decision_values = data.frame(
      kinase = kinase,
      ensemble = ensemble_idx,
      site = rownames(pred),
      decision_value_class_1 = pred_decision,
      stringsAsFactors = FALSE
    ),
    final_top = data.frame(
      kinase = kinase,
      ensemble = ensemble_idx,
      rank = seq_along(top_sites),
      site = top_sites,
      prob_class_1 = as.numeric(pred_col_1[match(top_sites, rownames(pred))]),
      stringsAsFactors = FALSE
    )
  )
}

trace_kinase_substrate_pred <- function(
  phosScoringMatrices,
  ensembleSize = 10,
  top = 50,
  cs = 0.8,
  inclusion = 20,
  iter = 5,
  trace_kinases = c("PRKAA1", "MAPK1", "MAPK9", "IRAK1", "TBK1", "LCK"),
  trace_top_n = 10
) {
  substrate_trace <- trace_substrate_list(
    phosScoringMatrices = phosScoringMatrices,
    top = top,
    cs = cs,
    inclusion = inclusion,
    trace_kinases = trace_kinases
  )
  substrate.list <- substrate_trace$substrate_list

  print("Predicting kinases for phosphosites:")
  featureMat <- phosScoringMatrices$combinedScoreMatrix
  predMatrix <- matrix(0, nrow = nrow(featureMat), ncol = length(substrate.list))
  colnames(predMatrix) <- names(substrate.list)
  rownames(predMatrix) <- rownames(featureMat)

  initial_negative_rows <- list()
  iteration_probability_rows <- list()
  iteration_decision_rows <- list()
  iteration_weight_rows <- list()
  iteration_sample_rows <- list()
  final_prediction_rows <- list()
  final_decision_rows <- list()
  final_top_rows <- list()
  initial_count <- 0
  probability_count <- 0
  decision_count <- 0
  weight_count <- 0
  sample_count <- 0
  final_prediction_count <- 0
  final_decision_count <- 0
  final_top_count <- 0

  for (i in seq_len(length(substrate.list))) {
    kinase <- names(substrate.list)[i]
    positive.train <- featureMat[substrate.list[[i]], , drop = FALSE]
    positive.cls <- rep(1, length(substrate.list[[i]]))
    negative.pool <- featureMat[!(rownames(featureMat) %in% substrate.list[[i]]), , drop = FALSE]
    cat(paste(i, ".", sep = ""))

    for (e in seq_len(ensembleSize)) {
      negativeSize <- length(substrate.list[[i]])
      idx <- sample(seq_len(nrow(negative.pool)), size = negativeSize, replace = FALSE)
      negative.samples <- rownames(negative.pool)[idx]
      negative.train <- featureMat[negative.samples, , drop = FALSE]
      negative.cls <- rep(2, length(negative.samples))
      train.mat <- rbind(positive.train, negative.train)
      cls <- as.factor(c(positive.cls, negative.cls))
      names(cls) <- rownames(train.mat)

      if (kinase %in% trace_kinases) {
        initial_count <- initial_count + 1
        initial_negative_rows[[initial_count]] <- data.frame(
          kinase = kinase,
          ensemble = e,
          draw = seq_along(negative.samples),
          site = negative.samples,
          stringsAsFactors = FALSE
        )

        trace_result <- trace_multi_ada_sampling(
          train.mat = train.mat,
          test.mat = featureMat,
          label = cls,
          kernelType = "radial",
          iter = iter,
          kinase = kinase,
          ensemble_idx = e,
          trace_top_n = trace_top_n
        )
        pred <- trace_result$pred
        probability_count <- probability_count + 1
        iteration_probability_rows[[probability_count]] <- trace_result$iteration_probabilities
        decision_count <- decision_count + 1
        iteration_decision_rows[[decision_count]] <- trace_result$iteration_decision_values
        weight_count <- weight_count + 1
        iteration_weight_rows[[weight_count]] <- trace_result$iteration_resampling_weights
        sample_count <- sample_count + 1
        iteration_sample_rows[[sample_count]] <- trace_result$iteration_samples
        final_prediction_count <- final_prediction_count + 1
        final_prediction_rows[[final_prediction_count]] <- trace_result$final_predictions
        final_decision_count <- final_decision_count + 1
        final_decision_rows[[final_decision_count]] <- trace_result$final_decision_values
        final_top_count <- final_top_count + 1
        final_top_rows[[final_top_count]] <- trace_result$final_top
      } else {
        pred <- multiAdaSampling_fn(
          train.mat,
          test.mat = featureMat,
          label = cls,
          kernelType = "radial",
          iter = iter
        )
      }

      predMatrix[rownames(pred), i] <- predMatrix[rownames(pred), i] + pred[, 1]
    }
  }

  predMatrix <- predMatrix / ensembleSize
  print("done")

  list(
    pred_matrix = predMatrix,
    substrate_list = substrate.list,
    trace_candidates = substrate_trace$candidates,
    trace_initial_negatives = if (length(initial_negative_rows) > 0) do.call(rbind, initial_negative_rows) else data.frame(),
    trace_iteration_probabilities = if (length(iteration_probability_rows) > 0) do.call(rbind, iteration_probability_rows) else data.frame(),
    trace_iteration_decision_values = if (length(iteration_decision_rows) > 0) do.call(rbind, iteration_decision_rows) else data.frame(),
    trace_iteration_resampling_weights = if (length(iteration_weight_rows) > 0) do.call(rbind, iteration_weight_rows) else data.frame(),
    trace_iteration_samples = if (length(iteration_sample_rows) > 0) do.call(rbind, iteration_sample_rows) else data.frame(),
    trace_final_ensemble_predictions = if (length(final_prediction_rows) > 0) do.call(rbind, final_prediction_rows) else data.frame(),
    trace_final_ensemble_decision_values = if (length(final_decision_rows) > 0) do.call(rbind, final_decision_rows) else data.frame(),
    trace_final_ensemble_top = if (length(final_top_rows) > 0) do.call(rbind, final_top_rows) else data.frame()
  )
}

main <- function() {
  args <- read_args()
  outdir <- args$outdir
  trace_kinases <- parse_csv_values(args$trace_kinases)
  trace_top_n <- as.integer(args$trace_top_n)
  if (is.na(trace_top_n) || trace_top_n < 1) {
    stop("trace_top_n must be a positive integer.", call. = FALSE)
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  trace_dir <- file.path(outdir, "prediction_trace")
  dir.create(trace_dir, recursive = TRUE, showWarnings = FALSE)

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
  native_num_motif <- 5L
  native_num_sub <- 1L
  native_top <- 30L
  native_score_threshold <- 0.6
  native_inclusion <- 5L

  L6.matrices <- kinase_substrate_score_compat(
    substrate.list = substrate_list,
    mat = L6.phos.std,
    seqs = L6.phos.seq,
    motif.mouse.list = motif.mouse.list,
    numMotif = native_num_motif,
    numSub = native_num_sub
  )

  native_substrate_map <- build_quantified_substrate_map(
    substrate.list = substrate_list,
    observed_sites = rownames(L6.phos.std),
    min_substrates = native_num_sub
  )
  native_profile_list <- PhosR::kinaseSubstrateProfile(substrate_list, L6.phos.std)
  native_profile_filtered <- native_profile_list[
    which(native_profile_list$NumSub >= native_num_sub)
  ]
  native_profile_matrix <- do.call(rbind, native_profile_filtered)
  rownames(native_profile_matrix) <- names(native_profile_filtered)
  native_profile_scores <- score_phosphosite_profile_compat(
    L6.phos.std,
    native_profile_filtered
  )
  native_motif_filtered <- motif.mouse.list[
    which(motif.mouse.list$NumInputSeq >= native_num_motif)
  ]
  native_motif_scores <- score_phosphosites_motifs_compat(
    mat = L6.phos.std,
    motif.mouse.list.filtered = native_motif_filtered,
    seqs = L6.phos.seq
  )
  native_motif_sizes <- data.frame(
    kinase = names(motif.mouse.list$NumInputSeq),
    motif_size = as.numeric(motif.mouse.list$NumInputSeq),
    stringsAsFactors = FALSE
  )
  native_combined_weights <- build_combined_weight_table(
    overlap_kinases = colnames(L6.matrices$combinedScoreMatrix),
    motif_sizes = motif.mouse.list$NumInputSeq,
    profile_sizes = native_profile_list$NumSub
  )
  native_candidate_trace <- trace_substrate_list(
    phosScoringMatrices = L6.matrices,
    top = native_top,
    cs = native_score_threshold,
    inclusion = native_inclusion,
    trace_kinases = character(0)
  )

  set.seed(1)
  L6.prediction <- trace_kinase_substrate_pred(
    L6.matrices,
    top = native_top,
    cs = native_score_threshold,
    inclusion = native_inclusion,
    trace_kinases = trace_kinases,
    trace_top_n = trace_top_n
  )
  L6.predMat <- L6.prediction$pred_matrix
  native_prediction_top30 <- build_prediction_top_table(L6.predMat, top_n = native_top)

  kinase_activity <- build_weighted_kinase_activity(L6.predMat, L6.phos.std)
  ksea <- build_ksea_scores(L6.predMat, L6.phos.std)

  write.csv(L6.phos.std, file.path(outdir, "l6_phospho_matrix.csv"), row.names = TRUE)
  write.csv(
    data.frame(site_id = names(L6.phos.seq), centralized_sequence = unname(L6.phos.seq)),
    file.path(outdir, "l6_site_sequences.csv"),
    row.names = FALSE
  )
  write.csv(
    flatten_grouped_mapping(native_substrate_map),
    file.path(outdir, "native_substrate_map.csv"),
    row.names = FALSE
  )
  write.csv(native_profile_matrix, file.path(outdir, "native_profile_matrix.csv"), row.names = TRUE)
  write.csv(native_profile_scores, file.path(outdir, "native_profile_scores.csv"), row.names = TRUE)
  write.csv(native_motif_scores, file.path(outdir, "native_motif_scores.csv"), row.names = TRUE)
  write.csv(native_motif_sizes, file.path(outdir, "native_motif_sizes.csv"), row.names = FALSE)
  write.csv(
    L6.matrices$combinedScoreMatrix,
    file.path(outdir, "native_combined_scores.csv"),
    row.names = TRUE
  )
  write.csv(
    native_combined_weights,
    file.path(outdir, "native_combined_weights.csv"),
    row.names = FALSE
  )
  write.csv(
    flatten_grouped_mapping(native_candidate_trace$substrate_list),
    file.path(outdir, "native_candidate_substrates.csv"),
    row.names = FALSE
  )
  write.csv(
    native_prediction_top30,
    file.path(outdir, "native_prediction_top30.csv"),
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
      n_targets = vapply(
        colnames(L6.predMat),
        function(k) sum(L6.predMat[, k] > 0.6, na.rm = TRUE),
        integer(1)
      )
    ),
    file.path(outdir, "kinase_target_counts.csv"),
    row.names = FALSE
  )

  write.csv(L6.prediction$trace_candidates, file.path(trace_dir, "trace_candidates.csv"), row.names = FALSE)
  write.csv(L6.prediction$trace_initial_negatives, file.path(trace_dir, "trace_initial_negatives.csv"), row.names = FALSE)
  write.csv(L6.prediction$trace_iteration_probabilities, file.path(trace_dir, "trace_iteration_probabilities.csv"), row.names = FALSE)
  write.csv(L6.prediction$trace_iteration_resampling_weights, file.path(trace_dir, "trace_iteration_resampling_weights.csv"), row.names = FALSE)
  write.csv(L6.prediction$trace_iteration_samples, file.path(trace_dir, "trace_iteration_samples.csv"), row.names = FALSE)
  write.csv(L6.prediction$trace_final_ensemble_predictions, file.path(trace_dir, "trace_final_ensemble_predictions.csv"), row.names = FALSE)
  write.csv(L6.prediction$trace_final_ensemble_top, file.path(trace_dir, "trace_final_ensemble_top.csv"), row.names = FALSE)
  write_trace_readme(trace_dir, trace_kinases, trace_top_n)

  write_session_info(outdir)
  message("Done. R L6 reference fixtures written to: ", outdir)
}

main()