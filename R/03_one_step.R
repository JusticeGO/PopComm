#' Analyze Ligand-Receptor Pair Correlations and Projection Scores (Specified Sender and Receiver)
#'
#' @description
#' Performs integrated analysis of ligand-receptor (LR) pairs through two consecutive phases:
#' 1. Filters LR pairs and analyzes correlations between specified cell types.
#' 2. Calculates projection scores based on regression models for valid pairs.
#' Returns comprehensive results combining statistical metrics.
#'
#' @param rna A Seurat object containing single-cell RNA expression data.
#' @param sender Cell type designated as the ligand sender (character).
#' @param receiver Cell type designated as the receptor receiver (character).
#' @param lr_database Data frame of LR pairs with columns "ligand_gene_symbol" and "receptor_gene_symbol".
#' @param sample_col Column name in metadata indicating sample identifiers (character).
#' @param cell_type_col Column name in metadata indicating cell type classifications (character).
#' @param min_cells Minimum cells required per sample for both sender and receiver (numeric, default 50).
#' @param min_samples Minimum valid samples required of LR filter (numeric, default 10).
#' @param min_cell_ratio Minimum ratio of cells expressing ligand and receptor genes in sender or receiver cells (numeric, default 0.1).
#' @param min_sample_ratio Minimum ratio of samples in which both the ligand and receptor genes must be expressed (numeric, default 0.1).
#' @param cor_method Correlation method: "spearman" (default), "pearson", or "kendall".
#' @param adjust_method P-value adjustment method (default "BH" for Benjamini-Hochberg).
#'        Options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param min_adjust_p Adjusted p-value threshold for significance (numeric, default 0.05).
#' @param min_cor Minimum correlation coefficient threshold (numeric, default 0). Must be ≥ 0.
#' @param num_cores Number of CPU cores for parallel processing (numeric, default 10). Automatically capped at (system cores - 1).
#'
#' @return Two data frames with columns:
#'   \item{ligand, receptor}{Ligand and receptor gene symbols (res1/res2).}
#'   \item{cor}{Correlation coefficient (res1/res2).}
#'   \item{p_val}{Raw p-value (res1/res2).}
#'   \item{adjust.p}{Adjusted p-value (res1/res2).}
#'   \item{sender, receiver}{Sender and receiver cell types (res1/res2).}
#'   \item{slope}{Slope of the linear regression model (res1/res2).}
#'   \item{intercept}{Intercept of the linear regression model (res1/res2).}
#'   \item{sample}{Sample identifier (res2).}
#'   \item{score}{Projection score (raw co-expression intensity) (res2).}
#'   \item{normalized_score}{Normalized score scaled between 0-1 (res2).}
#'   Returns \code{NULL} if no valid pairs in either phase.
#'
#' @export
#'
#' @importFrom dplyr %>% bind_rows
#' @importFrom stats cor.test p.adjust lm sd coef na.omit
#' @importFrom utils head packageVersion
#' @importFrom Matrix rowSums
#'
#' @examples
#' seurat_object <- load_example_seurat()
#' data(lr_db)
#'
#' # Integrated analysis with Cardiac -> Perivascular
#' res_single <- one_step_single(
#'   rna = seurat_object,
#'   sender = "Cardiac",
#'   receiver = "Perivascular",
#'   lr_database = lr_db,
#'   sample_col = "sample",
#'   cell_type_col = "cell.type",
#'   min_cells = 20,
#'   min_samples = 10,
#'   min_adjust_p = 0.5,
#'   num_cores = 1
#' )
#'
#' head(res_single$res1)
#' head(res_single$res2)
one_step_single <- function(rna, sender, receiver, lr_database = PopComm::lr_db,
                            sample_col, cell_type_col,
                            min_cells = 50, min_samples = 10,
                            min_cell_ratio = 0.1, min_sample_ratio = 0.1,
                            cor_method = "spearman", adjust_method = "BH",
                            min_adjust_p = 0.05, min_cor = 0,
                            num_cores = 10) {

  message("One-step analysis of receptor-ligand interaction: ", sender, " -> ", receiver)

  # Check parameters
  max_cores <- parallel::detectCores()
  if (num_cores > max_cores) {
    message("Warning: Using more cores (", num_cores, ") than available (", max_cores, ").")
    message("Using ", max_cores - 1, " cores instead of requested ", num_cores)
    num_cores <- max_cores - 1
  }

  # Pre-process metadata
  rna$sample <- rna@meta.data[, sample_col]
  rna$cell.type <- rna@meta.data[, cell_type_col]
  cell_types <- unique(rna@meta.data[[cell_type_col]])
  if (length(cell_types) < 1) stop("No cell types found.")

  selected_types <- unique(c(sender, receiver))
  missing_types <- setdiff(selected_types, cell_types)
  if (length(missing_types) > 0) {
    stop("Missing cell types: ", paste(missing_types, collapse = ", "))
  }

  # Determine the subset of data
  if (!setequal(selected_types, cell_types)) {
    message("Subsetting ", sender, " (sender) and ", receiver, " (receiver)...")
    cells.to.keep <- rna$cell.type %in% selected_types
    rna.data <- rna[, cells.to.keep]
  } else {
    rna.data <- rna
  }

  # Step 1: 01_filter_interaction
  message("\nStep 1: Filter interaction")
  message("Analyzing ligand-receptor interactions: ", sender, " -> ", receiver)

  # Filter samples based on cell counts per sample
  message("Filtering samples with cell counts...")
  cell_counts <- table(rna.data$sample, rna.data$cell.type)
  if (sender != receiver) {
    valid_samples <- names(which(
      cell_counts[, sender] > min_cells &
        cell_counts[, receiver] > min_cells
      ))
  } else {
    valid_samples <- names(which(
      cell_counts[, sender] > min_cells
      ))
  }
  message("Remaining samples after filtering: ", length(valid_samples))
  if (length(valid_samples) < min_samples) {
    message("Insufficient valid samples (", length(valid_samples), " < ", min_samples, "). Analysis stopped.")
    return(NULL)
  }
  rna.data <- subset(rna.data, sample %in% valid_samples)

  # Filter LR pairs based on minimum cell ratio
  message("Filtering LR pairs based on minimum cell ratio in sender or receiver cells...")
  lr <- lr_database
  lr <- lr[which(lr$ligand_gene_symbol %in% rownames(rna.data)), ]
  lr <- lr[which(lr$receptor_gene_symbol %in% rownames(rna.data)), ]

  sender_cells <- colnames(rna.data)[rna.data$cell.type == sender]
  receiver_cells <- colnames(rna.data)[rna.data$cell.type == receiver]

  seurat_version <- as.numeric(strsplit(as.character(packageVersion("Seurat")), "\\.")[[1]][[1]])
  if (seurat_version >= 5) {
    expr <- Seurat::GetAssayData(rna.data, layer = "data")  # Seurat v5
  } else {
    expr <- Seurat::GetAssayData(rna.data, slot = "data")   # Seurat v4
  }
  # expr <- Seurat::GetAssayData(rna.data, slot = "data")
  expr <- as.matrix(expr)

  sender_ratio <- rowSums(expr[, sender_cells, drop = FALSE] > 0) / length(sender_cells)
  receiver_ratio <- rowSums(expr[, receiver_cells, drop = FALSE] > 0) / length(receiver_cells)

  lr <- lr[lr$ligand_gene_symbol %in% names(sender_ratio[sender_ratio > min_cell_ratio]), ]
  lr <- lr[lr$receptor_gene_symbol %in% names(receiver_ratio[receiver_ratio > min_cell_ratio]), ]

  # Compute average expression for each sample-cell type group
  rna.data$group <- paste0(rna.data$sample, "-lr-", rna.data$cell.type)
  message("Computing average expression for each sample-cell type group...")
  rna.avg <- Seurat::AverageExpression(rna.data, group.by = "group")[[1]]   # seurat v4/v5
  rna.avg <- round(rna.avg, 5)

  avg.s <- rna.avg[, grep(sender, colnames(rna.avg))]
  avg.r <- rna.avg[, grep(receiver, colnames(rna.avg))]

  colnames(avg.s) <- stringr::str_match(colnames(avg.s), "^(.*)-lr-")[,2]
  colnames(avg.r) <- stringr::str_match(colnames(avg.r), "^(.*)-lr-")[,2]

  avg.r <- avg.r[, colnames(avg.s), drop = FALSE]

  avg.s <- avg.s[lr$ligand_gene_symbol, , drop = FALSE]
  avg.r <- avg.r[lr$receptor_gene_symbol, , drop = FALSE]

  message("Starting correlation and filtering process for ligand-receptor pairs...")

  # Calculate correlation and linear models
  calc_correlation <- function(i) {
    x <- as.numeric(avg.s[lr$ligand_gene_symbol[i], ])
    y <- as.numeric(avg.r[lr$receptor_gene_symbol[i], ])

    # filter sample
    # data_df <- data.frame(x = x, y = y)
    # data_df <- remove_outlier(data_df)
    data_df <- data.frame(x = x, y = y) |> remove_outlier() |> na.omit()
    p <- data_df$x
    q <- data_df$y

    if (nrow(data_df) < min_samples || length(p) == 0 || length(q) == 0 || sum(p) == 0 || sum(q) == 0 || sd(p) == 0 || sd(q) == 0) {
      return(NULL)
    }

    # Fitting Linear Models
    model <- tryCatch(
      lm(q ~ p),
      error = function(e) NULL
    )
    if (is.null(model)) return(NULL)

    slope <- round(coef(model)[2], 5)
    intercept <- round(coef(model)[1], 5)

    # cor.test
    res_cor <- tryCatch(
      cor.test(p, q, method = cor_method),
      error = function(e) NULL
      )
    if (is.null(res_cor)) return(NULL)

    pct1 <- round(sum(p > 0) / length(p), 3)
    pct2 <- round(sum(q > 0) / length(q), 3)

    lr_name <- paste0(row.names(avg.s)[i], "_", row.names(avg.r)[i])

    return(c(
      round(res_cor$estimate, 5),
      round(res_cor$p.value, 15),
      lr_name,
      pct1, pct2,
      slope, intercept
    ))
  }

  res_list <- run_parallel(
    1:nrow(avg.r),
    calc_correlation,
    num_cores = num_cores,
    export_vars = c("avg.s", "avg.r", "lr", "min_samples", "cor_method", "remove_outlier")
    )

  # res_mat <- do.call(rbind, res_list)
  res_mat <- do.call(rbind, Filter(Negate(is.null), res_list))

  if (is.null(res_mat)) {
    message("No valid ligand-receptor pairs passed the initial filtering.")
    return(NULL)
  }
  res <- data.frame(res_mat, stringsAsFactors = FALSE)

  colnames(res) <- c("cor", "p_val", "lr", "pct1", "pct2", "slope", "intercept")
  num_cols <- c("cor", "p_val", "pct1", "pct2", "slope", "intercept")
  res[num_cols] <- lapply(res[num_cols], as.numeric)

  res$adjust.p <- round(p.adjust(res$p_val, method = adjust_method), 15)

  # Filter the results based on adjusted p-value, correlation, and ratio thresholds
  message("Filtering results based on adjusted p-value, correlation, and ratio thresholds...")
  res <- res[which(res$adjust.p < min_adjust_p &
                     res$cor > min_cor &
                     res$pct1 > min_sample_ratio &
                     res$pct2 > min_sample_ratio), ]
  res <- res[order(res$adjust.p, -res$cor), ]
  if (nrow(res) > 0) {
    row.names(res) <- 1:nrow(res)
  }

  if (nrow(res) == 0) {
    message("No results meet the filtering criteria. Returning NULL.")
    return(NULL)
  }

  res$sender <- sender
  res$receiver <- receiver

  res$ligand <- stringr::str_match(res$lr, "^(.*)_")[,2]
  res$receptor <- stringr::str_match(res$lr, "_(.*)$")[,2]

  selected_cols <- c("ligand", "receptor", "cor", "p_val", "adjust.p", "sender", "receiver", "slope", "intercept")
  res <- res %>% select(all_of(selected_cols))

  message("Filter and correlation process complete.")
  message("Head of results (", nrow(res), "):")
  print(head(res))


  # Step 2: 02_score_sample
  message("\nStep 2: Calculate the Projection Score")
  message("Analyzing ligand-receptor projection scores: ", sender, " -> ", receiver)

  # Compute average expression for each sample-cell type group
  avg.s_sub <- avg.s[res$ligand, , drop = FALSE]
  avg.r_sub <- avg.r[res$receptor, , drop = FALSE]

  # Calculating projection scores
  calc_projection <- function(i) {
    x <- as.numeric(avg.s_sub[i, ])
    y <- as.numeric(avg.r_sub[i, ])

    if (length(x) == 0 || length(y) == 0 || all(is.na(x)) || all(is.na(y)) || sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
      return(data.frame())
    }

    slope <- res$slope[i]
    intercept <- res$intercept[i]

    projections <- t(sapply(
      1:length(x), function(j) {
        project_to_line(x[j], y[j], slope, intercept)
      }))

    dx <- projections[, 1] - min(projections[, 1])
    dy <- projections[, 2] - min(projections[, 2])
    score <- round(sqrt(dx^2 + dy^2), 5)
    normalized_score <- round(score / max(score), 5)

    lr_metadata <- res[i, ]
    result <- data.frame(
      lr_metadata,
      sample = colnames(avg.s_sub),
      score = score,
      normalized_score = normalized_score,
      row.names = NULL
    )

    result <- result[order(-result$score), ]
    return(result)
  }

  score_list <- run_parallel(
    1:nrow(avg.s_sub),
    calc_projection,
    num_cores = num_cores,
    export_vars = c("avg.s_sub", "avg.r_sub", "res", "project_to_line")
    )

  score_list <- score_list[sapply(score_list, function(x) nrow(x) > 0)]
  score.df <- bind_rows(score_list) %>% na.omit()

  message("Analyzing ligand-receptor projection scores process complete.")
  message("Head of projection score results (", nrow(score.df), "):")
  print(head(score.df))

  message("\n", sender, " -> ", receiver, " analyzing process complete.\n")

  return(list(res1 = res, res2 = score.df))
}




#' Analyze Ligand-Receptor Pair Correlations and Projection Scores (Across All Cell Types)
#'
#' @description
#' Performs integrated analysis of ligand-receptor (LR) pairs through two consecutive phases:
#' 1. Filters LR pairs and analyzes correlations across all cell types.
#' 2. Calculates projection scores based on regression models for valid pairs.
#' Returns comprehensive results combining statistical metrics.
#'
#' @param rna A Seurat object containing single-cell RNA expression data.
#' @param lr_database Data frame of LR pairs with columns "ligand_gene_symbol" and "receptor_gene_symbol".
#' @param sample_col Column name in metadata indicating sample identifiers (character).
#' @param cell_type_col Column name in metadata indicating cell type classifications (character).
#' @param min_cells Minimum cells required per sample for both sender and receiver (numeric, default 50).
#' @param min_samples Minimum valid samples required for LR filter (numeric, default 10).
#' @param min_cell_ratio Minimum ratio of cells expressing ligand and receptor genes in sender or receiver cells (numeric, default 0.1).
#' @param min_sample_ratio Minimum ratio of samples in which both the ligand and receptor genes must be expressed (numeric, default 0.1).
#' @param cor_method Correlation method: "spearman" (default), "pearson", or "kendall".
#' @param adjust_method P-value adjustment method (default "BH" for Benjamini-Hochberg).
#'        Options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param min_adjust_p Adjusted p-value threshold for significance (numeric, default 0.05).
#' @param min_cor Minimum correlation coefficient threshold (numeric, default 0). Must be ≥ 0.
#' @param num_cores Number of CPU cores for parallel processing (numeric, default 10). Automatically capped at (system cores - 1).
#'
#' @return Two data frames with columns:
#'   \item{ligand, receptor}{Ligand and receptor gene symbols (res1/res2).}
#'   \item{cor}{Correlation coefficient (res1/res2).}
#'   \item{p_val}{Raw p-value (res1/res2).}
#'   \item{adjust.p}{Adjusted p-value (res1/res2).}
#'   \item{sender, receiver}{Sender and receiver cell types (res1/res2).}
#'   \item{slope}{Slope of the linear regression model (res1/res2).}
#'   \item{intercept}{Intercept of the linear regression model (res1/res2).}
#'   \item{sample}{Sample identifier (res2).}
#'   \item{score}{Projection score (raw co-expression intensity) (res2).}
#'   \item{normalized_score}{Normalized score scaled between 0-1 (res2).}
#'   Returns \code{NULL} if no valid pairs in either phase.
#'
#' @export
#'
#' @importFrom dplyr %>% bind_rows
#' @importFrom stats cor.test p.adjust lm sd coef na.omit
#' @importFrom utils head
#' @importFrom Matrix rowSums
#'
#' @examples
#' seurat_object <- load_example_seurat()
#' data(lr_db)
#'
#' # Integrated analysis across all cell types
#' res_all <- one_step_all(
#'   rna = seurat_object,
#'   lr_database = lr_db,
#'   sample_col = "sample",
#'   cell_type_col = "cell.type",
#'   min_cells = 20,
#'   min_samples = 10,
#'   min_adjust_p = 0.5,
#'   num_cores = 1
#' )
#'
#' head(res_all$res1)
#' head(res_all$res2)
one_step_all <- function(rna, lr_database,
                         sample_col, cell_type_col,
                         min_cells = 50, min_samples = 10,
                         min_cell_ratio = 0.1, min_sample_ratio = 0.1,
                         cor_method = "spearman", adjust_method = "BH",
                         min_adjust_p = 0.05, min_cor = 0,
                         num_cores = 10) {

  message("\nOne-step analysis of receptor-ligand interaction: For all possible cell type pairs")

  # Check parameters
  max_cores <- parallel::detectCores()
  if (num_cores > max_cores) {
    message("Warning: Using more cores (", num_cores, ") than available (", max_cores, ").")
    message("Using ", max_cores - 1, " cores instead of requested ", num_cores)
    num_cores <- max_cores - 1
  }

  # Pre-process metadata
  rna$sample <- rna@meta.data[, sample_col]
  rna$cell.type <- rna@meta.data[, cell_type_col]
  cell_types <- unique(rna@meta.data[[cell_type_col]])
  if (length(cell_types) < 1) stop("No cell types found.")
  message("Cell types: ", paste(cell_types, collapse = ", "))

  run_one_step <- function(rna_obj, sender, receiver) {
    res <- one_step_single(
      rna = rna_obj,
      sender = sender,
      receiver = receiver,
      lr_database = lr_database,
      sample_col = "sample",
      cell_type_col = "cell.type",
      min_cells = min_cells,
      min_samples = min_samples,
      cor_method = cor_method,
      adjust_method = adjust_method,
      min_adjust_p = min_adjust_p,
      min_cor = min_cor,
      min_sample_ratio = min_sample_ratio,
      min_cell_ratio = min_cell_ratio,
      num_cores = num_cores
    )
    if (!is.null(res) && is.list(res) && !is.null(res$res1) && nrow(res$res1) > 0) {
      return(res)
    }
    return(NULL)
  }

  # Handles pairing of the same cell type
  message("\nProcessing same cell type pairs...")
  same_ct_results <- lapply(cell_types, function(ct) {
    message("\n  Analyzing pair: ", ct, " <-> ", ct)
    run_one_step(rna, sender = ct, receiver = ct)
  })
  same_ct_results <- Filter(Negate(is.null), same_ct_results)

  # Processing the pairing of different cell types
  message("\nProcessing different cell type pairs...")
  diff_ct_results <- list()
  unique_pairs <- combn(cell_types, 2, simplify = FALSE)

  for (pair in unique_pairs) {
    ct1 <- pair[1]
    ct2 <- pair[2]
    message("\n  Analyzing pair: ", ct1, " <-> ", ct2)

    # subset data
    cells_to_keep <- rna$cell.type %in% c(ct1, ct2)
    rna_subset <- rna[, cells_to_keep]
    if (length(unique(rna_subset$cell.type)) < 2) {
      message("  Skipping pair ", ct1, " and ", ct2, " due to insufficient cell types.")
      next
    }

    # ct1 -> ct2
    res_forward <- run_one_step(rna_subset, sender = ct1, receiver = ct2)
    # ct2 -> ct1
    res_reverse <- run_one_step(rna_subset, sender = ct2, receiver = ct1)

    if (!is.null(res_forward)) {
      diff_ct_results <- c(diff_ct_results, list(res_forward))
    }
    if (!is.null(res_reverse)) {
      diff_ct_results <- c(diff_ct_results, list(res_reverse))
    }
  }

  # Helper functions: Combine results
  combine_results <- function(results_list, element) {
    if (length(results_list) == 0) return(data.frame())
    do.call(rbind, lapply(results_list, function(x) x[[element]]))
  }

  all_results_res1 <- bind_rows(
    combine_results(same_ct_results, "res1"),
    combine_results(diff_ct_results, "res1")
  )
  all_results_res2 <- bind_rows(
    combine_results(same_ct_results, "res2"),
    combine_results(diff_ct_results, "res2")
  )

  final_res <- list(res1 = all_results_res1, res2 = all_results_res2)

  # Check that the result is not empty
  if (nrow(final_res$res1) == 0 && nrow(final_res$res2) == 0) {
    message("\nNo results meet the filtering criteria. Returning NULL.")
    return(NULL)
  }

  message("\nFor all possible cell type pairs analyzing process is complete.")
  message("Head of filter and correlation results (", nrow(final_res$res1), "):")
  print(head(final_res$res1))

  message("Head of LR projection scores results (", nrow(final_res$res2), "):")
  print(head(final_res$res2))

  return(final_res)
}
