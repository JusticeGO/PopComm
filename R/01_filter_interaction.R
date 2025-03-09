#' Filter and Analyze Ligand-Receptor Pair Correlations (Specified Sender and Receiver)
#'
#' @description
#' Filters ligand-receptor (LR) pairs and analyzes their correlations for specified sender and receiver cell types,
#' and returning significant LR pairs based on user-defined thresholds.
#'
#' @param rna A Seurat object containing single-cell RNA expression data.
#' @param sender Cell type designated as the ligand sender (character).
#' @param receiver Cell type designated as the receptor receiver (character).
#' @param lr_database A data frame of ligand-receptor pairs with columns "ligand_gene_symbol" and "receptor_gene_symbol".
#' @param sample_col Column name in Seurat metadata indicating sample identifiers (character).
#' @param cell_type_col Column name in Seurat metadata indicating cell type classifications (character).
#' @param min_cells Minimum cells required per sample for both sender and receiver (default 50).
#' @param min_samples Minimum valid samples required to proceed (default 10).
#' @param cor_method Correlation method: "spearman" (default), "pearson", or "kendall".
#' @param adjust_method P-value adjustment method (default "BH" for Benjamini-Hochberg).
#'        Options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param min_adjust_p Adjusted p-value threshold for significance (default 0.05).
#' @param min_cor Minimum correlation coefficient threshold (default 0). Must be ≥ 0.
#' @param min_pct Minimum percentage of non-zero expression in both ligand and receptor (default 0.1).
#' @param mc_cores Number of CPU cores for parallel processing (default 10). Automatically capped at (system cores - 1).
#'
#' @return A data frame of filtered LR pairs with columns:
#'   \item{\code{cor_<method>}}{Correlation coefficient (e.g., \code{cor_spearman}).}
#'   \item{\code{p_<method>}}{Raw p-value (e.g., \code{p_spearman}).}
#'   \item{pct1, pct2}{Percent non-zero expression in sender and receiver.}
#'   \item{lr}{Ligand-receptor pair identifier (\code{"Ligand_Receptor"}).}
#'   \item{adjust.p}{Adjusted p-value.}
#'   \item{sender, receiver}{Sending/receiving cell types.}
#'   Rows ordered by ascending \code{adjust.p} and descending \code{cor_<method>}.
#'   Returns \code{NULL} if no pairs pass filters.
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom stats lm sd coef cor.test p.adjust
#' @importFrom utils head
#'
#' @examples
#' seurat_object <- load_example_seurat()
#' data(lr_db)
#'
#' # Analyzing ligand-receptor interactions: Cardiac -> Perivascular
#' result01s <- filter_lr_single(
#'   rna = seurat_object,
#'   sender = "Cardiac",
#'   receiver = "Perivascular",
#'   lr_database = lr_db,
#'   sample_col = "sample",
#'   cell_type_col = "cell.type",
#'   min_cells = 20,
#'   min_samples = 10,
#'   min_adjust_p = 0.5,
#'   mc_cores = 1
#' )
filter_lr_single <- function(rna, sender, receiver, lr_database = PopComm::lr_db,
                             sample_col, cell_type_col,
                             min_cells = 50, min_samples = 10,
                             cor_method = "spearman", adjust_method = "BH",
                             min_adjust_p = 0.05, min_cor = 0, min_pct = 0.1,
                             mc_cores = 10) {

  # Check parameters
  max_cores <- parallel::detectCores()
  if (mc_cores > max_cores) {
    message("Warning: Using more cores (", mc_cores, ") than available (", max_cores, ").")
    message("Using ", max_cores - 1, " cores instead of requested ", mc_cores)
    mc_cores <- max_cores - 1
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

  # Step 1: Add sample and cell type columns to the RNA data object and subset
  message("Analyzing ligand-receptor interactions: ", sender, " -> ", receiver)

  # Determine the subset of data
  if (!setequal(selected_types, cell_types)) {
    message("Subsetting ", sender, " (sender) and ", receiver, " (receiver)...")
    cells.to.keep <- rna$cell.type %in% selected_types
    rna.data <- rna[, cells.to.keep]
  } else {
    rna.data <- rna
  }

  # Step 2: Filter samples based on cell counts per sample
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

  # Step 3: Filter ligand-receptor pairs that exist in RNA data
  lr <- lr_database
  lr <- lr[which(lr$ligand_gene_symbol %in% rownames(rna.data)), ]
  lr <- lr[which(lr$receptor_gene_symbol %in% rownames(rna.data)), ]

  # Step 4: Compute average expression for each sample-cell type group
  rna.data$group <- paste0(rna.data$sample, "-lr-", rna.data$cell.type)
  message("Computing average expression for each sample-cell type group...")
  rna.avg <- Seurat::AverageExpression(rna.data, group.by = "group")$RNA
  rna.avg <- round(rna.avg, 5)

  avg.s <- rna.avg[, grep(sender, colnames(rna.avg))]
  avg.r <- rna.avg[, grep(receiver, colnames(rna.avg))]

  colnames(avg.s) <- stringr::str_match(colnames(avg.s), "^(.*)-lr-")[,2]
  colnames(avg.r) <- stringr::str_match(colnames(avg.r), "^(.*)-lr-")[,2]

  avg.r <- avg.r[, colnames(avg.s), drop = FALSE]

  avg.s <- avg.s[lr$ligand_gene_symbol, , drop = FALSE]
  avg.r <- avg.r[lr$receptor_gene_symbol, , drop = FALSE]

  # Step 5: Compute correlations between ligand-receptor pairs
  message("Starting correlation and filtering process for ligand-receptor pairs...")

  cor_colname <- paste0("cor_", cor_method)
  p_colname <- paste0("p_", cor_method)

  # Calculate correlation and linear models
  calc_correlation <- function(i) {
    x <- as.numeric(avg.s[lr$ligand_gene_symbol[i], ])
    y <- as.numeric(avg.r[lr$receptor_gene_symbol[i], ])

    # filter sample
    data_df <- data.frame(x = x, y = y)
    data_df <- remove_outlier(data_df)
    p <- data_df$x
    q <- data_df$y

    if (nrow(data_df) < min_samples || sum(p) == 0 || sum(q) == 0) {
      return(NULL)
    }

    if (sd(p) == 0 || sd(q) == 0) {
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
    pct1 <- round(sum(p > 0) / length(p), 3)
    pct2 <- round(sum(q > 0) / length(q), 3)

    res_cor <- tryCatch(cor.test(p, q, method = cor_method), error = function(e) NULL)
    if (is.null(res_cor)) return(NULL)

    lr_name <- paste0(row.names(avg.s)[i], "_", row.names(avg.r)[i])

    return(c(
      round(res_cor$estimate, 5),
      round(res_cor$p.value, 15),
      pct1, pct2,
      lr_name,
      slope, intercept
    ))
  }

  res_list <- run_parallel(1:nrow(avg.r), calc_correlation, mc_cores = mc_cores)
  res_mat <- do.call(rbind, res_list)
  if (is.null(res_mat)) {
    message("No valid ligand-receptor pairs passed the initial filtering.")
    return(NULL)
  }
  res <- data.frame(res_mat, stringsAsFactors = FALSE)

  colnames(res) <- c(cor_colname, p_colname, "pct1", "pct2", "lr", "slope", "intercept")
  num_cols <- c(cor_colname, p_colname, "pct1", "pct2", "slope", "intercept")
  res[num_cols] <- lapply(res[num_cols], as.numeric)

  res$adjust.p <- round(p.adjust(res[[p_colname]], method = adjust_method), 15)

  # Step 6: Filter the results based on adjusted p-value, correlation, and percentage thresholds
  message("Filtering results based on adjusted p-value, correlation, and percentage thresholds...")
  res <- res[which(res$adjust.p < min_adjust_p &
                     res[[cor_colname]] > min_cor &
                     res$pct1 > min_pct &
                     res$pct2 > min_pct), ]
  res <- res[order(res$adjust.p, -res[[cor_colname]]), ]
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

  message("Filter and correlation process complete.")
  message("Head of results (", nrow(res), "):")
  print(head(res))

  return(res)
}



#' Filter and Analyze Ligand-Receptor Pair Correlations (All Cell Types)
#'
#' @description
#' Filters ligand-receptor (LR) pairs and analyzes their correlations for all possible cell type pairs,
#' and returning significant LR pairs based on user-defined thresholds.
#'
#' @param rna A Seurat object containing single-cell RNA expression data.
#' @param lr_database A data frame of ligand-receptor pairs with columns "ligand_gene_symbol" and "receptor_gene_symbol".
#' @param sample_col Column name in Seurat metadata indicating sample identifiers (character).
#' @param cell_type_col Column name in Seurat metadata indicating cell type classifications (character).
#' @param min_cells Minimum cells required per sample for both sender and receiver (default 50).
#' @param min_samples Minimum valid samples required to proceed (default 10).
#' @param cor_method Correlation method: "spearman" (default), "pearson", or "kendall".
#' @param adjust_method P-value adjustment method (default "BH" for Benjamini-Hochberg).
#'        Options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param min_adjust_p Adjusted p-value threshold for significance (default 0.05).
#' @param min_cor Minimum correlation coefficient threshold (default 0). Must be ≥ 0.
#' @param min_pct Minimum percentage of non-zero expression in both ligand and receptor (default 0.1).
#' @param mc_cores Number of CPU cores for parallel processing (default 10). Automatically capped at (system cores - 1).
#'
#' @return A data frame of filtered LR pairs with columns:
#'   \item{\code{cor_<method>}}{Correlation coefficient (e.g., \code{cor_spearman}).}
#'   \item{\code{p_<method>}}{Raw p-value (e.g., \code{p_spearman}).}
#'   \item{pct1, pct2}{Percent non-zero expression in sender and receiver.}
#'   \item{lr}{Ligand-receptor pair identifier (\code{"Ligand_Receptor"}).}
#'   \item{adjust.p}{Adjusted p-value.}
#'   \item{sender, receiver}{Sending/receiving cell types.}
#'   Rows ordered by ascending \code{adjust.p} and descending \code{cor_<method>}.
#'   Returns \code{NULL} if no pairs pass filters.
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom utils combn head
#'
#' @examples
#' seurat_object <- load_example_seurat()
#' data(lr_db)
#'
#' # Analyzing ligand-receptor interactions between all cell types
#' result01a <- filter_lr_all(
#'   rna = seurat_object,
#'   lr_database = lr_db,
#'   sample_col = "sample",
#'   cell_type_col = "cell.type",
#'   min_cells = 20,
#'   min_samples = 10,
#'   min_adjust_p = 0.5,
#'   mc_cores = 1
#' )
filter_lr_all <- function(rna, lr_database = PopComm::lr_db,
                          sample_col, cell_type_col,
                          min_cells = 50, min_samples = 10,
                          cor_method = "spearman", adjust_method = "BH",
                          min_adjust_p = 0.05, min_cor = 0, min_pct = 0.1,
                          mc_cores = 10) {

  # Check parameters
  max_cores <- parallel::detectCores()
  if (mc_cores > max_cores) {
    message("Warning: Using more cores (", mc_cores, ") than available (", max_cores, ").")
    message("Using ", max_cores - 1, " cores instead of requested ", mc_cores)
    mc_cores <- max_cores - 1
  }

  # Pre-process metadata
  rna$sample <- rna@meta.data[, sample_col]
  rna$cell.type <- rna@meta.data[, cell_type_col]
  cell_types <- unique(rna@meta.data[[cell_type_col]])
  if (length(cell_types) < 1) stop("No cell types found.")

  message("Analyzing ligand-receptor interactions between all cell types...")
  message("Cell types: ", paste(cell_types, collapse = ", "))

  run_filter <- function(rna_obj, sender, receiver) {
    res <- filter_lr_single(
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
      min_pct = min_pct,
      mc_cores = mc_cores
    )
    if (!is.null(res) && nrow(res) > 0) {
      return(res)
    }
    return(NULL)
  }

  # Step 1: Filter for ligand-receptor interactions where the same cell type is both sender and receiver
  message("\nProcessing same cell type pairs...")
  same_ct_results <- lapply(cell_types, function(ct) {
    message("\n  Analyzing pair: ", ct, " <-> ", ct)
    run_filter(rna, sender = ct, receiver = ct)
  })
  same_ct_results <- Filter(Negate(is.null), same_ct_results)

  # Step 2: Filter for ligand-receptor interactions where sender and receiver are different cell types
  message("\n\n\nProcessing different cell type pairs...")
  diff_ct_results <- list()
  unique_pairs <- combn(cell_types, 2, simplify = FALSE)

  for (pair in unique_pairs) {
    ct1 <- pair[1]
    ct2 <- pair[2]
    message("\n  Analyzing pair: ", ct1, " <-> ", ct2)

    # subset data
    message("Subsetting data for selected cell types: ", ct1, " and ", ct2)
    cells_to_keep <- rna$cell.type %in% c(ct1, ct2)
    rna_subset <- rna[, cells_to_keep]
    if (length(unique(rna_subset$cell.type)) < 2) {
      message("  Skipping pair ", ct1, " and ", ct2, " due to insufficient cell types.")
      next
    }

    # ct1 -> ct2
    res_forward <- run_filter(rna_subset, sender = ct1, receiver = ct2)
    # ct2 -> ct1
    res_reverse <- run_filter(rna_subset, sender = ct2, receiver = ct1)

    if (!is.null(res_forward)) {
      diff_ct_results <- c(diff_ct_results, list(res_forward))
    }
    if (!is.null(res_reverse)) {
      diff_ct_results <- c(diff_ct_results, list(res_reverse))
    }
  }

  final_res <- dplyr::bind_rows(c(same_ct_results, diff_ct_results))

  # Check if the filtered result is empty
  if (nrow(final_res) == 0) {
    message("\n\nNo results meet the filtering criteria. Returning NULL.")
    return(NULL)
  }

  message("\n\nAll cell types filter and correlation process complete.")
  message("Head of final results (", nrow(final_res), "):")
  print(head(final_res))

  return(final_res)
}
