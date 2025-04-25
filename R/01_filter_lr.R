#' Filter and Analyze Ligand-Receptor Pair Correlations (Specified Sender and Receiver)
#'
#' @description
#' Filters ligand-receptor (LR) pairs and analyzes their correlations for specified sender and receiver cell types,
#' and returns significant LR pairs based on user-defined thresholds.
#'
#' @param rna A Seurat object containing single-cell RNA expression data.
#' @param sender Cell type designated as the ligand sender (character).
#' @param receiver Cell type designated as the receptor receiver (character).
#' @param lr_database A data frame of ligand-receptor pairs with columns "ligand_gene_symbol" and "receptor_gene_symbol".
#' @param sample_col Column name in Seurat metadata indicating sample identifiers (character).
#' @param cell_type_col Column name in Seurat metadata indicating cell type classifications (character).
#' @param min_cells Minimum cells required per sample for both sender and receiver (numeric, default 50).
#' @param min_samples Minimum valid samples required to proceed (numeric, default 10).
#' @param min_cell_ratio Minimum ratio of cells expressing ligand and receptor genes in sender or receiver cells (numeric, default 0.1).
#' @param min_sample_ratio Minimum ratio of samples in which both the ligand and receptor genes must be expressed (numeric, default 0.1).
#' @param cor_method Correlation method: "spearman" (default), "pearson", or "kendall".
#' @param adjust_method P-value adjustment method (default "BH" for Benjamini-Hochberg).
#'        Options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param min_adjust_p Adjusted p-value threshold for significance (numeric, default 0.05).
#' @param min_cor Minimum correlation coefficient threshold (numeric, default 0). Must be \eqn{\ge}{>=} 0.
#' @param num_cores Number of CPU cores for parallel processing (numeric, default 10). Automatically capped at (system cores - 1).
#' @param verbose Logical indicating whether to print progress messages (logical, default: TRUE).
#'
#' @return A data frame includes LR pairs with sufficient correlation and expression support across samples.
#'   \item{ligand, receptor}{Ligand and receptor gene symbols.}
#'   \item{cor}{Correlation coefficient.}
#'   \item{p_val}{Raw p-value.}
#'   \item{adjust.p}{Adjusted p-value.}
#'   \item{sender, receiver}{Sender and receiver cell types.}
#'   \item{slope}{Slope of the linear regression model.}
#'   \item{intercept}{Intercept of the linear regression model.}
#' Rows are ordered by ascending \code{adjust.p} and descending \code{cor}.
#'
#' Returns \code{NULL} if:
#' \itemize{
#'   \item No cell types are found in the metadata.
#'   \item Insufficient samples or cells remain after filtering.
#'   \item No ligand-receptor pairs pass the filtering thresholds.
#' }
#'
#' @export
#'
#' @importFrom dplyr %>% bind_rows select
#' @importFrom stats lm sd coef cor.test p.adjust
#' @importFrom utils head packageVersion
#' @importFrom Matrix rowSums
#' @importFrom stringr str_match
#'
#' @examples
#' \donttest{
#'   # Long-running example (may take >10s)
#'   seurat_object <- load_example_seurat()
#'   data(lr_db)
#'
#'   # Analyzing ligand-receptor interactions: Cardiac -> Perivascular
#'   result01s <- filter_lr_single(
#'     rna = seurat_object,
#'     sender = "Cardiac",
#'     receiver = "Perivascular",
#'     lr_database = lr_db,
#'     sample_col = "sample",
#'     cell_type_col = "cell.type",
#'     min_cells = 20,
#'     min_samples = 10,
#'     min_adjust_p = 0.5,
#'     num_cores = 1,
#'     verbose = TRUE
#'   )
#'
#'   if (!is.null(result01s)) {
#'   print(head(result01s))
#'   }
#' }
filter_lr_single <- function(rna, sender, receiver, lr_database = PopComm::lr_db,
                             sample_col, cell_type_col,
                             min_cells = 50, min_samples = 10,
                             min_cell_ratio = 0.1, min_sample_ratio = 0.1,
                             cor_method = "spearman", adjust_method = "BH",
                             min_adjust_p = 0.05, min_cor = 0,
                             num_cores = 10, verbose = TRUE) {

  # Check parameters
  max_cores <- parallel::detectCores()
  if (num_cores > max_cores) {
    warning(
      "Requested cores (", num_cores, ") exceed available (", max_cores, "). ",
      "Using ", max_cores - 1, " cores.",
      immediate. = TRUE
    )
    num_cores <- max_cores - 1
  }

  # Pre-process metadata
  rna$sample <- rna@meta.data[, sample_col]
  rna$cell.type <- rna@meta.data[, cell_type_col]
  cell_types <- unique(rna@meta.data[[cell_type_col]])
  if (length(cell_types) < 1) {
    stop("No cell types found in column '", cell_type_col, "'.", call. = FALSE)
  }

  selected_types <- unique(c(sender, receiver))
  missing_types <- setdiff(selected_types, cell_types)
  if (length(missing_types) > 0) {
    stop("Missing cell types: ", paste(missing_types, collapse = ", "))
  }

  # Step 1: Add sample and cell type columns to the RNA data object and subset
  if (verbose) {
    message("Analyzing ligand-receptor interactions: ", sender, " -> ", receiver)
  }

  # Determine the subset of data
  if (!setequal(selected_types, cell_types)) {
    if (verbose) {
      message("Subsetting ", sender, " (sender) and ", receiver, " (receiver)...")
    }
    cells.to.keep <- rna$cell.type %in% selected_types
    rna.data <- rna[, cells.to.keep]
  } else {
    rna.data <- rna
  }

  # Step 2: Filter samples based on cell counts per sample
  if (verbose) {
    message("Filtering samples with cell counts...")
  }
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

  if (verbose) {
    message("Remaining samples after filtering: ", length(valid_samples))
  }
  if (length(valid_samples) < min_samples) {
    warning(
      sender, " -> ", receiver, ": insufficient valid samples (", length(valid_samples), " < ", min_samples, ").\n",
      "Check: min_cells (current=", min_cells, ") or sample collection.",
      immediate. = TRUE
    )
    return(NULL)
  }
  rna.data <- subset(rna.data, sample %in% valid_samples)

  # Step 3: Filter LR pairs based on minimum cell ratio
  if (verbose) {
    message("Filtering LR pairs based on minimum cell ratio in sender or receiver cells...")
  }
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

  # Step 4: Compute average expression for each sample-cell type group
  rna.data$group <- paste0(rna.data$sample, "-lr-", rna.data$cell.type)
  if (verbose) {
    message("Computing average expression for each sample-cell type group...")
  }

  rna.avg <- Seurat::AverageExpression(rna.data, group.by = "group")[[1]]   # seurat v4/v5
  rna.avg <- round(rna.avg, 5)

  avg.s <- rna.avg[, grep(sender, colnames(rna.avg))]
  avg.r <- rna.avg[, grep(receiver, colnames(rna.avg))]

  colnames(avg.s) <- str_match(colnames(avg.s), "^(.*)-lr-")[,2]
  colnames(avg.r) <- str_match(colnames(avg.r), "^(.*)-lr-")[,2]

  avg.r <- avg.r[, colnames(avg.s), drop = FALSE]

  avg.s <- avg.s[lr$ligand_gene_symbol, , drop = FALSE]
  avg.r <- avg.r[lr$receptor_gene_symbol, , drop = FALSE]

  # Step 5: Compute correlations between ligand-receptor pairs
  if (verbose) {
    message("Starting correlation and filtering process for ligand-receptor pairs...")
  }

  # Calculate correlation and linear models
  calc_correlation <- function(i) {
    x <- as.numeric(avg.s[lr$ligand_gene_symbol[i], ])
    y <- as.numeric(avg.r[lr$receptor_gene_symbol[i], ])

    # filter sample
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
    warning(
      sender, " -> ", receiver, ": no ligand-receptor pairs survived initial filtering.",
      immediate. = TRUE
    )
    return(NULL)
  }
  res <- data.frame(res_mat, stringsAsFactors = FALSE)

  colnames(res) <- c("cor", "p_val", "lr", "pct1", "pct2", "slope", "intercept")
  num_cols <- c("cor", "p_val", "pct1", "pct2", "slope", "intercept")
  res[num_cols] <- lapply(res[num_cols], as.numeric)

  res$adjust.p <- round(p.adjust(res$p_val, method = adjust_method), 15)

  # Step 6: Filter the results based on adjusted p-value, correlation, and ratio thresholds
  if (verbose) {
    message("Filtering results based on adjusted p-value, correlation, and ratio thresholds...")
  }
  res <- res[which(res$adjust.p < min_adjust_p &
                     res$cor > min_cor &
                     res$pct1 > min_sample_ratio &
                     res$pct2 > min_sample_ratio), ]
  res <- res[order(res$adjust.p, -res$cor), ]

  if (nrow(res) == 0) {
    message(sender, " -> ", receiver, ": no significant LR pairs found.")
    return(NULL)
  }

  row.names(res) <- 1:nrow(res)
  res$sender <- sender
  res$receiver <- receiver
  res$ligand <- str_match(res$lr, "^(.*)_")[,2]
  res$receptor <- str_match(res$lr, "_(.*)$")[,2]

  selected_cols <- c("ligand", "receptor", "cor", "p_val", "adjust.p", "sender", "receiver", "slope", "intercept")
  res <- res %>% select(all_of(selected_cols))

  if (verbose) {
    message("Filtered LR analysis complete. Identified ", nrow(res), " significant ligand-receptor pairs.")
  }

  return(res)
}



#' Filter and Analyze Ligand-Receptor Pair Correlations (All Cell Types)
#'
#' @description
#' Filters ligand-receptor (LR) pairs and analyzes their correlations for all possible cell type pairs,
#' and returns significant LR pairs based on user-defined thresholds.
#'
#' @param rna A Seurat object containing single-cell RNA expression data.
#' @param lr_database A data frame of ligand-receptor pairs with columns "ligand_gene_symbol" and "receptor_gene_symbol".
#' @param sample_col Column name in Seurat metadata indicating sample identifiers (character).
#' @param cell_type_col Column name in Seurat metadata indicating cell type classifications (character).
#' @param min_cells Minimum cells required per sample for both sender and receiver (numeric, default 50).
#' @param min_samples Minimum valid samples required to proceed (numeric, default 10).
#' @param min_cell_ratio Minimum ratio of cells expressing ligand and receptor genes in sender or receiver cells (numeric, default 0.1).
#' @param min_sample_ratio Minimum ratio of samples in which both the ligand and receptor genes must be expressed (numeric, default 0.1).
#' @param cor_method Correlation method: "spearman" (default), "pearson", or "kendall".
#' @param adjust_method P-value adjustment method (default "BH" for Benjamini-Hochberg).
#'        Options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param min_adjust_p Adjusted p-value threshold for significance (numeric, default 0.05).
#' @param min_cor Minimum correlation coefficient threshold (numeric, default 0). Must be \eqn{\ge}{>=} 0.
#' @param num_cores Number of CPU cores for parallel processing (numeric, default 10). Automatically capped at (system cores - 1).
#' @param verbose Logical indicating whether to print progress messages (logical, default: TRUE).
#'
#' @return A data frame includes LR pairs with sufficient correlation and expression support across samples.
#'   \item{ligand, receptor}{Ligand and receptor gene symbols.}
#'   \item{cor}{Correlation coefficient.}
#'   \item{p_val}{Raw p-value.}
#'   \item{adjust.p}{Adjusted p-value.}
#'   \item{sender, receiver}{Sender and receiver cell types.}
#'   \item{slope}{Slope of the linear regression model.}
#'   \item{intercept}{Intercept of the linear regression model.}
#' Rows are ordered by ascending \code{adjust.p} and descending \code{cor}.
#'
#' Returns \code{NULL} if:
#' \itemize{
#'   \item No cell types are found in the metadata.
#'   \item Insufficient samples or cells remain after filtering.
#'   \item No ligand-receptor pairs pass the filtering thresholds.
#' }
#'
#' @export
#'
#' @importFrom dplyr %>% bind_rows
#' @importFrom utils combn head
#'
#' @examples
#' \donttest{
#'   # Long-running example (may take >10s)
#'   seurat_object <- load_example_seurat()
#'   data(lr_db)
#'
#'   # Analyzing ligand-receptor interactions between all cell types
#'   result01a <- filter_lr_all(
#'     rna = seurat_object,
#'     lr_database = lr_db,
#'     sample_col = "sample",
#'     cell_type_col = "cell.type",
#'     min_cells = 20,
#'     min_samples = 10,
#'     min_adjust_p = 0.5,
#'     num_cores = 1,
#'     verbose = TRUE
#'   )
#'
#'   if (!is.null(result01a)) {
#'   print(head(result01a))
#'   }
#' }
filter_lr_all <- function(rna, lr_database = PopComm::lr_db,
                          sample_col, cell_type_col,
                          min_cells = 50, min_samples = 10,
                          min_cell_ratio = 0.1, min_sample_ratio = 0.1,
                          cor_method = "spearman", adjust_method = "BH",
                          min_adjust_p = 0.05, min_cor = 0,
                          num_cores = 10, verbose = TRUE) {

  # Check parameters
  max_cores <- parallel::detectCores()
  if (num_cores > max_cores) {
    warning(
      "Requested cores (", num_cores, ") exceed available (", max_cores, "). ",
      "Using ", max_cores - 1, " cores.",
      immediate. = TRUE
    )
    num_cores <- max_cores - 1
  }

  # Pre-process metadata
  rna$sample <- rna@meta.data[, sample_col]
  rna$cell.type <- rna@meta.data[, cell_type_col]
  cell_types <- unique(rna@meta.data[[cell_type_col]])
  if (length(cell_types) < 1) {
    stop("No cell types found in column '", cell_type_col, "'.", call. = FALSE)
  }

  if (verbose) {
    message("Analyzing ligand-receptor interactions between all cell types...")
    message("Cell types: ", paste(cell_types, collapse = ", "))
  }

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
      min_cell_ratio = min_cell_ratio,
      min_sample_ratio = min_sample_ratio,
      cor_method = cor_method,
      adjust_method = adjust_method,
      min_adjust_p = min_adjust_p,
      min_cor = min_cor,
      num_cores = num_cores,
      verbose = verbose
    )
    if (!is.null(res) && nrow(res) > 0) {
      return(res)
    }
    return(NULL)
  }

  # Step 1: Filter for ligand-receptor interactions where the same cell type is both sender and receiver
  if (verbose) {
    message("\nProcessing same cell type pairs...")
  }
  same_ct_results <- lapply(cell_types, function(ct) {
    if (verbose) {
      message("\n  Analyzing pair: ", ct, " <-> ", ct)
    }
    run_filter(rna, sender = ct, receiver = ct)
  })
  same_ct_results <- Filter(Negate(is.null), same_ct_results)

  # Step 2: Filter for ligand-receptor interactions where sender and receiver are different cell types
  if (verbose) {
    message("\n\n\nProcessing different cell type pairs...")
  }
  diff_ct_results <- list()
  unique_pairs <- combn(cell_types, 2, simplify = FALSE)

  for (pair in unique_pairs) {
    ct1 <- pair[1]
    ct2 <- pair[2]
    if (verbose) {
      message("\n  Analyzing pair: ", ct1, " <-> ", ct2)
    }

    # subset data
    if (verbose) {
      message("Subsetting data for selected cell types: ", ct1, " and ", ct2)
    }
    cells_to_keep <- rna$cell.type %in% c(ct1, ct2)
    rna_subset <- rna[, cells_to_keep]
    if (length(unique(rna_subset$cell.type)) < 2) {
      if (verbose) {
        message("  Skipping pair ", ct1, " and ", ct2, " due to insufficient cell types.")
      }
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

  final_res <- bind_rows(c(same_ct_results, diff_ct_results))

  # Check if the filtered result is empty
  if (nrow(final_res) == 0) {
    message("\n\nNo significant ligand-receptor pairs were identified in the filtered LR analysis.")
    return(NULL)
  }

  if (verbose) {
    message("\n\nFiltered LR analysis across all cell types complete. Identified ", nrow(final_res), " significant ligand-receptor pairs.")
  }

  return(final_res)
}

