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
#' @importFrom stats cor.test p.adjust
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
  message("\nAnalyzing ligand-receptor interactions: ", sender, " -> ", receiver)

  # Determine the subset of data
  if (!setequal(selected_types, cell_types)) {
    message("Subsetting ", sender, " (sender) and ", receiver, " (receiver)...")
    cells.to.keep <- rna$cell.type %in% selected_types
    rna.data <- rna[, cells.to.keep]
  } else {
    rna.data <- rna
  }

  # Step 2: Filter samples based on thresholds for the number of cells and samples
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
  lr <- lr[which(lr$ligand_gene_symbol %in% rownames(rna.data)),]
  lr <- lr[which(lr$receptor_gene_symbol %in% rownames(rna.data)),]

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

  # Check the operating system and use appropriate parallel function
  if (Sys.info()["sysname"] == "Windows") {
    # Use parLapply on Windows
    cl <- parallel::makeCluster(mc_cores)
    res <- parallel::parLapply(cl, 1:nrow(avg.r), function(i) {
      data <- data.frame(x = avg.s[i,], y = avg.r[i,])
      data <- remove_outlier(data)

      p <- data$x
      q <- data$y

      if (dim(data)[1] < min_samples || sum(p) == 0 || sum(q) == 0) {
        return(NULL)
      }

      pct1 <- sum(p > 0) / length(p)
      pct2 <- sum(q > 0) / length(q)

      res.cor <- cor.test(p, q, method = cor_method)

      lr <- paste0(row.names(avg.s)[i], "_", row.names(avg.r)[i])
      return(c(round(res.cor$estimate, 5), round(res.cor$p.value, 15), round(pct1, 3), round(pct2, 3), lr))
    })
    parallel::stopCluster(cl)
  } else {
    # Use pbmclapply for Linux systems
    res <- pbmcapply::pbmclapply(1:nrow(avg.r), function(i) {
      data <- data.frame(x = avg.s[i,], y = avg.r[i,])
      data <- remove_outlier(data)

      p <- data$x
      q <- data$y

      if (dim(data)[1] < min_samples || sum(p) == 0 || sum(q) == 0) {
        return(NULL)
      }

      pct1 <- sum(p > 0) / length(p)
      pct2 <- sum(q > 0) / length(q)

      res.cor <- cor.test(p, q, method = cor_method)

      lr <- paste0(row.names(avg.s)[i], "_", row.names(avg.r)[i])
      return(c(round(res.cor$estimate, 5), round(res.cor$p.value, 15), round(pct1, 3), round(pct2, 3), lr))
    }, mc.cores = mc_cores)
  }

  res <- do.call(rbind, res)
  res <- data.frame(res)

  colnames(res) <- c(cor_colname, p_colname, "pct1", "pct2", "lr")
  res[[cor_colname]] <- as.numeric(res[[cor_colname]])
  res[[p_colname]] <- as.numeric(res[[p_colname]])
  res$pct1 <- as.numeric(res$pct1)
  res$pct2 <- as.numeric(res$pct2)

  res$adjust.p <- round(p.adjust(res[[p_colname]], method = adjust_method), 15)

  # Step 6: Filter the results based on adjusted p-value, correlation, and percentage thresholds
  message("Filtering results based on adjusted p-value, correlation, and percentage thresholds...")
  res <- res[which(res$adjust.p < min_adjust_p &
                     res[[cor_colname]] > min_cor &
                     res$pct1 > min_pct &
                     res$pct2 > min_pct),]
  res <- res[order(res$adjust.p, -res[[cor_colname]]),]
  if (nrow(res) > 0) {
    row.names(res) <- 1:nrow(res)
  }

  # Check if the filtered result is empty
  if (nrow(res) == 0) {
    message("No results meet the filtering criteria. Returning NULL.")
    return(NULL)
  }

  res$sender <- sender
  res$receiver <- receiver

  # Final message and print a summary of results
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
  message("\nCell types: ", paste(cell_types, collapse = ", "))

  all_results <- list()

  # Step 1: Filter for ligand-receptor interactions where the same cell type is both sender and receiver
  message("\nProcessing same cell type pairs...")
  same_ct_res <- lapply(cell_types, function(ct) {
    message("\n\n  Analyzing pair: ", ct, " <-> ", ct)
    res <- filter_lr_single(
      rna = rna,
      sender = ct,
      receiver = ct,
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
    } else {
      return(NULL)
    }
  })
  same_ct_res <- Filter(Negate(is.null), same_ct_res)

  # Step 2: Filter for ligand-receptor interactions where sender and receiver are different cell types
  message("\n\n\nProcessing different cell type pairs...")
  unique_pairs <- combn(cell_types, 2, simplify = FALSE)

  diff_ct_res <- lapply(unique_pairs, function(pair) {
    ct1 <- pair[1]
    ct2 <- pair[2]
    message("\n\n  Analyzing pair: ", ct1, " <-> ", ct2)

    # subset data
    message("Subsetting data for selected cell types: ", ct1, " and ", ct2)
    cells.to.keep <- rna$cell.type %in% c(ct1, ct2)
    rna_subset <- rna[, cells.to.keep]
    if (length(unique(rna_subset$cell.type)) < 2) return(NULL)

    # ct1 -> ct2
    res1 <- filter_lr_single(
      rna = rna_subset,
      sender = ct1,
      receiver = ct2,
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

    # ct2 -> ct1
    res2 <- filter_lr_single(
      rna = rna_subset,
      sender = ct2,
      receiver = ct1,
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

    filtered_list <- list(res1, res2) %>%
      Filter(f = function(x) !is.null(x) && nrow(x) > 0)

    if (length(filtered_list) > 0) {
      return(filtered_list)
    } else {
      return(NULL)
    }
  })
  diff_ct_res <- unlist(diff_ct_res, recursive = FALSE)

  final_res <- dplyr::bind_rows(c(same_ct_res, diff_ct_res))

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
