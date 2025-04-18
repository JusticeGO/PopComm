#' Analyze Ligand-Receptor Projection Scores (Specified Sender and Receiver)
#'
#' @description
#' This function calculates the projection scores for ligand-receptor (LR) pairs
#' between specified sender and receiver cell types. The projection score is computed
#' based on linear regression models, measuring the normalized distance of each sample's
#' LR expression from the origin of the regression line.
#'
#' @param rna A Seurat object containing single-cell RNA expression data.
#' @param sender Cell type designated as the ligand sender (character).
#' @param receiver Cell type designated as the receptor receiver (character).
#' @param filtered_lr A data frame of filtered ligand-receptor pairs from prior analysis (e.g., output of `filter_lr_single`).
#'                  Must contain an "lr" column with pair identifiers in "Ligand_Receptor" format.
#' @param sample_col Column name in Seurat metadata indicating sample identifiers (character).
#' @param cell_type_col Column name in Seurat metadata indicating cell type classifications (character).
#' @param min_cells Minimum cells required per sample for both sender and receiver (numeric, default 50).
#' @param num_cores Number of CPU cores for parallel processing (numeric, default 10). Automatically capped at (system cores - 1).
#'
#' @return A data frame with projection scores per sample and LR pair. Columns:
#'   \item{All input from \code{filtered_lr}}{Original columns provided by the user in \code{filtered_lr}.}
#'   \item{sample}{Sample identifier.}
#'   \item{score}{Projection score (raw co-expression intensity).}
#'   \item{normalized_score}{Normalized score scaled between 0-1.}
#'   Rows ordered by \code{filtered_lr} columns and descending \code{score}.
#'   Returns \code{NULL} if no valid pairs.
#'
#' @export
#'
#' @importFrom dplyr %>% bind_rows
#' @importFrom stats lm sd coef na.omit
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
#'   num_cores = 1
#' )
#'
#' # Analyzing ligand-receptor projection scores: Cardiac -> Perivascular
#' result02s <- score_lr_single(
#'   rna = seurat_object,
#'   sender = "Cardiac",
#'   receiver = "Perivascular",
#'   filtered_lr = result01s,
#'   sample_col = "sample",
#'   cell_type_col = "cell.type",
#'   min_cells = 20,
#'   num_cores = 1
#' )
score_lr_single <- function(rna, sender, receiver, filtered_lr,
                            sample_col, cell_type_col,
                            min_cells = 50, num_cores = 10) {

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

  # Step 1: Add sample and cell type columns to the RNA data object and subset
  message("Analyzing ligand-receptor projection scores: ", sender, " -> ", receiver)

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
  if (length(valid_samples) < 2) {
    message("Insufficient valid samples. Analysis stopped.")
    return(NULL)
  }

  rna.data <- subset(rna.data, sample %in% valid_samples)

  # Step 3: Load the ligand-receptor pairs after filtering for interactions
  lr <- filtered_lr

  # Step 4: Compute average expression for each sample-cell type group
  rna.data$group <- paste0(rna.data$sample, "-lr-", rna.data$cell.type)
  message("Computing average expression for each sample-cell type group...")
  rna.avg <- Seurat::AverageExpression(rna.data, group.by = "group")[[1]]   # seurat v4/v5
  rna.avg <- round(rna.avg, 5)

  avg.s <- rna.avg[, grep(sender, colnames(rna.avg))]
  avg.r <- rna.avg[, grep(receiver, colnames(rna.avg))]

  colnames(avg.s) <- stringr::str_match(colnames(avg.s), "^(.*)-lr-")[,2]
  colnames(avg.r) <- stringr::str_match(colnames(avg.r), "^(.*)-lr-")[,2]

  avg.r <- avg.r[, colnames(avg.s), drop = FALSE]

  avg.s <- avg.s[lr$ligand, , drop = FALSE]
  avg.r <- avg.r[lr$receptor, , drop = FALSE]

  # Step 5: Calculating projection scores
  message("Calculating projection scores...")

  calc_projection <- function(i) {
    x <- as.numeric(avg.s[lr$ligand[i], ])
    y <- as.numeric(avg.r[lr$receptor[i], ])

    if (length(x) == 0 || length(y) == 0 || all(is.na(x)) || all(is.na(y)) || sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
      return(data.frame())
    }

    slope <- lr$slope[i]
    intercept <- lr$intercept[i]

    projections <- t(sapply(
      1:length(x), function(j) {
        project_to_line(x[j], y[j], slope, intercept)
      }))

    dx <- projections[, 1] - min(projections[, 1])
    dy <- projections[, 2] - min(projections[, 2])
    score <- round(sqrt(dx^2 + dy^2), 5)
    normalized_score <- round(score / max(score), 5)

    lr_metadata <- lr[i, ]
    result <- data.frame(
      lr_metadata,
      sample = colnames(avg.s),
      score = score,
      normalized_score = normalized_score,
      row.names = NULL
    )

    result <- result[order(-result$score), ]
    return(result)
  }

  score_list <- run_parallel(
    1:nrow(avg.s),
    calc_projection,
    num_cores = num_cores,
    export_vars = c("avg.s", "avg.r", "lr", "project_to_line")
    )

  # Combine the results into a single data frame and remove NAs
  score_list <- score_list[sapply(score_list, function(x) nrow(x) > 0)]
  score.df <- bind_rows(score_list) %>% na.omit()

  message("Analyzing ligand-receptor projection scores process complete.")
  message("Head of results (", nrow(score.df), "):")
  print(head(score.df))

  return(score.df)
}




#' Analyze Ligand-Receptor Projection Scores (Across All Cell Types)
#'
#' @description
#' This function calculates the ligand-receptor (LR) projection scores between all combinations
#' of sender and receiver cell types. The projection score is computed based on linear regression models,
#' measuring the normalized distance of each sample's LR expression from the origin of the regression line.
#'
#' @param rna A Seurat object containing single-cell RNA expression data.
#' @param filtered_lr A data frame of ligand-receptor pairs from prior analysis (e.g., output of `filter_lr_single`).
#'                  Must contain an "lr" column with pair identifiers in "Ligand_Receptor" format.
#' @param sample_col Column name in Seurat metadata indicating sample identifiers (character).
#' @param cell_type_col Column name in Seurat metadata indicating cell type classifications (character).
#' @param min_cells Minimum cells required per sample for both sender and receiver (numeric, default 50).
#' @param num_cores Number of CPU cores for parallel processing (numeric, default 10). Automatically capped at (system cores - 1).
#'
#' @return A data frame with projection scores per sample and LR pair. Columns:
#'   \item{All input from \code{filtered_lr}}{Original columns provided by the user in \code{filtered_lr}.}
#'   \item{sample}{Sample identifier.}
#'   \item{score}{Projection score (raw co-expression intensity).}
#'   \item{normalized_score}{Normalized score scaled between 0-1.}
#'   Rows ordered by \code{filtered_lr} columns and descending \code{score}.
#'   Returns \code{NULL} if no valid pairs.
#'
#' @export
#'
#' @importFrom dplyr %>% bind_rows
#' @importFrom utils head
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
#'   num_cores = 1
#' )
#'
#' # Analyzing ligand-receptor projection scores between all cell types
#' result02a <- score_lr_all(
#'   rna = seurat_object,
#'   filtered_lr = result01a,
#'   sample_col = "sample",
#'   cell_type_col = "cell.type",
#'   min_cells = 20,
#'   num_cores = 1
#' )
score_lr_all <- function(rna, filtered_lr,
                         sample_col, cell_type_col,
                         min_cells = 50, num_cores = 10) {

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

  message("Analyzing ligand-receptor projection scores between all cell types...")
  message("Cell types: ", paste(cell_types, collapse = ", "), "\n")

  all_results <- list()

  # split filtered_lr data
  split_data <- list()
  for (i in 1:nrow(filtered_lr)) {
    sender <- filtered_lr$sender[i]
    receiver <- filtered_lr$receiver[i]
    combination_name <- paste(sender, receiver, sep = "_")
    if (!combination_name %in% names(split_data)) {
      split_data[[combination_name]] <- filtered_lr[filtered_lr$sender == sender & filtered_lr$receiver == receiver, ]
    }
  }

  # Calculate the projection scores for ligand-receptor (LR) pairs
  res_list <- lapply(split_data, function(lr_s) {
    sender <- lr_s$sender[1]
    receiver <- lr_s$receiver[1]
    res <- score_lr_single(
      rna = rna,
      sender = sender,
      receiver = receiver,
      filtered_lr = lr_s,
      sample_col = "sample",
      cell_type_col = "cell.type",
      min_cells = min_cells,
      num_cores = num_cores
    )

    if (!is.null(res) && nrow(res) > 0) {
      return(res)
    } else {
      return(NULL)
    }
  })
  res_list <- Filter(Negate(is.null), res_list)
  final_res <- bind_rows(res_list)

  message("\n\nAll cell types analyzing ligand-receptor projection scores process complete.")
  message("Head of results (", nrow(final_res), "):")
  print(head(final_res))

  return(final_res)
}
