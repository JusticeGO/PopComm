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
#' @param lr_custom A data frame of ligand-receptor pairs from prior analysis (e.g., output of `filter_lr_single`).
#'                  Must contain an "lr" column with pair identifiers in "Ligand_Receptor" format.
#' @param sample_col Column name in Seurat metadata indicating sample identifiers (character).
#' @param cell_type_col Column name in Seurat metadata indicating cell type classifications (character).
#' @param min_cells Minimum cells required per sample for both sender and receiver (default 50).
#' @param min_samples Minimum valid samples required to proceed (default 10).
#' @param mc_cores Number of CPU cores for parallel processing (default 10). Automatically capped at (system cores - 1).
#'
#' @return A data frame with projection scores per sample and LR pair. Columns:
#'   \item{All input from \code{lr_custom}}{Original columns provided by the user in \code{lr_custom}.}
#'   \item{\code{sample}}{Sample identifier.}
#'   \item{\code{score}}{Projection score (raw co-expression intensity).}
#'   \item{\code{normalized_score}}{Normalized score scaled between 0-1.}
#'   Rows ordered by \code{lr_custom} columns and descending \code{score}.
#'   Returns \code{NULL} if no valid pairs.
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom stats lm sd coef na.omit
#' @importFrom utils head
#'
#' @examples
#' seurat_object <- load_example_seurat()
#' data(lr_db)
#'
#' # Analyzing ligand-receptor interactions: Cardiac -> Perivascular
#' result01s <- filter_lr_single_v2(
#'   rna = seurat_object,
#'   sender = "Cardiac",
#'   receiver = "Perivascular",
#'   lr_database = lr_db,
#'   sample_col = "sample",
#'   cell_type_col = "cell.type",
#'   min_cells = 20,
#'   min_samples = 10,
#'   min_adjust_p = 0.05,
#'   mc_cores = 1
#' )
#'
#' # Analyzing ligand-receptor projection scores: Cardiac -> Perivascular
#' result02s <- score_lr_single_v2(
#'   rna = seurat_object,
#'   sender = "Cardiac",
#'   receiver = "Perivascular",
#'   lr_custom = result01s,
#'   sample_col = "sample",
#'   cell_type_col = "cell.type",
#'   min_cells = 20,
#'   min_samples = 10,
#'   mc_cores = 1
#' )
score_lr_single_v2 <- function(rna, sender, receiver, lr_custom,
                               sample_col, cell_type_col,
                               min_cells = 50, min_samples = 10,
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
  message("\nAnalyzing ligand-receptor projection scores: ", sender, " -> ", receiver)

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

  # Step 3: Load the ligand-receptor pairs after filtering for interactions
  lr <- lr_custom
  lr$ligand <- stringr::str_match(lr$lr, "^(.*)_")[,2]
  lr$receptor <- stringr::str_match(lr$lr, "_(.*)$")[,2]

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

  avg.s <- avg.s[lr$ligand, , drop = FALSE]
  avg.r <- avg.r[lr$receptor, , drop = FALSE]

  # Step 5: Calculating projection scores
  message("Calculating projection scores...")

  # Check the operating system and use appropriate parallel function
  if (Sys.info()["sysname"] == "Windows") {
    # Use parLapply on Windows
    cl <- parallel::makeCluster(mc_cores)
    score.df <- parallel::parLapply(cl, 1:nrow(avg.s), function(i) {
      x <- as.numeric(avg.s[lr$ligand[i], ])
      y <- as.numeric(avg.r[lr$receptor[i], ])

      if (sd(x) == 0 || sd(y) == 0) {
        return(data.frame())
      }
#
#       model <- tryCatch(
#         lm(y ~ x),
#         error = function(e) NULL
#       )
#
#       if (is.null(model)) return(data.frame())
#
#       slope <- round(coef(model)[2], 5)
#       intercept <- round(coef(model)[1], 5)
      slope <- lr$slope[i]
      intercept <- lr$intercept[i]

      projections <- t(sapply(
        1:length(x),
        function(j) project_to_line(x[j], y[j], slope, intercept)
      ))

      dx <- projections[, 1] - min(projections[, 1])
      dy <- projections[, 2] - min(projections[, 2])
      score <- round(sqrt(dx^2 + dy^2), 5)
      normalized.score <- score / max(score)

      lr_metadata <- lr[i, ]
      result <- data.frame(
        lr_metadata,
        sample = colnames(avg.s),
        score = score,
        normalized_score = round(normalized.score, 5),
        row.names = NULL
      )

      result <- result[order(-result$score), ]

      return(result)
    })
    parallel::stopCluster(cl)
  } else {
    # Use pbmclapply for Linux systems
    score.df <- pbmcapply::pbmclapply(1:nrow(avg.s), function(i) {
      x <- as.numeric(avg.s[lr$ligand[i], ])
      y <- as.numeric(avg.r[lr$receptor[i], ])

      if (sd(x) == 0 || sd(y) == 0) {
        return(data.frame())
      }

      # model <- tryCatch(
      #   lm(y ~ x),
      #   error = function(e) NULL
      # )
      #
      # if (is.null(model)) return(data.frame())
      #
      # slope <- round(coef(model)[2], 5)
      # intercept <- round(coef(model)[1], 5)
      slope <- lr$slope[i]
      intercept <- lr$intercept[i]

      projections <- t(sapply(
        1:length(x),
        function(j) project_to_line(x[j], y[j], slope, intercept)
      ))

      dx <- projections[, 1] - min(projections[, 1])
      dy <- projections[, 2] - min(projections[, 2])
      score <- round(sqrt(dx^2 + dy^2), 5)
      normalized.score <- score / max(score)

      lr_metadata <- lr[i, ]
      result <- data.frame(
        lr_metadata,
        sample = colnames(avg.s),
        score = score,
        normalized_score = round(normalized.score, 5),
        row.names = NULL
      )

      result <- result[order(-result$score), ]

      return(result)
    }, mc.cores = mc_cores)
  }

  # Combine the results into a single data frame and remove NAs
  score.df <- score.df %>%
    dplyr::bind_rows() %>%
    na.omit()

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
#' @param lr_custom A data frame of ligand-receptor pairs from prior analysis (e.g., output of `filter_lr_single`).
#'                  Must contain an "lr" column with pair identifiers in "Ligand_Receptor" format.
#' @param sample_col Column name in Seurat metadata indicating sample identifiers (character).
#' @param cell_type_col Column name in Seurat metadata indicating cell type classifications (character).
#' @param min_cells Minimum cells required per sample for both sender and receiver (default 50).
#' @param min_samples Minimum valid samples required to proceed (default 10).
#' @param mc_cores Number of CPU cores for parallel processing (default 10). Automatically capped at (system cores - 1).
#'
#' @return A data frame containing projection scores for each sample and ligand-receptor (LR) pair, with the following columns:
#'   \item{All columns from the input `lr_custom`}{}
#'   \item{sample}{Sample identifier.}
#'   \item{score}{Projection score indicating the co-expression intensity of the LR pair in the sample.}
#'   \item{normalized_score}{Normalized projection score (range 0-1) indicating the relative co-expression intensity.}
#'   Rows are first ordered by `lr_custom` and then by the `score` column in descending order.
#'   Returns NULL if no valid pairs are found.
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom utils head
#'
#' @examples
#' seurat_object <- load_example_seurat()
#' data(lr_db)
#'
#' #' # Analyzing ligand-receptor interactions between all cell types
#' result01a <- filter_lr_all_v2(
#'   rna = seurat_object,
#'   lr_database = lr_db,
#'   sample_col = "sample",
#'   cell_type_col = "cell.type",
#'   min_cells = 20,
#'   min_samples = 10,
#'   min_adjust_p = 0.5,
#'   mc_cores = 1
#' )
#'
#' # Analyzing ligand-receptor projection scores between all cell types
#' result02a <- score_lr_all_v2(
#'   rna = seurat_object,
#'   lr_custom = result01a,
#'   sample_col = "sample",
#'   cell_type_col = "cell.type",
#'   min_cells = 20,
#'   min_samples = 10,
#'   mc_cores = 1
#' )
score_lr_all_v2 <- function(rna, lr_custom,
                            sample_col, cell_type_col,
                            min_cells = 50, min_samples = 10,
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

  message("Analyzing ligand-receptor projection scores between all cell types...")
  message("\nCell types: ", paste(cell_types, collapse = ", "))

  all_results <- list()

  # split lr_custom data
  split_data <- list()
  for (i in 1:nrow(lr_custom)) {
    sender <- lr_custom$sender[i]
    receiver <- lr_custom$receiver[i]
    combination_name <- paste(sender, receiver, sep = "_")
    if (!combination_name %in% names(split_data)) {
      split_data[[combination_name]] <- lr_custom[lr_custom$sender == sender & lr_custom$receiver == receiver, ]
    }
  }

  # Calculate the projection scores for ligand-receptor (LR) pairs
  res_list <- lapply(split_data, function(lr_s) {
    sender <- lr_s$sender[1]
    receiver <- lr_s$receiver[1]
    res <- score_lr_single_v2(
      rna = rna,
      sender = sender,
      receiver = receiver,
      lr_custom = lr_s,
      sample_col = "sample",
      cell_type_col = "cell.type",
      min_cells = min_cells,
      min_samples = min_samples,
      mc_cores = mc_cores
    )

    if (!is.null(res) && nrow(res) > 0) {
      return(res)
    } else {
      return(NULL)
    }
  })
  res_list <- Filter(Negate(is.null), res_list)
  final_res <- dplyr::bind_rows(res_list)

  message("\n\nAll cell types analyzing ligand-receptor projection scores process complete.")
  message("Head of results (", nrow(final_res), "):")
  print(head(final_res))

  return(final_res)
}
