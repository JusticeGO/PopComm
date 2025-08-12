#' Analyze Ligand-Receptor Projection Scores (Specified Sender and Receiver)
#'
#' @description
#' This function calculates the projection scores for ligand-receptor (LR) pairs
#' between specified sender and receiver cell types, and it supports both Seurat objects and
#' average expression matrices (matrix of gene expression data with cell types and samples as column names).
#' The projection score is computed based on linear regression models,
#' measuring the normalized distance of each sample's LR expression from the origin of the regression line.
#'
#' @param rna A Seurat object or a matrix containing single-cell RNA expression data.
#' @param sender Cell type designated as the ligand sender (character).
#' @param receiver Cell type designated as the receptor receiver (character).
#' @param filtered_lr A data frame of filtered ligand-receptor pairs from prior analysis (e.g., output of `filter_lr_single`).
#'                  Must contain an "lr" column with pair identifiers in "Ligand_Receptor" format.
#' @param sample_col Metadata column name (character) for sample identifiers in Seurat mode; Matrix mode uses column index (numeric).
#' @param cell_type_col Metadata column name (character) for cell type in Seurat mode; Matrix mode uses column index (numeric).
#' @param id_sep Separator used in matrix column names to split sample and cell type (e.g., `--` for "Cardiac--sample1"). Only used in Matrix mode.
#' @param min_cells Minimum number of cells per sample for both sender and receiver (numeric, default 50). Only used in Seurat mode.
#' @param num_cores Number of CPU cores for parallel processing (numeric, default 10). Automatically capped at (system cores - 1).
#' @param verbose Logical indicating whether to print progress messages (logical, default: TRUE).
#'
#' @return A data frame with projection scores per sample and LR pair. Columns:
#'   \item{All input from \code{filtered_lr}}{Original columns provided by the user in \code{filtered_lr}.}
#'   \item{sample}{Sample identifier.}
#'   \item{score}{Projection score (raw co-expression intensity).}
#'   \item{normalized_score}{Normalized score scaled between 0-1.}
#' Rows are ordered by \code{filtered_lr} columns and descending \code{score}.
#'
#' Returns \code{NULL} if:
#' \itemize{
#'   \item No cell types are found in the metadata.
#'   \item One or both of the specified sender and receiver cell types are missing in the data.
#'   \item Fewer than two valid samples remain after filtering based on minimum cell number per sample.
#' }

#'
#' @export
#'
#' @importFrom dplyr %>% bind_rows
#' @importFrom stats lm sd coef na.omit
#' @importFrom utils head
#' @importFrom stringr str_match
#'
#' @examples
#' \donttest{
#'   data(matrix_object)
#'   data(lr_db)
#'
#'   # Analyzing ligand-receptor interactions: Perivascular -> Endothelial
#'   result01s <- filter_lr_single(
#'     rna = matrix_object,
#'     sender = "Perivascular",
#'     receiver = "Endothelial",
#'     lr_database = lr_db,
#'     sample_col = 2,
#'     cell_type_col = 1,
#'     id_sep = "--",
#'     min_samples = 10,
#'     min_sample_ratio = 0.1,
#'     min_adjust_p = 0.05,
#'     num_cores = 1,
#'     verbose = TRUE
#'     )
#'
#'   # Analyzing ligand-receptor projection scores: Perivascular -> Endothelial
#'   result02s <- score_lr_single(
#'     rna = matrix_object,
#'     sender = "Perivascular",
#'     receiver = "Endothelial",
#'     filtered_lr = result01s,
#'     sample_col = 2,
#'     cell_type_col = 1,
#'     id_sep = "--",
#'     num_cores = 1,
#'     verbose = TRUE
#'     )
#'
#'   if (!is.null(result02s)) {
#'   print(head(result02s))
#'   }
#' }
score_lr_single <- function(rna, sender, receiver, filtered_lr,
                            sample_col, cell_type_col,
                            id_sep,                            # v0.2.0.0
                            min_cells = 50, num_cores = 10, verbose = TRUE) {

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
  if (missing(sample_col) || missing(cell_type_col)) {
    stop("When using expression matrix input, both 'sample_col' and 'cell_type_col' must be specified as integers.")
  }


  # Check if input is Seurat or average expression matrix v0.2.0.0
  is_seurat <- inherits(rna, "Seurat")
  is_matrix <- is.matrix(rna)

  if (!is_seurat && !is_matrix) {
    stop("'rna' must be either a Seurat object or a matrix of average expression data.")
  }

  if (is_seurat) {
    input_type <- "Seurat"

    if (!is.character(sample_col)) {
      stop("Seurat mode: 'sample_col' must be a character string.")
    }
    if (!is.character(cell_type_col)) {
      stop("Seurat mode: 'cell_type_col' must be a character string.")
    }

    if (!sample_col %in% colnames(rna@meta.data)) {
      stop(paste("Column", sample_col, "not found in Seurat metadata"))
    }
    if (!cell_type_col %in% colnames(rna@meta.data)) {
      stop(paste("Column", cell_type_col, "not found in Seurat metadata"))
    }
    rna$sample <- rna@meta.data[, sample_col]
    rna$cell.type <- rna@meta.data[, cell_type_col]

  } else if (is_matrix && nrow(rna) > 1 && ncol(rna) > 1) {
    input_type <- "Matrix"

    if (!is.numeric(sample_col) || sample_col <= 0 || sample_col %% 1 != 0) {
      stop("Matrix mode: 'sample_col' must be a positive integer.")
    }
    if (!is.numeric(cell_type_col) || cell_type_col <= 0 || cell_type_col %% 1 != 0) {
      stop("Matrix mode: 'cell_type_col' must be a positive integer.")
    }
    colnames_split <- strsplit(colnames(rna), id_sep, fixed = TRUE)
    if (!all(lengths(colnames_split) >= max(sample_col, cell_type_col))) {
      stop("Column names do not contain enough fields with the specified separator.")
    }
  } else {
    stop("Input must be either a Seurat object or an average expression matrix.")
  }


  if (input_type == "Seurat") {
    # Handle cell type check
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
      message("(Seurat mode) Analyzing ligand-receptor projection scores: ", sender, " -> ", receiver)
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
    if (length(valid_samples) < 2) {
      warning(
        sender, " -> ", receiver, ": insufficient valid samples (", length(valid_samples), " < 2).\n",
        "Check: min_cells (current=", min_cells, ") or sample collection.",
        immediate. = TRUE
      )
      return(NULL)
    }
    rna.data <- subset(rna.data, sample %in% valid_samples)

    # Step 3: Load the ligand-receptor pairs after filtering for interactions
    lr <- filtered_lr

    # Step 4: Compute average expression for each sample-cell type group
    rna.data$group <- paste0(rna.data$sample, "-lr-", rna.data$cell.type)
    if (verbose) {
      message("Computing average expression for each sample-cell type group...")
    }
    rna.avg <- Seurat::AverageExpression(rna.data, group.by = "group")[[1]]   # seurat v4/v5


  } else if (input_type == "Matrix") {
    # Split the column names using the specified separator
    cell_types <- unique(sapply(colnames_split, `[`, cell_type_col))
    # samples <- unique(sapply(colnames_split, `[`, sample_col))

    if (length(cell_types) < 1) {
      stop("No cell types found in column '", cell_type_col, "'.", call. = FALSE)
    }

    # Ensure sender and receiver are valid cell types in the matrix
    selected_types <- unique(c(sender, receiver))
    missing_types <- setdiff(selected_types, cell_types)
    if (length(missing_types) > 0) {
      stop("Missing cell types: ", paste(missing_types, collapse = ", "))
    }

    # Step 1: Add sample and cell type columns to the RNA data object and subset
    if (verbose) {
      message("(Matrix mode) Analyzing ligand-receptor projection scores: ", sender, " -> ", receiver)
    }

    # Determine the subset of data
    cell_types_loc <- sapply(colnames_split, `[`, cell_type_col)
    samples_loc <- sapply(colnames_split, `[`, sample_col)

    if (!setequal(selected_types, cell_types)) {
      if (verbose) {
        message("Subsetting ", sender, " (sender) and ", receiver, " (receiver)...")
      }
      cells.to.keep <- cell_types_loc %in% selected_types
      rna.data <- rna[, cells.to.keep]
    } else {
      rna.data <- rna
    }

    # Step 2: Filter samples based on min_samples
    valid_samples <- intersect(
      unique(samples_loc[cell_types_loc == sender]),
      unique(samples_loc[cell_types_loc == receiver])
    )

    if (verbose) {
      message("Remaining samples after filtering: ", length(valid_samples))
    }
    if (length(valid_samples) < 2) {
      warning(
        sender, " -> ", receiver, ": insufficient valid samples (", length(valid_samples), " < 2).\n",
        "Check: min_samples or sample collection.",
        immediate. = TRUE
      )
      return(NULL)
    }

    colnames_split <- strsplit(colnames(rna.data), id_sep, fixed = TRUE)
    cell_types_loc <- sapply(colnames_split, `[`, cell_type_col)
    samples_loc <- sapply(colnames_split, `[`, sample_col)

    rna.data <- rna.data[, samples_loc %in% valid_samples]

    # Step 3: Load the ligand-receptor pairs after filtering for interactions
    lr <- filtered_lr

    # Step 4: Reorganization average expression for each sample-cell type group
    if (verbose) {
      message("Reorganization average expression for each sample-cell type group...")
    }
    colnames_split <- strsplit(colnames(rna.data), id_sep, fixed = TRUE)
    cell_types_loc <- sapply(colnames_split, `[`, cell_type_col)
    samples_loc <- sapply(colnames_split, `[`, sample_col)

    colnames(rna.data) <- paste0(samples_loc, "-lr-", cell_types_loc)
    rna.avg <- rna.data
  }


  # Unified steps
  rna.avg <- round(rna.avg, 5)

  avg.s <- rna.avg[, grep(sender, colnames(rna.avg))]
  avg.r <- rna.avg[, grep(receiver, colnames(rna.avg))]

  colnames(avg.s) <- str_match(colnames(avg.s), "^(.*)-lr-")[,2]
  colnames(avg.r) <- str_match(colnames(avg.r), "^(.*)-lr-")[,2]

  avg.r <- avg.r[, colnames(avg.s), drop = FALSE]

  avg.s <- avg.s[lr$ligand, , drop = FALSE]
  avg.r <- avg.r[lr$receptor, , drop = FALSE]

  # Step 5: Calculating projection scores
  if (verbose) {
    message("Calculating projection scores...")
  }
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

  if (verbose) {
    message("LR projection score analysis complete. Identified ", nrow(score.df), " significant ligand-receptor pairs.\n")
  }

  return(score.df)
}



#' Analyze Ligand-Receptor Projection Scores (Across All Cell Types)
#'
#' @description
#' This function calculates the ligand-receptor (LR) projection scores between all combinations
#' of sender and receiver cell types, and it supports both Seurat objects and average expression
#' matrices (matrix of gene expression data with cell types and samples as column names).
#' The projection score is computed based on linear regression models, measuring
#' the normalized distance of each sample's LR expression from the origin of the regression line.
#'
#' @param rna A Seurat object or a matrix containing single-cell RNA expression data.
#' @param filtered_lr A data frame of ligand-receptor pairs from prior analysis (e.g., output of `filter_lr_single`).
#'                  Must contain an "lr" column with pair identifiers in "Ligand_Receptor" format.
#' @param sample_col Metadata column name (character) for sample identifiers in Seurat mode; Matrix mode uses column index (numeric).
#' @param cell_type_col Metadata column name (character) for cell type in Seurat mode; Matrix mode uses column index (numeric).
#' @param id_sep Separator used in matrix column names to split sample and cell type (e.g., `--` for "Cardiac--sample1"). Only used in Matrix mode.
#' @param min_cells Minimum number of cells per sample for both sender and receiver (numeric, default 50). Only used in Seurat mode.
#' @param num_cores Number of CPU cores for parallel processing (numeric, default 10). Automatically capped at (system cores - 1).
#' @param verbose Logical indicating whether to print progress messages (logical, default: TRUE).
#'
#' @return A data frame with projection scores per sample and LR pair. Columns:
#'   \item{All input from \code{filtered_lr}}{Original columns provided by the user in \code{filtered_lr}.}
#'   \item{sample}{Sample identifier.}
#'   \item{score}{Projection score (raw co-expression intensity).}
#'   \item{normalized_score}{Normalized score scaled between 0-1.}
#' Rows are ordered by \code{filtered_lr} columns and descending \code{score}.
#'
#' Returns \code{NULL} if:
#' \itemize{
#'   \item No cell types are found in the metadata.
#'   \item One or both of the specified sender and receiver cell types are missing in the data.
#'   \item Fewer than two valid samples remain after filtering based on minimum cell number per sample.
#' }
#'
#' @export
#'
#' @importFrom dplyr %>% bind_rows
#' @importFrom utils head
#'
#' @examples
#' \donttest{
#'   data(matrix_object)
#'   data(lr_db)
#'
#'   # Analyzing ligand-receptor interactions between all cell types
#'   result01a <- filter_lr_all(
#'     rna = matrix_object,
#'     lr_database = lr_db,
#'     sample_col = 2,
#'     cell_type_col = 1,
#'     id_sep = "--",
#'     min_samples = 10,
#'     min_sample_ratio = 0.1,
#'     min_adjust_p = 0.05,
#'     num_cores = 1,
#'     verbose = TRUE
#'     )
#'
#'   # Analyzing ligand-receptor projection scores between all cell types
#'   result02a <- score_lr_all(
#'     rna = matrix_object,
#'     filtered_lr = result01a,
#'     sample_col = 2,
#'     cell_type_col = 1,
#'     id_sep = "--",
#'     num_cores = 1,
#'     verbose = TRUE
#'     )
#'
#'   if (!is.null(result02a)) {
#'   print(head(result02a))
#'   }
#' }
score_lr_all <- function(rna, filtered_lr,
                         sample_col, cell_type_col,
                         id_sep,                            # v0.2.0.0
                         min_cells = 50, num_cores = 10, verbose = TRUE) {

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
  if (missing(sample_col) || missing(cell_type_col)) {
    stop("When using expression matrix input, both 'sample_col' and 'cell_type_col' must be specified as integers.")
  }

  # Check if input is Seurat or average expression matrix v0.2.0.0
  is_seurat <- inherits(rna, "Seurat")
  is_matrix <- is.matrix(rna)

  if (!is_seurat && !is_matrix) {
    stop("'rna' must be either a Seurat object or a matrix of average expression data.")
  }


  if (is_seurat) {
    input_type <- "Seurat"
    if (!is.character(sample_col)) {
      stop("Seurat mode: 'sample_col' must be a character string.")
    }
    if (!is.character(cell_type_col)) {
      stop("Seurat mode: 'cell_type_col' must be a character string.")
    }

    if (!sample_col %in% colnames(rna@meta.data)) {
      stop(paste("Column", sample_col, "not found in Seurat metadata"))
    }
    if (!cell_type_col %in% colnames(rna@meta.data)) {
      stop(paste("Column", cell_type_col, "not found in Seurat metadata"))
    }
    rna$sample <- rna@meta.data[, sample_col]
    rna$cell.type <- rna@meta.data[, cell_type_col]

  } else if (is_matrix && nrow(rna) > 1 && ncol(rna) > 1) {
    input_type <- "Matrix"
    if (!is.numeric(sample_col) || sample_col <= 0 || sample_col %% 1 != 0) {
      stop("Matrix mode: 'sample_col' must be a positive integer.")
    }
    if (!is.numeric(cell_type_col) || cell_type_col <= 0 || cell_type_col %% 1 != 0) {
      stop("Matrix mode: 'cell_type_col' must be a positive integer.")
    }
    colnames_split <- strsplit(colnames(rna), id_sep, fixed = TRUE)
    if (!all(lengths(colnames_split) >= max(sample_col, cell_type_col))) {
      stop("Column names do not contain enough fields with the specified separator.")
    }
  } else {
    stop("Input must be either a Seurat object or an average expression matrix.")
  }


  # check and set
  if (input_type == "Seurat") {
    # Handle cell type check
    cell_types <- unique(rna@meta.data[[cell_type_col]])
    if (length(cell_types) < 1) {
      stop("No cell types found in column '", cell_type_col, "'.", call. = FALSE)
    }

    if (verbose) {
      message("(Seurat mode) Analyzing ligand-receptor projection scores between all cell types...")
      message("Cell types: ", paste(cell_types, collapse = ", "), "\n")
    }

    sample_col <- "sample"
    cell_type_col <- "cell.type"
    id_sep <- ""


  } else if (input_type == "Matrix") {
    # Split the column names using the specified separator
    cell_types <- unique(sapply(colnames_split, `[`, cell_type_col))
    samples <- unique(sapply(colnames_split, `[`, sample_col))
    if (length(cell_types) < 1) {
      stop("No cell types found in column '", cell_type_col, "'.", call. = FALSE)
    }

    if (verbose) {
      message("(Matrix mode) Analyzing ligand-receptor projection scores between all cell types...")
      message("Cell types: ", paste(cell_types, collapse = ", "), "\n")
    }

    sample_col <- sample_col
    cell_type_col <- cell_type_col
  }


  # main
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
      sample_col = sample_col,
      cell_type_col = cell_type_col,
      id_sep = id_sep,
      min_cells = min_cells,
      num_cores = num_cores,
      verbose = verbose
    )

    if (!is.null(res) && nrow(res) > 0) {
      return(res)
    } else {
      return(NULL)
    }
  })
  res_list <- Filter(Negate(is.null), res_list)
  final_res <- bind_rows(res_list)

  # Check if the projection result is empty
  if (nrow(final_res) == 0) {
    message("\n\nNo significant ligand-receptor pairs were identified in the projection score analysis.")
    return(NULL)
  }

  if (verbose) {
    message("\n\nLR projection score analysis across all cell type combinations complete.
            Identified ", nrow(final_res), " significant ligand-receptor pairs.")
  }

  return(final_res)
}
