# 02_analyze_projection
# This function calculates the projection scores for ligand-receptor (LR) pairs between specified sender and receiver cell types. 
# The projection score is computed based on linear regression models, measuring the normalized distance of each sample's LR expression 
# from the regression line.
# 
# Parameters:
#   rna: A Seurat object containing single-cell RNA expression data.
#   sender: Cell type designated as the ligand sender (character).
#   receiver: Cell type designated as the receptor receiver (character).
#   lr_custom: A filtered data frame of ligand-receptor pairs from prior analysis (e.g., output of `filter_lr_single`).
#              Must contain an "lr" column with pair identifiers in "Ligand_Receptor" format.
#   sample_col: Column name in Seurat metadata indicating sample identifiers (character).
#   cell_type_col: Column name in Seurat metadata indicating cell type classifications (character).
#   min_cells: Minimum cells required per sample for both sender and receiver (default 50).
#   min_samples: Minimum valid samples required to proceed (default 10).
#   mc_cores: Number of CPU cores for parallel processing (default 10). Automatically capped at (system cores - 1).
# 
# Returns:
#   A data frame containing projection scores for each sample and LR pair, with the following columns:
#    - All columns from the input `lr_custom`.
#    - sample: Sample identifier.
#    - normalized_score: Normalized projection score (range 0-1) indicating the relative co-expression intensity of the LR pair in the sample.
#   Rows are ordered by `lr_custom`. 
#   Returns NULL if no valid pairs are found.



analyze_lr_projection <- function(rna, sender, receiver, lr_custom, 
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
  rna$sample <- rna@meta.data[,sample_col]
  rna$cell.type <- rna@meta.data[,cell_type_col]
  cell_types <- unique(rna@meta.data[[cell_type_col]])
  
  if (length(cell_types) < 1) stop("No cell types found.")
  selected_types <- unique(c(sender, receiver))
  missing_types <- setdiff(selected_types, cell_types)
  if (length(missing_types) > 0) {
    stop("Missing cell types: ", paste(missing_types, collapse = ", "))
  }
  
  # Step 1: Add sample and cell type columns to the RNA data object and subset
  message("\nAnalyzing ligand-receptor projection score: ", sender, " -> ", receiver)
  # Determine the subset of data
  if (!setequal(selected_types, cell_types)) {
    message("Subsetting ", sender, " (sender) and ", receiver, " (receiver).")
    rna.data <- subset(rna, cell.type %in% selected_types)
  } else {
    message("Using full dataset.")
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
  lr$ligand <- str_match(lr$lr, "^(.*)_")[,2]
  lr$receptor <- str_match(lr$lr, "_(.*)$")[,2]
  
  # Step 4: Compute average expression for each sample-cell type group
  rna.data$group <- paste0(rna.data$sample, "-lr-", rna.data$cell.type)
  message("Computing average expression for each sample-cell type group...")
  rna.avg <- AverageExpression(rna.data, group.by = "group")$RNA
  
  avg.s <- rna.avg[, grep(sender, colnames(rna.avg))]
  avg.r <- rna.avg[, grep(receiver, colnames(rna.avg))]
  
  colnames(avg.s) <- str_match(colnames(avg.s), "^(.*)-lr-")[,2]
  colnames(avg.r) <- str_match(colnames(avg.r), "^(.*)-lr-")[,2]
  
  avg.r <- avg.r[,colnames(avg.s)]
  
  avg.s <- avg.s[lr$ligand,]
  avg.r <- avg.r[lr$receptor,]

  # Step 5: Calculating projection scores
  message("Calculating projection scores...")
  score.df <- pbmclapply(
    1:nrow(avg.s), function(i) {
      x <- as.numeric(avg.s[lr$ligand[i], ])
      y <- as.numeric(avg.r[lr$receptor[i], ])
      
      if (sd(x) == 0 || sd(y) == 0) {
        return(data.frame())
      }
      
      model <- tryCatch(
        lm(y ~ x),
        error = function(e) NULL
      )
      
      if (is.null(model)) return(data.frame())
      
      slope <- coef(model)[2]
      intercept <- coef(model)[1]
      
      projections <- t(sapply(
        1:length(x),
        function(j) project_to_line(x[j], y[j], slope, intercept)
      ))
      
      dx <- projections[, 1] - min(projections[, 1])
      dy <- projections[, 2] - min(projections[, 2])
      score <- sqrt(dx^2 + dy^2)
      normalized.score <- score / max(score)
      
      lr_metadata <- lr[i, ]
      data.frame(
        lr_metadata,
        sample = colnames(avg.s),
        normalized_score = normalized.score,
        row.names = NULL
      )
    }, 
    mc.cores = mc_cores
  ) %>% 
    bind_rows() %>%
    na.omit()
  
  message("\n\nAnalyzing ligand-receptor projection score process complete.")
  message("Head of final results (", nrow(score.df), "):")
  print(head(score.df))
  
  return(score.df)
}
