#' Extracts Metadata from Seurat Object
#'
#' @description
#' Extracts metadata from a Seurat object and removes duplicates by sample.
#'
#' @param rna A Seurat object containing single-cell RNA expression data.
#' @param sample_col Column name in Seurat metadata indicating sample identifiers (character).
#'
#' @return A data frame with unique metadata per sample.
#'
#' @export
#'
#' @importFrom dplyr %>% bind_rows
#' @importFrom stats lm sd coef cor.test p.adjust
#' @importFrom utils head
#'
#' @examples
#' seurat_object <- load_example_seurat()
#' metadata <- get_sample_metadata(rna = seurat_object, sample_col = "sample")
get_sample_metadata <- function(rna, sample_col) {
  metadata <- rna@meta.data
  unique_metadata <- metadata[!duplicated(metadata[[sample_col]]), ]
  rownames(unique_metadata) <- unique_metadata[[sample_col]]
  return(unique_metadata)
}



#' Correlation Analysis between Ligand-Receptor Interaction Scores and Continuous Metadata
#'
#' @description
#' This function performs a Spearman correlation analysis between ligand-receptor (LR) interaction scores
#' (using raw scores) and a continuous metadata variable. It computes the correlation coefficient, p-value,
#' and adjusts p-values for multiple testing using the Benjamini-Hochberg method.
#'
#' @param lr_scores Data frame containing LR interaction scores per sample (data frame).
#' @param metadata Data frame containing sample metadata (data frame).
#' @param correlate_with Column name in \code{metadata} to compute correlation with (character).
#'
#' @return A data frame containing correlation results (correlation coefficient and p-values).
#'
#' @export
#'
#' @importFrom dplyr %>% mutate group_by summarise left_join filter arrange select
#' @importFrom tidyr pivot_wider
#' @importFrom stats cor.test p.adjust
#' @importFrom tibble column_to_rownames
#' @importFrom tidyselect all_of
#' @importFrom rlang .data
#'
#' @examples
#' data(lr_scores_eg)
#' data(metadata_eg)
#' res <- corr_lr_interaction(lr_scores_eg, metadata_eg, correlate_with = "Age_imputation")
#' head(res)
corr_lr_interaction <- function(lr_scores, metadata, correlate_with) {

  score_col <- "score"
  # Check and ensure metadata has 'sample' column; if not, use rownames as sample
  if (!"sample" %in% colnames(metadata)) {
    metadata <- data.frame(sample = rownames(metadata), metadata, row.names = NULL)
  }

  # Ensure sample names are character
  metadata$sample <- as.character(metadata$sample)
  lr_scores$sample <- as.character(lr_scores$sample)

  # Merge metadata with lr_scores for the correlate_with column
  metadata_sub <- metadata %>%
    select(sample, all_of(correlate_with))
  merged_data <- lr_scores %>%
    left_join(metadata_sub, by = "sample")

  # Create a new identifier for LR-CellType pairs
  merged_data <- merged_data %>%
    mutate(LRSR = paste(.data[["ligand"]], .data[["receptor"]], .data[["sender"]], .data[["receiver"]], sep = "_"))

  # Pivot table: calculate average raw score for each LRSR-sample pair
  heatmap_data <- merged_data %>%
    group_by(.data[["LRSR"]], sample) %>%
    summarise(mean_score = mean(.data[["score"]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = mean_score, values_fill = list(mean_score = 0))

  # Extract sample names (columns) and ensure continuous values are ordered accordingly
  sample_names <- colnames(heatmap_data)[-1]  # first column is LRSR


  if ("sample" %in% colnames(metadata)) {
    metadata_sub <- metadata %>% filter(sample %in% sample_names)
    missing_samples <- setdiff(sample_names, metadata_sub$sample)
    if (length(missing_samples) > 0) {
      warning("The following samples were not found in metadata: ", paste(missing_samples, collapse = ", "))
    }
    metadata_sub <- column_to_rownames(metadata_sub, var = "sample")
    metadata_sub <- metadata_sub[sample_names, , drop = FALSE]
  } else if (all(sample_names %in% rownames(metadata))) {
    metadata_sub <- metadata[rownames(metadata) %in% sample_names, , drop = FALSE]
    missing_samples <- setdiff(sample_names, rownames(metadata_sub))
    if (length(missing_samples) > 0) {
      warning("The following samples were not found in metadata: ", paste(missing_samples, collapse = ", "))
    }
    metadata_sub <- metadata_sub[sample_names, , drop = FALSE]
  } else {
    stop("No sample information corresponding to the columns was found in metadata.")
  }

  continuous_values <- metadata_sub[[correlate_with]]

  # Initialize vectors to store results
  correlations <- numeric(nrow(heatmap_data))
  p_values <- numeric(nrow(heatmap_data))

  # Compute Spearman correlation for each LRSR
  for (i in seq_len(nrow(heatmap_data))) {
    scores <- as.numeric(heatmap_data[i, -1])
    ct <- tryCatch({
      test <- cor.test(scores, continuous_values, method = "spearman")
      c(test$estimate, test$p.value)
    }, error = function(e) {
      warning(sprintf("Error in correlation test for LRSR: %s", heatmap_data$LRSR[i]))
      c(NA, NA)
    })
    correlations[i] <- ct[1]
    p_values[i] <- ct[2]
  }

  # Adjust p-values using Benjamini-Hochberg correction
  adjusted_p_values <- p.adjust(p_values, method = "BH")

  # Split LRSR into separate columns
  lr_split <- strsplit(heatmap_data$LRSR, "_")
  ligand <- sapply(lr_split, `[`, 1)
  receptor <- sapply(lr_split, `[`, 2)
  sender <- sapply(lr_split, `[`, 3)
  receiver <- sapply(lr_split, `[`, 4)

  # Construct result data frame
  result_df <- data.frame(
    ligand = ligand,
    receptor = receptor,
    sender = sender,
    receiver = receiver,
    correlation = correlations,
    p_value = p_values,
    adjusted_p_value = adjusted_p_values,
    stringsAsFactors = FALSE
  )

  # Sort by adjusted p-value
  result_df <- result_df[order(result_df$adjusted_p_value), ]

  return(result_df)
}



#' Differential Ligand-Receptor Interaction Analysis between Two Groups
#'
#' @description
#' This function performs a differential ligand-receptor (LR) interaction analysis between two sample groups.
#' It computes the average interaction scores per ligand-receptor-cell type (LRSR) pair for each group,
#' calculates the log fold change (logFC) between the groups, performs a Mann–Whitney U test (Wilcoxon rank-sum test)
#' for each LRSR pair, and adjusts the resulting p-values for multiple testing using the Benjamini-Hochberg method.
#'
#' @param lr_scores Data frame containing LR interaction scores per sample (data frame).
#' @param metadata Data frame containing sample metadata (data frame).
#' @param group_by Column name in \code{metadata} to define sample groups (character).
#' @param ident1 Value in \code{group_by} that defines the first group.
#' @param ident2 Value in \code{group_by} that defines the second group.
#'
#' @return A data frame containing differential interaction results (mean scores, logFC, p-values).
#'
#' @export
#'
#' @importFrom dplyr %>% mutate group_by summarise left_join filter select
#' @importFrom tidyr pivot_wider
#' @importFrom stats wilcox.test p.adjust
#' @importFrom tidyselect all_of
#' @importFrom rlang .data
#'
#' @examples
#' data(lr_scores_eg)
#' data(metadata_eg)
#' res <- diff_lr_interaction(lr_scores_eg, metadata_eg, group_by = "Age_group",
#'   ident1 = "Young", ident2 = "Old")
#' head(res)
diff_lr_interaction <- function(lr_scores, metadata, group_by, ident1, ident2) {

  score_col <- "score"
  # Check and ensure metadata has 'sample' column; if not, use rownames as sample
  if (!"sample" %in% colnames(metadata)) {
    metadata <- data.frame(sample = rownames(metadata), metadata, row.names = NULL)
  }

  # Ensure sample names are character
  metadata$sample <- as.character(metadata$sample)
  lr_scores$sample <- as.character(lr_scores$sample)

  metadata_sub_temp <- metadata %>% select(sample, all_of(group_by))
  merged_data <- lr_scores %>% left_join(metadata_sub_temp, by = "sample")

  # Create a new identifier for LR-CellType pairs
  merged_data <- merged_data %>%
    mutate(LRSR = paste(.data[["ligand"]], .data[["receptor"]], .data[["sender"]], .data[["receiver"]], sep = "_"))

  # Pivot table: calculate average raw score for each LRSR-sample pair
  heatmap_data <- merged_data %>%
    group_by(.data[["LRSR"]], sample) %>%
    summarise(mean_score = mean(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = mean_score, values_fill = list(mean_score = 0))

  # Retrieve sample names for each group from metadata
  sample_names <- colnames(heatmap_data)[-1]
  if ("sample" %in% colnames(metadata)) {
    metadata_sub <- metadata %>% filter(sample %in% sample_names)
    missing_samples <- setdiff(sample_names, metadata_sub$sample)
    if (length(missing_samples) > 0) {
      warning("The following samples were not found in metadata: ", paste(missing_samples, collapse = ", "))
    }
    metadata_sub <- column_to_rownames(metadata_sub, var = "sample")
    metadata_sub <- metadata_sub[sample_names, , drop = FALSE]
  } else if (all(sample_names %in% rownames(metadata))) {
    metadata_sub <- metadata[rownames(metadata) %in% sample_names, , drop = FALSE]
    missing_samples <- setdiff(sample_names, rownames(metadata_sub))
    if (length(missing_samples) > 0) {
      warning("The following samples were not found in metadata: ", paste(missing_samples, collapse = ", "))
    }
    metadata_sub <- metadata_sub[sample_names, , drop = FALSE]
  } else {
    stop("No sample information corresponding to the heatmap columns was found in metadata.")
  }

  samples_ident1 <- rownames(metadata_sub)[metadata_sub[[group_by]] == ident1]
  samples_ident2 <- rownames(metadata_sub)[metadata_sub[[group_by]] == ident2]

  # Calculate group-wise mean scores per LRSR pair
  group1_data <- rowMeans(as.data.frame(select(heatmap_data, all_of(samples_ident1))), na.rm = TRUE)
  group2_data <- rowMeans(as.data.frame(select(heatmap_data, all_of(samples_ident2))), na.rm = TRUE)


  # Compute log fold change, adding a small constant to avoid log(0)
  logFC <- log2(group1_data + 1e-6) - log2(group2_data + 1e-6)

  # Perform Mann–Whitney U test (Wilcoxon rank-sum test) for each LRSR pair
  p_values <- sapply(seq_len(nrow(heatmap_data)), function(i) {
    values1 <- as.numeric(heatmap_data[i, samples_ident1])
    values2 <- as.numeric(heatmap_data[i, samples_ident2])
    test <- tryCatch(
      wilcox.test(values1, values2, alternative = "two.sided"),
      error = function(e) {
        warning(sprintf("Error in Wilcoxon test for LRSR: %s", heatmap_data$LRSR[i]))
        return(NULL)
      }
    )
    if (is.null(test)) {
      return(NA)
    } else {
      return(test$p.value)
    }
  })

  # Adjust p-values using Benjamini-Hochberg correction
  adjusted_p_values <- p.adjust(p_values, method = "BH")

  # Split LRSR into separate columns
  lr_split <- strsplit(heatmap_data$LRSR, "_")
  ligand <- sapply(lr_split, `[`, 1)
  receptor <- sapply(lr_split, `[`, 2)
  sender <- sapply(lr_split, `[`, 3)
  receiver <- sapply(lr_split, `[`, 4)

  # Construct result data frame
  result_df <- data.frame(
    ligand = ligand,
    receptor = receptor,
    sender = sender,
    receiver = receiver,
    mean_group1 = group1_data,
    mean_group2 = group2_data,
    logFC = logFC,
    p_value = p_values,
    adjusted_p_value = adjusted_p_values,
    stringsAsFactors = FALSE
  )

  # Sort by adjusted p-value
  result_df <- result_df[order(result_df$adjusted_p_value), ]

  return(result_df)
}
