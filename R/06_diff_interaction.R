#' Extracts Metadata from Seurat Object
#'
#' @description
#' Extracts metadata from a Seurat object and removes duplicates by sample.
#'
#' @param rna A Seurat object containing single-cell RNA expression data.
#' @param sample_col Column name in Seurat metadata indicating sample identifiers (character).
#'
#' @return A data frame with unique metadata per sample.
#'#'
#' @importFrom dplyr %>% bind_rows
#' @importFrom stats lm sd coef cor.test p.adjust
#' @importFrom utils head
#'
#' @examples
#' seurat_object <- load_example_seurat()
#' metadata <- get_sample_metadata(rna = seurat_object, sample_col = "sample")
#'
#' @keywords internal
#' @noRd
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
#'
#' @keywords internal
#' @noRd
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
#'
#' @keywords internal
#' @noRd
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



#' Compare Ligand-Receptor Interaction Scores with Group Variable using Linear Regression
#'
#' @description
#' Perform linear regression analysis to compare ligand-receptor (LR) interaction scores
#' across groups, handling both continuous and binary group variables (ident1 vs ident2 or all others).
#'
#' @param lr_scores Data frame containing LR interaction scores per sample (data frame).
#' @param metadata Data frame containing sample metadata (data frame).
#' @param group_variable Column name in \code{metadata} to compare groups (categorical or continuous) (character).
#' @param ident1 If categorical, group to compare (coded as 1) (character).
#' @param ident2 Reference group or list of groups (coded as 0). If None, uses all others (character).
#' @param covariates Optional list of covariate column names (character vector).
#' @param fdr_threshold Significance cutoff for adjusted p-values (numeric, default: 0.05).
#'
#' @return Data frame with ligand, receptor, sender, receiver, coef (coefficient, logFC), p-values, and adjusted p-values.
#'
#' @export
#'
#' @importFrom dplyr filter mutate inner_join group_by summarise select across all_of where arrange
#' @importFrom tidyr pivot_wider separate
#' @importFrom broom tidy
#' @importFrom stats lm p.adjust reformulate
#' @importFrom tibble column_to_rownames
#' @importFrom purrr map_dfr
#'
#' @examples
#' \donttest{
#'   # Long-running example (may take >10s)
#'   data(lr_scores_eg)
#'   data(metadata_eg)
#'
#'   result <- lr_linear_model_discrete(
#'     lr_scores_eg, metadata_eg,
#'     group_variable = "IFN_type",
#'     ident1 = "high",
#'     covariates = c("Age_group", "Sex")
#'   )
#'   head(result)
#' }
lr_linear_model_discrete <- function(lr_scores,
                                     metadata,
                                     group_variable,
                                     ident1,
                                     ident2 = NULL,
                                     covariates = NULL,
                                     fdr_threshold = 0.05) {

  score_col <- "score"

  # Parameter validation
  if (!is.character(group_variable) || length(group_variable) != 1) {
    stop("group_variable must be a single character string.")
  }

  if (!group_variable %in% colnames(metadata)) {
    stop(paste0("Group variable '", group_variable, "' not found in metadata."))
  }

  if (!"sample" %in% colnames(metadata) || !"sample" %in% colnames(lr_scores)) {
    stop("Both lr_scores and metadata must contain a 'sample' column.")
  }

  if (!is.null(covariates)) {
    missing_covs <- setdiff(covariates, colnames(metadata))
    if (length(missing_covs) > 0) {
      stop("Covariates not found in metadata: ", paste(missing_covs, collapse = ", "))
    }
  }

  metadata <- metadata %>% dplyr::mutate(sample = as.character(sample))
  lr_scores <- lr_scores %>% dplyr::mutate(sample = as.character(sample))

  # Group treatment: Continuous or categorical variables
  if (is.numeric(metadata[[group_variable]])) {
    metadata$group_dummy <- metadata[[group_variable]]
    keep_samples <- metadata
  } else {
    if (!ident1 %in% metadata[[group_variable]]) {
      stop(paste0("Group ident1 '", ident1, "' not found in group_variable."))
    }

    ident2 <- if (is.null(ident2)) {
      setdiff(unique(metadata[[group_variable]]), ident1)
    } else {
      intersect(as.character(ident2), unique(metadata[[group_variable]]))
    }

    if (length(ident2) == 0) {
      stop("No valid ident2 groups found in metadata.")
    }

    keep_samples <- metadata %>%
      dplyr::filter(.data[[group_variable]] %in% c(ident1, ident2)) %>%
      dplyr::mutate(group_dummy = as.integer(.data[[group_variable]] == ident1))
  }

  # Merge with LR scores and Create LR interaction ID
  merged_data <- lr_scores %>%
    dplyr::inner_join(keep_samples, by = "sample") %>%
    dplyr::mutate(LRSR = paste(.data[["ligand"]], .data[["receptor"]], .data[["sender"]], .data[["receiver"]], sep = "_"))

  # Create LRSR × sample score matrix
  score_matrix <- merged_data %>%
    dplyr::group_by(.data[["LRSR"]], sample) %>%
    dplyr::summarise(mean_score = mean(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = sample, values_from = mean_score, values_fill = list(mean_score = 0))

  samples_order <- colnames(score_matrix)[-1]

  # Preparing model data
  model_data <- keep_samples %>%
    dplyr::filter(sample %in% colnames(score_matrix)) %>%
    dplyr::select(sample, .data[["group_dummy"]], dplyr::all_of(covariates)) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor)) %>%
    tibble::column_to_rownames("sample")

  model_data <- model_data[samples_order, , drop = FALSE]

  # linear model
  result_df <- score_matrix %>%
    split(~ LRSR) %>%
    purrr::map_dfr(function(row) {
      lrsr <- row$LRSR[1]
      y <- as.numeric(unlist(row[,-1]))

      tryCatch({
        model_frame <- model_data
        model_frame$y <- y

        model_formula <- if (is.null(covariates)) {
          y ~ group_dummy
        } else {
          reformulate(termlabels = c("group_dummy", covariates), response = "y")
        }

        # Fit model
        fit <- stats::lm(model_formula, data = model_frame)
        coef_row <- broom::tidy(fit) %>%
          dplyr::filter(.data[["term"]] == "group_dummy")

        data.frame(
          LRSR = lrsr,
          coef = coef_row$estimate,
          p_value = coef_row$p.value
        )
      }, error = function(e) NULL)
    })

  # Sort out data
  result_df <- result_df %>%
    dplyr::mutate(adjusted_p_value = stats::p.adjust(.data[["p_value"]], method = "BH")) %>%
    tidyr::separate(.data[["LRSR"]], into = c("ligand", "receptor", "sender", "receiver"), sep = "_") %>%
    dplyr::filter(.data[["adjusted_p_value"]] < fdr_threshold) %>%
    dplyr::arrange(.data[["adjusted_p_value"]]) %>%
    dplyr::select(dplyr::all_of(c("ligand", "receptor", "sender", "receiver", "coef", "p_value", "adjusted_p_value"))) %>%
    dplyr::mutate(coef = round(.data[["coef"]], 5),
                  p_value = round(.data[["p_value"]], 15),
                  adjusted_p_value = round(.data[["adjusted_p_value"]], 15))


  return(result_df)
}

