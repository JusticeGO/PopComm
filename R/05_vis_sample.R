#' Generate Heatmap of LR Interaction Scores
#'
#' @description
#' This function generates a heatmap to visualize the ligand-receptor (LR) interaction scores across samples.
#' Rows represent LR pairs and columns represent samples. Optionally, sample metadata can be used to annotate the columns.
#'
#' @param lr_scores Data frame containing LR interaction scores per sample (data frame).
#' @param metadata Data frame containing sample metadata (data frame).
#' @param score Character string indicating which score to use: "normalized" (default) or "raw" .
#' @param selected_sender Specific sender cell type to filter, default is None (use all) (character).
#' @param selected_receiver Specific receiver cell type to filter, default is None (use all) (character).
#' @param selected_metadata List of column names in \code{metadata} to annotate samples (default: None, use all)(character vector).
#'
#' @return Heatmap of average LR interaction scores per sample.
#'
#' @export
#'
#' @importFrom dplyr %>% filter mutate group_by summarise
#' @importFrom tidyr pivot_wider
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tibble column_to_rownames
#' @importFrom stats setNames
#' @importFrom grDevices colorRampPalette
#' @importFrom rlang .data
#'
#' @examples
#' # Heatmap of LR Interaction Scores
#' data(lr_scores_eg)
#' data(metadata_eg)
#' p <- heatmap_sample(lr_scores_eg, metadata_eg, score = "normalized", selected_sender = "Cardiac",
#'   selected_receiver = "Perivascular", selected_metadata = c("Sex", "Age_group"))
#' print(p)
heatmap_sample <- function(lr_scores, metadata, score = "normalized",
                           selected_sender = NULL,
                           selected_receiver = NULL,
                           selected_metadata = NULL) {

  # Parameter validation
  if (!score %in% c("normalized", "raw"))
    stop("score must be either 'normalized' or 'raw'")
  score_col <- ifelse(score == "normalized", "normalized_score", "score")

  if (!is.null(selected_sender)) {
    lr_scores <- lr_scores %>% filter(.data[["sender"]] == selected_sender)
  }
  if (!is.null(selected_receiver)) {
    lr_scores <- lr_scores %>% filter(.data[["receiver"]] == selected_receiver)
  }

  # Create a new identifier by concatenating ligand, receptor, sender, and receiver
  lr_scores <- lr_scores %>%
    mutate(LRSR = paste(.data[["ligand"]], .data[["receptor"]], .data[["sender"]], .data[["receiver"]], sep = "_"))

  # Pivot data: calculate the average score for each LRSR-sample pair
  heatmap_data <- lr_scores %>%
    group_by(.data[["LRSR"]], sample) %>%
    summarise(mean_score = mean(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = mean_score, values_fill = list(mean_score = 0))

  # Convert to matrix with LRSR as row names
  heatmap_matrix <- as.data.frame(heatmap_data)
  rownames(heatmap_matrix) <- heatmap_matrix$LRSR
  heatmap_matrix <- heatmap_matrix[,-1]

  # Prepare annotation for columns using metadata.
  annotation_col <- NULL

  if ("sample" %in% colnames(metadata)) {
    metadata_sub <- metadata %>%
      filter(sample %in% colnames(heatmap_matrix))
    missing_samples <- setdiff(colnames(heatmap_matrix), metadata_sub$sample)
    if (length(missing_samples) > 0) {
      warning("The following sample was not found in metadata: ", paste(missing_samples, collapse = ", "))
    }
    metadata_sub <- column_to_rownames(metadata_sub, var = "sample")
    metadata_sub <- metadata_sub[colnames(heatmap_matrix), , drop = FALSE]
  } else if (all(colnames(heatmap_matrix) %in% rownames(metadata))) {
    metadata_sub <- metadata[rownames(metadata) %in% colnames(heatmap_matrix), , drop = FALSE]
    missing_samples <- setdiff(colnames(heatmap_matrix), rownames(metadata_sub))
    if (length(missing_samples) > 0) {
      warning("The following sample was not found in metadata: ", paste(missing_samples, collapse = ", "))
    }
    metadata_sub <- metadata_sub[colnames(heatmap_matrix), , drop = FALSE]
  } else {
    stop("No sample information corresponding to the heatmap column was found in metadata.")
  }

  if (is.null(selected_metadata)) {
    selected_metadata <- colnames(metadata_sub)
    if ("sample" %in% selected_metadata) {
      selected_metadata <- setdiff(selected_metadata, "sample")
    }
  } else {
    missing_cols <- setdiff(selected_metadata, colnames(metadata_sub))
    if (length(missing_cols) > 0) {
      stop("The following selected metadata columns do not exist in metadata: ", paste(missing_cols, collapse = ", "))
    }
  }

  annotation_col <- metadata_sub[, selected_metadata, drop = FALSE]
  annotation_col[] <- lapply(annotation_col, as.factor)

  # set color
  auto_annotation_colors <- list()
  for(colname in colnames(annotation_col)){
    levels_current <- levels(annotation_col[[colname]])
    n_levels <- length(levels_current)
    if(n_levels < 3){
      base_colors <- brewer.pal(3, "Paired")[1:n_levels]
    } else {
      base_colors <- brewer.pal(min(9, n_levels), "Paired")
    }
    auto_colors <- colorRampPalette(base_colors)(n_levels)
    auto_annotation_colors[[colname]] <- setNames(auto_colors, levels_current)
  }

  # Plot the heatmap using pheatmap
  ph <- pheatmap(heatmap_matrix,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                 annotation_col = annotation_col,
                 annotation_colors = auto_annotation_colors,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 show_rownames = FALSE,
                 show_colnames = FALSE,
                 fontsize = 8,
                 border_color = "gray")

  return(ph)
}



#' Generate PCA of LR Interaction Scores
#'
#' @description
#' This function performs principal component analysis (PCA) on ligand-receptor (LR) interaction scores across samples,
#' and generates a scatter plot of the first two principal components. Optionally, sample metadata can be used to color the points.
#'
#' @param lr_scores Data frame containing LR interaction scores per sample (data frame).
#' @param metadata Data frame containing sample metadata (data frame).
#' @param score Character string indicating which score to use: "normalized" (default) or "raw" .
#' @param selected_sender Specific sender cell type to filter, default is None (use all) (character).
#' @param selected_receiver Specific receiver cell type to filter, default is None (use all) (character).
#' @param color_by \code{metadata} column name to color points in PCA plot (character).
#' @param n_components Number of principal components to extract (numeric, default: 2).
#'
#' @return A list with two elements: the first is a ggplot2 PCA scatter plot and the second is the PCA results data frame.
#'
#' @export
#'
#' @importFrom dplyr %>% filter mutate group_by summarise
#' @importFrom tidyr pivot_wider
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes geom_hline geom_point labs theme_minimal scale_color_brewer geom_vline element_rect
#' @importFrom rlang .data
#'
#' @examples
#' # PCA of LR Interaction Scores
#' data(lr_scores_eg)
#' data(metadata_eg)
#' res <- pca_sample(lr_scores_eg, metadata_eg, score = "normalized",
#'   selected_sender = "Cardiac", selected_receiver = "Perivascular", color_by = "Age_group")
#' print(res$pca_plot)
pca_sample <- function(lr_scores, metadata, score = "normalized",
                       selected_sender = NULL, selected_receiver = NULL,
                       color_by = NULL, n_components = 2) {

  # Parameter validation
  if (!score %in% c("normalized", "raw"))
    stop("score must be either 'normalized' or 'raw'")
  score_col <- ifelse(score == "normalized", "normalized_score", "score")

  # Filter by sender and receiver cell types if specified
  if (!is.null(selected_sender)) {
    lr_scores <- lr_scores %>% filter(.data[["sender"]] == selected_sender)
  }
  if (!is.null(selected_receiver)) {
    lr_scores <- lr_scores %>% filter(.data[["receiver"]] == selected_receiver)
  }

  # Create a new identifier by concatenating ligand, receptor, sender, and receiver
  lr_scores <- lr_scores %>%
    mutate(LRSR = paste(.data[["ligand"]], .data[["receptor"]], .data[["sender"]], .data[["receiver"]], sep = "_"))

  # Pivot table: calculate the average score for each LRSR-sample pair
  pca_data <- lr_scores %>%
    group_by(.data[["LRSR"]], sample) %>%
    summarise(mean_score = mean(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = mean_score, values_fill = list(mean_score = 0))

  # Convert to a matrix with LRSR as row names, then transpose so that rows are samples and columns are features
  pca_matrix <- as.data.frame(pca_data)
  rownames(pca_matrix) <- pca_matrix$LRSR
  pca_matrix <- as.matrix(pca_matrix[,-1])
  pca_matrix <- t(pca_matrix)

  # Standardize the data (center and scale)
  standardized_data <- scale(pca_matrix)

  # Perform PCA using prcomp
  pca_result <- prcomp(standardized_data, center = FALSE, scale. = FALSE)
  # Select the required number of principal components
  pca_scores <- as.data.frame(pca_result$x[, 1:n_components, drop = FALSE])

  # Add sample names as a column for merging
  pca_scores$sample <- rownames(pca_scores)

  # Merge PCA results with metadata if color_by is provided and exists in metadata
  if (!is.null(color_by) && color_by %in% colnames(metadata)) {
    if ("sample" %in% colnames(metadata)) {
      pca_scores <- merge(pca_scores, metadata[, c("sample", color_by), drop = FALSE],
                          by = "sample", all.x = TRUE)
    } else {
      warning("Metadata does not contain a 'sample' column; PCA result will not be annotated by metadata.")
    }
  }

  # Create PCA scatter plot using ggplot2
  if (!is.null(color_by) && color_by %in% colnames(pca_scores)) {
    p <- ggplot(pca_scores, aes(x = .data[["PC1"]], y = .data[["PC2"]], color = .data[[color_by]])) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
      geom_point(size = 3) +
      labs(x = "Principal Component 1", y = "Principal Component 2",
           title = "PCA of Ligand-Receptor Interaction Scores", color = color_by) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(size = 12, hjust = 0.5, color = "black"),
        legend.position = "right",
        axis.line = element_line(),
        axis.ticks = element_line(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      scale_color_brewer(palette = "Set1")
  } else {
    p <- ggplot(pca_scores, aes(x = .data[["PC1"]], y = .data[["PC2"]])) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
      geom_point(size = 3, color = "black") +
      labs(x = "Principal Component 1", y = "Principal Component 2",
           title = "PCA of Ligand-Receptor Interaction Scores") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(size = 12, hjust = 0.5, color = "black"),
        axis.line = element_line(),
        axis.ticks = element_line(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  }

  # Return a list containing the plot and the PCA results data frame
  return(list(pca_plot = p, pca_result = pca_scores))
}
