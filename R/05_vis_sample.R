#' Generate Heatmap of LR Interaction Scores
#'
#' @description
#' This function generates a heatmap to visualize the average ligand-receptor (LR) interaction scores across samples.
#' Rows represent LR pairs and columns represent samples. Optionally, sample metadata can be used to annotate the columns.
#'
#' @param lr_scores Data frame containing LR interaction scores (data frame).
#' @param metadata Data frame with sample metadata; when a "sample" column exists, it is used to match sample identifiers for annotation (data frame).
#' @param score Character string indicating which score to use: "normalized" (default) or "raw" .
#' @param selected_sender Specific sender cell type to filter, default is None (use all) (character/character vector).
#' @param selected_receiver Specific receiver cell type to filter, default is None (use all) (character/character vector).
#' @param selected_metadata List of column names in metadata to annotate samples (default: None, use all)(character vector).
#' @param export_file Optional character string specifying the file name to save the plot (character).
#' @param export_format Character string defining the export format, with options "pdf" (default) and "png" (character).
#'
#' @return Heatmap of average LR interaction scores per sample.
#'
#' @export
#'
#' @importFrom dplyr %>% filter mutate group_by summarise
#' @importFrom tidyr pivot_wider
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices pdf png dev.off
#'
#' @examples
#' \dontrun{
#'   # Example data
#'   lr_scores <- data.frame(
#'     ligand = c("L1", "L2", "L1", "L2"),
#'     receptor = c("R1", "R2", "R1", "R2"),
#'     sender = c("A", "A", "B", "B"),
#'     receiver = c("X", "X", "Y", "Y"),
#'     sample = rep(c("S1", "S2"), 2),
#'     normalized_score = runif(4)
#'   )
#'
#'   metadata <- data.frame(
#'     sample = c("S1", "S2"),
#'     condition = c("control", "treatment")
#'   )
#'
#'   # Generate heatmap with metadata annotation
#'   p <- heatmap_sample(lr_scores, metadata, selected_metadata = "condition")
#' }
heatmap_sample <- function(lr_scores, metadata, score = "normalized",
                           selected_sender = NULL,
                           selected_receiver = NULL,
                           selected_metadata = NULL) {

  # Parameter validation
  if (!score %in% c("normalized", "raw"))
    stop("score must be either 'normalized' or 'raw'")
  score_col <- ifelse(score == "normalized", "normalized_score", "score")

  if (!is.null(selected_sender)) {
    lr_scores <- lr_scores %>% filter(sender == selected_sender)
  }
  if (!is.null(selected_receiver)) {
    lr_scores <- lr_scores %>% filter(receiver == selected_receiver)
  }

  # Create a new identifier by concatenating ligand, receptor, sender, and receiver
  lr_scores <- lr_scores %>%
    mutate(LRSR = paste(ligand, receptor, sender, receiver, sep = "_"))

  # Pivot data: calculate the average score for each LRSR-sample pair
  heatmap_data <- lr_scores %>%
    group_by(LRSR, sample) %>%
    summarise(mean_score = mean(get(score_col), na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = mean_score, values_fill = list(mean_score = 0))

  # Convert to matrix with LRSR as row names
  heatmap_matrix <- as.data.frame(heatmap_data)
  rownames(heatmap_matrix) <- heatmap_matrix$LRSR
  heatmap_matrix <- heatmap_matrix[,-1]

  # Prepare annotation for columns if metadata is provided
  annotation_col <- NULL
  if (!is.null(selected_metadata)) {
    if ("sample" %in% colnames(metadata)) {
      metadata_sub <- metadata %>%
        filter(sample %in% colnames(heatmap_matrix)) %>%
        tibble::column_to_rownames(var = "sample")
    } else {
      metadata_sub <- metadata[rownames(metadata) %in% colnames(heatmap_matrix), ]
    }
    annotation_col <- metadata_sub[, selected_metadata, drop = FALSE]
    annotation_col[] <- lapply(annotation_col, as.factor)
  }

  # Plot the heatmap using pheatmap
  ph <- pheatmap(heatmap_matrix,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                 annotation_col = annotation_col,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 show_rownames = FALSE,
                 show_colnames = FALSE,
                 fontsize = 10,
                 cellwidth = 15,
                 cellheight = 15,
                 border_color = "gray")

  # Export the plot if export_file is specified
  if (!is.null(export_file)) {
    if (export_format == "pdf") {
      pdf(export_file, width = figsize[1], height = figsize[2])
      pheatmap(heatmap_matrix,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
               annotation_col = annotation_col,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = FALSE,
               show_colnames = FALSE,
               fontsize = 10,
               cellwidth = 15,
               cellheight = 15,
               border_color = "gray")
      dev.off()
    } else if (export_format == "png") {
      png(export_file, width = figsize[1] * 100, height = figsize[2] * 100)
      pheatmap(heatmap_matrix,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
               annotation_col = annotation_col,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = FALSE,
               show_colnames = FALSE,
               fontsize = 10,
               cellwidth = 15,
               cellheight = 15,
               border_color = "gray")
      dev.off()
    }
  }

  return(ph)
}

#' Perform PCA on LR Interaction Scores and Plot the Results
#'
#' @description
#' This function performs principal component analysis (PCA) on ligand-receptor (LR) interaction scores across samples,
#' and generates a scatter plot of the first two principal components. Optionally, sample metadata can be used to color the points.
#'
#' @param lr_scores A data frame containing LR interaction scores. Must include columns: \code{ligand}, \code{receptor}, \code{sender},
#'   \code{receiver}, \code{sample}, and score columns (\code{score} and/or \code{normalized_score}).
#' @param metadata A data frame containing sample metadata. If a \code{sample} column exists, it is used for merging.
#' @param score A character string indicating which score to use: either \code{"normalized"} (default) or \code{"raw"}.
#' @param selected_sender Optional character string. If provided, only rows with this sender cell type are retained.
#' @param selected_receiver Optional character string. If provided, only rows with this receiver cell type are retained.
#' @param color_by Optional character string. A column name from \code{metadata} used to color the PCA plot.
#' @param n_components Integer specifying the number of PCA components to compute (default: 2).
#' @param figsize A numeric vector of length 2 indicating the figure size in inches (default: \code{c(8, 6)}).
#'
#' @return A data frame containing the PCA results (first \code{n_components} principal components) for each sample.
#'
#' @export
#'
#' @importFrom dplyr %>% filter mutate group_by summarise arrange desc ungroup
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal scale_color_brewer
#'
#' @examples
#' \dontrun{
#'   # Example data
#'   lr_scores <- data.frame(
#'     ligand = c("L1", "L2", "L1", "L2"),
#'     receptor = c("R1", "R2", "R1", "R2"),
#'     sender = c("A", "A", "B", "B"),
#'     receiver = c("X", "X", "Y", "Y"),
#'     sample = rep(c("S1", "S2"), 2),
#'     normalized_score = runif(4)
#'   )
#'
#'   metadata <- data.frame(
#'     sample = c("S1", "S2"),
#'     group = c("control", "treatment")
#'   )
#'
#'   # Perform PCA and plot, coloring by 'group'
#'   pca_res <- pca_sample(lr_scores, metadata, color_by = "group")
#' }
pca_sample <- function(lr_scores, metadata, score = "normalized",
                       selected_sender = NULL, selected_receiver = NULL,
                       color_by = NULL, n_components = 2, figsize = c(8, 6)) {
  score_col <- ifelse(score == "normalized", "normalized_score", "score")

  if (!is.null(selected_sender)) {
    lr_scores <- lr_scores %>% filter(sender == selected_sender)
  }
  if (!is.null(selected_receiver)) {
    lr_scores <- lr_scores %>% filter(receiver == selected_receiver)
  }

  # Create a new identifier by concatenating ligand, receptor, sender, and receiver
  lr_scores <- lr_scores %>%
    mutate(LRSR = paste(ligand, receptor, sender, receiver, sep = "_"))

  # Pivot data: calculate the average score for each LRSR-sample pair
  pca_data <- lr_scores %>%
    group_by(LRSR, sample) %>%
    summarise(mean_score = mean(get(score_col), na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = mean_score, values_fill = list(mean_score = 0))

  # Convert to matrix with LRSR as row names
  pca_matrix <- as.data.frame(pca_data)
  rownames(pca_matrix) <- pca_matrix$LRSR
  pca_matrix <- pca_matrix[,-1]

  # Transpose the matrix so that rows are samples and columns are features
  pca_matrix_t <- t(pca_matrix)

  # Standardize the data
  standardized_data <- scale(pca_matrix_t)

  # Perform PCA
  pca_result <- prcomp(standardized_data, center = TRUE, scale. = TRUE)

  # Extract the first n_components
  pcs <- as.data.frame(pca_result$x[, 1:n_components])
  pcs$sample <- rownames(pcs)

  # Merge metadata for coloring if color_by is specified
  if (!is.null(color_by)) {
    if ("sample" %in% colnames(metadata)) {
      meta_sub <- metadata %>% select(sample, !!rlang::sym(color_by))
    } else {
      meta_sub <- metadata %>%
        tibble::rownames_to_column(var = "sample") %>%
        select(sample, !!rlang::sym(color_by))
    }
    pcs <- merge(pcs, meta_sub, by = "sample", all.x = TRUE)
  }

  # Plot PCA using ggplot2
  p <- ggplot(pcs, aes(x = PC1, y = PC2)) +
    geom_point(aes_string(color = if (!is.null(color_by)) color_by else "NULL"),
               size = 3) +
    labs(x = "Principal Component 1", y = "Principal Component 2",
         title = "PCA of Ligand-Receptor Interaction Scores") +
    theme_minimal()

  if (!is.null(color_by)) {
    p <- p + scale_color_brewer(palette = "Set1")
  }

  print(p)

  return(pcs)
}
