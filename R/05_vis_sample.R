#' Generate Heatmap of Ligand-Receptor Interaction Scores
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
#'   selected_receiver = "Perivascular", selected_metadata = c("Sex", "Age_group", "IFN_type"))
#' print(p)
heatmap_sample <- function(lr_scores,
                           metadata,
                           score = c("normalized", "raw"),
                           selected_sender = NULL,
                           selected_receiver = NULL,
                           selected_metadata = NULL) {

  # Parameter validation
  score <- match.arg(score)
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
    stop("No sample information corresponding to the column was found in metadata.")
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



#' Generate PCA of Ligand-Receptor Interaction Scores
#'
#' @description
#' This function performs principal component analysis (PCA) on ligand-receptor (LR) interaction scores across samples,
#' and generates a scatter plot of the first two principal components. Optionally, sample metadata can be used to color the points.
#'
#' @param lr_scores Data frame containing LR interaction scores per sample (data frame).
#' @param metadata Data frame containing sample metadata (data frame).
#' @param selected_sender Specific sender cell type to filter, default is None (use all) (character).
#' @param selected_receiver Specific receiver cell type to filter, default is None (use all) (character).
#' @param color_by \code{metadata} column name to color points in PCA plot (character).
#' @param n_components Number of principal components to extract (default: 2) (numeric).
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
#' res <- pca_sample(lr_scores_eg, metadata_eg, selected_sender = "Cardiac",
#'   selected_receiver = "Perivascular", color_by = "IFN_type")
#' print(res$pca_plot)
pca_sample <- function(lr_scores,
                       metadata,
                       selected_sender = NULL,
                       selected_receiver = NULL,
                       color_by = NULL,
                       n_components = 2) {

  # Parameter validation
  if ("normalized_score" %in% colnames(lr_scores)) {
    score_col <- "normalized_score"
  } else if ("score" %in% colnames(lr_scores)) {
    score_col <- "score"
  } else {
    stop("lr_scores must contain either a 'score' or a 'normalized_score' column.")
  }

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
  print(p)
  return(list(plot = p, df = pca_scores))
}



#' Boxplot Comparison of Ligand-Receptor Interaction Scores Across Groups
#'
#' @description
#' Generates  a boxplot comparing LR (ligand-receptor) interaction scores across sample groups.
#' with optional significance testing (t-test or Wilcoxon).
#'
#' @param lr_scores Data frame containing LR interaction scores per sample (data frame).
#' @param metadata Data frame containing sample metadata (data frame).
#' @param ligand Ligand gene name to filter (character).
#' @param receptor Receptor gene name to filter (character).
#' @param sender Sender cell type to filter (character).
#' @param receiver Receiver cell type to filter (character).
#' @param group_by Column name in \code{metadata} to group samples (character).
#' @param score Use 'normalized' or 'raw' score (default: "normalized") (character).
#' @param test Whether to add a statistical test annotation (default: TRUE) (logical).
#' @param paired Whether to treat the comparison as paired (default: FALSE) (logical).
#' @param test_method Statistical test to use: "t.test" or "wilcox.test" (default = "wilcox.test") (character).
#' @param colors Vector of colors for groups (default: c("#5fa9d1", "#154778")).
#' @param title Custom plot title (optional).
#'
#' @return A list containing:
#' \itemize{
#'   \item plot - ggplot object of the boxplot
#'   \item df - data frame used for plotting
#' }
#'
#' @export
#'
#' @importFrom dplyr %>% filter select left_join
#' @importFrom tibble rownames_to_column
#' @importFrom ggpubr ggboxplot stat_compare_means
#' @importFrom ggplot2 ggtitle xlab ylab theme_bw theme element_text element_blank element_line
#'
#' @examples
#' # Boxplot of LR Score by group
#' data(lr_scores_eg)
#' data(metadata_eg)
#' result <- boxplot_lr_group_comparison(
#'   lr_scores_eg, metadata_eg,
#'   ligand = "TAC4", receptor = "TACR1",
#'   sender = "Perivascular", receiver = "Cardiac",
#'   group_by = "IFN_type"
#' )
#' head(result$df)
boxplot_lr_group_comparison <- function(lr_scores, metadata,
                                        ligand, receptor,
                                        sender, receiver,
                                        group_by,
                                        score = c("normalized", "raw"),
                                        test = TRUE,
                                        paired = FALSE,
                                        test_method = c("wilcox.test", "t.test"),
                                        colors = c("#5fa9d1", "#154778"),
                                        title = NULL) {

  # Parameter validation
  test_method <- match.arg(test_method)
  score <- match.arg(score)
  score_col <- ifelse(score == "normalized", "normalized_score", "score")

  # Check required columns
  if (!score_col %in% colnames(lr_scores)) {
    stop(paste("Column", score_col, "not found in lr_scores"))
  }

  if (!group_by %in% colnames(metadata)) {
    stop(paste("Grouping variable", group_by, "not found in metadata"))
  }

  # Filter for LR pair and sender-receiver
  df <- lr_scores %>%
    dplyr::filter(
      ligand == !!ligand,
      receptor == !!receptor,
      sender == !!sender,
      receiver == !!receiver
    ) %>%
    dplyr::select(sample, .data[[score_col]])
    # dplyr::select(sample, all_of(score_col))

  if (nrow(df) == 0) {
    message("No data found for the specified LR pair and sender/receiver.")
    return(NULL)
  }

  # Join metadata
  if (!"sample" %in% colnames(metadata)) {
    metadata <- tibble::rownames_to_column(metadata, var = "sample")
  }

  df <- dplyr::left_join(df, metadata[, c("sample", group_by)], by = "sample")

  # Check group levels
  if (length(unique(df[[group_by]])) < 2 && test) {
    warning("Grouping variable has less than 2 levels, skipping statistical test")
    test <- FALSE
  }

  # Default title
  if (is.null(title)) {
    title <- paste0("LR Score Comparison: ", ligand, "-", receptor,
                    " (", sender, "\u2192", receiver, ")")
  }

  # Plot
  p <- ggpubr::ggboxplot(df,
                 x = group_by, y = score_col,
                 color = group_by,
                 palette = colors,
                 add = "jitter",
                 add.params = list(shape = 16, alpha = 0.6)) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(group_by) +
    ggplot2::ylab("Interaction Score") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title.y = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(color = "black", angle = 90,
                                          vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(color = "black"),
      legend.position = "none"
    )

  # Add significance test if enabled and valid
  if (test && length(unique(df[[group_by]])) >= 2) {
    p <- p + ggpubr::stat_compare_means(
      method = test_method,
      paired = paired,
      label.x.npc = 0.5,
      label.y.npc = 0.9
    )
  }

  print(p)
  return(list(plot = p, df = df))
}



#' Dotplot of Ligand-Receptor Interaction Scores Against Continuous Group Variable
#'
#' @description
#' Creates a dotplot (scatter plot) of LR interaction scores against a continuous variable
#' with optional regression line.
#'
#' @param lr_scores Data frame containing LR interaction scores per sample (data frame).
#' @param metadata Data frame containing sample metadata (data frame).
#' @param ligand Ligand gene name to filter (character).
#' @param receptor Receptor gene name to filter (character).
#' @param sender Sender cell type to filter (character).
#' @param receiver Receiver cell type to filter (character).
#' @param group_by Continuous variable column in \code{metadata} (e.g., age, severity score) (character).
#' @param score Use 'normalized' or 'raw' score (default: "normalized") (character).
#' @param point_size Size of the points in the plot (default: 3) (numeric).
#' @param point_color Color of the points in the plot (default: "dodgerblue4").
#' @param add_regression Whether to add regression line (default: TRUE) (logical).
#' @param title Custom plot title (optional).
#'
#' @return A list containing:
#' \itemize{
#'   \item plot - ggplot object of the dotplot
#'   \item df - data frame used for plotting
#' }
#'
#' @export
#'
#' @importFrom dplyr filter select left_join
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_text ggtitle geom_smooth
#' @importFrom ggpubr stat_regline_equation stat_cor
#' @importFrom stats na.omit
#'
#' @examples
#' # Dotplot of LR Score Against Continuous Group Variable
#' data(lr_scores_eg)
#' data(metadata_eg)
#' result <- dotplot_lr_continuous_group(
#'   lr_scores_eg, metadata_eg,
#'   ligand = "TAC4", receptor = "TACR1",
#'   sender = "Perivascular", receiver = "Cardiac",
#'   group_by = "IFN_type"
#' )
#' head(result$df)
dotplot_lr_continuous_group <- function(lr_scores, metadata,
                                        ligand, receptor,
                                        sender, receiver,
                                        group_by,
                                        score = c("normalized", "raw"),
                                        point_size = 3,
                                        point_color = "dodgerblue4",
                                        add_regression = TRUE,
                                        title = NULL) {

  # Parameter validation
  score <- match.arg(score)
  score_col <- ifelse(score == "normalized", "normalized_score", "score")

  # Check required columns
  if (!score_col %in% colnames(lr_scores)) {
    stop(paste("Column", score_col, "not found in lr_scores"))
  }

  if (!group_by %in% colnames(metadata)) {
    stop(paste("Grouping variable", group_by, "not found in metadata"))
  }

  # Filter for LR pair and sender-receiver
  df <- lr_scores %>%
    dplyr::filter(
      ligand == !!ligand,
      receptor == !!receptor,
      sender == !!sender,
      receiver == !!receiver
    ) %>%
    dplyr::select(sample, .data[[score_col]])

  if (nrow(df) == 0) {
    message("No data found for the specified LR pair and sender/receiver.")
    return(NULL)
  }

  # Join metadata
  if (!"sample" %in% colnames(metadata)) {
    metadata <- tibble::rownames_to_column(metadata, var = "sample")
  }

  df <- dplyr::left_join(df, metadata[, c("sample", group_by)], by = "sample")

  # Remove NA values
  df <- stats::na.omit(df)

  if (nrow(df) == 0) {
    message("No data remaining after removing NA values.")
    return(NULL)
  }

  if (!is.numeric(df[[group_by]])) {
    warning(paste(group_by, "is not numeric. Converting to numeric."))
    df[[group_by]] <- as.numeric(df[[group_by]])
  }

  # Default title
  if (is.null(title)) {
    title <- paste0("Dot Plot of LR Score vs ", group_by, "\n",
                    ligand, "-", receptor, " (", sender, "\u2192", receiver, ")")
  }

  # plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[group_by]], y = .data[[score_col]])) +
    # ggplot2::ggplot(df, ggplot2::aes(x = !!rlang::sym(group_by), y = !!rlang::sym(score_col))) +
    ggplot2::geom_point(size = point_size, color = point_color,
                        fill = point_color, alpha = 0.7, shape = 21, stroke = 0.8) +
    ggplot2::labs(
      x = group_by,
      y = "Interaction Score",
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      # panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black")
    )

  # Add regression line if requested
  if (add_regression) {
    p <- p +
      ggplot2::geom_smooth(method = "lm", se = TRUE,
                           color = "#c53929", fill = "gray80", linetype = "dashed", linewidth = 0.8) +
      ggpubr::stat_cor(
        aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~~")),
        label.x.npc = 0.05,
        label.y.npc = 0.95,
            size = 4,
            color = "black"
        ) +
      ggpubr::stat_regline_equation(
        aes(label = paste(after_stat(rr.label), after_stat(eq.label), sep = "~~")),
        label.x.npc = 0.05,
        label.y.npc = 0.85,
        size = 4,
        color = "black")
  }

  print(p)
  return(list(plot = p, df = df))
}


