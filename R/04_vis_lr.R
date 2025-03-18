#' Plot Circular Cell-Cell Interaction Network
#'
#' @description
#' Plots a circular cell-cell interaction network with curved directed edges.
#' Nodes are arranged in a circle, and edge widths and colors represent interaction strengths.
#'
#' @param filtered_lr A data frame of ligand-receptor pairs from prior analysis (e.g., output of `filter_lr_all`),
#'        containing at least the columns "sender", "receiver", and "cor".
#' @param edge_width Determines edge weights, either "cor" (correlation) or "count" (interaction count) (default: "count").
#' @param node_colors Named vector mapping cell types to colors. Example: c("Tcell" = "#E41A1C", "Macrophage" = "#377EB8"). If NULL, uses default palette.
#' @param show_self_interactions Logical indicating whether to display self-interactions (default: TRUE).
#'
#' @return A ggplot object representing the network plot.
#'
#' @export
#'
#' @importFrom scales hue_pal rescale
#' @importFrom dplyr %>% group_by summarise left_join rename
#' @importFrom igraph graph_from_data_frame V layout_in_circle
#' @importFrom ggraph ggraph geom_edge_arc
#' @importFrom ggplot2 geom_node_point geom_node_text scale_edge_width scale_edge_color_gradient2 theme_void ggtitle geom_curve arrow
#' @importFrom grid unit
#'
#' @examples
#' # Plot Circular Cell-Cell Interaction Network
#'
circle_plot <- function(filtered_lr,
                        edge_width = "count",
                        node_colors = NULL,
                        show_self_interactions = TRUE) {

  # Parameter validation
  if (!edge_width %in% c("cor", "count")) {
    stop("edge_width must be either 'cor' or 'count'")
  }
  required_cols <- c("sender", "receiver", "cor")
  if (!all(required_cols %in% colnames(filtered_lr))) {
    stop("filtered_lr must contain columns: sender, receiver, cor")
  }

  # Add nodes: extract unique cell types from sender and receiver
  cell_types <- unique(c(filtered_lr$sender, filtered_lr$receiver))

  # Set node colors: auto generate if not provided
  if (is.null(node_colors)) {
    default_colors <- hue_pal()(length(cell_types))
    names(default_colors) <- cell_types
    node_colors <- default_colors
  }

  # Calculate edge weights: use 'cor' if specified, else default to 1
  filtered_lr$weight <- if (edge_width == "cor") filtered_lr$cor else 1

  # Aggregate weights for the same sender-receiver pair
  edge_df <- filtered_lr %>%
    group_by(sender, receiver) %>%
    summarise(weight = sum(weight), .groups = "drop")

  # Build a directed graph using igraph
  g <- graph_from_data_frame(edge_df, vertices = cell_types, directed = TRUE)

  # Calculate node weights based on connected edge weights
  node_weights <- sapply(V(g)$name, function(n) {
    sum(edge_df$weight[edge_df$sender == n | edge_df$receiver == n])
  })
  # Normalize node sizes to range 300 to 2300
  node_sizes <- rescale(node_weights, to = c(300, 2300))

  # Compute circular layout for nodes
  layout <- layout_in_circle(g)
  layout_df <- as.data.frame(layout)
  colnames(layout_df) <- c("x", "y")
  layout_df$name <- V(g)$name

  # Merge node positions into edge data for plotting
  edge_df <- edge_df %>%
    left_join(layout_df, by = c("sender" = "name")) %>%
    rename(x_start = x, y_start = y) %>%
    left_join(layout_df, by = c("receiver" = "name")) %>%
    rename(x_end = x, y_end = y)

  # Normalize edge weights for edge width (range 0.5 to 5)
  edge_df$width_norm <- rescale(edge_df$weight, to = c(0.5, 5))

  # Define color limits for edge color gradient
  color_limits <- range(edge_df$weight)

  # Create the network plot using ggraph and ggplot2
  p <- ggraph(g, layout = 'circle') +
    geom_edge_arc(aes(edge_width = weight, edge_color = weight),
                  arrow = arrow(length = unit(3, 'mm'), type = "closed"),
                  curvature = 0.3, show.legend = TRUE) +
    geom_node_point(aes(x = x, y = y),
                    size = node_sizes / 100,
                    color = "black",
                    fill = node_colors[V(g)$name],
                    shape = 21, stroke = 0.5) +
    geom_node_text(aes(x = x, y = y, label = name),
                   repel = TRUE, bg.color = "white", size = 4) +
    scale_edge_width(range = c(0.5, 5)) +
    scale_edge_color_gradient2(low = "blue", mid = "white", high = "red",
                               midpoint = mean(color_limits), limits = color_limits) +
    theme_void() +
    ggtitle("Circular Cell-Cell Interaction Network")

  # Add self-interaction curves if required
  if (show_self_interactions) {
    self_edges <- edge_df[edge_df$sender == edge_df$receiver, ]
    if (nrow(self_edges) > 0) {
      p <- p + geom_curve(data = self_edges,
                          aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = weight),
                          curvature = 0.5, size = self_edges$width_norm,
                          arrow = arrow(length = unit(3, "mm")))
    }
  }

  return(p)
}



#' Create Ligand-Receptor Interaction Dot Plot
#'
#' @description
#' Generates a dot plot to visualize ligand-receptor interactions. Dot sizes are scaled by the correlation
#' coefficient and dot colors represent -log10(p-value). The function supports plotting the top interactions
#' per sender-receiver pair or user-specified ligand-receptor pairs.
#'
#' @param filtered_lr A data frame containing columns "ligand", "receptor", "sender", "receiver", "cor", and "p_val".
#' @param top_n Integer. Number of top interactions to select per sender-receiver pair (default is 5).
#' @param axis Character. Determines the configuration of rows and columns in the plot.
#'   Options: "LR-SR" (default, rows = ligand-receptor, columns = sender-receiver) or "SR-LR".
#' @param size_scale Numeric. Scaling factor for dot sizes (default is 1500).
#' @param min_size Numeric. Minimum dot size to ensure visibility (default is 5).
#' @param selected_LR Optional character vector of ligand-receptor pairs (e.g., c("EGF_EGFR", "TGFB1_TGFBR1")).
#'   If NULL, the top_n interactions per sender-receiver pair are used.
#' @param figsize Numeric vector of length 2 specifying figure width and height (default is c(12, 8)).
#'
#' @return A ggplot object representing the dot plot.
#'
#' @export
#'
#' @importFrom dplyr %>% mutate group_by arrange slice_head ungroup filter
#' @importFrom ggplot2 ggplot aes_string geom_point scale_color_gradient2 scale_size_continuous theme_bw theme element_text labs
#' @importFrom scales rescale
dot_plot <- function(filtered_lr, top_n = 5, axis = "LR-SR", size_scale = 1500, min_size = 5, selected_LR = NULL, figsize = c(12, 8)) {
  # Ensure required packages are available
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Package 'scales' is required.")

  # Construct ligand-receptor (LR) and sender-receiver (SR) identifiers
  filtered_lr <- filtered_lr %>%
    mutate(LR_pair = paste(ligand, receptor, sep = "_"),
           SR_pair = paste(sender, receiver, sep = "_"))

  # Calculate -log10(p_val) for color mapping
  filtered_lr <- filtered_lr %>%
    mutate(log_pval = -log10(p_val))

  # Normalize correlation (cor) to compute dot sizes and apply a square transformation for contrast
  min_cor <- min(filtered_lr$cor, na.rm = TRUE)
  max_cor <- max(filtered_lr$cor, na.rm = TRUE)
  filtered_lr <- filtered_lr %>%
    mutate(cor_size = ((cor - min_cor) / (max_cor - min_cor + 1e-6))^2 * size_scale + min_size)

  # Filter data: if selected_LR is provided, keep only those rows; otherwise, choose top_n per SR_pair
  if (!is.null(selected_LR)) {
    filtered_lr <- filtered_lr %>% filter(LR_pair %in% selected_LR)
  } else {
    filtered_lr <- filtered_lr %>%
      group_by(SR_pair) %>%
      arrange(desc(cor)) %>%
      slice_head(n = top_n) %>%
      ungroup()
  }

  # Decide row and column labels based on the axis parameter
  if (axis == "LR-SR") {
    row_label <- "LR_pair"
    col_label <- "SR_pair"
  } else if (axis == "SR-LR") {
    row_label <- "SR_pair"
    col_label <- "LR_pair"
  } else {
    stop("Invalid axis value. Choose 'LR-SR' or 'SR-LR'.")
  }

  # Convert rows and columns to factors to preserve the ordering
  filtered_lr[[row_label]] <- factor(filtered_lr[[row_label]], levels = unique(filtered_lr[[row_label]]))
  filtered_lr[[col_label]] <- factor(filtered_lr[[col_label]], levels = unique(filtered_lr[[col_label]]))

  # Create the dot plot using ggplot2
  p <- ggplot(filtered_lr, aes_string(x = col_label, y = row_label)) +
    geom_point(aes(size = cor_size, color = log_pval),
               shape = 21, stroke = 0.5, fill = "white", alpha = 0.8) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = median(filtered_lr$log_pval)) +
    scale_size_continuous(range = c(3, 15)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = if (is.null(selected_LR)) {
      paste0("Ligand-Receptor Interaction Dot Plot (Top ", top_n, " per Sender-Receiver)")
    } else {
      "Ligand-Receptor Interaction Dot Plot (Selected LR)"
    },
    x = "", y = "", color = expression(-log[10](p)), size = "Correlation")

  return(p)
}
