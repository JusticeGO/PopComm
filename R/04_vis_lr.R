#' Plot Circular Ligand-Receptor Interaction Network
#'
#' @description
#' Plots a circular ligand-receptor (LR) interaction network with curved directed edges.
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
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr %>% group_by summarise left_join rename
#' @importFrom igraph graph_from_data_frame V layout_in_circle
#' @importFrom ggraph ggraph geom_edge_arc geom_node_point geom_node_text scale_edge_color_gradient2 scale_edge_width create_layout
#' @importFrom ggplot2 aes theme_void ggtitle geom_curve arrow
#' @importFrom grid unit
#'
#' @examples
#' # Plot Circular Cell-Cell Interaction Network
#' data(filtered_lr_eg)
#' p <- circle_plot(filtered_lr_eg, edge_width = "count", show_self_interactions = FALSE)
circle_plot <- function(filtered_lr,
                        edge_width = "count",
                        node_colors = NULL,
                        show_self_interactions = FALSE) {

  # Parameter validation
  if (!edge_width %in% c("cor", "count"))
    stop("edge_width must be 'cor' or 'count'")
  if (!all(c("sender", "receiver", "cor") %in% colnames(filtered_lr)))
    stop("filtered_lr must contain sender, receiver and cor columns")

  # Add nodes: extract unique cell types from sender and receiver
  cell_types <- unique(c(filtered_lr$sender, filtered_lr$receiver))

  # Set node colors: auto generate if not provided
  if (is.null(node_colors)) {
    node_colors <- hue_pal()(length(cell_types))
    names(node_colors) <- cell_types
  } else {
    missing_colors <- setdiff(cell_types, names(node_colors))
    if (length(missing_colors) > 0) {
      warning("Some nodes do not provide colors, default colors will be assigned.")
      extra_colors <- hue_pal()(length(missing_colors))
      names(extra_colors) <- missing_colors
      node_colors <- c(node_colors, extra_colors)
    }
  }

  # Calculate edge weights: use 'cor' if specified, else default to 1
  filtered_lr$weight <- if (edge_width == "cor") filtered_lr$cor else 1

  # Aggregate weights for the same sender-receiver pair
  edge_df <- filtered_lr %>%
    group_by(sender, receiver) %>%
    summarise(weight = sum(weight), .groups = "drop")

  # Determine the data range for normalized calculations based on the show_self_interactions parameter
  if (show_self_interactions) {
    norm_data <- edge_df  # 全部边
  } else {
    norm_data <- edge_df %>% filter(sender != receiver)
  }

  edge_df$width_norm <- rescale(edge_df$weight,
                                to = c(0.5, 5),
                                from = range(norm_data$weight))

  non_self_edges <- edge_df %>% filter(sender != receiver)
  self_edges <- edge_df %>% filter(sender == receiver)

  # Build a directed graph using igraph
  g <- graph_from_data_frame(non_self_edges, vertices = data.frame(name = cell_types), directed = TRUE)

  p <- g

  # p <- ggraph(g, layout = 'circle') +
  #   geom_edge_arc(aes(edge_width = weight, edge_color = weight),
  #                 arrow = arrow(length = unit(3, 'mm'), type = "closed"),
  #                 curvature = 0.3, show.legend = TRUE) +
  #   geom_node_point(aes(x = x, y = y),
  #                   size = {
  #                     node_weights <- sapply(cell_types, function(n) {
  #                       sum(edge_df$weight[edge_df$sender == n | edge_df$receiver == n])
  #                     })
  #                     rescale(node_weights, to = c(300, 2300)) / 100
  #                   },
  #                   color = "black",
  #                   fill = node_colors[V(g)$name],
  #                   shape = 21, stroke = 0.5) +
  #   geom_node_text(aes(x = x, y = y, label = name),
  #                  size = 4, color = "black",
  #                  fontface = "bold", show.legend = FALSE) +
  #   scale_edge_width(range = c(0.5, 5)) +
  #   scale_edge_color_gradient2(low = "#377EB8", mid = "grey80", high = "#E41A1C",
  #                              midpoint = mean(norm_data$weight),
  #                              limits = range(norm_data$weight)) +
  #   theme_void() +
  #   ggtitle("Circular Cell-Cell Interaction Network") +
  #   coord_fixed()
  #
  # # 如果需要显示自交互边，则单独添加
  # if (show_self_interactions && nrow(self_edges) > 0) {
  #   # 使用 create_layout 获取节点坐标
  #   layout_df <- create_layout(g, layout = "circle")
  #   # 将自交互边与对应节点坐标合并
  #   self_edges <- self_edges %>%
  #     left_join(layout_df, by = c("sender" = "name")) %>%
  #     rename(x = x, y = y) %>%
  #     mutate(x_start_new = x,
  #            y_start_new = y,
  #            x_end_new = x + 0.1,  # 根据需要调整偏移量
  #            y_end_new = y)
  #
  #   p <- p + geom_curve(data = self_edges,
  #                       aes(x = x_start_new, y = y_start_new,
  #                           xend = x_end_new, yend = y_end_new,
  #                           color = weight, size = width_norm),
  #                       curvature = 0.5,
  #                       arrow = arrow(length = unit(3, "mm")))
  # }

  return(p)
}



#' Create Ligand-Receptor Interaction Dot Plot
#'
#' @description
#' Generates a dot plot to visualize ligand-receptor (LR) interaction. Dot sizes are scaled by the correlation
#' coefficient and dot colors represent -log10(adjust.p). The function supports plotting the top interactions
#' per sender-receiver pair or user-specified ligand-receptor pairs.
#'
#' @param filtered_lr A data frame containing ligand-receptor interaction data.
#' @param top_n Integer specifying the number of top interactions to select per sender-receiver pair (numeric, default: 5).
#' @param axis Character indicating the configuration of rows and columns in the plot.
#'        Options: "LR-SR" (default, rows = ligand-receptor pairs, columns = sender-receiver pairs) or "SR-LR".
#' @param type_scale Character indicating the scaling method for the plot.
#'        Options: "size" (default, uses `scale_size()` for point scaling) or "radius" (uses `scale_radius()` for point scaling).
#' @param selected_LR Optional character vector of ligand-receptor pair identifiers (e.g., c("TIMP1_CD63", "DSCAM_DCC")).
#'        If NULL, the top_n interactions per sender-receiver pair are used.
#'
#' @return A ggplot object representing the dot plot.
#'
#' @export
#'
#' @importFrom dplyr %>% mutate filter group_by arrange desc slice_head ungroup
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn theme_minimal theme element_text element_line element_blank labs scale_size scale_radius coord_fixed
#' @importFrom rlang .data
#'
#' @examples
#' # Plot LR Interaction Dot Plot
#' data(filtered_lr_eg)
#' p <- dot_plot(filtered_lr_eg, axis = "LR-SR", type_scale = "size")
#' print(p)
dot_plot <- function(filtered_lr,
                     top_n = 5,
                     axis = "LR-SR",
                     type_scale = "size",
                     selected_LR = NULL) {

  # Parameter validation
  if (!all(c("ligand", "receptor", "sender", "receiver", "cor", "adjust.p") %in% colnames(filtered_lr)))
    stop("filtered_lr must contain columns: ligand, receptor, sender, receiver, cor, and adjust.p")

  if (!is.numeric(top_n) || length(top_n) != 1 || top_n <= 0)
    stop("top_n must be a single positive numeric value")

  if (!axis %in% c("LR-SR", "SR-LR"))
    stop("axis must be either 'LR-SR' or 'SR-LR'")

  if (!type_scale %in% c("size", "radius"))
    stop("type_scale must be either 'size' or 'radius'")

  if (!is.null(selected_LR) && !is.character(selected_LR))
    stop("selected_LR must be a character vector if provided")

  # Construct ligand-receptor (LR) and sender-receiver (SR) identifiers
  filtered_lr <- filtered_lr %>%
    mutate(LR_pair = paste(.data[["ligand"]], .data[["receptor"]], sep = "_"),
           SR_pair = paste(.data[["sender"]], .data[["receiver"]], sep = "_"))

  # adjust.p capped at 1e-20
  # Calculate -log10(adjust.p) for color mapping
  filtered_lr$adjust.p <- pmax(filtered_lr$adjust.p, 1e-20)
  filtered_lr$p.scaled <- -log10(filtered_lr$adjust.p)

  # Filter data: if selected_LR is provided, keep only those rows; otherwise, choose top_n per SR_pair
  if (!is.null(selected_LR)) {
    filtered_lr <- filtered_lr %>% filter(.data[["LR_pair"]] %in% selected_LR)
  } else {
    filtered_lr <- filtered_lr %>%
      group_by(.data[["SR_pair"]]) %>%
      arrange(desc(.data[["cor"]])) %>%
      slice_head(n = top_n) %>%
      ungroup()
  }

  if (axis == "LR-SR") {
    row_label <- "LR_pair"
    col_label <- "SR_pair"
  } else if (axis == "SR-LR") {
    row_label <- "SR_pair"
    col_label <- "LR_pair"
  } else {
    stop("Invalid axis value. Choose 'LR-SR' or 'SR-LR'.")
  }

  scale_param <- if (type_scale == "size") {
    scale_size(name = "Correlation (cor)")
  } else if (type_scale == "radius") {
    scale_radius(name = "Correlation (cor)")
  } else {
    stop("Invalid type_scale value. Choose 'size' or 'radius'.")
  }

  # Convert rows and columns to factors to preserve the ordering
  filtered_lr[[row_label]] <- factor(filtered_lr[[row_label]], levels = unique(filtered_lr[[row_label]]))
  filtered_lr[[col_label]] <- factor(filtered_lr[[col_label]], levels = unique(filtered_lr[[col_label]]))

  # Create the dot plot using ggplot2
  p <- ggplot(filtered_lr, aes(x = .data[[col_label]], y = .data[[row_label]], size = .data[["cor"]], color = .data[["p.scaled"]])) +
    geom_point() +
    scale_color_gradientn(
      colors = c("#c8ebf6", "#5fa9d1", "#154778"),
      name = "Adjusted P-value (-log10)",
      limits = range(filtered_lr$p.scaled, na.rm = TRUE)) +
    scale_param +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      plot.title = element_text(size = 12, hjust = 0.5, color = "black"),
      legend.position = "right",
      axis.line = element_line(),
      axis.ticks = element_line(),
      panel.grid.major = element_line(color = "gray95"),
      panel.grid.minor = element_blank()
    ) +
    labs(title = if (is.null(selected_LR)) {
      paste0("Ligand-Receptor Interaction Dot Plot\n(Top ", top_n, " per Sender-Receiver)")
    } else {
      "Ligand-Receptor Interaction Dot Plot\n(Selected LR)"
    },
    x = "", y = "") +
    coord_fixed(ratio = 1)

  return(p)
}
