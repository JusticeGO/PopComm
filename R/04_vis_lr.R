#' Plot Circular Ligand-Receptor Interaction Network
#'
#' @description
#' Plots a circular ligand-receptor (LR) interaction network with curved directed edges.
#' Nodes are arranged in a circle, and edge widths and colors represent interaction strengths.
#'
#' @param filtered_lr A data frame of ligand-receptor pairs from prior analysis (e.g., output of `filter_lr_all`),
#'        containing at least the columns "sender", "receiver", and "cor".
#' @param edge_width Determines edge weights, either "cor" (correlation) or "count" (interaction count) (default: "count").
#' @param node_colors Named vector mapping cell types to colors. Example: c("Cardiac" = "#E41A1C", "Fibroblast" = "#377EB8"). If NULL, uses default palette.
#' @param show_self_interactions Logical indicating whether to display self-interactions (logical, default: TRUE).
#' @param cutoff Minimum edge weight to display (numeric, default: 0).
#'
#' @return A ggplot object representing the network plot.
#'
#' @export
#'
#' @importFrom scales hue_pal rescale alpha col_numeric
#' @importFrom dplyr %>% count rename left_join
#' @importFrom stats aggregate
#' @importFrom igraph graph_from_adjacency_matrix V E layout_in_circle ends delete_edges plot.igraph
#' @importFrom reshape2 dcast
#' @importFrom rlang .data
#' @importFrom grDevices recordPlot
#' @importFrom ggplot2 ggplot ggplotGrob geom_tile scale_fill_gradientn guides guide_colorbar
#' @importFrom cowplot plot_grid
#' @importFrom utils packageVersion
#'
#' @examples
#' # Plot Circular Cell-Cell Interaction Network
#' data(filtered_lr_eg)
#'
#' p <- circle_plot(
#'   filtered_lr = filtered_lr_eg,
#'   edge_width = "count",
#'   show_self_interactions = TRUE
#' )
#'
#' print(p)
circle_plot <- function(filtered_lr,
                        edge_width = c("count", "cor"),
                        node_colors = NULL,
                        show_self_interactions = TRUE,
                        cutoff = 0) {

  # Parameter validation
  edge_width <- match.arg(edge_width)
  if (!all(c("sender", "receiver", "cor") %in% colnames(filtered_lr))) {
    stop("Input must contain 'sender', 'receiver', and 'cor' columns.")
  }

    # Build interaction matrix
  filtered_lr$weight <- if (edge_width == "cor") filtered_lr$cor else 1

  net <- filtered_lr %>%
    dplyr::count(.data[["sender"]], .data[["receiver"]]) %>%
    dplyr::rename(count = .data[["n"]])

  # Add correlation if needed
  if (edge_width == "cor") {
    cor_data <- stats::aggregate(weight ~ sender + receiver, data = filtered_lr, FUN = sum)
    # net <- merge(net, cor_data, by = c("sender", "receiver"), all.x = TRUE)
    net <- dplyr::left_join(net, cor_data, by = c("sender", "receiver"))
    net$count <- net$weight
  }

  # Filter edges by cutoff
  net <- net[net$count >= cutoff, ]

  mat <- reshape2::dcast(net, sender ~ receiver, value.var = "count", fill = 0)
  rownames(mat) <- mat[,1]
  mat <- mat[,-1]

  g <- igraph::graph_from_adjacency_matrix(as.matrix(mat), mode = "directed", weighted = TRUE)

  # Node sizes by degree
  vertex.weight <- rowSums(mat) + colSums(mat)
  vertex.size <- scales::rescale(vertex.weight, to = c(10, 30))

  # Node colors
  cell_types <- igraph::V(g)$name
  if (is.null(node_colors)) {
    color.use <- setNames(scales::alpha(scales::hue_pal()(length(cell_types)), 0.8), cell_types)
  } else {
    missing <- setdiff(cell_types, names(node_colors))
    if (length(missing) > 0) {
      auto_colors <- setNames(scales::hue_pal()(length(missing)), missing)
      node_colors <- c(node_colors, auto_colors)
    }
    color.use <- scales::alpha(node_colors[cell_types], 0.8)
  }

  # Plotting layout
  layout <- igraph::layout_in_circle(g)

  # # Self-loops
  # if (show_self_interactions) {
  #   loop.idx <- which(igraph::ends(g, igraph::E(g))[,1] == igraph::ends(g, igraph::E(g))[,2])
  #   if (length(loop.idx) > 0) {
  #     loop.nodes <- igraph::ends(g, igraph::E(g))[loop.idx, 1]
  #     loop.idx.vertices <- match(loop.nodes, igraph::V(g)$name)
  #     x <- layout[loop.idx.vertices, 1]
  #     y <- layout[loop.idx.vertices, 2]
  #     angles <- ifelse(x > 0, -atan(y / x), pi - atan(y / x))
  #     igraph::E(g)$loop.angle <- NA
  #     igraph::E(g)$loop.angle[loop.idx] <- angles
  #   }
  # } else {
  #   g <- igraph::delete_edges(g, which(igraph::ends(g, igraph::E(g))[,1] == igraph::ends(g, igraph::E(g))[,2]))
  # }

  # Handle Self-loops based on igraph version
  if (utils::packageVersion("igraph") <= "2.1.4") {
    # For igraph 2.1.4 and earlier: manual loop angle calculation
    if (show_self_interactions) {
      loop.idx <- which(igraph::ends(g, igraph::E(g))[,1] == igraph::ends(g, igraph::E(g))[,2])
      if (length(loop.idx) > 0) {
        loop.nodes <- igraph::ends(g, igraph::E(g))[loop.idx, 1]
        loop.idx.vertices <- match(loop.nodes, igraph::V(g)$name)
        x <- layout[loop.idx.vertices, 1]
        y <- layout[loop.idx.vertices, 2]
        angles <- ifelse(x > 0, -atan(y / x), pi - atan(y / x))
        igraph::E(g)$loop.angle <- NA
        igraph::E(g)$loop.angle[loop.idx] <- angles
      }
    } else {
      g <- igraph::delete_edges(g, which(igraph::ends(g, igraph::E(g))[,1] == igraph::ends(g, igraph::E(g))[,2]))
    }

  } else {
    # For igraph 2.1.5 and later: automatic loop angle handling
    if (!show_self_interactions) {
      loop.idx <- which(igraph::ends(g, igraph::E(g))[,1] == igraph::ends(g, igraph::E(g))[,2])
      if (length(loop.idx) > 0) {
        g <- igraph::delete_edges(g, loop.idx)
      }
    }
  }

  edge.width <- scales::rescale(igraph::E(g)$weight, to = c(1, 10))

  w <- igraph::E(g)$weight
  igraph::E(g)$color <- scales::col_numeric(
    # palette = c("#377EB8", "#E41A1C"),
    palette = c("#377EB8", "grey80", "#E41A1C"),
    domain = range(w)
  )(w)

  # Title selection
  plot_title <- ifelse(edge_width == "count", "Interaction Count", "Interaction Correlation")

  # Conditional margin setup
  if (show_self_interactions) {
    margin_value <- 0.2  # Add margin only when self-interactions are shown
  } else {
    margin_value <- 0  # No extra margin when no self-interactions
  }

  plot(g,
       layout = layout,
       edge.arrow.size = 0.6,
       edge.arrow.width = 1.2,
       edge.curved = 0.2,
       edge.width = edge.width,
       edge.color = igraph::E(g)$color,
       vertex.color = color.use,
       vertex.label.family = "sans",
       vertex.label.color = "black",
       vertex.label.cex = 1.2,
       vertex.size = vertex.size,
       vertex.frame.color = "black",
       vertex.frame.width = 1.5,
       margin = margin_value,
       main = plot_title
  )

  p_net <- grDevices::recordPlot()

  # legend
  w_range <- range(w)
  legend_df <- data.frame(
    x = 1,
    y = seq(w_range[1], w_range[2], length.out = 100),
    val = seq(w_range[1], w_range[2], length.out = 100)
  )

  legend_plot <- ggplot2::ggplot(legend_df, aes(x = x, y = y, fill = .data[["val"]])) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colours = c("#377EB8", "grey80", "#E41A1C"),
      name = plot_title,
      limits = w_range,
      breaks = pretty(w_range, n = 6)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      legend.title.position = "left",
      legend.title = ggplot2::element_text(size = 12, angle = 90, hjust = 0.5),
      legend.text = ggplot2::element_text(size = 10)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        barwidth = 0.8,
        barheight = 10
      )
    )

  get_legend <- function(plot, legend = NULL) {
    gt <- ggplot2::ggplotGrob(plot)
    pattern <- "guide-box"
    if (!is.null(legend)) {
      pattern <- paste0(pattern, "-", legend)
    }
    indices <- grep(pattern, gt$layout$name)

    not_empty <- !vapply(
      gt$grobs[indices],
      inherits, what = "zeroGrob",
      FUN.VALUE = logical(1)
    )
    indices <- indices[not_empty]

    if (length(indices) > 0) {
      return(gt$grobs[[indices[1]]])
    }
    return(NULL)
  }

  legend_only <- get_legend(legend_plot)

  combined_plot <- cowplot::plot_grid(
    p_net, legend_only,
    rel_widths = c(6, 1),
    nrow = 1
  )

  return(combined_plot)
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
#'
#' p <- dot_plot(
#'   filtered_lr = filtered_lr_eg,
#'   top_n = 3,
#'   axis = "LR-SR",
#'   type_scale = "size",
#'   )
#'
#' print(p)
dot_plot <- function(filtered_lr,
                     top_n = 5,
                     axis = c("LR-SR", "SR-LR"),
                     type_scale = c("size", "radius"),
                     selected_LR = NULL) {

  # Parameter validation
  axis <- match.arg(axis)
  type_scale <- match.arg(type_scale)
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
