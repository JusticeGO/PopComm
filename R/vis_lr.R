#' Plot Circular Cell-Cell Interaction Network
#'
#' @description
#' Plots a circular cell-cell interaction network with curved directed edges.
#' Nodes are arranged in a circle, and edge widths and colors represent interaction strengths.
#'
#' @param filtered_lr A data frame of ligand-receptor pairs from prior analysis (e.g., output of `filter_lr_all`),
#'        containing at least the columns "sender", "receiver", and "cor".
#' @param edge_width Determines edge weights, either "cor" (correlation) or "count" (interaction count) (default: "count").
#' @param node_colors Named vector mapping cell types to colors. If NULL, uses default palette.
#' @param show_self_interactions Logical indicating whether to display self-interactions (default: TRUE).
#'
#' @return A ggplot object representing the network plot.
#'
#' @export
#'
#' @importFrom dplyr %>% group_by summarise
#' @importFrom igraph graph_from_data_frame V E
#' @import ggraph
#' @import ggplot2
#'
#' @examples
#' # Plot Circular Cell-Cell Interaction Network
#'

circle_plot <- function(filtered_lr,
                        edge_width = "count",
                        node_colors = NULL,
                        show_self_interactions = TRUE) {

  # Step 1: Set up the function framework and parameter check
  if (!edge_width %in% c("cor", "count")) {
    stop("edge_width must be either 'cor' or 'count'")
  }
  required_cols <- c("sender", "receiver", "cor")
  if (!all(required_cols %in% colnames(filtered_lr))) {
    stop("filtered_lr must contain columns: sender, receiver, cor")
  }


  # Step 2: Create the network graph and calculate the edge weights
  # Create edge data
  edges <- filtered_lr %>%
    group_by(sender, receiver) %>%
    summarise(
      weight = if (edge_width == "cor") sum(cor) else n(),
      .groups = "drop"
    )

  # Create node data
  nodes <- data.frame(
    name = unique(c(edges$sender, edges$receiver))
  )

  # Create igraph object
  g <- graph_from_data_frame(edges, vertices = nodes)


  # Step 3: Calculate node size and color
  # 计算节点权重
  node_weights <- sapply(V(g), function(v) {
    sum(edges$weight[edges$sender == v$name | edges$receiver == v$name])
  })

  # 标准化节点大小
  node_sizes <- scales::rescale(node_weights, to = c(3, 20))

  # 设置节点颜色
  if (is.null(node_colors)) {
    node_colors <- scales::hue_pal()(nrow(nodes))
    names(node_colors) <- nodes$name
  }


  # Step 4: Create a basic network diagram
  # 创建圆形布局
  layout <- create_layout(g, layout = "linear", circular = TRUE)

  # 初始化绘图
  p <- ggraph(layout) +
    theme_void()


  # Step 5: Add edges and arrows
  # 添加边
  p <- p +
    geom_edge_arc(
      aes(
        width = weight,
        color = weight,
        filter = if (!show_self_interactions) sender != receiver
      ),
      arrow = arrow(
        type = "closed",
        length = unit(3, "mm")
      ),
      strength = 0.3
    ) +
    scale_edge_width_continuous(range = c(0.5, 3)) +
    scale_edge_color_gradientn(
      colors = RColorBrewer::brewer.pal(9, "Cool")
    )


  # Step 6: Add nodes and labels
  p <- p +
    geom_node_point(
      aes(size = node_sizes, fill = name),
      shape = 21,
      stroke = 0.5
    ) +
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      size = 3
    ) +
    scale_size_identity() +
    scale_fill_manual(values = node_colors)


  # Step 7: Handling self crossing edges
  if (show_self_interactions) {
    self_edges <- edges %>% filter(sender == receiver)
    if (nrow(self_edges) > 0) {
      p <- p +
        geom_edge_loop(
          aes(
            width = weight,
            color = weight,
            direction = 45,
            span = 90
          )
        )
    }
  }


  # Step 8: Add legends and adjust themes
  p <- p +
    guides(
      edge_width = guide_legend(title = if (edge_width == "cor") "Correlation" else "Count"),
      edge_color = "none",
      fill = "none",
      size = "none"
    ) +
    theme(
      legend.position = "bottom",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )

  return(p)
}






