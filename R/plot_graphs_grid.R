# R/plot_graphs_grid.R

#' Plot Multiple Graphs in a Grid Layout
#'
#' This function plots multiple graphs (adjacency matrices) in a grid layout
#' within the same image.
#'
#' @param graph_list A list of square adjacency matrices representing graphs.
#' @param communities_list A list of community vectors for each graph.
#'   If NULL, communities will be detected for each graph. Default is NULL.
#' @param nrow Number of rows in the grid layout. Default is determined automatically.
#' @param ncol Number of columns in the grid layout. Default is determined automatically.
#' @param layout A layout function or list of layouts for arranging the graph nodes.
#'   Default is \code{layout_with_fr}.
#' @param vertex_size Numeric value indicating the size of the vertices. Default is 5.
#' @param vertex_label Logical indicating whether to display vertex labels. Default is FALSE.
#' @param edge_width Numeric value indicating the width of the edges. Default is 1.
#' @param titles A character vector of titles for each graph. Default is NULL.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A grid of plots is produced. No return value.
#' @export
#'
#' @examples
#' # Generate a list of graphs
#' graphs <- list(
#'   generate_community_graph(n_nodes = 30, n_communities = 3, intra_prob = 0.8, inter_prob = 0.2),
#'   generate_community_graph(n_nodes = 30, n_communities = 3, intra_prob = 0.7, inter_prob = 0.3)
#' )
#'
#' # Plot the graphs in a grid
#' plot_graphs_grid(graphs)
plot_graphs_grid <- function(graph_list, communities_list = NULL, nrow = NULL, ncol = NULL,
                             layout = igraph::layout_with_fr, vertex_size = 5, vertex_label = FALSE,
                             edge_width = 1, titles = NULL, ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package is required for this function. Please install it using install.packages('igraph').")
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("The 'gridExtra' package is required for this function. Please install it using install.packages('gridExtra').")
  }

  num_graphs <- length(graph_list)

  if (is.null(nrow) && is.null(ncol)) {
    ncol <- ceiling(sqrt(num_graphs))
    nrow <- ceiling(num_graphs / ncol)
  } else if (is.null(nrow)) {
    nrow <- ceiling(num_graphs / ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(num_graphs / nrow)
  }

  if (is.null(communities_list)) {
    communities_list <- vector("list", num_graphs)
  }

  if (length(communities_list) != num_graphs) {
    stop("Length of communities_list must match length of graph_list.")
  }

  # Prepare individual plots
  plot_list <- vector("list", num_graphs)
  for (i in seq_along(graph_list)) {
    G <- graph_list[[i]]

    if (!is.matrix(G) || nrow(G) != ncol(G)) {
      stop("Each graph in graph_list must be a square adjacency matrix.")
    }

    graph <- igraph::graph_from_adjacency_matrix(G, mode = "undirected", diag = FALSE)

    # Detect communities if not provided
    if (is.null(communities_list[[i]])) {
      community_detection <- igraph::cluster_louvain(graph)
      communities <- igraph::membership(community_detection)
    } else {
      communities <- communities_list[[i]]
    }

    # Assign colors to communities
    unique_communities <- sort(unique(communities))
    color_palette <- rainbow(length(unique_communities))
    community_colors <- color_palette[match(communities, unique_communities)]

    # Set vertex attributes
    igraph::V(graph)$color <- community_colors
    igraph::V(graph)$size <- vertex_size
    igraph::V(graph)$label <- if (vertex_label) igraph::V(graph)$name else NA

    # Set edge attributes
    igraph::E(graph)$width <- edge_width

    # Generate the plot
    p <- function() {
      plot(graph,
           layout = layout,
           vertex.color = igraph::V(graph)$color,
           vertex.size = igraph::V(graph)$size,
           vertex.label = igraph::V(graph)$label,
           edge.width = igraph::E(graph)$width,
           main = if (!is.null(titles)) titles[i] else NULL,
           ...)
    }

    # Capture the plot as a grob
    plot_list[[i]] <- ggplotify::as.grob(p)
  }

  # Arrange plots in a grid
  gridExtra::grid.arrange(grobs = plot_list, nrow = nrow, ncol = ncol)
}
