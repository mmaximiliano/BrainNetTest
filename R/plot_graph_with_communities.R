# R/plot_graph_with_communities.R

#' Plot a Graph with Highlighted Communities
#'
#' This function plots a graph (adjacency matrix), highlighting the communities
#' in different colors and displaying relevant information in a legend.
#'
#' @param G A square adjacency matrix representing the graph.
#' @param communities A vector indicating the community membership of each node.
#'   If NULL, communities will be detected using a community detection algorithm (Louvain method).
#'   Default is NULL.
#' @param layout A layout function or matrix for arranging the graph nodes.
#'   Default is \code{layout_with_fr} (Fruchterman-Reingold force-directed algorithm).
#' @param vertex_size Numeric value indicating the size of the vertices. Default is 5.
#' @param vertex_label Logical indicating whether to display vertex labels. Default is FALSE.
#' @param edge_width Numeric value indicating the width of the edges. Default is 1.
#' @param main A character string for the main title of the plot. Default is NULL.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A plot is produced. No return value.
#' @export
#'
#' @examples
#' # Generate a graph with community structure
#' G <- generate_community_graph(n_nodes = 50, n_communities = 3, intra_prob = 0.8, inter_prob = 0.2)
#'
#' # Plot the graph with communities highlighted
#' plot_graph_with_communities(G)
plot_graph_with_communities <- function(G, communities = NULL, layout = igraph::layout_with_fr,
                                        vertex_size = 5, vertex_label = FALSE,
                                        edge_width = 1, main = NULL, ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package is required for this function. Please install it using install.packages('igraph').")
  }

  if (!is.matrix(G) || nrow(G) != ncol(G)) {
    stop("G must be a square adjacency matrix.")
  }

  # Create igraph object
  graph <- igraph::graph_from_adjacency_matrix(G, mode = "undirected", diag = FALSE)

  # Detect communities if not provided
  if (is.null(communities)) {
    community_detection <- igraph::cluster_louvain(graph)
    communities <- igraph::membership(community_detection)
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

  # Plot the graph
  plot(graph,
       layout = layout,
       vertex.color = igraph::V(graph)$color,
       vertex.size = igraph::V(graph)$size,
       vertex.label = igraph::V(graph)$label,
       edge.width = igraph::E(graph)$width,
       main = main,
       ...)

  # Add legend
  legend("topright",
         legend = paste("Community", unique_communities),
         col = color_palette,
         pch = 19,
         pt.cex = vertex_size / 5,
         bty = "n")
}
