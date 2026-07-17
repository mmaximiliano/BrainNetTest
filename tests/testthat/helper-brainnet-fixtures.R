graph_from_edges <- function(n_nodes, edges = matrix(integer(), ncol = 2L)) {
  graph <- matrix(0L, nrow = n_nodes, ncol = n_nodes)
  if (length(edges) > 0L) {
    edges <- matrix(edges, ncol = 2L)
    graph[edges] <- 1L
    graph[cbind(edges[, 2L], edges[, 1L])] <- 1L
  }
  graph
}

constant_group <- function(graph, n = 3L) {
  replicate(n, graph, simplify = FALSE)
}

tiny_brainnet_data <- function() {
  empty <- graph_from_edges(3L)
  edge_12 <- graph_from_edges(3L, c(1L, 2L))
  edge_13 <- graph_from_edges(3L, c(1L, 3L))

  brainnet_data(list(
    GroupA = list(empty, edge_12, edge_13),
    GroupB = list(edge_12, edge_12, edge_13)
  ))
}
