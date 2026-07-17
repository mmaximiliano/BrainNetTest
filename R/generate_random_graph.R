#' Generate a Random Symmetric Adjacency Matrix
#'
#' This function generates a random symmetric adjacency matrix representing
#' a brain network. The adjacency matrix is binary, with edges present based
#' on a specified probability.
#'
#' @param n_nodes An integer specifying the number of nodes (brain regions).
#' @param edge_prob A numeric value between 0 and 1 specifying the probability
#'   of an edge existing between any two nodes.
#'
#' @return A symmetric binary adjacency matrix with no self-loops.
#' @export
#'
#' @importFrom stats rbinom
#' @examples
#' graph <- generate_random_graph(n_nodes = 10, edge_prob = 0.1)
generate_random_graph <- function(n_nodes, edge_prob = 0.1) {
  n_nodes <- .assert_whole_number(n_nodes, "n_nodes", minimum = 2L)
  .assert_scalar_number(edge_prob, "edge_prob", lower = 0, upper = 1)

  graph <- matrix(0L, nrow = n_nodes, ncol = n_nodes)
  upper_indices <- which(upper.tri(graph), arr.ind = TRUE)
  graph[upper_indices] <- stats::rbinom(
    nrow(upper_indices),
    1L,
    edge_prob
  )
  graph <- graph + t(graph)
  diag(graph) <- 0L
  graph
}
