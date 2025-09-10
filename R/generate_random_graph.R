# R/generate_random_graph.R

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
#' G <- generate_random_graph(n_nodes = 10, edge_prob = 0.1)
generate_random_graph <- function(n_nodes, edge_prob = 0.1) {
  if (!is.numeric(n_nodes) || length(n_nodes) != 1 || n_nodes <= 0) {
    stop("n_nodes must be a positive integer.")
  }
  
  if (!is.numeric(edge_prob) || edge_prob < 0 || edge_prob > 1) {
    stop("edge_prob must be a numeric value between 0 and 1.")
  }
  
  # Initialize a matrix of zeros
  G <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  
  # Get indices of upper triangle excluding the diagonal
  upper_indices <- which(upper.tri(G), arr.ind = TRUE)
  
  # Assign random edges to upper triangle
  G[upper_indices] <- rbinom(nrow(upper_indices), 1, edge_prob)
  
  # Mirror the upper triangle to the lower triangle to make it symmetric
  G <- G + t(G)
  
  # Ensure diagonal is zero (no self-loops)
  diag(G) <- 0
  
  return(G)
}
