# R/compute_distance.R

#' Compute the Manhattan Norm Distance Between Two Graphs
#'
#' This function calculates the Manhattan norm (L1 norm) distance between two
#' adjacency matrices representing brain networks.
#'
#' @param G A square adjacency matrix representing a brain network.
#' @param M A square adjacency matrix representing the central brain network.
#'
#' @return A numeric value representing the Manhattan distance between G and M.
#' @export
#'
#' @importFrom stats rbinom
#' @examples
#' # Generate synthetic data
#' G1 <- generate_random_graph(n_nodes = 5, edge_prob = 0.1)
#' G2 <- generate_random_graph(n_nodes = 5, edge_prob = 0.1)
#' central_graph <- compute_central_graph(list(G1, G2))
#' distance <- compute_distance(G1, central_graph)
compute_distance <- function(G, M) {
  if (!is.matrix(G) || !is.matrix(M)) {
    stop("Both G and M must be matrices.")
  }
  
  if (!all(dim(G) == dim(M))) {
    stop("G and M must have the same dimensions.")
  }
  
  return(sum(abs(G - M)))
}
