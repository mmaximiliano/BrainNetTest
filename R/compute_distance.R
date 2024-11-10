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
#' @examples
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
