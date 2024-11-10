# R/compute_central_graph.R

#' Compute the Central (Representative) Graph for a Population
#'
#' This function computes the central graph for a given population by averaging
#' the adjacency matrices of all graphs within that population.
#'
#' @param graph_list A list of adjacency matrices representing brain networks.
#'
#' @return A single adjacency matrix representing the central graph of the population.
#' @export
#'
#' @examples
#' central_graph <- compute_central_graph(Control)
compute_central_graph <- function(graph_list) {
  if (length(graph_list) == 0) {
    stop("graph_list is empty.")
  }
  
  # Check if all matrices are square and have the same dimensions
  n_nodes <- nrow(graph_list[[1]])
  for (G in graph_list) {
    if (!is.matrix(G)) {
      stop("All elements of graph_list must be matrices.")
    }
    if (nrow(G) != n_nodes || ncol(G) != n_nodes) {
      stop("All adjacency matrices must be square and of the same dimensions.")
    }
  }
  
  # Initialize a matrix of zeros
  central_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  
  # Sum all adjacency matrices
  for (G in graph_list) {
    central_matrix <- central_matrix + G
  }
  
  # Average the summed matrix
  central_matrix <- central_matrix / length(graph_list)
  
  return(central_matrix)
}
