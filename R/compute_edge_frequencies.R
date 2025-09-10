# R/compute_edge_frequencies.R

#' Compute Edge Frequencies in Populations
#'
#' This function computes the frequency of each edge (connection between nodes)
#' across all graphs within each population.
#'
#' @param populations A list where each element is a population containing a list of graphs (adjacency matrices).
#' @return A list containing edge counts and edge proportions for each population.
#' @export
#'
#' @examples
#' # Generate synthetic populations
#' control_graphs <- generate_category_graphs(
#'   n_graphs = 5,
#'   n_nodes = 10,
#'   n_communities = 2,
#'   base_intra_prob = 0.8,
#'   base_inter_prob = 0.2,
#'   intra_prob_variation = 0.05,
#'   inter_prob_variation = 0.05,
#'   seed = 1
#' )
#' disease_graphs <- generate_category_graphs(
#'   n_graphs = 5,
#'   n_nodes = 10,
#'   n_communities = 2,
#'   base_intra_prob = 0.6,
#'   base_inter_prob = 0.4,
#'   intra_prob_variation = 0.05,
#'   inter_prob_variation = 0.05,
#'   seed = 2
#' )
#' populations <- list(Control = control_graphs, Disease = disease_graphs)
#'
#' # Compute edge frequencies
#' frequencies <- compute_edge_frequencies(populations)
#' # View edge counts for the first population
#' edge_counts_control <- frequencies$edge_counts[,,1]
#' print(edge_counts_control)
compute_edge_frequencies <- function(populations) {
  num_populations <- length(populations)
  N <- sapply(populations, length)
  
  # Assuming all graphs have the same dimensions
  n_nodes <- nrow(populations[[1]][[1]])
  
  # Initialize arrays to store edge counts and frequencies
  edge_counts <- array(0, dim = c(n_nodes, n_nodes, num_populations))
  
  # For each population
  for (k in seq_along(populations)) {
    # Sum adjacency matrices across all graphs in the population
    edge_counts[,,k] <- Reduce("+", populations[[k]])
  }
  
  # Calculate edge presence proportions
  edge_proportions <- edge_counts
  for (k in seq_along(populations)) {
    edge_proportions[,,k] <- edge_counts[,,k] / N[k]
  }
  
  return(list(edge_counts = edge_counts, edge_proportions = edge_proportions))
}
