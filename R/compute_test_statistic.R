# R/compute_test_statistic.R

#' Compute the Test Statistic T for Brain Network Populations
#'
#' This function computes the test statistic T to assess whether different
#' populations of brain networks originate from the same distribution.
#'
#' @param populations A named list where each element is a list of adjacency
#'   matrices for a population. Example:
#'   \code{list(Control = list(G1, G2, ...), Schizophrenia = list(G1, G2, ...), Alzheimer = list(G1, G2, ...))}
#' @param a A normalization constant. Default is 1.
#'
#' @return A numeric value representing the test statistic T.
#' @export
#'
#' @importFrom stats rbinom
#' @examples
#' # Generate synthetic populations data
#' Control <- list(
#'   generate_random_graph(n_nodes = 5, edge_prob = 0.1),
#'   generate_random_graph(n_nodes = 5, edge_prob = 0.1)
#' )
#' Schizophrenia <- list(
#'   generate_random_graph(n_nodes = 5, edge_prob = 0.15),
#'   generate_random_graph(n_nodes = 5, edge_prob = 0.15)
#' )
#' Alzheimer <- list(
#'   generate_random_graph(n_nodes = 5, edge_prob = 0.2),
#'   generate_random_graph(n_nodes = 5, edge_prob = 0.2)
#' )
#' 
#' populations <- list(Control = Control, Schizophrenia = Schizophrenia, Alzheimer = Alzheimer)
#' 
#' # Compute the test statistic T
#' T_value <- compute_test_statistic(populations, a = 1)
#' print(T_value)
compute_test_statistic <- function(populations, a = 1) {
  if (!is.list(populations) || length(populations) == 0) {
    stop("populations must be a non-empty list.")
  }
  
  m <- length(populations)  # Number of populations
  n_i <- sapply(populations, length)  # Number of graphs in each population
  n <- sum(n_i)  # Total number of graphs
  
  if (any(n_i < 2)) {
    stop("Each population must have at least two graphs.")
  }
  
  # Compute central graphs for each population
  central_graphs <- lapply(populations, compute_central_graph)
  
  # Compute average distance of each population to its central graph
  avg_d_Gi_Mi <- numeric(m)
  names(avg_d_Gi_Mi) <- names(populations)
  
  for (i in seq_along(populations)) {
    distances <- sapply(populations[[i]], compute_distance, M = central_graphs[[i]])
    avg_d_Gi_Mi[i] <- mean(distances)
  }
  
  # Compute average distance of all graphs to each central graph
  avg_d_G_Mi <- numeric(m)
  names(avg_d_G_Mi) <- names(populations)
  
  for (i in seq_along(central_graphs)) {
    all_distances <- numeric(n)
    idx <- 1
    for (j in seq_along(populations)) {
      distances <- sapply(populations[[j]], compute_distance, M = central_graphs[[i]])
      all_distances[idx:(idx + length(distances) - 1)] <- distances
      idx <- idx + length(distances)
    }
    avg_d_G_Mi[i] <- mean(all_distances)
  }
  
  # Compute the sum part of the test statistic
  sum_term <- 0
  for (i in seq_along(populations)) {
    term1 <- (n_i[i] / (n_i[i] - 1)) * avg_d_Gi_Mi[i]
    term2 <- (n / (n - 1)) * avg_d_G_Mi[i]
    sum_term <- sum_term + sqrt(n_i[i]) * (term1 - term2)
  }
  
  # Compute the test statistic T
  T <- (sqrt(m) / a) * sum_term
  
  return(T)
}
