# R/compute_edge_pvalues.R

#' Compute P-Values for Edge Differences Between Populations
#'
#' This function computes p-values for the differences in edge presence between populations
#' using statistical tests (e.g., Fisher's exact test, chi-squared test).
#'
#' @param edge_counts An array of edge counts from \code{compute_edge_frequencies}.
#' @param N A vector of sample sizes for each population.
#' @param method The statistical test method to use: "fisher", "chi.squared", or "prop". Default is "fisher".
#' @param adjust_method The method for p-value adjustment for multiple testing. Default is "none".
#' @importFrom stats fisher.test chisq.test prop.test p.adjust
#' @return A matrix of p-values for each edge.
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
#' edge_counts <- frequencies$edge_counts
#' N <- sapply(populations, length)
#'
#' # Compute p-values for edge differences
#' edge_pvalues <- compute_edge_pvalues(edge_counts, N)
#' # View p-values for the first few edges
#' print(edge_pvalues[1:5, 1:5])
compute_edge_pvalues <- function(edge_counts, N, method = "fisher", adjust_method = "none") {
  n_nodes <- dim(edge_counts)[1]
  num_populations <- dim(edge_counts)[3]
  edge_pvalues <- matrix(1, nrow = n_nodes, ncol = n_nodes)
  
  # For each pair of nodes (edge)
  for (i in 1:(n_nodes - 1)) {
    for (j in (i + 1):n_nodes) {
      # Build contingency table
      counts <- matrix(0, nrow = num_populations, ncol = 2)
      for (k in 1:num_populations) {
        counts[k, 1] <- edge_counts[i, j, k]                # Edge present
        counts[k, 2] <- N[k] - edge_counts[i, j, k]         # Edge absent
      }
      
      # Perform the chosen statistical test
      if (method == "fisher") {
        test_result <- fisher.test(counts)
        p_value <- test_result$p.value
      } else if (method == "chi.squared") {
        test_result <- chisq.test(counts)
        p_value <- test_result$p.value
      } else if (method == "prop") {
        # Use prop.test for two populations
        if (num_populations != 2) stop("prop.test is only implemented for two populations.")
        x <- c(edge_counts[i, j, 1], edge_counts[i, j, 2])
        n <- c(N[1], N[2])
        test_result <- prop.test(x, n)
        p_value <- test_result$p.value
      } else {
        stop("Invalid method specified.")
      }
      
      edge_pvalues[i, j] <- p_value
      edge_pvalues[j, i] <- p_value  # Symmetric matrix
    }
  }
  
  # Adjust p-values for multiple testing if needed
  if (adjust_method != "none") {
    # Convert upper triangle of p-value matrix to vector
    p_values_vec <- edge_pvalues[upper.tri(edge_pvalues)]
    adjusted_p_values <- p.adjust(p_values_vec, method = adjust_method)
    # Map back to matrix
    edge_pvalues[upper.tri(edge_pvalues)] <- adjusted_p_values
    edge_pvalues[lower.tri(edge_pvalues)] <- adjusted_p_values
  }
  
  return(edge_pvalues)
}
