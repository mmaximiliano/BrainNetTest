# R/identify_critical_links.R

#' Identify Critical Links in Graphs Between Populations
#'
#' This function identifies the critical links (edges) in graphs that contribute
#' to significant differences observed between populations in a hypothesis test.
#' It iteratively removes edges based on their statistical significance and
#' re-tests until the difference is no longer significant.
#'
#' @param populations A list of populations, each containing a list of graphs (adjacency matrices).
#' @param alpha The significance level for the hypothesis test. Default is 0.05.
#' @param method The statistical test method to use for edge comparison. Options are "fisher", "chi.squared", "prop".
#'   Default is "fisher".
#' @param adjust_method The method for p-value adjustment for multiple testing. See ?p.adjust. Default is "none".
#' @param batch_size Number of edges to remove in each iteration. Default is 1.
#' @param n_bootstrap Number of bootstrap resamples for computing p-value. Default is 1000.
#' @importFrom stats pnorm sample
#' @return A list containing critical edges, edges removed, and modified populations.
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
#' # Identify critical links
#' result <- identify_critical_links(populations, alpha = 0.05)
#' # View critical edges
#' critical_edges <- result$critical_edges
#' print(critical_edges)
identify_critical_links <- function(populations, alpha = 0.05, method = "fisher",
                                    adjust_method = "none", batch_size = 1,
                                    n_bootstrap = 1000) {
  N <- sapply(populations, length)
  num_populations <- length(populations)

  # Step 1 & 2: Compute edge frequencies and p-values
  frequencies <- compute_edge_frequencies(populations)
  edge_pvalues <- compute_edge_pvalues(frequencies$edge_counts, N, method, adjust_method)

  # Step 3: Rank edges (ascending order of p-values)
  edge_df <- rank_edges(edge_pvalues)

  # Initialize variables
  edges_removed <- list()
  modified_populations <- populations

  # Compute initial test statistic
  initial_T <- compute_test_statistic(modified_populations)

  # Define a function to compute p-value from T using bootstrap
  compute_p_value_from_T <- function(T_observed, populations, n_bootstrap = 1000) {
    # Flatten all graphs from all populations into a single list
    all_graphs <- unlist(populations, recursive = FALSE)
    n_total <- length(all_graphs)

    # Store bootstrap T values
    t_bootstrap <- numeric(n_bootstrap)

    # Original population sizes
    pop_sizes <- sapply(populations, length)

    # Perform bootstrap resampling
    for(b in 1:n_bootstrap) {
      # For each bootstrap iteration, create random populations by resampling
      # from the pooled set of graphs (under the null hypothesis that all graphs
      # come from the same distribution)
      bootstrap_pops <- list()
      start_idx <- 1
      
      for(p in 1:length(populations)) {
        # Sample with replacement from all graphs
        sampled_indices <- sample(1:n_total, pop_sizes[p], replace = TRUE)
        bootstrap_pops[[p]] <- all_graphs[sampled_indices]
      }

      names(bootstrap_pops) <- names(populations)
      
      # Calculate test statistic for this bootstrap sample
      t_bootstrap[b] <- compute_test_statistic(bootstrap_pops)
    }
    
    # Calculate p-value as the proportion of bootstrap samples with T statistics
    # more extreme than or equal to the observed value
    p_value <- mean(t_bootstrap <= T_observed)
    
    # For numerical stability, ensure p_value is never exactly 0
    p_value <- max(p_value, 1/n_bootstrap)
    
    return(p_value)
  }

  # Compute initial p-value using bootstrap
  initial_p_value <- compute_p_value_from_T(initial_T, modified_populations, n_bootstrap)

  # Check if initial test is significant
  if (initial_p_value > alpha) {
    warning("Initial test is not significant. No critical links to identify.")
    return(list(
      critical_edges = NULL,
      edges_removed = list(),
      modified_populations = modified_populations
    ))
  }

  significant <- TRUE
  idx <- 1
  while (significant && idx <= nrow(edge_df)) {
    # Remove edges in batches
    batch_indices <- idx:min(idx + batch_size - 1, nrow(edge_df))
    batch_edges <- edge_df[batch_indices, ]

    # Remove edges from all graphs
    for (edge_row in seq_len(nrow(batch_edges))) {
      i <- batch_edges$node1[edge_row]
      j <- batch_edges$node2[edge_row]
      for (k in 1:num_populations) {
        for (g in 1:N[k]) {
          modified_populations[[k]][[g]][i, j] <- 0
          modified_populations[[k]][[g]][j, i] <- 0  # Symmetric
        }
      }
      edges_removed[[length(edges_removed) + 1]] <- c(i, j)
    }
    # Recompute test statistic
    T_value <- compute_test_statistic(modified_populations)
    
    # Compute p-value from T_value using bootstrap
    p_value <- compute_p_value_from_T(T_value, modified_populations, n_bootstrap)
    
    # Check if test is no longer significant
    if (p_value > alpha) {
      significant <- FALSE
      break  # Stop removing edges
    }
    
    idx <- idx + batch_size
  }
  
  # Critical edges are the edges removed
  critical_edges <- edge_df[1:(idx - 1), c("node1", "node2", "p_value")]
  
  return(list(
    critical_edges = critical_edges,
    edges_removed = edges_removed,
    modified_populations = modified_populations
  ))
}
