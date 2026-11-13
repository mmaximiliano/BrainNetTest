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
#'   seed = 1
#' )
#' patient_graphs <- generate_category_graphs(
#'   n_graphs = 5,
#'   n_nodes = 10,
#'   n_communities = 2,
#'   base_intra_prob = 0.6,
#'   base_inter_prob = 0.4,
#'   seed = 2
#' )
#' populations <- list(Control = control_graphs, Patient = patient_graphs)
#'
#' # Compute edge frequencies
#' frequencies <- compute_edge_frequencies(populations)
#' edge_counts <- frequencies$edge_counts
#' N <- lengths(populations)
#'
#' # Compute p-values for edge differences
#' edge_pvalues <- compute_edge_pvalues(edge_counts, N)
#' # View p-values for the first few edges
#' print(edge_pvalues[1:5, 1:5])
compute_edge_pvalues <- function(edge_counts, N, method = "fisher", adjust_method = "none") {
  n_nodes <- dim(edge_counts)[1]
  num_populations <- dim(edge_counts)[3]
  edge_pvalues <- matrix(1, nrow = n_nodes, ncol = n_nodes)
  
  # For each pair of nodes (edge).  seq_len() rather than 1:(n_nodes - 1) so
  # that a single-node network yields no edges instead of a reversed sequence.
  for (i in seq_len(n_nodes - 1L)) {
    for (j in (i + 1L):n_nodes) {
      # Build contingency table
      counts <- matrix(0, nrow = num_populations, ncol = 2)
      for (k in seq_len(num_populations)) {
        counts[k, 1] <- edge_counts[i, j, k]                # Edge present
        counts[k, 2] <- N[k] - edge_counts[i, j, k]         # Edge absent
      }
      
      # Perform the chosen statistical test
      if (method == "fisher") {
        test_result <- fisher.test(counts)
        p_value <- test_result$p.value
      } else if (method == "chi.squared") {
        # Check if any row or column sum is zero
        row_sums <- rowSums(counts)
        col_sums <- colSums(counts)
        
        if (any(row_sums == 0) || any(col_sums == 0)) {
          # If zero marginals, fall back to Fisher's exact test
          test_result <- fisher.test(counts)
          p_value <- test_result$p.value
        } else {
          # Check if expected counts might be too small for chi-square approximation
          total <- sum(counts)
          has_small_expected <- FALSE
          
          # Calculate expected counts and check if any are < 5
          for (r in seq_len(nrow(counts))) {
            for (c in seq_len(ncol(counts))) {
              expected <- row_sums[r] * col_sums[c] / total
              if (expected < 5) {
                has_small_expected <- TRUE
                break
              }
            }
            if (has_small_expected) break
          }
          
          # Use standard chi-square test (simulation can fail with certain data patterns)
          test_result <- chisq.test(counts, correct = TRUE)
          p_value <- test_result$p.value
        }
      } else if (method == "prop") {
        # Use prop.test for two populations
        if (num_populations != 2) stop("prop.test is only implemented for two populations.")
        
        # Check for zero values that might cause issues
        x <- c(edge_counts[i, j, 1], edge_counts[i, j, 2])
        n <- c(N[1], N[2])
        
        # Handle edge cases
        if (all(x == 0) || all(x == n)) {
          p_value <- 1.0  # No difference when all are identical
        } else {
          test_result <- prop.test(x, n)
          p_value <- test_result$p.value
        }
      } else {
        stop("Invalid method specified.")
      }
      
      # Exact tests (notably fisher.test) can return p-values fractionally
      # outside [0, 1] because of floating-point rounding in their internal
      # summations. The extra slack is hidden by extended (long double)
      # precision on most platforms but surfaces on builds without it, such as
      # CRAN's noLD check. Clamp to a valid probability so downstream code and
      # tests always see p in [0, 1].
      p_value <- min(max(p_value, 0), 1)
      edge_pvalues[i, j] <- p_value
      edge_pvalues[j, i] <- p_value  # Symmetric matrix
    }
  }
  
  # Adjust p-values for multiple testing if needed
  edge_pvalues <- .adjust_edge_pvalues(edge_pvalues, adjust_method)

  return(edge_pvalues)
}

#' Multiplicity-adjust a symmetric matrix of edge p-values
#'
#' @param edge_pvalues Symmetric matrix of unadjusted edge p-values.
#' @param adjust_method Method passed to `stats::p.adjust()`.
#' @return The matrix with its off-diagonal entries adjusted.
#' @noRd
#' @keywords internal
.adjust_edge_pvalues <- function(edge_pvalues, adjust_method) {
  if (adjust_method != "none") {
    # Adjust using only the upper triangle, then mirror to the lower triangle.
    # Note: assigning the same vector to both upper.tri() and lower.tri() would
    # pair elements incorrectly (R traverses each triangle in column-major
    # order over different index sets), producing an asymmetric matrix. We
    # therefore write the upper triangle first and then transpose.
    ut <- upper.tri(edge_pvalues)
    edge_pvalues[ut] <- p.adjust(edge_pvalues[ut], method = adjust_method)
    edge_pvalues[lower.tri(edge_pvalues)] <- t(edge_pvalues)[lower.tri(edge_pvalues)]
  }
  edge_pvalues
}
