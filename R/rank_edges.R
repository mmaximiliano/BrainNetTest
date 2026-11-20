# R/rank_edges.R

#' Rank Edges Based on P-Values
#'
#' This function ranks edges based on their p-values obtained from statistical tests,
#' ordering them from lowest to highest p-value (most to least significant).
#' The ranking compares the stored p-values directly. For equal values it
#' uses the optional secondary values, then upper-triangle order (second
#' endpoint first, then first endpoint).
#'
#' @param edge_pvalues A square matrix of p-values for each edge, typically obtained from \code{compute_edge_pvalues}.
#' @param tie_values Optional square numeric matrix with the same dimensions as
#'   \code{edge_pvalues} (such as unadjusted marginal p-values when
#'   \code{edge_pvalues} has been multiplicity-adjusted with
#'   \code{\link[stats]{p.adjust}}), used to order edges with exactly equal
#'   primary p-values. Default \code{NULL}.
#' @return A data frame with edges and their corresponding p-values, ordered from most significant.
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
#' 
#' # Rank edges based on p-values
#' edge_df <- rank_edges(edge_pvalues)
#' # View the top ranked edges
#' head(edge_df)
rank_edges <- function(edge_pvalues, tie_values = NULL) {
  if (!is.matrix(edge_pvalues) || !is.numeric(edge_pvalues) ||
      nrow(edge_pvalues) != ncol(edge_pvalues)) {
    stop("`edge_pvalues` must be a square numeric matrix.", call. = FALSE)
  }
  if (!is.null(tie_values) &&
      (!is.matrix(tie_values) || !is.numeric(tie_values) ||
       !identical(dim(tie_values), dim(edge_pvalues)))) {
    stop("`tie_values` must be NULL or a numeric matrix with the same dimensions as `edge_pvalues`.",
         call. = FALSE)
  }

  # Convert upper triangle of p-value matrix to a vector
  edge_indices <- which(upper.tri(edge_pvalues), arr.ind = TRUE)
  p_values <- edge_pvalues[edge_indices]
  
  # Create a data frame
  edge_df <- data.frame(
    node1 = edge_indices[, 1],
    node2 = edge_indices[, 2],
    p_value = p_values
  )

  # The upper triangle is already in endpoint order. Stable sorting preserves
  # that order when the primary and any supplied secondary values are equal.
  ord <- if (is.null(tie_values)) order(p_values) else
    order(p_values, tie_values[edge_indices])
  edge_df <- edge_df[ord, ]
  rownames(edge_df) <- NULL
  edge_df
}
