# R/rank_edges.R

#' Rank Edges Based on P-Values
#'
#' This function ranks edges based on their p-values obtained from statistical tests,
#' ordering them from lowest to highest p-value (most to least significant).
#' Values differing by at most 64 machine epsilons relative to the smallest
#' value in a tie group are treated as numerical ties. Tied edges are ordered
#' by their second endpoint, then their first endpoint (upper-triangle order).
#' This removes floating-point noise from the ordering without rounding small
#' p-values to zero. The returned p-values are unchanged; this tolerance is
#' used only for ranking, not for significance thresholds.
#'
#' @param edge_pvalues A square matrix of p-values for each edge, typically obtained from \code{compute_edge_pvalues}.
#' @param tie_values Optional square numeric matrix with the same dimensions as
#'   \code{edge_pvalues} (such as unadjusted marginal p-values when
#'   \code{edge_pvalues} has been multiplicity-adjusted with
#'   \code{\link[stats]{p.adjust}}), used to break ties within each primary
#'   numerical tie group before falling back to endpoint order. Default
#'   \code{NULL}.
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

  # Anchor every tie group at its smallest p-value, so a chain of close
  # neighbours cannot join values whose total separation exceeds tolerance.
  # Optional secondary values (e.g. raw p-values under p.adjust() plateaus)
  # refine ties before falling back to upper-triangle endpoint order.
  groups <- .pvalue_tie_groups(p_values)
  secondary <- if (is.null(tie_values)) rep(0L, length(p_values)) else
    .pvalue_tie_groups(tie_values[edge_indices])
  edge_df <- edge_df[order(groups, secondary, edge_df$node2, edge_df$node1), ]
  rownames(edge_df) <- NULL

  return(edge_df)
}

#' @noRd
#' @keywords internal
.pvalue_tie_groups <- function(p_values) {
  ord <- order(p_values)
  groups <- integer(length(p_values))
  group <- 0L
  anchor <- NA_real_
  for (i in ord) {
    value <- p_values[i]
    if (group == 0L || !is.finite(value) || !is.finite(anchor) ||
        value - anchor > .PVALUE_REL_TOL * abs(anchor)) {
      group <- group + 1L
      anchor <- value
    }
    groups[i] <- group
  }
  groups
}
