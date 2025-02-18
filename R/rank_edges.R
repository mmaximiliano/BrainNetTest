# R/rank_edges.R

#' Rank Edges Based on P-Values
#'
#' This function ranks edges based on their p-values obtained from statistical tests,
#' ordering them from lowest to highest p-value (most to least significant).
#'
#' @param edge_pvalues A square matrix of p-values for each edge, typically obtained from \code{compute_edge_pvalues}.
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
#' disease_graphs <- generate_category_graphs(
#'   n_graphs = 5,
#'   n_nodes = 10,
#'   n_communities = 2,
#'   base_intra_prob = 0.6,
#'   base_inter_prob = 0.4,
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
#' 
#' # Rank edges based on p-values
#' edge_df <- rank_edges(edge_pvalues)
#' # View the top ranked edges
#' head(edge_df)
rank_edges <- function(edge_pvalues) {
  # Convert upper triangle of p-value matrix to a vector
  edge_indices <- which(upper.tri(edge_pvalues), arr.ind = TRUE)
  p_values <- edge_pvalues[edge_indices]
  
  # Create a data frame
  edge_df <- data.frame(
    node1 = edge_indices[, 1],
    node2 = edge_indices[, 2],
    p_value = p_values
  )
  
  # Order by p-value ascending (most significant first)
  edge_df <- edge_df[order(edge_df$p_value, decreasing = FALSE), ]
  
  return(edge_df)
}
