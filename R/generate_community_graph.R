# R/generate_community_graph.R

#' Generate a Random Symmetric Adjacency Matrix with Community Structure
#'
#' This function generates a random symmetric adjacency matrix representing
#' a brain network with community structure. Nodes within the same community
#' have a higher probability of being connected compared to nodes from different
#' communities.
#'
#' @param n_nodes An integer specifying the total number of nodes (brain regions).
#'   Default is 100.
#' @param n_communities An integer specifying the number of communities. Default is 4.
#' @param community_sizes An integer vector specifying the sizes of each community.
#'   If NULL, communities are of equal size. Default is NULL.
#' @param intra_prob A numeric value between 0 and 1 specifying the probability
#'   of an edge existing between nodes within the same community. Default is 0.8.
#' @param inter_prob A numeric value between 0 and 1 specifying the probability
#'   of an edge existing between nodes from different communities. Default is 0.2.
#' @param seed An optional integer for setting the random seed to ensure reproducibility.
#'   Default is NULL.
#'
#' @return A symmetric binary adjacency matrix with no self-loops, representing
#'   a brain network with community structure.
#' @export
#'
#' @importFrom stats rbinom
#' @examples
#' # Generate a brain network with community structure
#' G <- generate_community_graph(n_nodes = 100, n_communities = 4, intra_prob = 0.8, inter_prob = 0.2)
generate_community_graph <- function(n_nodes = 100, n_communities = 4, community_sizes = NULL,
                                     intra_prob = 0.8, inter_prob = 0.2, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.numeric(n_nodes) || length(n_nodes) != 1 || n_nodes <= 0) {
    stop("n_nodes must be a positive integer.")
  }
  n_nodes <- as.integer(n_nodes)
  
  if (!is.numeric(n_communities) || length(n_communities) != 1 || n_communities <= 0) {
    stop("n_communities must be a positive integer.")
  }
  n_communities <- as.integer(n_communities)
  
  if (is.null(community_sizes)) {
    # If community sizes are not provided, divide nodes equally
    base_size <- n_nodes %/% n_communities
    remainder <- n_nodes %% n_communities
    community_sizes <- rep(base_size, n_communities)
    if (remainder > 0) {
      community_sizes[1:remainder] <- community_sizes[1:remainder] + 1
    }
  } else {
    if (length(community_sizes) != n_communities) {
      stop("Length of community_sizes must equal n_communities.")
    }
    if (sum(community_sizes) != n_nodes) {
      stop("Sum of community_sizes must equal n_nodes.")
    }
  }

  if (!is.numeric(intra_prob) || intra_prob < 0 || intra_prob > 1) {
    stop("intra_prob must be a numeric value between 0 and 1.")
  }
  
  if (!is.numeric(inter_prob) || inter_prob < 0 || inter_prob > 1) {
    stop("inter_prob must be a numeric value between 0 and 1.")
  }

  # Initialize adjacency matrix
  G <- matrix(0, nrow = n_nodes, ncol = n_nodes)

  # Assign nodes to communities
  node_indices <- 1:n_nodes
  community_assignments <- rep(1:n_communities, times = community_sizes)

  # Create community blocks
  for (i in 1:n_communities) {
    nodes_in_i <- node_indices[community_assignments == i]
    # Edges within community i
    if (length(nodes_in_i) > 1) {
      G_intra <- rbinom(length(nodes_in_i) * (length(nodes_in_i) - 1) / 2, 1, intra_prob)
      G_intra_matrix <- matrix(0, nrow = length(nodes_in_i), ncol = length(nodes_in_i))
      G_intra_matrix[upper.tri(G_intra_matrix)] <- G_intra
      G_intra_matrix <- G_intra_matrix + t(G_intra_matrix)
      G[nodes_in_i, nodes_in_i] <- G_intra_matrix
    }
  }

  # Create inter-community edges
  for (i in 1:(n_communities - 1)) {
    for (j in (i + 1):n_communities) {
      nodes_in_i <- node_indices[community_assignments == i]
      nodes_in_j <- node_indices[community_assignments == j]
      # Edges between community i and community j
      G_inter <- matrix(rbinom(length(nodes_in_i) * length(nodes_in_j), 1, inter_prob),
                        nrow = length(nodes_in_i), ncol = length(nodes_in_j))
      G[nodes_in_i, nodes_in_j] <- G_inter
      G[nodes_in_j, nodes_in_i] <- t(G_inter)
    }
  }

  # Remove self-loops
  diag(G) <- 0

  return(G)
}
