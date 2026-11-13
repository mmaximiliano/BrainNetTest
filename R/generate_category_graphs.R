# R/generate_category_graphs.R

#' Generate a Set of Graphs with Similar Community Structures for a Category
#'
#' This function generates a list of adjacency matrices representing brain networks
#' belonging to the same category (e.g., Control group). Every graph has the same
#' community structure, and each edge is an independent Bernoulli variable with the
#' within- or between-community probability. By default all graphs share the same
#' probabilities; optionally, the probabilities of each graph can be shifted at random
#' (see \code{intra_prob_variation} and \code{inter_prob_variation}).
#'
#' @param n_graphs An integer specifying the number of graphs to generate. Default is 10.
#' @param n_nodes An integer specifying the total number of nodes (brain regions). Default is 100.
#' @param n_communities An integer specifying the number of communities. Default is 4.
#' @param community_sizes An integer vector specifying the sizes of each community.
#'   If NULL, communities are of equal size. Default is NULL.
#' @param base_intra_prob A numeric value between 0 and 1, or a numeric vector
#'   of length \code{n_communities}, specifying the base probability of an edge
#'   existing between nodes within the same community. If a vector, each element
#'   sets the base probability for the corresponding community. Default is 0.8.
#' @param base_inter_prob A numeric value between 0 and 1 specifying the base probability
#'   of an edge existing between nodes from different communities. Default is 0.2.
#' @param intra_prob_variation A numeric value \eqn{v \ge 0}. If positive, the
#'   within-community probability of each community is shifted, in every graph, by an
#'   amount drawn uniformly from \eqn{[-v, v]}, independently for each community and
#'   each graph (and truncated to \eqn{[0, 1]}). Default is 0: no variation, all
#'   graphs share the same probabilities.
#' @param inter_prob_variation A numeric value \eqn{v \ge 0}. If positive, the
#'   between-community probability of every graph is shifted by an amount drawn
#'   uniformly from \eqn{[-v, v]}, independently for each graph. Default is 0.
#' @param seed An optional integer for setting the random seed to ensure reproducibility.
#'   Default is NULL.
#'
#' @return A list of symmetric binary adjacency matrices with no self-loops, representing
#'   brain networks with similar community structures.
#' @export
#'
#' @importFrom stats runif
#' @examples
#' # Generate a set of 5 graphs for the Control category
#' control_graphs <- generate_category_graphs(n_graphs = 5, n_nodes = 100, n_communities = 4,
#'                                            base_intra_prob = 0.8, base_inter_prob = 0.2)
generate_category_graphs <- function(n_graphs = 10, n_nodes = 100, n_communities = 4, community_sizes = NULL,
                                     base_intra_prob = 0.8, base_inter_prob = 0.2,
                                     intra_prob_variation = 0, inter_prob_variation = 0, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n_graphs      <- .check_count(n_graphs, "n_graphs")
  n_nodes       <- .check_count(n_nodes, "n_nodes")
  n_communities <- .check_count(n_communities, "n_communities")
  if (n_communities > n_nodes) {
    stop("n_communities (", n_communities, ") must not exceed n_nodes (",
         n_nodes, ").")
  }
  
  # Validate probabilities and variations
  if (!is.numeric(base_intra_prob) || any(base_intra_prob < 0) || any(base_intra_prob > 1)) {
    stop("base_intra_prob must be between 0 and 1.")
  }
  if (!is.numeric(base_inter_prob) || base_inter_prob < 0 || base_inter_prob > 1) {
    stop("base_inter_prob must be between 0 and 1.")
  }
  if (!is.numeric(intra_prob_variation) || intra_prob_variation < 0 || intra_prob_variation > 1) {
    stop("intra_prob_variation must be between 0 and 1.")
  }
  if (!is.numeric(inter_prob_variation) || inter_prob_variation < 0 || inter_prob_variation > 1) {
    stop("inter_prob_variation must be between 0 and 1.")
  }
  
  # Generate a base community assignment
  if (is.null(community_sizes)) {
    # Equal community sizes
    base_size <- n_nodes %/% n_communities
    remainder <- n_nodes %% n_communities
    community_sizes <- rep(base_size, n_communities)
    if (remainder > 0) {
      community_sizes[seq_len(remainder)] <-
        community_sizes[seq_len(remainder)] + 1
    }
  } else {
    community_sizes <- .check_counts(community_sizes, "community_sizes")
    if (length(community_sizes) != n_communities) {
      stop("Length of community_sizes must equal n_communities.")
    }
    if (sum(community_sizes) != n_nodes) {
      stop("Sum of community_sizes must equal n_nodes.")
    }
  }
  
  # List to store generated graphs
  graph_list <- vector("list", n_graphs)
  
  for (g in seq_len(n_graphs)) {
    # Optional per-graph variation of the probabilities (none by default)
    intra_prob <- base_intra_prob + runif(length(base_intra_prob),
                                          -intra_prob_variation, intra_prob_variation)
    inter_prob <- base_inter_prob + runif(1, -inter_prob_variation, inter_prob_variation)
    
    # Ensure probabilities remain within [0, 1]
    intra_prob <- pmin(pmax(intra_prob, 0), 1)
    inter_prob <- min(max(inter_prob, 0), 1)
    
    # Generate the graph using the base community assignments
    G <- generate_community_graph(n_nodes = n_nodes, n_communities = n_communities,
                                  community_sizes = community_sizes,
                                  intra_prob = intra_prob, inter_prob = inter_prob)
    
    graph_list[[g]] <- G
  }
  
  return(graph_list)
}
