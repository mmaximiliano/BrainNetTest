# R/generate_category_graphs.R

#' Generate a Set of Graphs with Similar Community Structures for a Category
#'
#' This function generates a list of adjacency matrices representing brain networks
#' belonging to the same category (e.g., Control group). The generated graphs have
#' similar community structures but vary slightly in their intra-community and
#' inter-community connection probabilities to reflect natural variability.
#'
#' @param n_graphs An integer specifying the number of graphs to generate. Default is 10.
#' @param n_nodes An integer specifying the total number of nodes (brain regions). Default is 100.
#' @param n_communities An integer specifying the number of communities. Default is 4.
#' @param community_sizes An integer vector specifying the sizes of each community.
#'   If NULL, communities are of equal size. Default is NULL.
#' @param base_intra_prob A numeric value between 0 and 1 specifying the base probability
#'   of an edge existing between nodes within the same community. Default is 0.8.
#' @param base_inter_prob A numeric value between 0 and 1 specifying the base probability
#'   of an edge existing between nodes from different communities. Default is 0.2.
#' @param intra_prob_variation A numeric value specifying the maximum variation to apply
#'   to the intra-community probability for each graph. Default is 0.05.
#' @param inter_prob_variation A numeric value specifying the maximum variation to apply
#'   to the inter-community probability for each graph. Default is 0.05.
#' @param seed An optional integer for setting the random seed to ensure reproducibility.
#'   Default is NULL.
#'
#' @return A list of symmetric binary adjacency matrices with no self-loops, representing
#'   brain networks with similar community structures.
#' @export
#'
#' @examples
#' # Generate a set of 5 graphs for the Control category
#' control_graphs <- generate_category_graphs(n_graphs = 5, n_nodes = 100, n_communities = 4,
#'                                            base_intra_prob = 0.8, base_inter_prob = 0.2)
generate_category_graphs <- function(n_graphs = 10, n_nodes = 100, n_communities = 4, community_sizes = NULL,
                                     base_intra_prob = 0.8, base_inter_prob = 0.2,
                                     intra_prob_variation = 0.05, inter_prob_variation = 0.05, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (!is.numeric(n_graphs) || length(n_graphs) != 1 || n_graphs <= 0) {
    stop("n_graphs must be a positive integer.")
  }
  n_graphs <- as.integer(n_graphs)
  
  # Validate probabilities and variations
  if (!is.numeric(base_intra_prob) || base_intra_prob < 0 || base_intra_prob > 1) {
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
  
  # List to store generated graphs
  graph_list <- vector("list", n_graphs)
  
  for (g in 1:n_graphs) {
    # Apply slight variations to probabilities
    intra_prob <- base_intra_prob + runif(1, -intra_prob_variation, intra_prob_variation)
    inter_prob <- base_inter_prob + runif(1, -inter_prob_variation, inter_prob_variation)
    
    # Ensure probabilities remain within [0, 1]
    intra_prob <- min(max(intra_prob, 0), 1)
    inter_prob <- min(max(inter_prob, 0), 1)
    
    # Generate the graph using the base community assignments
    G <- generate_community_graph(n_nodes = n_nodes, n_communities = n_communities,
                                  community_sizes = community_sizes,
                                  intra_prob = intra_prob, inter_prob = inter_prob)
    
    graph_list[[g]] <- G
  }
  
  return(graph_list)
}
