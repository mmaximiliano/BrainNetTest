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
#'   If `NULL`, nodes are divided as evenly as possible. Default is `NULL`.
#' @param base_intra_prob A numeric value between 0 and 1, or a numeric vector
#'   of length \code{n_communities}, specifying the base probability of an edge
#'   existing between nodes within the same community. If a vector, each element
#'   sets the base probability for the corresponding community. Default is 0.8.
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
#' @importFrom stats runif
#' @examples
#' # Generate a set of 5 graphs for the Control category
#' control_graphs <- generate_category_graphs(
#'   n_graphs = 5, n_nodes = 100, n_communities = 4,
#'   base_intra_prob = 0.8, base_inter_prob = 0.2
#' )
generate_category_graphs <- function(n_graphs = 10,
                                     n_nodes = 100,
                                     n_communities = 4,
                                     community_sizes = NULL,
                                     base_intra_prob = 0.8,
                                     base_inter_prob = 0.2,
                                     intra_prob_variation = 0.05,
                                     inter_prob_variation = 0.05,
                                     seed = NULL) {
  n_graphs <- .assert_whole_number(n_graphs, "n_graphs")
  n_nodes <- .assert_whole_number(n_nodes, "n_nodes", minimum = 2L)
  n_communities <- .assert_whole_number(
    n_communities,
    "n_communities"
  )
  if (n_communities > n_nodes) {
    stop("`n_communities` cannot exceed `n_nodes`.", call. = FALSE)
  }
  community_sizes <- .normalize_community_sizes(
    n_nodes,
    n_communities,
    community_sizes
  )
  base_intra_prob <- .validate_probability_vector(
    base_intra_prob,
    "base_intra_prob",
    n_communities
  )
  .assert_scalar_number(
    base_inter_prob,
    "base_inter_prob",
    lower = 0,
    upper = 1
  )
  .assert_scalar_number(
    intra_prob_variation,
    "intra_prob_variation",
    lower = 0,
    upper = 1
  )
  .assert_scalar_number(
    inter_prob_variation,
    "inter_prob_variation",
    lower = 0,
    upper = 1
  )
  if (!is.null(seed)) {
    .assert_whole_number(seed, "seed", minimum = 0L)
  }

  .with_preserved_seed(seed, function() {
    graphs <- vector("list", n_graphs)
    for (graph_index in seq_len(n_graphs)) {
      intra_prob <- base_intra_prob + stats::runif(
        n_communities,
        -intra_prob_variation,
        intra_prob_variation
      )
      inter_prob <- base_inter_prob + stats::runif(
        1L,
        -inter_prob_variation,
        inter_prob_variation
      )

      graphs[[graph_index]] <- generate_community_graph(
        n_nodes = n_nodes,
        n_communities = n_communities,
        community_sizes = community_sizes,
        intra_prob = pmin(pmax(intra_prob, 0), 1),
        inter_prob = min(max(inter_prob, 0), 1)
      )
    }
    graphs
  })
}
