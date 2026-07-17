#' Generate a Random Symmetric Adjacency Matrix with Community Structure
#'
#' This function generates a random symmetric adjacency matrix representing
#' a network with community structure. Users specify within-community and
#' between-community probabilities independently.
#'
#' @param n_nodes An integer specifying the total number of nodes (brain regions).
#'   Default is 100.
#' @param n_communities An integer specifying the number of communities. Default is 4.
#' @param community_sizes An integer vector specifying the sizes of each community.
#'   If `NULL`, nodes are divided as evenly as possible. Default is `NULL`.
#' @param intra_prob A numeric value between 0 and 1, or a numeric vector of
#'   length \code{n_communities}, specifying the probability of an edge existing
#'   between nodes within the same community. If a scalar, the same probability
#'   is used for all communities. Default is 0.8.
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
#' graph <- generate_community_graph(
#'   n_nodes = 20,
#'   n_communities = 2,
#'   intra_prob = 0.8,
#'   inter_prob = 0.2
#' )
generate_community_graph <- function(n_nodes = 100,
                                     n_communities = 4,
                                     community_sizes = NULL,
                                     intra_prob = 0.8,
                                     inter_prob = 0.2,
                                     seed = NULL) {
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
  intra_prob <- .validate_probability_vector(
    intra_prob,
    "intra_prob",
    n_communities
  )
  .assert_scalar_number(inter_prob, "inter_prob", lower = 0, upper = 1)
  if (!is.null(seed)) {
    .assert_whole_number(seed, "seed", minimum = 0L)
  }

  .with_preserved_seed(seed, function() {
    graph <- matrix(0L, nrow = n_nodes, ncol = n_nodes)
    node_indices <- seq_len(n_nodes)
    assignments <- rep.int(seq_len(n_communities), community_sizes)

    for (community in seq_len(n_communities)) {
      nodes <- node_indices[assignments == community]
      if (length(nodes) > 1L) {
        block <- matrix(0L, nrow = length(nodes), ncol = length(nodes))
        block[upper.tri(block)] <- stats::rbinom(
          sum(upper.tri(block)),
          1L,
          intra_prob[[community]]
        )
        graph[nodes, nodes] <- block + t(block)
      }
    }

    if (n_communities >= 2L) {
      for (first in seq_len(n_communities - 1L)) {
        for (second in seq.int(first + 1L, n_communities)) {
          first_nodes <- node_indices[assignments == first]
          second_nodes <- node_indices[assignments == second]
          block <- matrix(
            stats::rbinom(
              length(first_nodes) * length(second_nodes),
              1L,
              inter_prob
            ),
            nrow = length(first_nodes),
            ncol = length(second_nodes)
          )
          graph[first_nodes, second_nodes] <- block
          graph[second_nodes, first_nodes] <- t(block)
        }
      }
    }

    graph
  })
}
