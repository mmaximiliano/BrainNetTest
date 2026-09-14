# R/permutation-null.R
# Internal building blocks of the permutation test: the edge-wise
# decomposition of T and the vectorised generation of the permutation null.
# Shared by global_test() and identify_critical_links(), so that the two
# evaluate exactly the same statistic on exactly the same permutations.
# None of these are exported.

#' Upper-triangle edge index of a graph on `n_nodes` nodes
#'
#' @return An |E| x 2 integer matrix of (i, j) pairs with i < j, in the
#'   column-major order that `upper.tri()` produces.
#' @noRd
#' @keywords internal
.edge_index <- function(n_nodes) {
  which(upper.tri(matrix(0, n_nodes, n_nodes)), arr.ind = TRUE)
}

#' Edge-wise contributions Delta_e of the observed data
#'
#' Evaluates the decomposition T = sum_e Delta_e for the observed edge
#' counts, so that sum(deltas) is the test statistic of
#' `compute_test_statistic()` (for binary graphs) without any distance
#' computation.
#'
#' @param edge_counts `n_nodes x n_nodes x m` array of per-population edge
#'   counts, as returned in `compute_edge_frequencies()$edge_counts`.
#' @param group_sizes Integer vector of the m population sizes.
#' @param a Normalisation constant of T.
#' @return A list with `indices`, the |E| x 2 edge index, and `deltas`, the
#'   numeric vector of Delta_e in the same order.
#' @noRd
#' @keywords internal
.observed_edge_deltas <- function(edge_counts, group_sizes, a = 1) {
  m   <- length(group_sizes)
  n   <- sum(group_sizes)
  idx <- .edge_index(dim(edge_counts)[1])
  if (nrow(idx) == 0L)
    return(list(indices = idx, deltas = numeric(0)))

  counts <- do.call(cbind, lapply(seq_len(m),
                   function(k) edge_counts[ , , k][idx]))
  if (!is.matrix(counts))
    counts <- matrix(counts, nrow = 1L)

  p_mat <- sweep(counts, 2, group_sizes, "/")
  p_tot <- rowSums(counts) / n
  Ptot  <- matrix(p_tot, nrow = nrow(counts), ncol = m)

  d_mat <- 2 * p_mat * (1 - p_mat)
  D_mat <- p_mat + Ptot - 2 * p_mat * Ptot

  coef_d <- sqrt(group_sizes) * (group_sizes / (group_sizes - 1))
  coef_D <- sqrt(group_sizes) * (n           / (n           - 1))

  delta  <- (sqrt(m) / a) *
            (d_mat %*% coef_d - D_mat %*% coef_D)[, 1]

  list(indices = idx, deltas = delta)
}

#' Edge-wise contributions Delta_e for every permutation replicate
#'
#' The same formulas as `.observed_edge_deltas()`, evaluated column-wise on
#' the |E| x B count matrices of all replicates at once.
#'
#' @param counts_list List of m matrices, each |E| x B, giving the edge
#'   counts of population k under each of the B permutations.
#' @inheritParams .observed_edge_deltas
#' @return An |E| x B matrix whose column b holds the Delta_e of replicate b.
#' @noRd
#' @keywords internal
.permuted_edge_deltas <- function(counts_list, group_sizes, a = 1) {
  m <- length(group_sizes)
  n <- sum(group_sizes)
  E <- nrow(counts_list[[1]])
  B <- ncol(counts_list[[1]])

  p_list  <- lapply(seq_len(m), function(k) counts_list[[k]] / group_sizes[k])
  C_total <- Reduce("+", counts_list)
  p_tot   <- C_total / n

  coef_d <- sqrt(group_sizes) * (group_sizes / (group_sizes - 1))
  coef_D <- sqrt(group_sizes) * (n           / (n           - 1))

  delta_mat <- matrix(0, nrow = E, ncol = B)
  for (k in seq_len(m)) {
    pk  <- p_list[[k]]
    d_k <- 2 * pk * (1 - pk)
    D_k <- pk + p_tot - 2 * pk * p_tot
    delta_mat <- delta_mat + coef_d[k] * d_k - coef_D[k] * D_k
  }
  (sqrt(m) / a) * delta_mat
}

#' Graph-by-edge matrix X of the pooled sample
#'
#' Row g holds the upper-triangle entries of the g-th graph, populations
#' concatenated in order, so that `X[g, e] = A_e(G_g)`.
#'
#' @param populations The populations list.
#' @param indices The |E| x 2 edge index from `.edge_index()`.
#' @return An n x |E| numeric (double) matrix.
#' @noRd
#' @keywords internal
.graph_edge_matrix <- function(populations, indices) {
  all_graphs <- unlist(populations, recursive = FALSE)
  X <- do.call(rbind, lapply(all_graphs, function(A) A[indices]))
  storage.mode(X) <- "double"
  X
}

#' Edge counts of every population under B random label permutations
#'
#' Draws B permutations of the pooled graph labels, preserving the group
#' sizes, and returns the permuted edge counts as the matrix products
#' C_k = X' W^(k), where W^(k) is the n x B indicator matrix of the graphs
#' assigned to population k in each replicate. This is the only place in the
#' permutation phase that consumes random numbers: one `sample.int()` draw
#' per replicate.
#'
#' @param X The n x |E| matrix from `.graph_edge_matrix()`.
#' @param group_sizes Integer vector of the m population sizes.
#' @param n_permutations Number of replicates B.
#' @return A list of m matrices, each |E| x B.
#' @noRd
#' @keywords internal
.permuted_edge_counts <- function(X, group_sizes, n_permutations) {
  m <- length(group_sizes)
  n <- sum(group_sizes)

  perm_mat <- replicate(n_permutations, sample.int(n))
  cum_sizes <- c(0L, cumsum(group_sizes))
  W_list <- lapply(seq_len(m), function(k) {
    rows    <- (cum_sizes[k] + 1L):cum_sizes[k + 1L]
    idx_all <- as.vector(perm_mat[rows, , drop = FALSE])
    rep_id  <- rep(seq_len(n_permutations), each = group_sizes[k])
    lin_idx <- (rep_id - 1L) * n + idx_all
    W <- matrix(tabulate(lin_idx, nbins = n * n_permutations),
                nrow = n, ncol = n_permutations)
    storage.mode(W) <- "double"
    W
  })

  lapply(W_list, function(W) crossprod(X, W))
}
