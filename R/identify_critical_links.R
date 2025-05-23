# R/identify_critical_links.R
# Fast, mathematically exact implementation
# ──────────────────────────────────────────────────────────────────────────────

#' Identify critical edges that explain population differences
#' @inheritParams compute_edge_pvalues
#' @param populations Named list; each element is itself a list of adjacency
#'   matrices (binary, symmetric, zero diagonal).
#' @param alpha Significance level (global + iterative tests).
#' @param batch_size How many edges to drop per iteration (≥ 1).
#' @param n_bootstrap Size of the bootstrap null distribution.
#' @param a Normalisation constant used in the original statistic.
#' @param seed Optional integer for reproducibility.
#' @return List with `critical_edges`, `edges_removed`, `modified_populations`.
#' @export
identify_critical_links <- function(populations,
                                    alpha         = 0.05,
                                    method        = "fisher",
                                    adjust_method = "none",
                                    batch_size    = 1,
                                    n_bootstrap   = 1000,
                                    a             = 1,
                                    seed          = NULL) {

  if (!is.null(seed)) set.seed(seed)
  if (!is.list(populations) || length(populations) < 2)
    stop("`populations` must be a list with ≥ 2 groups.")
  if (batch_size  < 1L) stop("`batch_size` must be ≥ 1.")
  if (n_bootstrap < 1L) stop("`n_bootstrap` must be ≥ 1.")

  ## 0 ── metadata ------------------------------------------------------------
  Npop <- sapply(populations, length)      # n_k
  m    <- length(populations)
  n    <- sum(Npop)

  make_key <- function(i, j) paste(i, j, sep = "-")

  ## 1 ── rank edges by marginal p-value --------------------------------------
  freq       <- compute_edge_frequencies(populations)
  edge_pvals <- compute_edge_pvalues(freq$edge_counts, Npop,
                                     method        = method,
                                     adjust_method = adjust_method)
  edge_df    <- rank_edges(edge_pvals)               # ascending p
  n_edges    <- nrow(edge_df)
  if (n_edges == 0L)
    return(list(critical_edges       = NULL,
                edges_removed        = list(),
                modified_populations = populations))

  ## 2 ── Δₑ helper -----------------------------------------------------------
  .edge_deltas <- function(edge_counts) {
    n_nodes <- dim(edge_counts)[1]
    idx     <- which(upper.tri(matrix(0, n_nodes, n_nodes)), arr.ind = TRUE)
    if (nrow(idx) == 0L)
      return(list(indices = idx, deltas = numeric(0)))

    ## counts[#edges × m] – number of 1's per edge & population
    counts <- do.call(cbind, lapply(seq_len(m),
                     function(k) edge_counts[ , , k][idx]))
    # for a single edge counts is 1 × m (matrix, not vector)
    if (!is.matrix(counts))
      counts <- matrix(counts, nrow = 1L)

    ## probabilities ----------------------------------------------------------
    p_mat <- sweep(counts, 2, Npop, "/")            # p_{k,e}
    p_tot <- rowSums(counts) / n                    # p_e
    Ptot  <- matrix(p_tot, nrow = nrow(counts), ncol = m)

    ## within- and cross-group deviations -------------------------------------
    d_mat <- 2 * p_mat * (1 - p_mat)                # d_{k,e}
    D_mat <- p_mat + Ptot - 2 * p_mat * Ptot        # D_{k,e}

    coef_d <- sqrt(Npop) * (Npop / (Npop - 1))
    coef_D <- sqrt(Npop) * (n    / (n    - 1))

    delta  <- (sqrt(m) / a) *
              (d_mat %*% coef_d - D_mat %*% coef_D)[, 1]

    list(indices = idx, deltas = delta)
  }

  obs_delta_info <- .edge_deltas(freq$edge_counts)

  key_all   <- make_key(obs_delta_info$indices[, 1], obs_delta_info$indices[, 2])
  map_idx   <- match(make_key(edge_df$node1, edge_df$node2), key_all)
  delta_ord <- obs_delta_info$deltas[map_idx]

  prefix_obs <- cumsum(delta_ord)          # T after k removals (obs)
  T0         <- sum(obs_delta_info$deltas)

  ## 3 ── bootstrap -----------------------------------------------------------
  all_graphs  <- unlist(populations, recursive = FALSE)
  boot_deltas <- matrix(0, nrow = n_bootstrap, ncol = n_edges)

  for (b in seq_len(n_bootstrap)) {
    boot_pops <- lapply(Npop, function(sz)
      all_graphs[sample.int(n, sz, TRUE)])

    boot_counts       <- compute_edge_frequencies(boot_pops)$edge_counts
    boot_deltas[b, ]  <- .edge_deltas(boot_counts)$deltas[map_idx]
  }
  T_boot0     <- rowSums(boot_deltas)
  prefix_boot <- t(apply(boot_deltas, 1, cumsum))   # B × n_edges

  ## 4 ── global test ---------------------------------------------------------
  initial_p <- mean(abs(T_boot0) >= abs(T0))        # two-tailed
  if (initial_p > alpha) {
    warning("Initial test is not significant (p = ",
            round(initial_p, 4),
            "). Populations may be identical or differences too small to detect.")
    return(list(critical_edges       = NULL,
                edges_removed        = list(),
                modified_populations = populations))
  }

  ## 5 ── iterative edge removal ---------------------------------------------
  edges_removed <- list()
  k_removed     <- 0L
  continue      <- TRUE

  while (continue && k_removed < n_edges) {
    batch_end <- min(k_removed + batch_size, n_edges)

    T_obs_k  <- T0 - prefix_obs[batch_end]
    T_boot_k <- T_boot0 - prefix_boot[ , batch_end]
    p_k      <- mean(abs(T_boot_k) >= abs(T_obs_k))   # two-tailed

    ## physically remove those edges -----------------------------------------
    idx_batch <- (k_removed + 1L):batch_end
    for (idx in idx_batch) {
      i <- edge_df$node1[idx]; j <- edge_df$node2[idx]
      for (pop in seq_along(populations))
        for (g in seq_len(Npop[pop])) {
          populations[[pop]][[g]][i, j] <- 0
          populations[[pop]][[g]][j, i] <- 0
        }
      edges_removed[[length(edges_removed) + 1L]] <- c(i, j)
    }

    k_removed <- batch_end
    continue  <- (p_k <= alpha)
  }

  critical_edges <- if (k_removed > 0)
                      edge_df[seq_len(k_removed), ]
                    else
                      NULL

  list(critical_edges       = critical_edges,
       edges_removed        = edges_removed,
       modified_populations = populations)
}
