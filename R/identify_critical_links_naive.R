# R/identify_critical_links_naive.R
# Naive (brute-force) reference implementation of identify_critical_links.
# ---------------------------------------------------------------------------
# This function recomputes the test statistic T from scratch after every edge
# removal and regenerates the full permutation null at each step.  Its
# complexity is O(K * B * |E| * m) -- the "naive" baseline described in the
# paper.  It is kept internally for benchmarking and as a correctness
# reference for the optimised `identify_critical_links()`; it is NOT exported.
#
#' @noRd
#' @keywords internal
identify_critical_links_naive <- function(populations,
                                          alpha          = 0.05,
                                          method         = "fisher",
                                          adjust_method  = "none",
                                          batch_size     = 1,
                                          n_permutations = 1000,
                                          a              = 1,
                                          seed           = NULL) {

  if (!is.null(seed)) set.seed(seed)
  if (!is.list(populations) || length(populations) < 2)
    stop("`populations` must be a list with at least 2 groups.")
  if (batch_size  < 1L) stop("`batch_size` must be >= 1.")
  if (n_permutations < 1L) stop("`n_permutations` must be >= 1.")

  ## 0. Metadata ------------------------------------------------------------
  Npop <- sapply(populations, length)
  m    <- length(populations)
  n    <- sum(Npop)

  ## 1. Rank edges by marginal p-value (same as optimised) -----------------
  freq       <- compute_edge_frequencies(populations)
  edge_pvals <- compute_edge_pvalues(freq$edge_counts, Npop,
                                     method        = method,
                                     adjust_method = adjust_method)
  edge_df    <- rank_edges(edge_pvals)
  n_edges    <- nrow(edge_df)
  if (n_edges == 0L)
    return(list(critical_edges       = NULL,
                edges_removed        = list(),
                modified_populations = populations))

  ## Helper: permutation p-value for the current populations ----------------
  .permutation_pvalue <- function(pops, Npop_vec, n_total, n_perm, a_val) {
    T_obs      <- compute_test_statistic(pops, a = a_val)
    all_graphs <- unlist(pops, recursive = FALSE)
    T_perm     <- numeric(n_perm)
    for (b in seq_len(n_perm)) {
      perm_idx  <- sample.int(n_total)
      cumN      <- c(0L, cumsum(Npop_vec))
      perm_pops <- lapply(seq_along(Npop_vec), function(k)
        all_graphs[perm_idx[(cumN[k] + 1L):cumN[k + 1L]]])
      names(perm_pops) <- names(pops)
      T_perm[b] <- compute_test_statistic(perm_pops, a = a_val)
    }
    mean(T_perm < T_obs)
  }

  ## 2. Global test ---------------------------------------------------------
  initial_p <- .permutation_pvalue(populations, Npop, n, n_permutations, a)
  if (initial_p > alpha) {
    warning("Initial test is not significant (p = ",
            round(initial_p, 4),
            "). Populations may be identical or differences too small to detect.")
    return(list(critical_edges       = NULL,
                edges_removed        = list(),
                modified_populations = populations))
  }

  ## 3. Iterative edge removal (naive: full recomputation each step) ------
  edges_removed <- list()
  k_removed     <- 0L
  continue      <- TRUE

  while (continue && k_removed < n_edges) {
    batch_end <- min(k_removed + batch_size, n_edges)
    idx_batch <- (k_removed + 1L):batch_end

    ## Physically remove edges from all graphs -----------------------------
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

    ## Recompute T and permutation null FROM SCRATCH -----------------------
    p_k      <- .permutation_pvalue(populations, Npop, n, n_permutations, a)
    continue <- (p_k <= alpha)
  }

  critical_edges <- if (k_removed > 0)
                      edge_df[seq_len(k_removed), ]
                    else
                      NULL

  list(critical_edges       = critical_edges,
       edges_removed        = edges_removed,
       modified_populations = populations)
}
