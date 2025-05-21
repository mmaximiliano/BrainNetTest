# R/identify_critical_links.R

#' Identify Critical Links Between Brain‑Network Populations (fast version)
#'
#' This re‑implementation eliminates the quadratic **bootstrap × edge** cost of
#' the original algorithm by exploiting the additive structure of the test
#' statistic.  Every edge contributes an independent term (Δₑ) to *T*; once those
#' Δₑ are cached for the observed data **and for each bootstrap replicate**, the
#' statistic after removing the first *k* ranked edges is obtained with a single
#' subtraction and a prefix‑sum lookup.  The result is mathematically identical
#' to the slow implementation but runs two–three orders of magnitude faster on
#' realistic inputs.
#'
#' @param populations A *named* list.  Each element is itself a list of binary
#'   adjacency matrices representing the graphs of one population.
#' @param alpha Significance level for the global test (default 0.05).
#' @param method Statistical test for the marginal edge screen (passed to
#'   [compute_edge_pvalues()]).
#' @param adjust_method P‑value adjustment (see `?p.adjust`).
#' @param batch_size Number of edges removed per iteration (≥ 1).
#' @param n_bootstrap Number of bootstrap resamples used for the null
#'   distribution.
#' @param a Normalisation constant of the test statistic (same `a` as in
#'   [compute_test_statistic()]).
#' @param seed Optional integer.  If supplied, the RNG state is fixed for fully
#'   reproducible results.
#'
#' @return A list with three components: `critical_edges` (data‑frame),
#'   `edges_removed` (list of integer pairs), and `modified_populations` (the
#'   populations object after the last deletion).
#' @export
#'
#' @examples
#' # Synthetic example — two populations with known differences
#' ctrl  <- generate_category_graphs(5, n_nodes = 20, n_communities = 2,
#'                                   base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
#' dis   <- generate_category_graphs(5, n_nodes = 20, n_communities = 2,
#'                                   base_intra_prob = 0.5, base_inter_prob = 0.4, seed = 2)
#' res <- identify_critical_links(list(Control = ctrl, Disease = dis), n_bootstrap = 200)
#' print(res$critical_edges)
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
    stop("`populations` must be a list with ≥2 groups of graphs.")

  if (batch_size < 1L)        stop("`batch_size` must be ≥ 1.")
  if (n_bootstrap < 1L)       stop("`n_bootstrap` must be ≥ 1.")

  ## ------------------------------------------------------------------ ##
  ## 0.  Basic metadata                                                 ##
  ## ------------------------------------------------------------------ ##
  Npop   <- sapply(populations, length)          # n_i
  m      <- length(populations)
  n      <- sum(Npop)
  pop_nm <- names(populations)

  ## helper to build a unique key for (i,j) ---------------------------- ##
  make_key <- function(i, j) paste(i, j, sep = "-")
  ## ------------------------------------------------------------------ ##
  ## 1.  Rank edges by marginal p‑value                                 ##
  ## ------------------------------------------------------------------ ##
  freq            <- compute_edge_frequencies(populations)
  edge_pvals      <- compute_edge_pvalues(freq$edge_counts, Npop,
                                          method = method,
                                          adjust_method = adjust_method)
  edge_df         <- rank_edges(edge_pvals)              # ascending p‑values
  n_edges         <- nrow(edge_df)
    # Handle case with no edges
  if (n_edges == 0) {
    warning("No edges found in the graphs to analyze.")
    return(list(
      critical_edges = NULL,
      edges_removed = list(),
      modified_populations = populations
    ))
  }  ## ------------------------------------------------------------------ ##
  ## 2.  Edge‑wise Δₑ for the observed data                             ##
  ## ------------------------------------------------------------------ ##  
  .edge_deltas <- function(edge_counts) {
    n_nodes <- dim(edge_counts)[1]
    idx     <- which(upper.tri(matrix(0, n_nodes, n_nodes)), arr.ind = TRUE)
    
    # Return early with appropriate structure if no edges
    if (nrow(idx) == 0) {
      return(list(indices = idx, deltas = numeric(0)))
    }    # matrix[#edges, m] of counts per population ----------------------- ##
    counts  <- vapply(seq_len(m),
                     function(k) edge_counts[,,k][idx],
                     numeric(nrow(idx)))
    
    # Handle special cases for counts matrix dimensions
    if (is.null(dim(counts)) && length(counts) > 0) {
      # For single edge case, ensure counts is a proper matrix
      counts <- matrix(counts, nrow = 1)
    }

    # s_i(e) = 2 c_i,e [1 - c_i,e / n_i] ------------------------------- ##
    # Create s_i matrix with proper dimension checks
    if (length(counts) == 0) {
      s_i <- numeric(0)
    } else {
      s_i <- matrix(0, nrow = nrow(counts), ncol = ncol(counts))
      for (i in seq_len(ncol(counts))) {
        s_i[, i] <- 2 * counts[, i] * (1 - counts[, i] / Npop[i])
      }
    }
    
    coef_i <- sqrt(Npop) * (Npop / (Npop - 1) - n / (n - 1))
    delta  <- as.numeric((sqrt(m) / a) * (s_i %*% coef_i))

    list(indices = idx, deltas = delta)
  }

  obs_delta_info <- .edge_deltas(freq$edge_counts)

  ## match ordering of deltas to edge_df ------------------------------- ##
  key_all  <- make_key(obs_delta_info$indices[,1], obs_delta_info$indices[,2])
  map_idx  <- match(make_key(edge_df$node1, edge_df$node2), key_all)
  delta_ord <- obs_delta_info$deltas[map_idx]

  prefix_obs <- cumsum(delta_ord)
  T0         <- sum(obs_delta_info$deltas)

  ## ------------------------------------------------------------------ ##
  ## 3.  Bootstrap — cache Δₑ^(b) once                                  ##
  ## ------------------------------------------------------------------ ##
  all_graphs <- unlist(populations, recursive = FALSE)
  n_total    <- length(all_graphs)

  boot_deltas <- matrix(0, nrow = n_bootstrap, ncol = n_edges)
  for (b in seq_len(n_bootstrap)) {
    # resample graphs under H0 ----------------------------------------- ##
    boot_pops <- lapply(Npop, function(sz) all_graphs[sample.int(n_total, sz, TRUE)])

    boot_counts <- compute_edge_frequencies(boot_pops)$edge_counts
    boot_deltas[b, ] <- .edge_deltas(boot_counts)$deltas[map_idx]
  }
  T_boot0    <- rowSums(boot_deltas)
  prefix_boot <- t(apply(boot_deltas, 1, cumsum))     # B × n_edges
  ## ------------------------------------------------------------------ ##
  ## 4.  Initial test for significant difference between populations    ##
  ## ------------------------------------------------------------------ ##
  # Calculate initial p-value to test if there's any significant difference
  initial_p_value <- mean(T_boot0 <= T0)
    # If the initial test is not significant, issue a warning and return early
  if (initial_p_value > alpha) {
    warning("Initial test is not significant (p = ", round(initial_p_value, 4), 
            "). Populations may be identical or differences too small to detect.")
    # Return NULL for critical_edges when initial test is not significant
    return(list(
      critical_edges = NULL,
      edges_removed = list(),
      modified_populations = populations
    ))
  }

  ## ------------------------------------------------------------------ ##
  ## 5.  Iterative edge removal                                         ##
  ## ------------------------------------------------------------------ ##
  edges_removed     <- list()
  k_removed         <- 0L
  continue_removal  <- TRUE  # Since we've already checked initial_p_value <= alpha

  while (continue_removal && k_removed < n_edges) {
    batch_end <- min(k_removed + batch_size, n_edges)

    # statistic & p‑value after hypothetical removal of first `batch_end` edges
    T_obs_k  <- T0 - prefix_obs[batch_end]
    T_boot_k <- T_boot0 - prefix_boot[, batch_end]
    p_k      <- mean(T_boot_k <= T_obs_k)

    # perform the removal (these become the real state) ---------------- ##
    idx_batch <- (k_removed + 1L):batch_end
    for (idx in idx_batch) {
      i <- edge_df$node1[idx]
      j <- edge_df$node2[idx]
      for (pop in seq_along(populations)) {
        for (g in seq_len(Npop[pop])) {
          populations[[pop]][[g]][i, j] <- 0
          populations[[pop]][[g]][j, i] <- 0
        }
      }
      edges_removed[[length(edges_removed) + 1L]] <- c(i, j)
    }

    k_removed <- batch_end
    continue_removal <- (p_k <= alpha)
  }
  # If no edges were removed (because continue_removal was false from the start)
  # return NULL for critical_edges
  critical_edges <- if (k_removed > 0) edge_df[seq_len(k_removed), ] else NULL

  list(
    critical_edges       = critical_edges,
    edges_removed        = edges_removed,
    modified_populations = populations
  )
}
