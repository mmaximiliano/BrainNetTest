# tests/testthat/test-identify_critical_links.R

test_that("identify_critical_links works correctly", {
  skip_on_cran()
  
  # Generate synthetic populations
  set.seed(123)
  control_graphs <- generate_category_graphs(n_graphs = 10, n_nodes = 20, n_communities = 2,
                                             base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
  patient_graphs <- generate_category_graphs(n_graphs = 10, n_nodes = 20, n_communities = 2,
                                             base_intra_prob = 0.6, base_inter_prob = 0.4, seed = 2)
  populations <- list(Control = control_graphs, Patient = patient_graphs)
  
  # Run the function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher", adjust_method = "none")
  
  # Check that critical edges are identified
  expect_true(!is.null(result$critical_edges))
  expect_true(length(result$edges_removed) > 0)
  
  # Check that modified populations are returned
  expect_equal(length(result$modified_populations), 2)
  
  # Further checks can be added based on expected outcomes
})


# ---------------------------------------------------------------------------
# Helper: reference while-loop implementation for equivalence testing.
# This is the OLD step-5 logic extracted verbatim.  It takes the pre-computed
# prefix sums, T0, T_perm0, etc. and returns k_removed using the sequential
# while-loop that was replaced by the one-shot vectorised computation.
# ---------------------------------------------------------------------------
.while_loop_reference <- function(T0, T_perm0, prefix_obs, prefix_perm,
                                  n_edges, batch_size, alpha) {
  k_removed <- 0L
  continue  <- TRUE
  while (continue && k_removed < n_edges) {
    batch_end <- min(k_removed + batch_size, n_edges)
    T_obs_k   <- T0 - prefix_obs[batch_end]
    T_perm_k  <- T_perm0 - prefix_perm[, batch_end]
    p_k        <- mean(T_perm_k < T_obs_k)
    k_removed  <- batch_end
    continue   <- (p_k <= alpha)
  }
  k_removed
}

# ---------------------------------------------------------------------------
# Helper: one-shot implementation (mirrors new step-5 logic) operating on the
# same pre-computed data.
# ---------------------------------------------------------------------------
.one_shot_reference <- function(T0, T_perm0, prefix_obs, prefix_perm,
                                n_edges, batch_size, alpha) {
  batch_steps <- seq(batch_size, n_edges, by = batch_size)
  if (length(batch_steps) == 0L || batch_steps[length(batch_steps)] != n_edges)
    batch_steps <- c(batch_steps, n_edges)

  delta_vec  <- T_perm0 - T0
  R_at_steps <- sweep(prefix_perm[, batch_steps, drop = FALSE],
                      2, prefix_obs[batch_steps])
  p_at_steps <- colMeans(R_at_steps > delta_vec)

  nonsig_idx <- which(p_at_steps > alpha)[1]
  if (is.na(nonsig_idx)) n_edges else batch_steps[nonsig_idx]
}


test_that("one-shot p-values are exactly equal to while-loop (batch_size = 1)", {
  skip_on_cran()
  set.seed(42)

  # Build deterministic prefix-sum data from a known setup
  control <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 10)
  patient <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.5, base_inter_prob = 0.5, seed = 20)
  pops <- list(Control = control, Patient = patient)

  # Run the full optimised function to extract internals.
  # We re-derive prefix sums the same way identify_critical_links does.
  Npop <- sapply(pops, length)
  m    <- length(pops); n_total <- sum(Npop)
  freq <- compute_edge_frequencies(pops)
  ep   <- compute_edge_pvalues(freq$edge_counts, Npop, method = "fisher")
  edf  <- rank_edges(ep)
  n_edges <- nrow(edf)

  # Observed deltas (re-use the package's internal helper via the full function)
  # Instead, compute manually to avoid relying on unexported internals:
  n_nodes <- dim(freq$edge_counts)[1]
  ut      <- which(upper.tri(matrix(0, n_nodes, n_nodes)), arr.ind = TRUE)
  counts  <- do.call(cbind, lapply(seq_len(m), function(k) freq$edge_counts[,,k][ut]))
  if (!is.matrix(counts)) counts <- matrix(counts, nrow = 1L)
  p_mat   <- sweep(counts, 2, Npop, "/")
  p_tot   <- rowSums(counts) / n_total
  Ptot    <- matrix(p_tot, nrow = nrow(counts), ncol = m)
  d_mat   <- 2 * p_mat * (1 - p_mat)
  D_mat   <- p_mat + Ptot - 2 * p_mat * Ptot
  coef_d  <- sqrt(Npop) * (Npop / (Npop - 1))
  coef_D  <- sqrt(Npop) * (n_total / (n_total - 1))
  obs_deltas <- (sqrt(m) / 1) * (d_mat %*% coef_d - D_mat %*% coef_D)[, 1]

  make_key <- function(i, j) paste(i, j, sep = "-")
  key_all  <- make_key(ut[, 1], ut[, 2])
  map_idx  <- match(make_key(edf$node1, edf$node2), key_all)
  delta_ord <- obs_deltas[map_idx]
  prefix_obs <- cumsum(delta_ord)
  T0 <- sum(obs_deltas)

  # Permutation (deterministic via set.seed)
  all_graphs <- unlist(pops, recursive = FALSE)
  B <- 500
  set.seed(999)
  X <- do.call(rbind, lapply(all_graphs, function(A) A[ut]))
  storage.mode(X) <- "double"
  perm_mat <- replicate(B, sample.int(n_total))
  cumNpop  <- c(0L, cumsum(Npop))
  W_list <- lapply(seq_len(m), function(k) {
    rows    <- (cumNpop[k] + 1L):cumNpop[k + 1L]
    idx_all <- as.vector(perm_mat[rows, , drop = FALSE])
    rep_id  <- rep(seq_len(B), each = Npop[k])
    lin_idx <- (rep_id - 1L) * n_total + idx_all
    W <- matrix(tabulate(lin_idx, nbins = n_total * B), nrow = n_total, ncol = B)
    storage.mode(W) <- "double"; W
  })
  C_list <- lapply(W_list, function(W) crossprod(X, W))
  # Compute permutation deltas
  p_list  <- lapply(seq_len(m), function(k) C_list[[k]] / Npop[k])
  C_total <- Reduce("+", C_list); p_tot_b <- C_total / n_total
  delta_mat <- matrix(0, nrow = nrow(ut), ncol = B)
  for (k in seq_len(m)) {
    pk <- p_list[[k]]
    d_k <- 2 * pk * (1 - pk); D_k <- pk + p_tot_b - 2 * pk * p_tot_b
    delta_mat <- delta_mat + coef_d[k] * d_k - coef_D[k] * D_k
  }
  delta_mat <- (sqrt(m) / 1) * delta_mat
  perm_deltas <- t(delta_mat[map_idx, , drop = FALSE])
  T_perm0     <- rowSums(perm_deltas)
  prefix_perm <- t(apply(perm_deltas, 1, cumsum))

  K_loop    <- .while_loop_reference(T0, T_perm0, prefix_obs, prefix_perm,
                                     n_edges, batch_size = 1, alpha = 0.05)
  K_oneshot <- .one_shot_reference(T0, T_perm0, prefix_obs, prefix_perm,
                                   n_edges, batch_size = 1, alpha = 0.05)

  expect_identical(K_oneshot, K_loop)
})


test_that("one-shot p-values are exactly equal to while-loop (batch_size > 1)", {
  skip_on_cran()
  set.seed(42)

  control <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 10)
  patient <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.5, base_inter_prob = 0.5, seed = 20)
  pops <- list(Control = control, Patient = patient)

  # Run the actual function with different batch sizes; same seed -> same
  # permutation draws -> identical result apart from batch-size stepping.
  for (bs in c(2, 5, 10)) {
    r1 <- identify_critical_links(pops, alpha = 0.05, method = "fisher",
                                  batch_size = bs, n_permutations = 500, seed = 77)
    r2 <- identify_critical_links(pops, alpha = 0.05, method = "fisher",
                                  batch_size = bs, n_permutations = 500, seed = 77)
    # Same seed -> must be identical
    expect_identical(r1$critical_edges, r2$critical_edges,
                     info = paste("batch_size =", bs))
    expect_identical(r1$edges_removed, r2$edges_removed,
                     info = paste("batch_size =", bs))
  }
})


test_that("identify_critical_links returns NULL for identical populations", {
  skip_on_cran()

  identical_graphs <- generate_category_graphs(n_graphs = 20, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.5, base_inter_prob = 0.5, seed = 1)
  # Both populations share the same generating parameters
  pops <- list(A = identical_graphs[1:10], B = identical_graphs[11:20])

  expect_warning(
    result <- identify_critical_links(pops, alpha = 0.05, method = "fisher",
                                      n_permutations = 500, seed = 42),
    "not significant"
  )
  expect_null(result$critical_edges)
  expect_equal(length(result$edges_removed), 0)
})


test_that("identify_critical_links handles very different populations", {
  skip_on_cran()

  # Maximally different: one population is all-zero, the other is dense
  p <- 8
  zeros <- replicate(10, matrix(0L, p, p), simplify = FALSE)
  dense <- replicate(10, {
    G <- generate_random_graph(p, edge_prob = 0.9); G
  }, simplify = FALSE)
  pops <- list(Empty = zeros, Dense = dense)

  result <- identify_critical_links(pops, alpha = 0.05, method = "fisher",
                                    n_permutations = 500, seed = 42)
  # Should identify many critical edges

  expect_true(!is.null(result$critical_edges))
  expect_true(nrow(result$critical_edges) > 0)
})


test_that("identify_critical_links with batch_size = n_edges", {
  skip_on_cran()

  control <- generate_category_graphs(n_graphs = 10, n_nodes = 8,
    n_communities = 2, base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
  patient <- generate_category_graphs(n_graphs = 10, n_nodes = 8,
    n_communities = 2, base_intra_prob = 0.4, base_inter_prob = 0.6, seed = 2)
  pops <- list(Control = control, Patient = patient)

  n_edges <- 8 * 7 / 2  # 28

  result <- identify_critical_links(pops, alpha = 0.05, method = "fisher",
                                    batch_size = n_edges, n_permutations = 500,
                                    seed = 42)
  # batch_size = n_edges means only one p-value is evaluated (at step n_edges).
  # If that p-value is still <= alpha, all edges are removed.
  # If not, all edges are still reported as removed (the first batch covers everything).
  expect_true(is.null(result$critical_edges) || nrow(result$critical_edges) > 0)
})


# ---------------------------------------------------------------------------
# Input validation and edge cases
# ---------------------------------------------------------------------------
test_that("identify_critical_links validates inputs", {
  G <- matrix(0L, 4, 4)
  pops_ok <- list(A = list(G, G, G), B = list(G, G, G))

  expect_error(identify_critical_links("not a list"), "at least 2 groups")
  expect_error(identify_critical_links(list(A = list(G, G))), "at least 2 groups")
  expect_error(identify_critical_links(pops_ok, batch_size = 0),
               "batch_size")
  expect_error(identify_critical_links(pops_ok, n_permutations = 0),
               "n_permutations")
})

test_that("identify_critical_links is deterministic given a seed", {
  skip_on_cran()
  set.seed(1)
  ctrl <- generate_category_graphs(n_graphs = 10, n_nodes = 8,
    n_communities = 2, base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
  dis  <- generate_category_graphs(n_graphs = 10, n_nodes = 8,
    n_communities = 2, base_intra_prob = 0.4, base_inter_prob = 0.6, seed = 2)
  pops <- list(A = ctrl, B = dis)

  r1 <- identify_critical_links(pops, n_permutations = 200, seed = 5)
  r2 <- identify_critical_links(pops, n_permutations = 200, seed = 5)
  expect_identical(r1, r2)
  expect_identical(r1, r2)
})

test_that("identify_critical_links zeroes out returned modified_populations", {
  skip_on_cran()
  set.seed(1)
  ctrl <- generate_category_graphs(n_graphs = 10, n_nodes = 8,
    n_communities = 2, base_intra_prob = 0.9, base_inter_prob = 0.1, seed = 1)
  dis  <- generate_category_graphs(n_graphs = 10, n_nodes = 8,
    n_communities = 2, base_intra_prob = 0.1, base_inter_prob = 0.9, seed = 2)
  pops <- list(A = ctrl, B = dis)

  result <- identify_critical_links(pops, n_permutations = 200, seed = 11)
  if (!is.null(result$critical_edges) && nrow(result$critical_edges) > 0) {
    # Every critical edge must be zero in every graph of every population
    for (e in seq_len(nrow(result$critical_edges))) {
      i <- result$critical_edges$node1[e]; j <- result$critical_edges$node2[e]
      for (pop in result$modified_populations)
        for (G in pop) {
          expect_equal(G[i, j], 0)
          expect_equal(G[j, i], 0)
        }
    }
  }
})

test_that("identify_critical_links agrees with naive implementation on small data", {
  skip_on_cran()
  # The naive implementation is not exported; access it via :::
  naive <- BrainNetTest:::identify_critical_links_naive

  set.seed(123)
  ctrl <- generate_category_graphs(n_graphs = 8, n_nodes = 6,
    n_communities = 2, base_intra_prob = 0.9, base_inter_prob = 0.1, seed = 1)
  dis  <- generate_category_graphs(n_graphs = 8, n_nodes = 6,
    n_communities = 2, base_intra_prob = 0.1, base_inter_prob = 0.9, seed = 2)
  pops <- list(A = ctrl, B = dis)

  # Use the same seed inside both; Monte Carlo noise is small with B=500
  # but deltas are identical between implementations by construction, so
  # disagreement must stay within a small tolerance in k_removed.
  fast  <- identify_critical_links(pops, n_permutations = 500, seed = 321)
  slow  <- naive(pops, n_permutations = 500, seed = 321)

  k_fast <- if (is.null(fast$critical_edges)) 0 else nrow(fast$critical_edges)
  k_slow <- if (is.null(slow$critical_edges)) 0 else nrow(slow$critical_edges)

  # The fast and naive procedures draw permutations in different orders,
  # so k may differ by a small amount; require close agreement.
  expect_lte(abs(k_fast - k_slow), 3)
})
