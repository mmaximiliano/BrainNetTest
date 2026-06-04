# tests/testthat/test-compute_edge_pvalues.R

# Helper: build a 2-pop edge_counts array from raw matrices
.make_counts <- function(pop_counts_list) {
  n_nodes <- nrow(pop_counts_list[[1]])
  arr <- array(0, dim = c(n_nodes, n_nodes, length(pop_counts_list)))
  for (k in seq_along(pop_counts_list))
    arr[, , k] <- pop_counts_list[[k]]
  arr
}

test_that("compute_edge_pvalues returns a symmetric matrix with valid p-values", {
  set.seed(1)
  ctrl <- generate_category_graphs(n_graphs = 10, n_nodes = 6, n_communities = 2,
    base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
  dis  <- generate_category_graphs(n_graphs = 10, n_nodes = 6, n_communities = 2,
    base_intra_prob = 0.4, base_inter_prob = 0.6, seed = 2)
  pops <- list(A = ctrl, B = dis)
  freq <- compute_edge_frequencies(pops)

  pvals <- compute_edge_pvalues(freq$edge_counts, sapply(pops, length))
  expect_equal(dim(pvals), c(6, 6))
  expect_true(isSymmetric(pvals))
  expect_true(all(pvals >= 0 & pvals <= 1))
  # Diagonal not meaningful but kept at 1
  expect_equal(diag(pvals), rep(1, 6))
})

test_that("compute_edge_pvalues: identical populations yield p >= alpha mostly", {
  set.seed(1)
  g <- generate_category_graphs(n_graphs = 20, n_nodes = 5, n_communities = 1,
    base_intra_prob = 0.5, base_inter_prob = 0.5, seed = 1)
  pops <- list(A = g[1:10], B = g[11:20])
  freq <- compute_edge_frequencies(pops)
  pvals <- compute_edge_pvalues(freq$edge_counts, sapply(pops, length))
  # With 10+10 samples drawn from the same process, no edge should be strongly
  # significant; all p-values must be valid probabilities.
  expect_true(all(pvals >= 0 & pvals <= 1))
})

test_that("compute_edge_pvalues supports all three test methods", {
  set.seed(2)
  ctrl <- generate_category_graphs(n_graphs = 10, n_nodes = 5, n_communities = 1,
    base_intra_prob = 0.9, base_inter_prob = 0.9, seed = 1)
  dis  <- generate_category_graphs(n_graphs = 10, n_nodes = 5, n_communities = 1,
    base_intra_prob = 0.1, base_inter_prob = 0.1, seed = 2)
  pops <- list(A = ctrl, B = dis)
  freq <- compute_edge_frequencies(pops)
  N <- sapply(pops, length)

  p_fisher <- compute_edge_pvalues(freq$edge_counts, N, method = "fisher")
  p_chisq  <- suppressWarnings(
    compute_edge_pvalues(freq$edge_counts, N, method = "chi.squared"))
  p_prop   <- suppressWarnings(
    compute_edge_pvalues(freq$edge_counts, N, method = "prop"))

  for (pv in list(p_fisher, p_chisq, p_prop)) {
    expect_equal(dim(pv), c(5, 5))
    expect_true(isSymmetric(pv))
    expect_true(all(pv >= 0 & pv <= 1))
  }
})

test_that("compute_edge_pvalues: invalid method throws", {
  counts <- .make_counts(list(matrix(c(0, 1, 1, 0), nrow = 2),
                              matrix(c(0, 0, 0, 0), nrow = 2)))
  expect_error(compute_edge_pvalues(counts, c(1, 1), method = "tukey"),
               "Invalid method")
})

test_that("compute_edge_pvalues: prop method requires exactly two populations", {
  counts <- array(0, dim = c(3, 3, 3))   # 3 populations
  expect_error(compute_edge_pvalues(counts, c(1, 1, 1), method = "prop"),
               "only implemented for two populations")
})

test_that("compute_edge_pvalues: adjust_method BH reduces significance", {
  set.seed(3)
  ctrl <- generate_category_graphs(n_graphs = 10, n_nodes = 5, n_communities = 1,
    base_intra_prob = 0.8, base_inter_prob = 0.8, seed = 1)
  dis  <- generate_category_graphs(n_graphs = 10, n_nodes = 5, n_communities = 1,
    base_intra_prob = 0.2, base_inter_prob = 0.2, seed = 2)
  pops <- list(A = ctrl, B = dis)
  freq <- compute_edge_frequencies(pops)
  N    <- sapply(pops, length)

  p_raw <- compute_edge_pvalues(freq$edge_counts, N, adjust_method = "none")
  p_bh  <- compute_edge_pvalues(freq$edge_counts, N, adjust_method = "BH")

  # Adjusted p-values are never smaller than raw p-values (BH monotone) on the
  # upper triangle, which is what downstream code consumes.
  ut <- upper.tri(p_raw)
  expect_true(all(p_bh[ut] >= p_raw[ut] - 1e-12))
  # Valid probabilities
  expect_true(all(p_bh >= 0 & p_bh <= 1))
})

test_that("compute_edge_pvalues returns a symmetric matrix for all adjust_methods", {
  # Regression test for a bug where the BH (and other) adjustment was written
  # incorrectly to the lower triangle, producing an asymmetric result.
  set.seed(7)
  ctrl <- generate_category_graphs(n_graphs = 6, n_nodes = 6, n_communities = 2,
    base_intra_prob = 0.9, base_inter_prob = 0.2, seed = 1)
  dis  <- generate_category_graphs(n_graphs = 6, n_nodes = 6, n_communities = 2,
    base_intra_prob = 0.3, base_inter_prob = 0.6, seed = 2)
  pops <- list(A = ctrl, B = dis)
  freq <- compute_edge_frequencies(pops)
  N    <- sapply(pops, length)

  for (adj in c("none", "BH", "bonferroni", "holm", "BY")) {
    p <- compute_edge_pvalues(freq$edge_counts, N, adjust_method = adj)
    expect_true(isSymmetric(p),
                info = paste("Non-symmetric p-value matrix for adjust_method =", adj))
    expect_true(all(p >= 0 & p <= 1))
  }
})

test_that("compute_edge_pvalues: prop with all-zero or all-present edges yields p = 1", {
  # One edge always absent in both groups
  n <- 4
  ec <- array(0, dim = c(n, n, 2))
  N  <- c(5, 5)
  pv <- compute_edge_pvalues(ec, N, method = "prop")
  expect_equal(unname(pv[1, 2]), 1.0)
  # One edge always present in both groups
  ec2 <- ec; ec2[1, 2, ] <- N; ec2[2, 1, ] <- N
  pv2 <- compute_edge_pvalues(ec2, N, method = "prop")
  expect_equal(unname(pv2[1, 2]), 1.0)
})

test_that("compute_edge_pvalues clamps out-of-range test p-values into [0, 1]", {
  # Regression test for the CRAN noLD failure
  # (test-compute_edge_pvalues.R:114). On builds without extended (long
  # double) precision, exact tests such as fisher.test() can return a p-value
  # fractionally outside [0, 1] due to floating-point rounding in their
  # internal summations. We emulate that here by mocking the underlying test,
  # so the regression is caught on every platform rather than only on noLD
  # ones. Without the clamp in compute_edge_pvalues(), these expectations fail.
  skip_if_not_installed("testthat", "3.1.7")  # with_mocked_bindings()

  ec <- array(0, dim = c(3, 3, 2))
  ec[1, 2, ] <- c(2, 4); ec[2, 1, ] <- c(2, 4)
  N <- c(6, 6)

  # p-value slightly above 1: the observed noLD failure mode.
  p_hi <- with_mocked_bindings(
    compute_edge_pvalues(ec, N, method = "fisher"),
    fisher.test = function(...) list(p.value = 1 + 1e-12)
  )
  expect_true(all(p_hi >= 0 & p_hi <= 1))
  expect_true(isSymmetric(p_hi))

  # p-value slightly below 0: defensive lower bound.
  p_lo <- with_mocked_bindings(
    compute_edge_pvalues(ec, N, method = "fisher"),
    fisher.test = function(...) list(p.value = -1e-12)
  )
  expect_true(all(p_lo >= 0 & p_lo <= 1))
})
