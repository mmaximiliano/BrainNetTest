## Tests of the combined stopping criterion (T together with the smallest
## marginal p-value, in one permutation test).

## Background probability 0.2 everywhere, except one edge whose probability in
## the second group is p1: the situation in which T alone is weakest.
single_edge_pops <- function(n_nodes = 30, n_graphs = 60, p1 = 0.9, seed = 1) {
  set.seed(seed)
  make <- function(p_edge) lapply(seq_len(n_graphs), function(g) {
    P <- matrix(0.2, n_nodes, n_nodes); P[1, 2] <- P[2, 1] <- p_edge
    A <- matrix(0L, n_nodes, n_nodes); ut <- upper.tri(A)
    A[ut] <- stats::rbinom(sum(ut), 1, P[ut]); A + t(A) })
  list(A = make(0.2), B = make(p1))
}

test_that("the combined criterion detects a single strongly altered edge", {
  pops <- single_edge_pops()
  cb <- identify_critical_links(pops, n_permutations = 300, seed = 2)
  expect_identical(cb$settings$criterion, "combined")
  expect_false(is.null(cb$critical_edges))
  expect_identical(as.integer(unlist(cb$critical_edges[1, c("node1", "node2")])),
                   c(1L, 2L))
  # with T alone the single edge is diluted among the 435 edges
  gt <- global_test(pops, n_permutations = 300, seed = 2, criterion = "original")
  expect_gt(gt$p_value, cb$p_value)
})

test_that("combined p-values lie on the grid of B + 1 values", {
  pops <- single_edge_pops(p1 = 0.5, seed = 3)
  B <- 199
  gt <- global_test(pops, n_permutations = B, seed = 4)
  for (p in c(gt$p_value, gt$p_value_global, gt$p_value_local)) {
    expect_equal(p * (B + 1), round(p * (B + 1)))
    expect_gte(p, 1 / (B + 1))
    expect_lte(p, 1)
  }
  expect_equal(gt$local_statistic,
               min(compute_edge_pvalues(compute_edge_frequencies(pops)$edge_counts,
                                        lengths(pops))[upper.tri(diag(30))]),
               tolerance = 1e-12)
})

test_that("the Fisher table reproduces fisher.test()", {
  H <- BrainNetTest:::.fisher_pvalue_table(7, 9)
  for (ct in c(0, 3, 8, 16)) for (c1 in max(0, ct - 9):min(7, ct)) {
    tab <- matrix(c(c1, ct - c1, 7 - c1, 9 - (ct - c1)), 2)
    expect_equal(H[ct + 1, c1 + 1], fisher.test(tab)$p.value, tolerance = 1e-12)
  }
})

test_that("the combined test keeps its level under the null hypothesis", {
  skip_on_cran()
  rej <- vapply(1:60, function(s) {
    same <- generate_category_graphs(n_graphs = 30, n_nodes = 8,
      n_communities = 2, base_intra_prob = 0.5, base_inter_prob = 0.3, seed = s)
    global_test(list(A = same[1:15], B = same[16:30]),
                n_permutations = 99, seed = s)$p_value <= 0.05
  }, logical(1))
  expect_lte(mean(rej), 0.15)   # 5% nominal; loose bound for 60 data sets
})

test_that("other marginal tests and more than two groups are supported", {
  pops <- single_edge_pops(n_nodes = 8, n_graphs = 20, seed = 5)
  three <- list(A = pops$A[1:10], B = pops$A[11:20], C = pops$B[1:10])
  for (m in c("chi.squared", "prop")) {
    res <- suppressWarnings(identify_critical_links(pops, method = m,
                                                    n_permutations = 99, seed = 1))
    expect_s3_class(res, "critical_links")
  }
  res3 <- suppressWarnings(identify_critical_links(three, n_permutations = 99, seed = 1))
  expect_s3_class(res3, "critical_links")
  expect_true(res3$p_value >= 1 / 100 && res3$p_value <= 1)
})
