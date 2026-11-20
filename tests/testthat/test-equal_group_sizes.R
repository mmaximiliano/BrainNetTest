test_that("T-based inference rejects unequal groups in either order", {
  empty <- matrix(0, 3, 3)
  full <- matrix(1, 3, 3) - diag(3)
  groups <- list(A = rep(list(empty), 4), B = rep(list(full), 8))
  for (pops in list(groups, rev(groups), c(groups, list(C = groups$A)))) {
    expect_error(compute_test_statistic(pops), "requires equal group sizes")
    for (criterion in c("original", "combined")) {
      expect_error(global_test(pops, n_permutations = 9,
                               criterion = criterion), "requires equal group sizes")
      expect_error(identify_critical_links(pops, n_permutations = 9,
                                           criterion = criterion),
                   "requires equal group sizes")
    }
    expect_error(BrainNetTest:::identify_critical_links_naive(
      pops, n_permutations = 9), "requires equal group sizes")
  }
})

test_that("unequal groups cannot bypass the guard through an empty edge set", {
  pops <- list(A = rep(list(matrix(0, 1, 1)), 2),
               B = rep(list(matrix(0, 1, 1)), 3))
  expect_error(global_test(pops, 9), "requires equal group sizes")
  expect_error(identify_critical_links(pops, n_permutations = 9),
               "requires equal group sizes")
  expect_error(BrainNetTest:::identify_critical_links_naive(
    pops, n_permutations = 9), "requires equal group sizes")
})

test_that("balanced designs and descriptive summaries remain available", {
  g <- matrix(0, 3, 3)
  for (m in 2:4) {
    pops <- rep(list(rep(list(g), 3)), m)
    expect_equal(compute_test_statistic(pops), 0)
    expect_equal(global_test(pops, 9, criterion = "original")$p_value, 1)
  }
  pops <- list(A = rep(list(g), 2), B = rep(list(g), 3))
  expect_equal(compute_edge_frequencies(pops)$edge_counts, array(0, c(3, 3, 2)))
})

test_that("compute_test_statistic rejects single-group inputs and invalid a", {
  g <- matrix(0, 3, 3)
  expect_error(compute_test_statistic(list(A = list(g, g))),
               "at least 2 groups")
  pops <- list(A = list(g, g), B = list(g, g))
  for (bad_a in list(0, -1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(compute_test_statistic(pops, a = bad_a), "`a`")
  }
})
