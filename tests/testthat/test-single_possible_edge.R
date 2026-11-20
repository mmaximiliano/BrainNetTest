test_that("combined identification agrees with the global test for one edge", {
  empty <- matrix(0L, 2, 2)
  full <- matrix(c(0L, 1L, 1L, 0L), 2, 2)
  pops <- list(A = rep(list(empty), 10), B = rep(list(full), 10))
  for (B in c(1L, 2L, 3L, 98L, 99L, 100L)) {
    gt <- global_test(pops, n_permutations = B, seed = 1)
    id <- suppressWarnings(identify_critical_links(pops,
      n_permutations = B, seed = 1))
    expect_identical(id$p_value, gt$p_value, info = paste("B =", B))
    expect_equal(id$n_edges, 1L)
    if (gt$p_value <= 0.05) {
      expect_equal(nrow(id$critical_edges), 1L)
      expect_equal(unname(unlist(id$critical_edges[1, c("node1", "node2")])),
                   c(1, 2))
      expect_true(all(vapply(unlist(id$modified_populations, recursive = FALSE),
                             function(g) all(g == 0), logical(1))))
    } else {
      expect_null(id$critical_edges)
    }
  }
})

test_that("one-edge combined inference preserves the null case", {
  empty <- matrix(0L, 2, 2)
  pops <- list(A = rep(list(empty), 4), B = rep(list(empty), 4))
  for (B in c(1L, 2L, 99L, 100L)) {
    gt <- global_test(pops, n_permutations = B, seed = 2)
    id <- suppressWarnings(identify_critical_links(pops,
      n_permutations = B, seed = 2))
    expect_identical(id$p_value, gt$p_value)
    expect_equal(id$p_value, 1)
    expect_null(id$critical_edges)
  }
})
