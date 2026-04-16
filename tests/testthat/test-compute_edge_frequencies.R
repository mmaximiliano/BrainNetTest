# tests/testthat/test-compute_edge_frequencies.R

test_that("compute_edge_frequencies returns correct shape and values", {
  # Build two deterministic populations (2 pops x 2 graphs x 3 nodes)
  G_full  <- matrix(c(0, 1, 1,
                      1, 0, 1,
                      1, 1, 0), nrow = 3, byrow = TRUE)
  G_empty <- matrix(0, 3, 3)

  pops <- list(A = list(G_full, G_full),
               B = list(G_empty, G_full))

  freq <- compute_edge_frequencies(pops)

  # Array dims: n_nodes x n_nodes x num_populations
  expect_equal(dim(freq$edge_counts), c(3, 3, 2))
  expect_equal(dim(freq$edge_proportions), c(3, 3, 2))

  # Pop A: both graphs full -> count = 2 on off-diagonal, proportion = 1
  expect_equal(freq$edge_counts[, , 1], G_full * 2)
  expect_equal(freq$edge_proportions[, , 1], G_full)

  # Pop B: one empty + one full -> count = 1 on off-diagonal, proportion = 0.5
  expect_equal(freq$edge_counts[, , 2], G_full * 1)
  expect_equal(freq$edge_proportions[, , 2], G_full * 0.5)
})

test_that("compute_edge_frequencies handles unequal group sizes", {
  G <- matrix(c(0, 1, 1, 0), nrow = 2)
  pops <- list(A = list(G, G, G),    # n_A = 3
               B = list(G))          # n_B = 1
  freq <- compute_edge_frequencies(pops)
  expect_equal(freq$edge_counts[, , 1], G * 3)
  expect_equal(freq$edge_counts[, , 2], G * 1)
  expect_equal(freq$edge_proportions[, , 1], G)
  expect_equal(freq$edge_proportions[, , 2], G)
})
