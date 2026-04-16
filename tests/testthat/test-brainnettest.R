# tests/testthat/test-brainnettest.R

library(testthat)
library(BrainNetTest)

# ---------------------------------------------------------------------------
# compute_central_graph
# ---------------------------------------------------------------------------
test_that("compute_central_graph averages adjacency matrices", {
  G1 <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
  G2 <- matrix(c(0, 0, 1, 0, 0, 1, 1, 0, 0), nrow = 3)
  central <- compute_central_graph(list(G1, G2))
  expect_equal(central, (G1 + G2) / 2)
  # Entries are the fraction of graphs containing each edge
  expect_true(all(central >= 0 & central <= 1))
})

test_that("compute_central_graph validates inputs", {
  expect_error(compute_central_graph(list()), "empty")
  expect_error(compute_central_graph(list(1, 2)), "must be matrices")
  M1 <- matrix(0, 3, 3); M2 <- matrix(0, 4, 4)
  expect_error(compute_central_graph(list(M1, M2)),
               "same dimensions")
})

test_that("compute_central_graph of a singleton returns that graph", {
  G <- matrix(c(0, 1, 1, 0), 2, 2)
  expect_equal(compute_central_graph(list(G)), G)
})

# ---------------------------------------------------------------------------
# compute_distance
# ---------------------------------------------------------------------------
test_that("compute_distance returns the L1 norm", {
  G <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
  M <- matrix(c(0, 0, 1, 0, 0, 1, 1, 0, 0), nrow = 3)
  expect_equal(compute_distance(G, M), sum(abs(G - M)))
})

test_that("compute_distance(G, G) == 0", {
  set.seed(1)
  G <- generate_random_graph(n_nodes = 6, edge_prob = 0.3)
  expect_equal(compute_distance(G, G), 0)
})

test_that("compute_distance is symmetric", {
  set.seed(1)
  G1 <- generate_random_graph(6, 0.3)
  G2 <- generate_random_graph(6, 0.4)
  expect_equal(compute_distance(G1, G2), compute_distance(G2, G1))
})

test_that("compute_distance validates inputs", {
  expect_error(compute_distance("a", matrix(0, 2, 2)), "matrices")
  expect_error(compute_distance(matrix(0, 2, 2), list()), "matrices")
  expect_error(compute_distance(matrix(0, 2, 2), matrix(0, 3, 3)),
               "same dimensions")
})

# ---------------------------------------------------------------------------
# compute_test_statistic
# ---------------------------------------------------------------------------
test_that("compute_test_statistic returns a finite numeric", {
  set.seed(1)
  ctrl <- replicate(5, generate_random_graph(6, 0.3), simplify = FALSE)
  dis  <- replicate(5, generate_random_graph(6, 0.6), simplify = FALSE)
  T_value <- compute_test_statistic(list(A = ctrl, B = dis), a = 1)
  expect_type(T_value, "double")
  expect_length(T_value, 1L)
  expect_true(is.finite(T_value))
})

test_that("compute_test_statistic is invariant to population order", {
  set.seed(2)
  A <- replicate(6, generate_random_graph(6, 0.3), simplify = FALSE)
  B <- replicate(6, generate_random_graph(6, 0.5), simplify = FALSE)
  T1 <- compute_test_statistic(list(A = A, B = B), a = 1)
  T2 <- compute_test_statistic(list(B = B, A = A), a = 1)
  expect_equal(unname(T1), unname(T2))
})

test_that("compute_test_statistic scales as 1 / a", {
  set.seed(3)
  A <- replicate(5, generate_random_graph(5, 0.3), simplify = FALSE)
  B <- replicate(5, generate_random_graph(5, 0.5), simplify = FALSE)
  T_a1 <- compute_test_statistic(list(A = A, B = B), a = 1)
  T_a2 <- compute_test_statistic(list(A = A, B = B), a = 2)
  expect_equal(T_a2, T_a1 / 2)
})

test_that("compute_test_statistic validates inputs", {
  expect_error(compute_test_statistic(list()), "non-empty list")
  expect_error(compute_test_statistic("bad"),  "non-empty list")
  # Each population must have at least 2 graphs
  G <- matrix(0, 3, 3)
  expect_error(
    compute_test_statistic(list(A = list(G), B = list(G, G))),
    "at least two graphs")
})

