# tests/testthat/test-generate_random_graph.R

test_that("generate_random_graph returns a valid adjacency matrix", {
  set.seed(1)
  G <- generate_random_graph(n_nodes = 20, edge_prob = 0.3)

  expect_true(is.matrix(G))
  expect_equal(dim(G), c(20, 20))
  expect_true(isSymmetric(G))
  expect_equal(diag(G), rep(0, 20))
  expect_true(all(G %in% c(0, 1)))
})

test_that("generate_random_graph respects edge_prob boundaries", {
  # edge_prob = 0 -> empty graph
  G0 <- generate_random_graph(n_nodes = 10, edge_prob = 0)
  expect_equal(sum(G0), 0)

  # edge_prob = 1 -> complete graph (off-diagonal)
  G1 <- generate_random_graph(n_nodes = 10, edge_prob = 1)
  expect_equal(sum(G1), 10 * 9)   # all off-diagonal entries
  expect_equal(diag(G1), rep(0, 10))
})

test_that("generate_random_graph validates inputs", {
  expect_error(generate_random_graph(n_nodes = 0, edge_prob = 0.1),
               "positive integer")
  expect_error(generate_random_graph(n_nodes = -5, edge_prob = 0.1),
               "positive integer")
  expect_error(generate_random_graph(n_nodes = c(5, 10), edge_prob = 0.1),
               "positive integer")
  expect_error(generate_random_graph(n_nodes = 10, edge_prob = -0.1),
               "between 0 and 1")
  expect_error(generate_random_graph(n_nodes = 10, edge_prob = 1.1),
               "between 0 and 1")
  expect_error(generate_random_graph(n_nodes = 10, edge_prob = "a"),
               "numeric")
})
