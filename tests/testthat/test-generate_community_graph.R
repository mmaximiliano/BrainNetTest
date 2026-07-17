# tests/testthat/test-generate_community_graph.R

test_that("generate_community_graph returns a valid adjacency matrix", {
  set.seed(123)
  G <- generate_community_graph(
    n_nodes = 10, n_communities = 2,
    intra_prob = 0.9, inter_prob = 0.1
  )

  expect_equal(dim(G), c(10, 10))
  expect_true(isSymmetric(G))
  expect_equal(diag(G), rep(0, 10))
  expect_true(all(G %in% c(0, 1)))
})

test_that("generate_community_graph respects probability boundaries", {
  graph <- generate_community_graph(
    n_nodes = 6L,
    n_communities = 2L,
    intra_prob = 1,
    inter_prob = 0,
    seed = 1L
  )
  communities <- rep(1:2, each = 3L)

  expect_true(all(graph[outer(communities, communities, `==`) &
    row(graph) != col(graph)] == 1))
  expect_true(all(graph[outer(communities, communities, `!=`)] == 0))
})

test_that("generate_community_graph is reproducible with seed", {
  set.seed(7)
  before <- .Random.seed
  G1 <- generate_community_graph(20, 2, intra_prob = 0.7, inter_prob = 0.2, seed = 99)
  expect_identical(.Random.seed, before)
  G2 <- generate_community_graph(20, 2, intra_prob = 0.7, inter_prob = 0.2, seed = 99)
  expect_identical(G1, G2)
})

test_that("generate_community_graph respects custom community_sizes", {
  G <- generate_community_graph(
    n_nodes = 12, n_communities = 3,
    community_sizes = c(2, 4, 6),
    intra_prob = 0.8, inter_prob = 0.2, seed = 1
  )
  expect_equal(dim(G), c(12, 12))
  expect_true(isSymmetric(G))
})

test_that("generate_community_graph accepts a vector intra_prob", {
  set.seed(1)
  # Community 1 is dense, community 2 is sparse
  G <- generate_community_graph(
    n_nodes = 40, n_communities = 2,
    intra_prob = c(0.95, 0.05),
    inter_prob = 0.1
  )
  comm <- rep(1:2, each = 20)
  dens1 <- mean(G[comm == 1, comm == 1][upper.tri(G[comm == 1, comm == 1])])
  dens2 <- mean(G[comm == 2, comm == 2][upper.tri(G[comm == 2, comm == 2])])
  expect_true(dens1 > dens2)
})

test_that("generate_community_graph validates inputs", {
  expect_error(
    generate_community_graph(n_nodes = 0, n_communities = 2),
    "n_nodes"
  )
  expect_error(
    generate_community_graph(n_nodes = 10.5, n_communities = 2),
    "whole number"
  )
  expect_error(
    generate_community_graph(n_nodes = 10, n_communities = 0),
    "n_communities"
  )
  expect_error(
    generate_community_graph(n_nodes = 3, n_communities = 4),
    "cannot exceed"
  )
  expect_error(
    generate_community_graph(
      n_nodes = 10, n_communities = 2,
      community_sizes = c(3, 3)
    ),
    "sum"
  )
  expect_error(
    generate_community_graph(
      n_nodes = 10, n_communities = 3,
      community_sizes = c(5, 5)
    ),
    "community_sizes"
  )
  expect_error(
    generate_community_graph(
      n_nodes = 10, n_communities = 2,
      intra_prob = 1.5
    ),
    "intra_prob"
  )
  expect_error(
    generate_community_graph(
      n_nodes = 10, n_communities = 2,
      inter_prob = -0.1
    ),
    "inter_prob"
  )
  # Wrong-length intra_prob vector
  expect_error(
    generate_community_graph(
      n_nodes = 10, n_communities = 2,
      intra_prob = c(0.5, 0.5, 0.5)
    ),
    "intra_prob"
  )
})

test_that("generate_community_graph works with n_communities = 1", {
  # Regression: n_communities = 1 previously triggered `1:(n_communities-1)`
  # which is `1:0`, corrupting the inter-community loop.
  set.seed(1)
  G <- generate_community_graph(
    n_nodes = 10, n_communities = 1,
    intra_prob = 0.7, inter_prob = 0.0
  )
  expect_equal(dim(G), c(10, 10))
  expect_true(isSymmetric(G))
  expect_equal(diag(G), rep(0, 10))
  expect_true(all(G %in% c(0, 1)))
})
