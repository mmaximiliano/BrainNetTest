# tests/testthat/test-generate_category_graphs.R

test_that("generate_category_graphs works correctly", {
  set.seed(123)
  n_graphs <- 5
  n_nodes <- 50
  n_communities <- 3
  graphs <- generate_category_graphs(
    n_graphs = n_graphs, n_nodes = n_nodes, n_communities = n_communities,
    base_intra_prob = 0.8, base_inter_prob = 0.2,
    intra_prob_variation = 0.05, inter_prob_variation = 0.05
  )

  expect_equal(length(graphs), n_graphs)

  for (G in graphs) {
    expect_equal(dim(G), c(n_nodes, n_nodes))
    expect_true(isSymmetric(G))
    expect_equal(diag(G), rep(0, n_nodes))
    expect_true(all(G %in% c(0, 1)))
  }
})

test_that("generate_category_graphs respects deterministic probability bounds", {
  graphs <- generate_category_graphs(
    n_graphs = 3L,
    n_nodes = 6L,
    n_communities = 2L,
    base_intra_prob = 1,
    base_inter_prob = 0,
    intra_prob_variation = 0,
    inter_prob_variation = 0,
    seed = 1L
  )
  communities <- rep(1:2, each = 3L)

  for (graph in graphs) {
    expect_true(all(graph[outer(communities, communities, `==`) &
      row(graph) != col(graph)] == 1))
    expect_true(all(graph[outer(communities, communities, `!=`)] == 0))
  }
})

test_that("generate_category_graphs is reproducible with seed", {
  set.seed(7)
  before <- .Random.seed
  g1 <- generate_category_graphs(
    n_graphs = 5, n_nodes = 10, n_communities = 2,
    base_intra_prob = 0.7, base_inter_prob = 0.2, seed = 42
  )
  expect_identical(.Random.seed, before)
  g2 <- generate_category_graphs(
    n_graphs = 5, n_nodes = 10, n_communities = 2,
    base_intra_prob = 0.7, base_inter_prob = 0.2, seed = 42
  )
  expect_identical(g1, g2)
})

test_that("generate_category_graphs validates inputs", {
  expect_error(
    generate_category_graphs(n_graphs = 0, n_nodes = 10, n_communities = 2),
    "n_graphs"
  )
  expect_error(
    generate_category_graphs(n_graphs = 0.7, n_nodes = 10, n_communities = 2),
    "whole number"
  )
  expect_error(
    generate_category_graphs(n_graphs = 3, n_nodes = 0.7, n_communities = 2),
    "whole number"
  )
  expect_error(
    generate_category_graphs(n_graphs = 3, n_nodes = 3, n_communities = 4),
    "cannot exceed"
  )
  expect_error(
    generate_category_graphs(
      n_graphs = 3, n_nodes = 10, n_communities = 2,
      base_intra_prob = 1.2
    ),
    "base_intra_prob"
  )
  expect_error(
    generate_category_graphs(
      n_graphs = 3, n_nodes = 10, n_communities = 2,
      base_inter_prob = -0.1
    ),
    "base_inter_prob"
  )
  expect_error(
    generate_category_graphs(
      n_graphs = 3, n_nodes = 10, n_communities = 2,
      intra_prob_variation = 2
    ),
    "intra_prob_variation"
  )
  expect_error(
    generate_category_graphs(
      n_graphs = 3, n_nodes = 10, n_communities = 2,
      inter_prob_variation = -0.5
    ),
    "inter_prob_variation"
  )
  expect_error(
    generate_category_graphs(
      n_graphs = 3, n_nodes = 10, n_communities = 3,
      community_sizes = c(5, 5)
    ),
    "community_sizes"
  )
  expect_error(
    generate_category_graphs(
      n_graphs = 3, n_nodes = 10, n_communities = 2,
      community_sizes = c(3, 3)
    ),
    "sum"
  )
})

test_that("generate_category_graphs supports vector base_intra_prob", {
  set.seed(1)
  g <- generate_category_graphs(
    n_graphs = 4, n_nodes = 20, n_communities = 2,
    base_intra_prob = c(0.9, 0.1),
    base_inter_prob = 0.2,
    intra_prob_variation = 0.02,
    inter_prob_variation = 0.02,
    seed = 1
  )
  expect_equal(length(g), 4)
  for (G in g) {
    expect_equal(dim(G), c(20, 20))
    expect_true(isSymmetric(G))
    expect_equal(diag(G), rep(0, 20))
  }
})
