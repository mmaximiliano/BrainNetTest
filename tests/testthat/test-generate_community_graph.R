# tests/testthat/test-generate_community_graph.R

test_that("generate_community_graph works correctly", {
  set.seed(123)
  G <- generate_community_graph(n_nodes = 10, n_communities = 2, intra_prob = 0.9, inter_prob = 0.1)

  expect_equal(dim(G), c(10, 10))
  expect_true(isSymmetric(G))
  expect_equal(diag(G), rep(0, 10))

  # Check that intra-community edges are more frequent than inter-community edges
  community_assignments <- rep(1:2, each = 5)
  intra_edges <- sum(G[community_assignments == 1, community_assignments == 1]) +
                 sum(G[community_assignments == 2, community_assignments == 2])
  inter_edges <- sum(G[community_assignments == 1, community_assignments == 2])

  expect_true(intra_edges > inter_edges)
})
