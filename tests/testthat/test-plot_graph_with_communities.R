# tests/testthat/test-plot_graph_with_communities.R

test_that("plot_graph_with_communities works without errors", {
  skip_on_cran()
  G <- generate_community_graph(n_nodes = 20, n_communities = 2, intra_prob = 0.8, inter_prob = 0.2)
  expect_silent(plot_graph_with_communities(G))
})
