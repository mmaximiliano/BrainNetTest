# tests/testthat/test-plot_graph_with_communities.R

test_that("plot_graph_with_communities works without errors", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  G <- generate_community_graph(n_nodes = 20, n_communities = 2,
                                intra_prob = 0.8, inter_prob = 0.2)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot_graph_with_communities(G))
})

test_that("plot_graph_with_communities accepts a user-supplied community vector", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  G <- generate_community_graph(n_nodes = 12, n_communities = 3,
                                intra_prob = 0.9, inter_prob = 0.1, seed = 1)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(
    plot_graph_with_communities(G, communities = rep(1:3, each = 4),
                                main = "custom"))
})

test_that("plot_graph_with_communities validates input", {
  expect_error(plot_graph_with_communities("not a matrix"),
               "square adjacency matrix")
  expect_error(plot_graph_with_communities(matrix(0, 3, 4)),
               "square adjacency matrix")
})

