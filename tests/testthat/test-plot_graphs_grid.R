# tests/testthat/test-plot_graphs_grid.R

test_that("plot_graphs_grid works without errors", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("ggplotify")
  graphs <- list(
    generate_community_graph(n_nodes = 20, n_communities = 2, intra_prob = 0.8, inter_prob = 0.2),
    generate_community_graph(n_nodes = 20, n_communities = 2, intra_prob = 0.7, inter_prob = 0.3)
  )
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(plot_graphs_grid(graphs))
})

test_that("plot_graphs_grid validates communities_list length", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("gridExtra")
  graphs <- list(
    generate_random_graph(5, 0.3),
    generate_random_graph(5, 0.3)
  )
  expect_error(
    plot_graphs_grid(graphs, communities_list = list(rep(1, 5))),
    "communities_list"
  )
})

