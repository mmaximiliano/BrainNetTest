# tests/testthat/test-plot_graphs_grid.R

test_that("plot_graphs_grid works without errors", {
  skip_on_cran()
  graphs <- list(
    generate_community_graph(n_nodes = 20, n_communities = 2, intra_prob = 0.8, inter_prob = 0.2),
    generate_community_graph(n_nodes = 20, n_communities = 2, intra_prob = 0.7, inter_prob = 0.3)
  )
  expect_silent(plot_graphs_grid(graphs))
})
