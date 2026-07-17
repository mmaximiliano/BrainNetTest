test_that("version 1.0 exposes only the frozen public API", {
  expected <- c(
    "ablation_path",
    "brainnet_data",
    "brainnet_test",
    "central_graph",
    "generate_category_graphs",
    "generate_community_graph",
    "generate_random_graph",
    "graph_distance",
    "selected_edges",
    "selected_nodes"
  )

  expect_setequal(getNamespaceExports("BrainNetTest"), expected)

  removed <- c(
    "compute_central_graph",
    "compute_distance",
    "compute_edge_frequencies",
    "compute_edge_pvalues",
    "compute_test_statistic",
    "get_critical_nodes",
    "identify_critical_links",
    "plot_critical_edges",
    "rank_edges"
  )
  expect_false(any(removed %in% getNamespaceExports("BrainNetTest")))
})
