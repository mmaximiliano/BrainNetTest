# tests/testthat/test-identify_critical_links.R

test_that("identify_critical_links works correctly", {
  skip_on_cran()
  
  # Generate synthetic populations
  set.seed(123)
  control_graphs <- generate_category_graphs(n_graphs = 10, n_nodes = 20, n_communities = 2,
                                             base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
  disease_graphs <- generate_category_graphs(n_graphs = 10, n_nodes = 20, n_communities = 2,
                                             base_intra_prob = 0.6, base_inter_prob = 0.4, seed = 2)
  populations <- list(Control = control_graphs, Disease = disease_graphs)
  
  # Run the function
  result <- identify_critical_links(populations, alpha = 0.05, method = "fisher", adjust_method = "none")
  
  # Check that critical edges are identified
  expect_true(!is.null(result$critical_edges))
  expect_true(length(result$edges_removed) > 0)
  
  # Check that modified populations are returned
  expect_equal(length(result$modified_populations), 2)
  
  # Further checks can be added based on expected outcomes
})
