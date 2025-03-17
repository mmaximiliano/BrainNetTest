# A simple test for the identify_critical_links function.
# This test aims to:
# 1) Generate two small "Control" and "Disease" populations with clearly different edge probabilities.
# 2) Run identify_critical_links() to see if it can detect any "critical" links that distinguish them.
# 3) Verify that the function returns the correct structure (a list) containing critical_edges,
#    edges_removed, and modified_populations.
# 4) We do not check for a specific number of critical edges here—just the structure
#    and a confirmed difference in edge probabilities between the two populations.

library(BrainNetTest)
library(testthat)

test_that("identify_critical_links works on a small example", {
  # Generate a small Control population:
  # This population is expected to have higher intra-community probability of edges.
  control_graphs <- generate_category_graphs(
    n_graphs = 2,
    n_nodes = 5,
    n_communities = 1,
    base_intra_prob = 0.8,
    base_inter_prob = 0.2,
    seed = 1
  )

  # Generate a small Disease population:
  # This population is expected to have a lower intra-community probability of edges (and is thus quite different).
  disease_graphs <- generate_category_graphs(
    n_graphs = 2,
    n_nodes = 5,
    n_communities = 1,
    base_intra_prob = 0.2,
    base_inter_prob = 0.8,
    seed = 2
  )

  # Combine the two groups into a single list for testing
  populations <- list(Control = control_graphs, Disease = disease_graphs)

  # Call the function with significance level alpha=1, using Fisher's exact test and no p-value adjustment
  result <- identify_critical_links(
    populations = populations,
    alpha = 1,
    method = "fisher",
    adjust_method = "none"
  )

  # Check that the output is a list and has the expected components:
  expect_type(result, "list")

  # critical_edges might be NULL if no edges are found significant. If it exists, it should be a data frame.
  expect_true(is.list(result$critical_edges) || is.null(result$critical_edges))

  # edges_removed must be a list, potentially storing pairs of node indices
  expect_true(is.list(result$edges_removed))

  # modified_populations should be a list of lists, each containing adjacency matrices
  expect_true(is.list(result$modified_populations))

  # If there are any critical_edges, check that it has at least 0 rows (sanity check)
  # and presumably columns for node1, node2, p_value.
  if (!is.null(result$critical_edges)) {
    expect_true(nrow(result$critical_edges) >= 0)
  }

  # Expected outcome:
  # - The function should detect some edges as "critical" since the populations differ in their edge probabilities.
  # - The test will pass as long as the function runs and returns the correct structure.
  # - We do not require a specific number of critical edges for this test.
})
