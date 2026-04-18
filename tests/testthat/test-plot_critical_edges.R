# tests/testthat/test-plot_critical_edges.R

make_populations <- function(seed = 1) {
  set.seed(seed)
  control <- generate_category_graphs(
    n_graphs = 10, n_nodes = 12, n_communities = 3,
    community_sizes = c(4, 4, 4),
    base_intra_prob = rep(0.8, 3), base_inter_prob = 0.1,
    intra_prob_variation = 0.02, inter_prob_variation = 0.02, seed = seed)
  patient <- generate_category_graphs(
    n_graphs = 10, n_nodes = 12, n_communities = 3,
    community_sizes = c(4, 4, 4),
    base_intra_prob = c(0.4, 0.8, 0.8), base_inter_prob = 0.1,
    intra_prob_variation = 0.02, inter_prob_variation = 0.02,
    seed = seed + 1)
  list(Control = control, Patient = patient)
}

test_that("plot_critical_edges runs without errors on a typical workflow", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  populations <- make_populations()
  result <- identify_critical_links(populations, alpha = 0.05,
    method = "fisher", n_permutations = 100, seed = 1)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(
    plot_critical_edges(populations, result,
                        communities = rep(1:3, each = 4)))
})

test_that("plot_critical_edges accepts pre-computed central graphs and a layout matrix", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  populations <- make_populations(2)
  result <- identify_critical_links(populations, alpha = 0.05,
    method = "fisher", n_permutations = 100, seed = 2)
  centrals <- lapply(populations, compute_central_graph)
  lay <- igraph::layout_in_circle(
    igraph::graph_from_adjacency_matrix(centrals[[1]] > 0.5,
                                        mode = "undirected"))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(
    plot_critical_edges(populations, result,
                        central_graphs = centrals,
                        layout = lay,
                        reference = "Patient",
                        vertex_label = FALSE))
})

test_that("plot_critical_edges validates input", {
  expect_error(plot_critical_edges(list(), list(critical_edges = NULL)),
               "non-empty list")
  expect_error(
    plot_critical_edges(list(A = list()), list(not_critical = NULL)),
    "identify_critical_links")
  populations <- make_populations(3)
  fake_result <- list(critical_edges = NULL)
  expect_error(
    plot_critical_edges(populations, fake_result,
                        communities = 1:5),
    "communities")
  expect_error(
    plot_critical_edges(populations, fake_result, reference = 99L),
    "out of range")
  expect_error(
    plot_critical_edges(populations, fake_result, reference = "Nope"),
    "not found")
})
