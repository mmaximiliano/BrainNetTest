# tests/testthat/test-critical_links_methods.R

make_result <- function(seed = 42, n_permutations = 200) {
  control <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
  patient <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.5, base_inter_prob = 0.5, seed = 2)
  populations <- list(Control = control, Patient = patient)
  list(populations = populations,
       result = identify_critical_links(populations,
                                        n_permutations = n_permutations,
                                        seed = seed))
}

test_that("identify_critical_links returns a classed list", {
  skip_on_cran()
  res <- make_result()$result

  expect_s3_class(res, "critical_links")
  expect_true(is.list(res))
  # The components documented before 0.2.2 keep their names, positions and
  # contents, so existing code that extracts them is unaffected.
  expect_identical(names(res)[1:3],
                   c("critical_edges", "edges_removed",
                     "modified_populations"))
  expect_s3_class(res$critical_edges, "data.frame")
  expect_type(res$edges_removed, "list")
  expect_length(res$modified_populations, 2L)
})

test_that("the object records the global test and the settings", {
  skip_on_cran()
  res <- make_result()$result

  expect_true(is.numeric(res$p_value) && res$p_value >= 0 && res$p_value <= 1)
  expect_identical(res$n_edges, 45L)   # 10 nodes -> 45 candidate edges
  expect_identical(res$settings$method, "fisher")
  expect_identical(res$settings$n_permutations, 200L)
  expect_identical(res$settings$seed, 42)
  expect_true(is.call(res$call))
})

test_that("print returns its argument invisibly and reports the counts", {
  skip_on_cran()
  res <- make_result()$result

  out <- capture.output(printed <- withVisible(print(res)))
  expect_false(printed$visible)
  expect_identical(printed$value, res)

  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Control, Patient")
  expect_match(txt, "Candidate edges: 45")
  expect_match(txt, sprintf("Critical edges:  %d of 45", nrow(res$critical_edges)))
  expect_match(txt, "Fisher's exact test")
})

test_that("print reports an empty critical set without failing", {
  skip_on_cran()
  identical_graphs <- generate_category_graphs(n_graphs = 20, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.5, base_inter_prob = 0.5, seed = 1)
  pops <- list(A = identical_graphs[1:10], B = identical_graphs[11:20])

  suppressWarnings(
    res <- identify_critical_links(pops, n_permutations = 500, seed = 42))
  expect_null(res$critical_edges)

  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "none")
})

test_that("summary carries the node-level view and prints", {
  skip_on_cran()
  res <- make_result()$result
  s <- summary(res)

  expect_s3_class(s, "summary.critical_links")
  expect_identical(s$n_critical, nrow(res$critical_edges))
  expect_identical(s$critical_nodes, get_critical_nodes(res))

  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "Settings:")
  expect_match(out, "permutations")
  expect_match(out, "critical degree")
})

test_that("summary passes node labels through", {
  skip_on_cran()
  res <- make_result()$result
  labels <- paste0("R", 1:10)
  s <- summary(res, node_labels = labels)

  expect_true("label" %in% names(s$critical_nodes))
  expect_identical(s$critical_nodes$label, labels[s$critical_nodes$node])
})

test_that("plot method delegates to plot_critical_edges", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  fixture <- make_result()

  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(fixture$result, fixture$populations,
                       communities = rep(1:2, each = 5)))
})

test_that("plot requires the original populations", {
  skip_on_cran()
  res <- make_result()$result
  # modified_populations has the critical edges removed, so it cannot stand in
  # for the input; the method says so rather than drawing something wrong.
  expect_error(plot(res), "populations")
})

test_that("the naive implementation returns the same structure", {
  skip_on_cran()
  naive <- BrainNetTest:::identify_critical_links_naive

  ctrl <- generate_category_graphs(n_graphs = 8, n_nodes = 6,
    n_communities = 2, base_intra_prob = 0.9, base_inter_prob = 0.1, seed = 1)
  dis  <- generate_category_graphs(n_graphs = 8, n_nodes = 6,
    n_communities = 2, base_intra_prob = 0.1, base_inter_prob = 0.9, seed = 2)
  pops <- list(A = ctrl, B = dis)

  slow <- naive(pops, n_permutations = 200, seed = 321)
  fast <- identify_critical_links(pops, n_permutations = 200, seed = 321)

  expect_s3_class(slow, "critical_links")
  expect_identical(names(slow), names(fast))
  expect_identical(slow$n_edges, fast$n_edges)
  expect_identical(slow$settings, fast$settings)
  expect_no_error(capture.output(print(slow)))
})
