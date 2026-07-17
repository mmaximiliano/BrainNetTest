test_that("edge inference tests the complete upper-triangle family", {
  empty <- graph_from_edges(3L)
  edge_12 <- graph_from_edges(3L, c(1L, 2L))
  data <- brainnet_data(list(
    GroupA = constant_group(empty, 8L),
    GroupB = constant_group(edge_12, 8L)
  ))

  result <- brainnet_test(
    data,
    exact = "never",
    n_permutations = 99L,
    seed = 1L
  )
  edges <- as.data.frame(result)

  expect_equal(nrow(edges), 3L)
  expect_identical(edges$node1, c(1L, 1L, 2L))
  expect_identical(edges$node2, c(2L, 3L, 3L))
  expect_identical(result$edge_inference$family_size, 3L)
  expect_identical(result$edge_inference$adjustment, "holm")

  changed <- edges[edges$node1 == 1L & edges$node2 == 2L, ]
  expect_equal(changed$risk_difference, 1)
  expect_true(changed$selected)
  expect_true(changed$p_adjusted >= changed$p_value)
})

test_that("selected edge and node extractors use adjusted decisions", {
  empty <- graph_from_edges(3L)
  edge_12 <- graph_from_edges(3L, c(1L, 2L))
  data <- brainnet_data(
    list(
      GroupA = constant_group(empty, 8L),
      GroupB = constant_group(edge_12, 8L)
    ),
    node_labels = c("Frontal", "Central", "Occipital")
  )
  result <- brainnet_test(
    data,
    exact = "never",
    n_permutations = 99L,
    seed = 1L
  )

  expect_equal(nrow(selected_edges(result)), 1L)
  nodes <- selected_nodes(result)
  expect_identical(nodes$node, c(1L, 2L))
  expect_identical(nodes$label, c("Frontal", "Central"))
  expect_true(all(nodes$selected_degree == 1L))
})

test_that("multi-group effects are omnibus ranges with named extrema", {
  empty <- graph_from_edges(3L)
  edge_12 <- graph_from_edges(3L, c(1L, 2L))
  mixed <- c(constant_group(empty, 2L), constant_group(edge_12, 2L))
  data <- brainnet_data(list(
    Low = constant_group(empty, 4L),
    Middle = mixed,
    High = constant_group(edge_12, 4L)
  ))

  result <- brainnet_test(
    data,
    exact = "never",
    n_permutations = 39L,
    seed = 3L
  )
  changed <- as.data.frame(result)
  changed <- changed[changed$node1 == 1L & changed$node2 == 2L, ]

  expect_equal(changed$effect_range, 1)
  expect_identical(changed$min_group, "Low")
  expect_identical(changed$max_group, "High")
  expect_true(is.na(changed$risk_difference))
})

test_that("colliding syntactic group names retain distinct proportion columns", {
  empty <- graph_from_edges(3L)
  edge_12 <- graph_from_edges(3L, c(1L, 2L))
  data <- brainnet_data(list(
    "A B" = constant_group(empty, 4L),
    "A.B" = constant_group(edge_12, 4L)
  ))

  result <- brainnet_test(
    data,
    exact = "never",
    n_permutations = 19L,
    seed = 1L
  )
  columns <- unname(result$edge_inference$proportion_columns)

  expect_length(unique(columns), 2L)
  expect_true(all(columns %in% names(result$edges)))
  expect_identical(
    names(result$edge_inference$proportion_columns),
    c("A B", "A.B")
  )
})

test_that("BH adjustment warns about its dependence assumption", {
  data <- tiny_brainnet_data()

  expect_warning(
    result <- brainnet_test(
      data,
      adjust = "BH",
      exact = "never",
      n_permutations = 19L,
      seed = 1L
    ),
    "positive dependence"
  )
  expect_identical(result$edge_inference$adjustment, "BH")
})

test_that("BY and chi-squared paths expose their inferential metadata", {
  data <- tiny_brainnet_data()
  by_result <- brainnet_test(
    data,
    adjust = "BY",
    exact = "never",
    n_permutations = 19L,
    seed = 1L
  )
  expect_identical(by_result$edge_inference$adjustment, "BY")

  expect_warning(
    chi_result <- brainnet_test(
      data,
      edge_method = "chisq",
      exact = "never",
      n_permutations = 19L,
      seed = 1L
    ),
    "expected counts"
  )
  expect_identical(chi_result$edge_inference$method, "chisq")
})

test_that("edge p-values are clamped to the probability range", {
  counts <- matrix(c(1, 2, 2, 1, 0, 3), nrow = 3L, byrow = TRUE)
  indices <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))

  result <- testthat::with_mocked_bindings(
    BrainNetTest:::.brainnet_edge_inference(
      counts = counts,
      group_sizes = c(3L, 3L),
      edge_indices = indices,
      group_names = c("A", "B"),
      alpha = 0.05,
      method = "fisher",
      adjust = "holm"
    ),
    .edge_test_result = function(...) {
      list(p_value = 1 + .Machine$double.eps, small_expected = FALSE)
    },
    .package = "BrainNetTest"
  )

  expect_true(all(result$edges$p_value <= 1))
  expect_true(all(result$edges$p_value >= 0))
})
