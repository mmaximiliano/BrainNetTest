test_that("brainnet_data has concise print and summary methods", {
  data <- tiny_brainnet_data()

  expect_output(print(data), "2 groups")
  expect_output(print(data), "3 nodes")

  summary <- summary(data)
  expect_s3_class(summary, "summary_brainnet_data")
  expect_identical(summary$group_sizes, data$group_sizes)
  expect_length(summary$mean_density, 2L)
})

test_that("brainnet_result methods expose complete and selected edge views", {
  data <- tiny_brainnet_data()
  result <- brainnet_test(
    data,
    exact = "never",
    n_permutations = 19L,
    seed = 1L
  )

  expect_output(print(result), "Global permutation test")
  expect_output(print(result), "Holm")
  expect_output(print(result$global), "T =")
  expect_s3_class(summary(result), "summary_brainnet_result")
  expect_identical(as.data.frame(result), result$edges)
  expect_true(all(selected_edges(result)$selected))
})

test_that("ablation is explicitly descriptive and deterministic", {
  data <- tiny_brainnet_data()
  without <- brainnet_test(
    data,
    exact = "never",
    n_permutations = 19L,
    seed = 11L
  )
  with_1 <- brainnet_test(
    data,
    exact = "never",
    n_permutations = 19L,
    seed = 11L,
    ablation = TRUE
  )
  with_2 <- brainnet_test(
    data,
    exact = "never",
    n_permutations = 19L,
    seed = 11L,
    ablation = TRUE
  )

  expect_error(ablation_path(without), "not computed")

  path <- ablation_path(with_1)
  expect_equal(nrow(path), 4L)
  expect_named(
    path,
    c("n_removed", "statistic", "descriptive_tail_fraction")
  )
  expect_false("p_value" %in% names(path))
  expect_equal(path$statistic[1L], unname(with_1$global$statistic))
  expect_identical(path, ablation_path(with_2))
})

test_that("plot method uses stored central graphs", {
  skip_if_not_installed("igraph")
  data <- tiny_brainnet_data()
  result <- brainnet_test(
    data,
    exact = "never",
    n_permutations = 19L,
    seed = 1L
  )

  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output)
  on.exit(
    {
      grDevices::dev.off()
      unlink(output)
    },
    add = TRUE
  )

  expect_invisible(plot(result))
})

test_that("brainnet_result methods reject corrupted result objects", {
  result <- brainnet_test(
    tiny_brainnet_data(),
    exact = "never",
    n_permutations = 19L,
    seed = 1L
  )

  corrupted_edges <- result
  corrupted_edges$edges <- corrupted_edges$edges[-1L, ]
  expect_error(print(corrupted_edges), "internally inconsistent")

  corrupted_statistic <- result
  names(corrupted_statistic$global$statistic) <- NULL
  expect_error(summary(corrupted_statistic), "internally inconsistent")
})
