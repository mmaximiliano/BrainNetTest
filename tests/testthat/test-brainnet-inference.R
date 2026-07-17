test_that("graph_distance counts each undirected edge once", {
  empty <- graph_from_edges(3L)
  one_edge <- graph_from_edges(3L, c(1L, 2L))

  expect_equal(graph_distance(empty, one_edge), 1)
  expect_equal(graph_distance(one_edge, empty), 1)
  expect_equal(graph_distance(one_edge, one_edge), 0)
})

test_that("central_graph is the element-wise graph mean", {
  empty <- graph_from_edges(3L)
  one_edge <- graph_from_edges(3L, c(1L, 2L))

  observed <- central_graph(list(empty, one_edge))
  expect_equal(observed, (empty + one_edge) / 2)
  expect_equal(observed[1L, 2L], 0.5)
})

test_that("two-node populations are supported end to end", {
  empty <- graph_from_edges(2L)
  connected <- graph_from_edges(2L, c(1L, 2L))
  data <- brainnet_data(list(
    GroupA = list(empty, empty, connected),
    GroupB = list(connected, connected, empty)
  ))

  result <- brainnet_test(data, exact = "always")

  expect_s3_class(result, "brainnet_result")
  expect_equal(nrow(result$edges), 1L)
  expect_identical(result$edges$node1, 1L)
  expect_identical(result$edges$node2, 2L)
})

test_that("the canonical score equals its edge contributions", {
  data <- tiny_brainnet_data()
  encoded <- BrainNetTest:::.brainnet_encode(data)
  counts <- BrainNetTest:::.brainnet_group_counts(
    encoded$edge_matrix,
    encoded$group
  )

  score <- BrainNetTest:::.brainnet_score(counts, data$group_sizes)
  contributions <- BrainNetTest:::.brainnet_edge_contributions(
    counts,
    data$group_sizes
  )

  expect_equal(score, sum(contributions), tolerance = 1e-12)
})

test_that("prefix subtraction equals physical edge removal", {
  data <- tiny_brainnet_data()
  encoded <- BrainNetTest:::.brainnet_encode(data)
  counts <- BrainNetTest:::.brainnet_group_counts(
    encoded$edge_matrix,
    encoded$group
  )
  contributions <- BrainNetTest:::.brainnet_edge_contributions(
    counts,
    data$group_sizes
  )

  removed_populations <- lapply(data$populations, function(group) {
    lapply(group, function(graph) {
      graph[1L, 2L] <- graph[2L, 1L] <- 0L
      graph
    })
  })
  removed_data <- brainnet_data(removed_populations)
  removed_encoded <- BrainNetTest:::.brainnet_encode(removed_data)
  removed_counts <- BrainNetTest:::.brainnet_group_counts(
    removed_encoded$edge_matrix,
    removed_encoded$group
  )

  expect_equal(
    BrainNetTest:::.brainnet_score(removed_counts, data$group_sizes),
    BrainNetTest:::.brainnet_score(counts, data$group_sizes) -
      contributions[1L],
    tolerance = 1e-12
  )
})

test_that("exact assignments are unique and preserve group sizes", {
  assignments <- BrainNetTest:::.enumerate_group_assignments(c(2L, 2L))

  expect_equal(nrow(assignments), 6L)
  expect_equal(nrow(unique(assignments)), 6L)
  expect_true(all(apply(assignments, 1L, function(x) {
    identical(tabulate(x, nbins = 2L), c(2L, 2L))
  })))
})

test_that("Monte Carlo assignments sample the complete fixed-size orbit", {
  set.seed(123)
  assignments <- BrainNetTest:::.sample_group_assignments(
    c(2L, 2L),
    5000L
  )
  frequencies <- table(apply(assignments, 1L, paste, collapse = ""))

  expect_length(frequencies, 6L)
  expect_lt(max(abs(as.numeric(frequencies) - 5000 / 6)), 100)
})

test_that("Monte Carlo randomization controls size on a finite null orbit", {
  group_sizes <- c(2L, 3L)
  edge_matrix <- matrix(c(0, 0, 1, 1, 1), ncol = 1L)
  observed_assignments <- BrainNetTest:::.enumerate_group_assignments(
    group_sizes
  )
  rejected <- logical(nrow(observed_assignments) * 500L)
  output_index <- 0L

  set.seed(20260713)
  for (observed_index in seq_len(nrow(observed_assignments))) {
    observed_labels <- observed_assignments[observed_index, ]
    observed <- BrainNetTest:::.brainnet_score(
      BrainNetTest:::.brainnet_group_counts(
        edge_matrix,
        observed_labels,
        n_groups = 2L
      ),
      group_sizes
    )

    for (monte_carlo_index in seq_len(500L)) {
      assignments <- BrainNetTest:::.sample_group_assignments(
        group_sizes,
        19L
      )
      scores <- BrainNetTest:::.brainnet_scores_for_assignments(
        edge_matrix,
        assignments,
        group_sizes
      )
      tolerance <- sqrt(.Machine$double.eps) *
        max(1, abs(observed), abs(scores))
      p_value <- (
        1 + sum(abs(scores) >= abs(observed) - tolerance)
      ) / 20
      output_index <- output_index + 1L
      rejected[[output_index]] <- p_value <= 0.05
    }
  }

  expect_lte(mean(rejected), 0.06)
})

test_that("exact assignment enumeration supports multiple groups", {
  assignments <- BrainNetTest:::.enumerate_group_assignments(c(2L, 2L, 2L))

  expect_equal(nrow(assignments), 90L)
  expect_equal(nrow(unique(assignments)), 90L)
  expect_true(all(apply(assignments, 1L, function(x) {
    identical(tabulate(x, nbins = 3L), c(2L, 2L, 2L))
  })))
})

test_that("vectorized assignment scores equal direct scores", {
  data <- tiny_brainnet_data()
  encoded <- BrainNetTest:::.brainnet_encode(data)
  assignments <- BrainNetTest:::.enumerate_group_assignments(data$group_sizes)

  vectorized <- BrainNetTest:::.brainnet_scores_for_assignments(
    encoded$edge_matrix,
    assignments,
    data$group_sizes
  )
  direct <- apply(assignments, 1L, function(labels) {
    BrainNetTest:::.brainnet_score(
      BrainNetTest:::.brainnet_group_counts(
        encoded$edge_matrix,
        labels,
        n_groups = length(data$group_sizes)
      ),
      data$group_sizes
    )
  })

  expect_equal(vectorized, direct, tolerance = 1e-12)
})

test_that("vectorized assignment scores equal direct scores for three groups", {
  empty <- graph_from_edges(3L)
  edge_12 <- graph_from_edges(3L, c(1L, 2L))
  data <- brainnet_data(list(
    GroupA = list(empty, edge_12),
    GroupB = list(empty, empty),
    GroupC = list(edge_12, edge_12)
  ))
  encoded <- BrainNetTest:::.brainnet_encode(data)
  assignments <- BrainNetTest:::.enumerate_group_assignments(data$group_sizes)

  vectorized <- BrainNetTest:::.brainnet_scores_for_assignments(
    encoded$edge_matrix,
    assignments,
    data$group_sizes
  )
  direct <- apply(assignments, 1L, function(labels) {
    BrainNetTest:::.brainnet_score(
      BrainNetTest:::.brainnet_group_counts(
        encoded$edge_matrix,
        labels,
        n_groups = length(data$group_sizes)
      ),
      data$group_sizes
    )
  })

  expect_equal(vectorized, direct, tolerance = 1e-12)
})

test_that("exact global p-value matches exhaustive enumeration", {
  data <- tiny_brainnet_data()
  encoded <- BrainNetTest:::.brainnet_encode(data)
  assignments <- BrainNetTest:::.enumerate_group_assignments(data$group_sizes)
  scores <- BrainNetTest:::.brainnet_scores_for_assignments(
    encoded$edge_matrix,
    assignments,
    data$group_sizes
  )
  observed <- BrainNetTest:::.brainnet_score(
    BrainNetTest:::.brainnet_group_counts(
      encoded$edge_matrix,
      encoded$group
    ),
    data$group_sizes
  )
  tolerance <- sqrt(.Machine$double.eps) *
    max(1, abs(observed), abs(scores))
  expected <- mean(abs(scores) >= abs(observed) - tolerance)

  expect_warning(
    result <- brainnet_test(data, exact = "always", alpha = 0.01),
    "cannot attain"
  )

  expect_s3_class(result, "brainnet_result")
  expect_identical(result$global$permutation_method, "exact")
  expect_equal(result$global$p_value, expected)
  expect_equal(unname(result$global$statistic), observed)
})

test_that("floating-point equivalents are included as permutation ties", {
  empty <- graph_from_edges(3L)
  edge_12 <- graph_from_edges(3L, c(1L, 2L))
  data <- brainnet_data(list(
    GroupA = list(empty, edge_12),
    GroupB = list(empty, empty),
    GroupC = list(edge_12, edge_12)
  ))
  encoded <- BrainNetTest:::.brainnet_encode(data)
  assignments <- BrainNetTest:::.enumerate_group_assignments(data$group_sizes)
  scores <- BrainNetTest:::.brainnet_scores_for_assignments(
    encoded$edge_matrix,
    assignments,
    data$group_sizes
  )
  observed <- BrainNetTest:::.brainnet_score(
    BrainNetTest:::.brainnet_group_counts(
      encoded$edge_matrix,
      encoded$group
    ),
    data$group_sizes
  )
  tolerance <- sqrt(.Machine$double.eps) *
    max(1, abs(observed), abs(scores))
  expected <- mean(abs(scores) >= abs(observed) - tolerance)

  result <- brainnet_test(data, exact = "always")

  expect_equal(result$global$p_value, expected)
  expect_gt(result$global$n_ties, 0L)
})

test_that("global and edge inference are invariant to group and node order", {
  data <- tiny_brainnet_data()
  original <- brainnet_test(data, exact = "always")

  swapped_data <- brainnet_data(data$populations[c("GroupB", "GroupA")])
  swapped <- brainnet_test(swapped_data, exact = "always")
  expect_equal(original$global$p_value, swapped$global$p_value)
  expect_equal(
    sort(original$edges$p_value),
    sort(swapped$edges$p_value)
  )

  permutation <- c(3L, 1L, 2L)
  permuted_populations <- lapply(data$populations, function(group) {
    lapply(group, function(graph) {
      graph[permutation, permutation]
    })
  })
  permuted <- brainnet_test(
    brainnet_data(permuted_populations),
    exact = "always"
  )
  expect_equal(original$global$p_value, permuted$global$p_value)
  expect_equal(
    sort(original$edges$p_value),
    sort(permuted$edges$p_value)
  )
})

test_that("degenerate tied data returns p = 1", {
  empty <- graph_from_edges(3L)
  data <- brainnet_data(list(
    GroupA = constant_group(empty, 3L),
    GroupB = constant_group(empty, 3L)
  ))

  exact <- brainnet_test(data, exact = "always")
  monte_carlo <- brainnet_test(
    data,
    exact = "never",
    n_permutations = 19L,
    seed = 42L
  )

  expect_equal(exact$global$p_value, 1)
  expect_equal(monte_carlo$global$p_value, 1)
  expect_equal(exact$global$n_ties, exact$global$n_assignments)
})

test_that("Monte Carlo inference is reproducible and preserves caller RNG", {
  data <- tiny_brainnet_data()
  set.seed(100)
  before <- .Random.seed

  result_1 <- brainnet_test(
    data,
    exact = "never",
    n_permutations = 49L,
    seed = 99L
  )
  after <- .Random.seed
  result_2 <- brainnet_test(
    data,
    exact = "never",
    n_permutations = 49L,
    seed = 99L
  )

  expect_identical(after, before)
  expect_identical(result_1$global$p_value, result_2$global$p_value)
  expect_equal(result_1$global$resolution, 1 / 50)
  expect_true(result_1$global$p_value >= result_1$global$resolution)
})

test_that("seeded inference restores an absent RNG state even after errors", {
  global <- globalenv()
  had_seed <- exists(".Random.seed", envir = global, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = global, inherits = FALSE)
  }
  on.exit(
    {
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = global)
      } else if (exists(".Random.seed", envir = global, inherits = FALSE)) {
        rm(".Random.seed", envir = global)
      }
    },
    add = TRUE
  )
  if (exists(".Random.seed", envir = global, inherits = FALSE)) {
    rm(".Random.seed", envir = global)
  }

  brainnet_test(
    tiny_brainnet_data(),
    exact = "never",
    n_permutations = 19L,
    seed = 1L
  )
  expect_false(exists(".Random.seed", envir = global, inherits = FALSE))

  expect_error(
    BrainNetTest:::.with_preserved_seed(1L, function() stop("planned")),
    "planned"
  )
  expect_false(exists(".Random.seed", envir = global, inherits = FALSE))
})

test_that("brainnet_test validates inferential controls", {
  data <- tiny_brainnet_data()

  expect_error(brainnet_test(list()), "brainnet_data")
  expect_error(brainnet_test(data, alpha = 0), "alpha")
  expect_error(brainnet_test(data, alpha = 1), "alpha")
  expect_error(
    brainnet_test(data, n_permutations = 10.5, exact = "never"),
    "whole number"
  )
  expect_error(brainnet_test(data, exact = "sometimes"), "exact")
  expect_error(
    brainnet_test(data, exact = "always", max_exact = 10L),
    "above"
  )
  expect_error(
    brainnet_test(data, max_exact = 1000001L),
    "cannot exceed"
  )
  expect_error(
    brainnet_test(
      data,
      exact = "never",
      n_permutations = .Machine$integer.max + 1
    ),
    "through"
  )
  expect_equal(
    BrainNetTest:::.number_group_assignments(c(100L, 100L)),
    Inf
  )
})
