test_that("all inference entry points reject invalid methods eagerly", {
  invalid_methods <- list(NULL, NA_character_, 1, c("fisher", "prop"),
                          "unknown", "fish", character())
  for (n_nodes in c(1L, 3L)) {
    graph <- matrix(0, n_nodes, n_nodes)
    populations <- list(A = rep(list(graph), 4), B = rep(list(graph), 4))
    freq <- compute_edge_frequencies(populations)
    for (method in invalid_methods) {
      for (criterion in c("original", "combined")) {
        expect_error(global_test(populations, n_permutations = 1,
                                 method = method, criterion = criterion),
                     "Invalid method specified")
        expect_error(identify_critical_links(populations, n_permutations = 1,
                                            method = method, criterion = criterion),
                     "Invalid method specified")
      }
      expect_error(compute_edge_pvalues(freq$edge_counts, c(4L, 4L), method),
                   "Invalid method specified")
      expect_error(BrainNetTest:::identify_critical_links_naive(
        populations, n_permutations = 1, method = method),
        "Invalid method specified")
    }
  }
})

test_that("prop consistently requires two populations even without edges", {
  for (n_nodes in c(1L, 3L)) {
    graph <- matrix(0, n_nodes, n_nodes)
    populations <- rep(list(rep(list(graph), 4)), 3)
    freq <- compute_edge_frequencies(populations)
    for (criterion in c("original", "combined")) {
      expect_error(global_test(populations, n_permutations = 1,
                               method = "prop", criterion = criterion),
                   "only implemented for two populations")
      expect_error(identify_critical_links(populations, n_permutations = 1,
                                          method = "prop", criterion = criterion),
                   "only implemented for two populations")
    }
    expect_error(compute_edge_pvalues(freq$edge_counts, c(4L, 4L, 4L), "prop"),
                 "only implemented for two populations")
  }
})

test_that("supported methods retain agreement across inference functions", {
  for (n_populations in c(2L, 3L)) {
    populations <- rep(list(rep(list(matrix(0, 3, 3)), 4)), n_populations)
    methods <- c("fisher", "chi.squared")
    if (n_populations == 2L) methods <- c(methods, "prop")
    for (method in methods) {
      for (criterion in c("original", "combined")) {
        global <- global_test(populations, n_permutations = 3, seed = 1,
                              method = method, criterion = criterion)
        identified <- suppressWarnings(identify_critical_links(
          populations, n_permutations = 3, seed = 1,
          method = method, criterion = criterion))
        expect_equal(global$p_value, 1)
        expect_identical(global$p_value, identified$p_value)
        expect_identical(identified$settings$method, method)
      }
    }
  }
})

test_that("adjust_method is validated eagerly across edge-ranking functions", {
  invalid_adj <- list(NULL, NA_character_, 1, c("BH", "holm"),
                      "unknown", "bh", "bonf", character())
  for (n_nodes in c(1L, 3L)) {
    graph <- matrix(0, n_nodes, n_nodes)
    populations <- list(A = rep(list(graph), 4), B = rep(list(graph), 4))
    freq <- compute_edge_frequencies(populations)
    for (adj in invalid_adj) {
      for (criterion in c("original", "combined")) {
        expect_error(identify_critical_links(
          populations, n_permutations = 1, adjust_method = adj,
          criterion = criterion), "Invalid `adjust_method` specified")
      }
      expect_error(compute_edge_pvalues(
        freq$edge_counts, c(4L, 4L), adjust_method = adj),
        "Invalid `adjust_method` specified")
      expect_error(BrainNetTest:::identify_critical_links_naive(
        populations, n_permutations = 1, adjust_method = adj),
        "Invalid `adjust_method` specified")
    }
    for (adj in stats::p.adjust.methods) {
      pv <- compute_edge_pvalues(freq$edge_counts, c(4L, 4L),
                                 adjust_method = adj)
      expect_equal(dim(pv), c(n_nodes, n_nodes))
    }
  }
})
