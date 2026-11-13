test_that("both inference entry points reject invalid graph data", {
  G <- matrix(0, 3, 3)
  bad <- list(
    asymmetric = { x <- G; x[2, 1] <- 1; x },
    self_loop = { x <- G; x[2, 2] <- 1; x },
    weighted = { x <- G; x[1, 2] <- x[2, 1] <- 0.5; x },
    missing = { x <- G; x[1, 2] <- NA_real_; x },
    infinite = { x <- G; x[1, 2] <- Inf; x },
    nan = { x <- G; x[1, 2] <- NaN; x },
    empty = matrix(numeric(), 0, 0)
  )
  errors <- c("symmetric", "zero diagonal", rep("not binary", 4),
              "at least one node")
  for (f in list(global_test, identify_critical_links)) {
    for (criterion in c("original", "combined")) {
      for (i in seq_along(bad)) {
        # Repeat the invalid graph to isolate its own input error.
        pops <- list(A = rep(list(bad[[i]]), 2), B = rep(list(bad[[i]]), 2))
        expect_error(f(pops, criterion = criterion, n_permutations = 9),
                     errors[i], info = names(bad)[i])
      }
    }
  }
})

test_that("every inference population needs two graphs, even without edges", {
  for (n_nodes in c(1L, 3L)) {
    G <- matrix(0, n_nodes, n_nodes)
    for (pops in list(list(A = list(G), B = list(G, G)),
                      list(A = list(G, G), B = list(G)))) {
      for (f in list(global_test, identify_critical_links)) {
        for (criterion in c("original", "combined")) {
          expect_error(f(pops, criterion = criterion, n_permutations = 9),
                       "at least two graphs")
        }
      }
      expect_error(BrainNetTest:::identify_critical_links_naive(pops),
                   "at least two graphs")
    }
  }
})

test_that("valid one-node networks keep their no-inference result", {
  G <- matrix(FALSE, 1, 1)
  pops <- list(A = list(G, G), B = list(G, G))
  for (criterion in c("original", "combined")) {
    gt <- global_test(pops, criterion = criterion, n_permutations = 9)
    id <- identify_critical_links(pops, criterion = criterion, n_permutations = 9)
    expect_equal(gt$statistic, 0)
    expect_true(is.na(gt$p_value))
    expect_true(is.na(id$p_value))
    expect_null(id$critical_edges)
  }
})

test_that("generic graph helpers retain finite weighted input support", {
  G <- matrix(c(0.1, 0.3, 0.7, 0.2), 2, 2)
  pops <- list(A = list(G, G), B = list(G, G))
  fr <- compute_edge_frequencies(pops)
  expect_equal(fr$edge_counts[, , 1], 2 * G)
  expect_equal(compute_central_graph(list(G, G)), G)
  expect_true(is.finite(compute_test_statistic(pops)))
  G[1, 1] <- Inf
  expect_error(compute_edge_frequencies(list(A = list(G))), "finite entries")
})
