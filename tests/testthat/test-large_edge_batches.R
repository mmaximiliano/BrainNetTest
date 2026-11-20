test_that("oversized batches remove the complete candidate set safely", {
  empty <- matrix(0L, 4, 4)
  full <- matrix(1L, 4, 4); diag(full) <- 0L
  pops <- list(A = rep(list(empty), 10), B = rep(list(full), 10))
  for (criterion in c("original", "combined")) {
    whole <- identify_critical_links(pops, criterion = criterion,
      batch_size = 6, n_permutations = 99, seed = 1)
    expect_lte(whole$p_value, 0.05)
    for (batch in c(7L, 100L, .Machine$integer.max)) {
      res <- identify_critical_links(pops, criterion = criterion,
        batch_size = batch, n_permutations = 99, seed = 1)
      expect_equal(res$critical_edges, whole$critical_edges)
      expect_equal(res$modified_populations, whole$modified_populations)
      expect_equal(length(res$edges_removed), 6L)
      expect_equal(res$settings$batch_size, batch)
    }
  }
})

test_that("a final partial batch removes the remaining edges", {
  empty <- matrix(0L, 4, 4)
  full <- matrix(1L, 4, 4); diag(full) <- 0L
  pops <- list(A = rep(list(empty), 10), B = rep(list(full), 10))
  for (criterion in c("original", "combined")) {
    for (batch in c(2L, 4L, 5L)) {
      res <- identify_critical_links(pops, criterion = criterion,
        batch_size = batch, n_permutations = 99, seed = 1)
      expect_equal(nrow(res$critical_edges), 6L)
      expect_true(all(vapply(unlist(res$modified_populations, recursive = FALSE),
                             function(g) all(g == 0), logical(1))))
    }
  }
})

test_that("oversized batches do not bypass a nonsignificant initial test", {
  G <- matrix(0L, 3, 3)
  pops <- list(A = rep(list(G), 4), B = rep(list(G), 4))
  for (criterion in c("original", "combined")) {
    res <- suppressWarnings(identify_critical_links(pops, criterion = criterion,
      batch_size = 100, n_permutations = 99, seed = 1))
    expect_equal(res$p_value, 1)
    expect_null(res$critical_edges)
    expect_equal(res$modified_populations, pops)
  }
})
