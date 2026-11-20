sparse_populations <- function(q, counts) {
  lapply(counts, function(count) lapply(seq_len(q), function(i) {
    A <- matrix(0, 3, 3)
    A[1, 2] <- A[2, 1] <- as.integer(i <= count)
    A
  }))
}
test_that("permutation p-values use the stored values and include exact ties", {
  T <- c(1, 1 + .Machine$double.eps, 1, 2, 3)
  U <- c(0.2, 0.5, 0.2, 0.1, 0.9)
  expect_identical(.plus_one_pvalue(T), 2 / 5)
  expect_identical(.combined_pvalue(T, U), 3 / 5)
  expect_identical(.combined_pvalue(rep(1, 5), rep(0.1, 5)), 1)
})

test_that("reported global p-values reproduce the manuscript rank on the same draws", {
  pops <- sparse_populations(20L, c(3L, 0L))
  original <- global_test(pops, n_permutations = 99, seed = 1, criterion = "original")
  observed_and_permuted <- c(original$statistic, original$null_distribution)
  expect_identical(original$p_value, mean(observed_and_permuted <= original$statistic))

  for (criterion in c("original", "combined")) {
    global <- global_test(pops, n_permutations = 99, seed = 1, criterion = criterion)
    selection <- suppressWarnings(identify_critical_links(pops, n_permutations = 99,
      seed = 1, criterion = criterion))
    expect_identical(selection$p_value, global$p_value)
  }
})

test_that("positive normalisation leaves edge selection and p-values unchanged", {
  pops <- sparse_populations(20, c(20, 0))
  pops[[1]][[1]][1, 3] <- pops[[1]][[1]][3, 1] <- 1
  for (criterion in c("original", "combined")) {
    result <- identify_critical_links(pops, n_permutations = 199, seed = 2,
                                      criterion = criterion)
    for (a in c(.Machine$double.xmin, .1, 1e10, .Machine$double.xmax)) {
      rescaled <- identify_critical_links(pops, n_permutations = 199, seed = 2,
                                          criterion = criterion, a = a)
      expect_identical(rescaled$p_value, result$p_value)
      expect_identical(rescaled$edges_removed, result$edges_removed)
    }
  }
})

test_that("constant all-zero and all-one edges leave the global test unchanged", {
  small <- sparse_populations(10, c(2, 0))
  padded <- lapply(small, function(pop) lapply(pop, function(G) {
    M <- matrix(0L, 12, 12)
    M[1:3, 1:3] <- G
    ## Add constant all-one edges across every graph of both populations.
    M[4, 5] <- M[5, 4] <- 1L
    M[6, 7] <- M[7, 6] <- 1L
    M
  }))
  g_small  <- global_test(small, n_permutations = 199, seed = 9,
                          criterion = "original")
  g_padded <- global_test(padded, n_permutations = 199, seed = 9,
                          criterion = "original")
  expect_identical(g_small$statistic, g_padded$statistic)
  expect_identical(g_small$p_value, g_padded$p_value)
})

