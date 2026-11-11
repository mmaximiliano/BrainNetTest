## Independent count oracle: for balanced groups of size q, T is a positive
## constant times this INTEGER score (not the package's floating formula).
integer_score <- function(counts, q) {
  n <- length(counts) * q
  total <- sum(counts)
  (n - q) * total - (n - 1) * sum(counts^2) + (q - 1) * total^2
}
sparse_populations <- function(q, counts) {
  lapply(counts, function(count) lapply(seq_len(q), function(i) {
    A <- matrix(0, 3, 3)
    A[1, 2] <- A[2, 1] <- as.integer(i <= count)
    A
  }))
}
test_that("exact balanced count assignments retain the correct T ties", {
  q <- 3L
  for (m in 2:4) {
    grid <- as.matrix(expand.grid(rep(list(0:q), m)))
    for (total in 0:(m * q)) {
      counts <- grid[rowSums(grid) == total, , drop = FALSE]
      copies <- apply(counts, 1, function(x) prod(choose(q, x)))
      oracle <- apply(counts, 1, integer_score, q = q)
      deltas <- .permuted_edge_deltas(lapply(seq_len(m), function(k)
        matrix(counts[, k], nrow = 1)), rep(q, m))[1, ]
      expanded <- rep(seq_len(nrow(counts)), copies)
      tolerance <- .statistic_tolerance(rep(q, m), 1L)
      expect_equal(.permutation_ranks(deltas[expanded], tolerance),
        rank(oracle[expanded], ties.method = "max"),
        info = paste("groups", m, "total", total))
      for (row in seq_len(nrow(counts))) {
        ec <- array(0, c(2, 2, m)); ec[1, 2, ] <- counts[row, ]
        expect_identical(.observed_edge_deltas(ec, rep(q, m))$deltas,
                         unname(deltas[row]))
      }
    }
  }
})
test_that("roundoff near zero cannot invent a significant sparse difference", {
  pops <- sparse_populations(3, c(0, 1))
  for (criterion in c("original", "combined")) {
    gt <- global_test(pops, n_permutations = 99, seed = 1,
                      criterion = criterion)
    expect_identical(gt$p_value, 1)
    result <- suppressWarnings(identify_critical_links(pops,
      n_permutations = 99, seed = 1, criterion = criterion))
    expect_identical(result$p_value, 1)
    expect_length(result$edges_removed, 0)
  }
  slow <- suppressWarnings(identify_critical_links_naive(pops,
    n_permutations = 99, seed = 1))
  expect_identical(slow$p_value, 1)
  expect_length(slow$edges_removed, 0)
})
test_that("20-per-group sparse results match an integer oracle in either order", {
  ## Exhaust all choose(40, 3) placements of the three present graphs.
  counts <- cbind(0:3, 3:0)
  copies <- apply(counts, 1, function(x) prod(choose(20, x)))
  expanded <- rep(seq_len(nrow(counts)), copies)
  deltas <- .permuted_edge_deltas(lapply(1:2, function(k)
    matrix(counts[, k], nrow = 1)), c(20, 20))[1, ]
  exact <- .plus_one_pvalue(deltas[expanded],
                            .statistic_tolerance(c(20, 20), 1))
  expect_equal(exact, 3 / 13)  # 0.230769..., not the old 0.001
  for (counts in list(c(3, 0), c(0, 3))) {
    pops <- sparse_populations(20, counts)
    observed <- integer_score(counts, 20)
    present <- unlist(lapply(counts, function(k) seq_len(20) <= k))
    set.seed(1)
    score <- replicate(999, {
      order <- sample.int(40)
      integer_score(c(sum(present[order[1:20]]),
                      sum(present[order[21:40]])), 20)
    })
    expected <- (1 + sum(score <= observed)) / 1000
    for (criterion in c("original", "combined")) {
      gt <- global_test(pops, n_permutations = 999, seed = 1,
                        criterion = criterion)
      expect_identical(gt$p_value, expected)
      expect_gt(gt$p_value, .2)
      fast <- suppressWarnings(identify_critical_links(pops,
        n_permutations = 999, seed = 1, criterion = criterion))
      expect_identical(fast$p_value, expected)
      expect_length(fast$edges_removed, 0)
    }
    slow <- suppressWarnings(identify_critical_links_naive(pops,
      n_permutations = 999, seed = 1))
    expect_identical(slow$p_value, expected)
    expect_length(slow$edges_removed, 0)
  }
})
test_that("small genuine differences and tiny marginal p-values stay distinct", {
  expect_identical(.plus_one_pvalue(c(1e-100, 2e-100, 3e-100)), 1 / 3)
  expect_equal(.permutation_ranks(c(1e-100, 2e-100, 3e-100)), 1:3)
  expect_equal(.permutation_ranks(c(-1e-100, 0, 1e-100)), 1:3)
  expect_identical(.combined_pvalue(c(1, 1, 1),
    c(1e-100, 2e-100, 3e-100)), 1 / 3)
  expect_equal(.permutation_ranks(c(0, 1e-12, 2e-12), 1e-14), 1:3)
  expect_equal(.permutation_ranks(c(-1e-16, 0, 1e-16), 1e-14), rep(3, 3))
  expect_error(.plus_one_pvalue(c(0, NA_real_)), "finite")
  expect_error(.combined_pvalue(c(0, Inf), c(1, 1)), "finite")
})
test_that("combined calibration is equivariant to the observed assignment", {
  T <- c(0, 1e-16, -1e-16, -1, -1 + 1e-16, 2)
  U <- c(1e-100, 2e-100, 3e-100, .5, .2, .8)
  p <- vapply(seq_along(T), function(j) {
    indices <- c(j, setdiff(seq_along(T), j))
    .combined_pvalue(T[indices], U[indices], tolerance = 1e-14)
  }, numeric(1))
  for (alpha in seq_along(T) / length(T))
    expect_lte(mean(p <= alpha), alpha)
  for (j in seq_along(T)) {
    indices <- c(j, rev(setdiff(seq_along(T), j)))
    expect_identical(.combined_pvalue(T[indices], U[indices], 1e-14), p[j])
  }
})
test_that("retained sums avoid subtracting a nearly equal removed prefix", {
  expect_identical(.suffix_sum(c(1e16, 1))[[2]], 1)
  expect_identical(.suffix_sum(c(1e16, -1e16, 1))[[3]], 1)
  expect_identical(.statistic_tolerance(c(20, 20), 0), 0)
  pops <- sparse_populations(20, c(20, 0))
  pops[[1]][[1]][1, 3] <- pops[[1]][[1]][3, 1] <- 1
  for (criterion in c("original", "combined")) {
    result <- identify_critical_links(pops, n_permutations = 199, seed = 2,
                                      criterion = criterion)
    expect_identical(result$edges_removed, list(c(1L, 2L)))
    for (a in c(.Machine$double.xmin, .1, 1e10, .Machine$double.xmax)) {
      rescaled <- identify_critical_links(pops, n_permutations = 199, seed = 2,
                                          criterion = criterion, a = a)
      expect_identical(rescaled$p_value, result$p_value)
      expect_identical(rescaled$edges_removed, result$edges_removed)
    }
  }
  slow <- identify_critical_links_naive(pops, n_permutations = 199, seed = 2)
  expect_identical(slow$edges_removed, list(c(1L, 2L)))
})

