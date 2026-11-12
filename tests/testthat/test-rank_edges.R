# tests/testthat/test-rank_edges.R

test_that("rank_edges orders by ascending p-value", {
  pv <- matrix(1, 4, 4)
  pv[1, 2] <- pv[2, 1] <- 0.01
  pv[1, 3] <- pv[3, 1] <- 0.50
  pv[1, 4] <- pv[4, 1] <- 0.20
  pv[2, 3] <- pv[3, 2] <- 0.05
  pv[2, 4] <- pv[4, 2] <- 0.80
  pv[3, 4] <- pv[4, 3] <- 0.001

  edf <- rank_edges(pv)

  expect_s3_class(edf, "data.frame")
  expect_equal(names(edf), c("node1", "node2", "p_value"))
  # 4 nodes -> 6 upper-triangle edges
  expect_equal(nrow(edf), 6L)
  # Sorted ascending
  expect_true(all(diff(edf$p_value) >= 0))
  # Most significant first: (3,4) with 0.001
  expect_equal(edf$node1[1], 3L)
  expect_equal(edf$node2[1], 4L)
  expect_equal(edf$p_value[1], 0.001)
})

test_that("rank_edges only reports upper-triangle edges (i < j)", {
  set.seed(1)
  pv <- matrix(runif(25), 5, 5)
  pv <- (pv + t(pv)) / 2            # symmetrise
  diag(pv) <- 1
  edf <- rank_edges(pv)
  expect_true(all(edf$node1 < edf$node2))
  expect_equal(nrow(edf), 5 * 4 / 2)
})

test_that("rank_edges returns consecutive row names", {
  pv <- matrix(1, 4, 4)
  pv[1, 2] <- pv[2, 1] <- 0.90
  pv[1, 3] <- pv[3, 1] <- 0.10
  pv[2, 4] <- pv[4, 2] <- 0.01
  edf <- rank_edges(pv)
  # Sorting used to leave the pre-sort positions as row names, so the ranking
  # printed as e.g. 31, 39, 41 instead of 1, 2, 3.
  expect_identical(rownames(edf), as.character(seq_len(nrow(edf))))
})

test_that("rank_edges handles the 2-node degenerate case", {
  pv <- matrix(c(1, 0.02, 0.02, 1), 2, 2)
  edf <- rank_edges(pv)
  expect_equal(nrow(edf), 1L)
  expect_equal(edf$node1, 1L)
  expect_equal(edf$node2, 2L)
  expect_equal(edf$p_value, 0.02)
})

test_that("machine-level Fisher ties use a stable endpoint order", {
  # These values exchanged edges at the stopping boundary of the manuscript's
  # 64-node distance example across R versions, despite the same random seed.
  tied <- c(0.048337915123163772979, 0.048337915123163759101)
  for (values in list(tied, rev(tied))) {
    pv <- matrix(1, 64, 64)
    pv[1, 2] <- values[1]
    pv[38, 44] <- values[2]
    edf <- rank_edges(pv)
    expect_identical(edf$node1[1:2], c(1L, 38L))
    expect_identical(edf$node2[1:2], c(2L, 44L))
    expect_identical(edf$p_value[1:2], values)
  }
})

test_that("relative ties preserve tiny p-values and distinguish zero", {
  pv <- matrix(1, 4, 4)
  pv[1, 2] <- 1e-300 * (1 + 4 * .Machine$double.eps)
  pv[1, 3] <- 1e-300
  pv[2, 3] <- 2e-300
  pv[1, 4] <- 0
  edf <- rank_edges(pv)
  expect_identical(edf$node1[1:4], c(1L, 1L, 1L, 2L))
  expect_identical(edf$node2[1:4], c(4L, 2L, 3L, 3L))
  expect_identical(edf$p_value[1:4], pv[cbind(c(1, 1, 1, 2),
                                           c(4, 2, 3, 3))])
})

test_that("meaningfully distinct p-values keep their original order", {
  pv <- matrix(1, 4, 4)
  pv[1, 2] <- 0.05 + 1e-12
  pv[1, 3] <- 0.05 - 1e-12
  pv[2, 3] <- 0.1 * (1 + 128 * .Machine$double.eps)
  pv[1, 4] <- 0.1
  edf <- rank_edges(pv)
  expect_identical(edf$node1[1:4], c(1L, 1L, 1L, 2L))
  expect_identical(edf$node2[1:4], c(3L, 2L, 4L, 3L))
  # The ranking policy must never round a p-value across a test threshold.
  expect_identical(edf$p_value < 0.05,
                   pv[cbind(edf$node1, edf$node2)] < 0.05)
})

test_that("numerical tie groups do not grow by chaining neighbours", {
  pv <- matrix(1, 3, 3)
  pv[1, 2] <- 0.1 * (1 + 96 * .Machine$double.eps)
  pv[1, 3] <- 0.1 * (1 + 48 * .Machine$double.eps)
  pv[2, 3] <- 0.1
  edf <- rank_edges(pv)
  expect_identical(edf$node1, c(1L, 2L, 1L))
  expect_identical(edf$node2, c(3L, 3L, 2L))
})

test_that("adjusted p-value plateaus break ties by raw marginal evidence", {
  # Under Bonferroni on 6 edges, raw p = 0.20 and raw p = 1.00 both clamp to
  # adjusted p = 1.00, while raw p = 0.001 becomes 0.006. Passing raw p-values
  # as tie_values ranks the stronger raw edge (3, 4) ahead of dead edge (1, 2)
  # while keeping the adjusted p-values in `p_value`.
  raw <- matrix(1, 4, 4)
  raw[1, 3] <- raw[3, 1] <- 0.001
  raw[3, 4] <- raw[4, 3] <- 0.20
  adj <- matrix(1, 4, 4)
  adj[1, 3] <- adj[3, 1] <- stats::p.adjust(raw[upper.tri(raw)], "bonferroni")[2]

  edf <- rank_edges(adj, tie_values = raw)
  expect_identical(edf$node1[1:2], c(1L, 3L))
  expect_identical(edf$node2[1:2], c(3L, 4L))
  expect_identical(edf$p_value[1:2], c(0.006, 1))
  expect_error(rank_edges(adj, tie_values = matrix(1, 3, 3)), "same dimensions")
})
