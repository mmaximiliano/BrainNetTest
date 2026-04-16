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

test_that("rank_edges handles the 2-node degenerate case", {
  pv <- matrix(c(1, 0.02, 0.02, 1), 2, 2)
  edf <- rank_edges(pv)
  expect_equal(nrow(edf), 1L)
  expect_equal(edf$node1, 1L)
  expect_equal(edf$node2, 2L)
  expect_equal(edf$p_value, 0.02)
})
