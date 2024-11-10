# tests/testthat/test-brainnettest.R

library(testthat)
library(BrainNetTest)

test_that("compute_central_graph works correctly", {
  G1 <- matrix(c(0,1,0,1,0,1,0,1,0), nrow = 3)
  G2 <- matrix(c(0,0,1,0,0,1,1,0,0), nrow = 3)
  central <- compute_central_graph(list(G1, G2))
  expected <- (G1 + G2) / 2
  expect_equal(central, expected)
})

test_that("compute_distance computes Manhattan norm correctly", {
  G <- matrix(c(0,1,0,1,0,1,0,1,0), nrow = 3)
  M <- matrix(c(0,0,1,0,0,1,1,0,0), nrow = 3)
  distance <- compute_distance(G, M)
  expected <- sum(abs(G - M))
  expect_equal(distance, expected)
})

test_that("compute_test_statistic computes T correctly", {
  # Generate simple synthetic data
  G1 <- matrix(c(0,1,0,1,0,1,0,1,0), nrow = 3)
  G2 <- matrix(c(0,0,1,0,0,1,1,0,0), nrow = 3)
  populations <- list(
    Control = list(G1, G1),
    Disease = list(G2, G2)
  )
  T_value <- compute_test_statistic(populations, a = 1)
  expect_type(T_value, "double")
})
