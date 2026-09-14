# tests/testthat/test-global_test.R

make_populations <- function(inter_patient = 0.4, n_graphs = 20) {
  control <- generate_category_graphs(n_graphs = n_graphs, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.6, base_inter_prob = 0.2, seed = 1)
  patient <- generate_category_graphs(n_graphs = n_graphs, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.6, base_inter_prob = inter_patient,
    seed = 2)
  list(Control = control, Patient = patient)
}

## Two populations that differ only slightly, so that the p value falls well
## inside (0, 1) and comparisons between implementations are informative.
make_null_like <- function() {
  a <- generate_category_graphs(n_graphs = 15, n_nodes = 10, n_communities = 2,
    base_intra_prob = 0.5, base_inter_prob = 0.30, seed = 11)
  b <- generate_category_graphs(n_graphs = 15, n_nodes = 10, n_communities = 2,
    base_intra_prob = 0.5, base_inter_prob = 0.36, seed = 12)
  list(A = a, B = b)
}

## The permutation test written out directly from its definition: permute the
## pooled labels, split by the original group sizes, recompute T from the
## graph distances. Draws sample.int() once per replicate, as global_test()
## does, so the same seed gives the same permutations.
direct_null <- function(populations, n_permutations, seed) {
  set.seed(seed)
  graphs <- unlist(populations, recursive = FALSE)
  sizes  <- lengths(populations)
  bounds <- c(0L, cumsum(sizes))
  vapply(seq_len(n_permutations), function(b) {
    idx  <- sample.int(length(graphs))
    pops <- lapply(seq_along(sizes), function(k)
      graphs[idx[(bounds[k] + 1L):bounds[k + 1L]]])
    names(pops) <- names(populations)
    compute_test_statistic(pops, a = 1)
  }, numeric(1))
}

test_that("global_test returns a classed list with the documented components", {
  pops <- make_populations()
  gt <- global_test(pops, n_permutations = 200, seed = 42)

  expect_s3_class(gt, "global_test")
  expect_named(gt, c("statistic", "p_value", "null_distribution",
                     "n_permutations", "group_sizes", "n_nodes", "call"))
  expect_length(gt$statistic, 1L)
  expect_null(names(gt$statistic))
  expect_true(gt$p_value >= 0 && gt$p_value <= 1)
  expect_length(gt$null_distribution, 200L)
  expect_identical(gt$n_permutations, 200L)
  expect_identical(gt$group_sizes, c(Control = 20L, Patient = 20L))
  expect_identical(gt$n_nodes, 10L)
  expect_true(is.call(gt$call))
})

test_that("the statistic is compute_test_statistic() with a = 1", {
  pops <- make_populations()
  gt <- global_test(pops, n_permutations = 50, seed = 1)
  expect_equal(gt$statistic, compute_test_statistic(pops, a = 1))

  three <- list(A = pops$Control[1:7], B = pops$Control[8:14],
                C = pops$Patient[1:7])
  expect_equal(global_test(three, n_permutations = 50, seed = 1)$statistic,
               compute_test_statistic(three, a = 1))
})

test_that("the null distribution matches the direct permutation loop", {
  pops <- make_null_like()
  B <- 100
  gt <- global_test(pops, n_permutations = B, seed = 7)
  direct <- direct_null(pops, B, seed = 7)

  # Same permutations, same statistic: the vectorised path must reproduce the
  # distance-based computation to floating-point precision.
  expect_equal(gt$null_distribution, direct, tolerance = 1e-10)
})

test_that("the p value is a proportion of the replicates", {
  gt <- global_test(make_null_like(), n_permutations = 400, seed = 7)
  expect_equal(gt$p_value * 400, round(gt$p_value * 400))
  expect_equal(gt$p_value, mean(gt$null_distribution < gt$statistic))
})

test_that("the p value is the one identify_critical_links() records", {
  pops <- make_null_like()
  for (seed in 1:3) {
    gt  <- global_test(pops, n_permutations = 300, seed = seed)
    icl <- suppressWarnings(
      identify_critical_links(pops, n_permutations = 300, seed = seed))
    expect_identical(gt$p_value, icl$p_value)
  }

  pops <- make_populations()
  gt  <- global_test(pops, n_permutations = 200, seed = 42)
  icl <- identify_critical_links(pops, n_permutations = 200, seed = 42)
  expect_identical(gt$p_value, icl$p_value)
})

test_that("a clear difference is detected and near-identical populations are not", {
  strong <- global_test(make_populations(), n_permutations = 500, seed = 42)
  expect_lt(strong$statistic, 0)
  expect_identical(strong$p_value, 0)

  same <- generate_category_graphs(n_graphs = 30, n_nodes = 10,
    n_communities = 2, base_intra_prob = 0.5, base_inter_prob = 0.5, seed = 3)
  null <- global_test(list(A = same[1:15], B = same[16:30]),
                      n_permutations = 500, seed = 42)
  expect_gt(null$p_value, 0.05)
})

test_that("the seed makes the replicates reproducible", {
  pops <- make_null_like()
  g1 <- global_test(pops, n_permutations = 100, seed = 99)
  g2 <- global_test(pops, n_permutations = 100, seed = 99)
  g3 <- global_test(pops, n_permutations = 100, seed = 100)
  expect_identical(g1$null_distribution, g2$null_distribution)
  expect_false(identical(g1$null_distribution, g3$null_distribution))
})

test_that("print reports the comparison and returns its argument invisibly", {
  gt <- global_test(make_populations(), n_permutations = 200, seed = 42)

  out <- capture.output(printed <- withVisible(print(gt)))
  expect_false(printed$visible)
  expect_identical(printed$value, gt)

  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Populations:  2 \\(Control, Patient\\), 20 / 20 graphs")
  expect_match(txt, "Nodes:        10 \\(45 possible edges\\)")
  expect_match(txt, "T = -13.2", fixed = TRUE)
  # p = 0 out of 200 replicates is reported at the resolution of the estimate
  expect_match(txt, "p value:      < 0.005 (200 permutations)", fixed = TRUE)

  mid <- global_test(make_null_like(), n_permutations = 400, seed = 7)
  expect_match(paste(capture.output(print(mid)), collapse = "\n"),
               sprintf("p value:      %.4f (400 permutations)", mid$p_value),
               fixed = TRUE)
})

test_that("the resolution bound is rounded to two significant digits", {
  gt <- global_test(make_populations(), n_permutations = 3000, seed = 1)
  expect_identical(gt$p_value, 0)
  expect_match(paste(capture.output(print(gt)), collapse = "\n"),
               "< 0.00033", fixed = TRUE)
})

test_that("input validation names the problem", {
  pops <- make_populations()
  expect_error(global_test(list(A = pops$Control), 10), "at least 2 groups")
  expect_error(global_test(list(A = pops$Control[1], B = pops$Patient), 10),
               "at least two graphs")
  expect_error(global_test(pops, n_permutations = 0.5), "n_permutations")
  expect_error(global_test(pops, n_permutations = 0), "n_permutations")

  weighted <- lapply(pops$Control, function(A) A * 0.5)
  expect_error(global_test(list(A = weighted, B = pops$Patient), 10),
               "not binary")

  smaller <- list(matrix(0, 5, 5), matrix(0, 5, 5))
  expect_error(global_test(list(A = pops$Control, B = smaller), 10),
               "same vertex set")
})

test_that("a single-node network gives an NA p value and still prints", {
  one <- replicate(3, matrix(0, 1, 1), simplify = FALSE)
  gt <- global_test(list(A = one, B = one), n_permutations = 10)
  expect_identical(gt$statistic, 0)
  expect_true(is.na(gt$p_value))
  expect_identical(gt$null_distribution, rep(0, 10))
  expect_match(paste(capture.output(print(gt)), collapse = "\n"),
               "not evaluated")
})

test_that("more than two populations are handled", {
  pops <- make_populations()
  three <- list(A = pops$Control[1:10], B = pops$Control[11:20],
                C = pops$Patient[1:10])
  gt <- global_test(three, n_permutations = 200, seed = 5)
  expect_identical(gt$group_sizes, c(A = 10L, B = 10L, C = 10L))
  expect_equal(gt$null_distribution, direct_null(three, 200, seed = 5),
               tolerance = 1e-10)
  expect_match(paste(capture.output(print(gt)), collapse = "\n"),
               "3 (A, B, C), 10 / 10 / 10 graphs", fixed = TRUE)
})
