# Independent small-sample oracles. These intentionally avoid the package's
# central-graph, edge-count and permutation helpers when calculating references.
reference_graphs <- function(codes, nodes = 4L) {
  lapply(codes, function(code) {
    graph <- matrix(0, nodes, nodes)
    graph[upper.tri(graph)] <- as.integer(intToBits(code))[
      seq_len(nodes * (nodes - 1L) / 2L)]
    graph + t(graph)
  })
}

reference_populations <- function(graphs, sizes) {
  split(graphs, rep(seq_along(sizes), sizes))
}

# For binary coordinates, distance to the sample mean equals the average
# pairwise distance. This formulation works directly on graph pairs, without
# using Bernoulli frequencies or the optimized edge-wise algebra.
reference_pairwise_T <- function(populations, a = 1) {
  graphs <- unlist(populations, recursive = FALSE)
  sizes <- lengths(populations)
  n <- length(graphs)
  distances <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    distances[i, j] <- sum(abs(graphs[[i]] - graphs[[j]]))
  }
  groups <- split(seq_len(n), rep(seq_along(sizes), sizes))
  terms <- vapply(groups, function(ids) {
    nk <- length(ids)
    sqrt(nk) * (sum(distances[ids, ids, drop = FALSE]) / (nk * (nk - 1)) -
                sum(distances[, ids, drop = FALSE]) / (nk * (n - 1)))
  }, numeric(1))
  unname(sqrt(length(sizes)) * sum(terms) / a)
}

test_that("observed and permuted statistics match a graph-pair oracle", {
  graphs <- reference_graphs(c(0, 1, 3, 7, 10, 17, 31, 47, 63))
  for (sizes in list(c(3L, 3L), c(3L, 3L, 3L))) {
    pops <- reference_populations(graphs[seq_len(sum(sizes))], sizes)
    ref <- reference_pairwise_T(pops)
    expect_equal(compute_test_statistic(pops), ref, tolerance = 1e-12)
    counts <- array(0, c(4, 4, length(sizes)))
    for (k in seq_along(pops)) for (graph in pops[[k]]) {
      counts[, , k] <- counts[, , k] + graph
    }
    contributions <- BrainNetTest:::.observed_edge_deltas(counts, sizes)
    expect_equal(2 * sum(contributions$deltas), ref, tolerance = 1e-12)
    expect_equal(compute_test_statistic(pops, a = 7), ref / 7,
                 tolerance = 1e-12)

    # Reconstruct the exact sampled label assignments without the count engine.
    seed <- 430L
    B <- 11L
    set.seed(seed)
    assignments <- replicate(B, sample.int(sum(sizes)), simplify = FALSE)
    pooled <- unlist(pops, recursive = FALSE)
    reference_null <- vapply(assignments, function(indices) {
      reference_pairwise_T(reference_populations(pooled[indices], sizes))
    }, numeric(1))
    got <- global_test(pops, n_permutations = B, seed = seed,
                       criterion = "original")
    expect_equal(got$statistic, ref, tolerance = 1e-12)
    expect_equal(got$null_distribution, reference_null, tolerance = 1e-12)
  }
})

test_that("edge counts preserve individual assignments and totals", {
  graphs <- reference_graphs(c(0, 1, 3, 7, 10, 17, 31, 47, 63))
  sizes <- c(2L, 3L, 4L)
  X <- t(vapply(graphs, function(graph) graph[upper.tri(graph)], numeric(6)))
  B <- 13L
  set.seed(312)
  assignments <- replicate(B, sample.int(nrow(X)), simplify = FALSE)
  set.seed(312)
  got <- BrainNetTest:::.permuted_edge_counts(X, sizes, B)
  positions <- split(seq_len(nrow(X)), rep(seq_along(sizes), sizes))
  expected <- lapply(positions, function(ids) {
    vapply(assignments, function(indices) {
      colSums(X[indices[ids], , drop = FALSE])
    }, numeric(ncol(X)))
  })
  expect_equal(unname(got), unname(expected), tolerance = 0)
  expect_equal(Reduce(`+`, got), matrix(colSums(X), ncol(X), B), tolerance = 0)

  # Unequal sizes are intentional here: this tests the low-level algebra,
  # independently of which designs the public inference API supports.
  deltas <- BrainNetTest:::.permuted_edge_deltas(got, sizes)
  reference <- vapply(assignments, function(indices) {
    reference_pairwise_T(reference_populations(graphs[indices], sizes))
  }, numeric(1))
  expect_equal(2 * colSums(deltas), reference, tolerance = 1e-12)
})

test_that("every small Fisher table agrees with exact combinatorial weights", {
  for (sizes in list(c(2L, 3L), c(4L, 4L), c(7L, 5L))) {
    n1 <- sizes[1]
    n2 <- sizes[2]
    total <- n1 + n2
    expected <- matrix(NA_real_, total + 1L, n1 + 1L)
    for (successes in 0:total) {
      feasible <- max(0L, successes - n2):min(n1, successes)
      # All weights are exactly represented small integers. No dhyper() or
      # fisher.test() is used, including for deciding probability ties.
      weights <- round(choose(n1, feasible) * choose(n2, successes - feasible))
      for (i in seq_along(feasible)) {
        expected[successes + 1L, feasible[i] + 1L] <-
          sum(weights[weights <= weights[i]]) / round(choose(total, successes))
      }
    }
    got <- BrainNetTest:::.fisher_pvalue_table(n1, n2)
    expect_equal(got, expected, tolerance = 1e-12)
  }
})

test_that("permutation p-values obey exhaustive finite-orbit level bounds", {
  # All 70 assignments of eight fixed subjects to two groups of four. The
  # two binary edges have integer-valued evidence scores: larger squared
  # differences are stronger. A second score isolates the strongest edge.
  X <- rbind(c(0, 0), c(0, 0), c(0, 1), c(0, 1),
             c(1, 0), c(1, 0), c(1, 0), c(1, 1))
  assignments <- combn(nrow(X), 4L, simplify = FALSE)
  scores <- vapply(assignments, function(ids) {
    delta <- colSums(X[ids, , drop = FALSE]) -
      colSums(X[-ids, , drop = FALSE])
    c(T = -sum(delta^2), U = -max(delta^2))
  }, numeric(2))
  N <- ncol(scores)
  # Exhaustive counting is the oracle, rather than a second call to rank().
  tails <- t(apply(scores, 1L, function(values) {
    vapply(values, function(value) sum(values <= value), integer(1))
  }))
  Q <- apply(tails, 2L, min)
  expected_combined <- vapply(Q, function(value) sum(Q <= value) / N, numeric(1))
  p_combined <- p_original <- numeric(N)
  for (i in seq_len(N)) {
    order_i <- c(i, setdiff(seq_len(N), i))
    p_combined[i] <- BrainNetTest:::.combined_pvalue(
      scores[1, order_i], scores[2, order_i])
    p_original[i] <- BrainNetTest:::.plus_one_pvalue(scores[1, order_i])
  }
  expect_equal(p_original, tails[1, ] / N, tolerance = 0)
  expect_equal(p_combined, expected_combined, tolerance = 0)
  for (p in list(p_original, p_combined)) {
    expect_true(all(p >= 1 / N & p <= 1))
    expect_equal(p * N, round(p * N), tolerance = 1e-12)
    # Verify every attainable significance level, not a noisy Monte Carlo
    # rejection rate with a permissive upper confidence limit.
    expect_true(all(vapply(seq_len(N), function(k) {
      sum(p <= k / N) <= k
    }, logical(1))))
  }
  expect_equal(BrainNetTest:::.combined_pvalue(rep(0, N), rep(1, N)), 1)
  expect_equal(BrainNetTest:::.plus_one_pvalue(rep(0, N)), 1)
})
