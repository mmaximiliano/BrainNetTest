# BrainNetTest 1.0.0

## Brain-network workflow

* Refocused the package interface and documentation on independent
  populations of subject-level binary brain networks, while keeping the
  mathematical data contract usable for other aligned network populations.
* Added a reproducible ABIDE resting-state fMRI application to the paper
  materials, including a sex-, age-, and motion-balanced sensitivity sample.

## Inferential contract

* Replaced the strict-tail permutation estimate with an inclusive,
  finite-sample randomization test. Exhaustive tests include the observed
  assignment; Monte Carlo tests use the plus-one correction and can never
  return zero.
* Canonicalized L1 distance to the upper triangle so each undirected edge is
  counted once.
* Replaced the direction-dependent global rule with two-sided extremeness
  `abs(T)`, supporting unbalanced and multi-group designs without assuming
  that every alternative shifts `T` in the same direction.
* Separated global inference from edge-level inference. Every possible edge is
  tested with a two-sided Fisher exact test, with Holm adjustment by default.
* Demoted adaptively ranked edge removal to an optional descriptive ablation
  path. It no longer returns a first-non-significant or “critical” set.
* Retained the exact prefix-sum computation for that path and added a
  numerical-identity and runtime benchmark against repeated tail summation.

## Breaking API migration

* Construct input with `brainnet_data()`; raw nested lists are no longer
  accepted by the main analysis.
* Use `brainnet_test()` instead of `identify_critical_links()`.
* Use `selected_edges()` instead of `result$critical_edges`.
* Use `selected_nodes()` instead of `get_critical_nodes()`.
* Use `plot(result)` instead of `plot_critical_edges(populations, result)`.
* Use `central_graph()` and `graph_distance()` instead of
  `compute_central_graph()` and `compute_distance()`.
* `compute_test_statistic()`, `compute_edge_frequencies()`,
  `compute_edge_pvalues()`, and `rank_edges()` are no longer public. Their
  responsibilities are represented in the classed `brainnet_result`.
* The unused normalization argument `a` and misleading `adjust_method`
  ranking argument were removed.

## R classes and validation

* Added validated `brainnet_data` and `brainnet_result` S3 classes.
* Added `print()`, `summary()`, `plot()`, and `as.data.frame()` methods plus
  `selected_edges()`, `selected_nodes()`, and `ablation_path()` extractors.
* Added complete checks for group structure, sample sizes, dimensions,
  binary values, finite entries, symmetry, zero diagonals, and node ordering.
* Random seeds supplied to simulation or inference functions no longer alter
  the caller's random-number state.

# BrainNetTest 0.2.1

* `compute_edge_pvalues()` now clamps every returned p-value to the valid
  probability range `[0, 1]`. On platforms built without extended (long
  double) precision, exact tests such as `fisher.test()` can return a value
  fractionally greater than 1 due to floating-point rounding, which caused a
  test failure under CRAN's noLD check. (Reported by the CRAN team.)
* Added a regression test that mocks the underlying test to return an
  out-of-range p-value, so the clamping is verified on every platform rather
  than only on noLD builds.

# BrainNetTest 0.2.0

* Removed `plot_graph_with_communities()` and `plot_graphs_grid()` to
  streamline the API. The recommended plotting function is
  `plot_critical_edges()`, which produces a multi-panel visualisation of
  the analysis results.
* Removed unused Suggests: `ggplotify`, `gridExtra`.
* Removed unused imports: `grDevices::rainbow`, `graphics::legend`.

# BrainNetTest 0.1.0

* Initial CRAN submission.
* Implements the L1-distance ANOVA test for populations of brain networks
  of Fraiman and Fraiman (2018) <doi:10.1038/s41598-018-23152-5>.
* Fast permutation procedure for identifying critical edges via a prefix-sum
  decomposition of the test statistic, reducing complexity from
  O(K * B * |E| * m) to O(B * |E| * m).
* Helpers to generate synthetic community-structured graphs and to visualise
  brain networks with communities.
