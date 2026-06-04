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
  of Fraiman and Fraiman (2018) <doi:10.1038/s41598-018-21688-0>.
* Fast permutation procedure for identifying critical edges via a prefix-sum
  decomposition of the test statistic, reducing complexity from
  O(K * B * |E| * m) to O(B * |E| * m).
* Helpers to generate synthetic community-structured graphs and to visualise
  brain networks with communities.
