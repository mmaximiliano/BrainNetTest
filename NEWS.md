# BrainNetTest 0.2.2

* `identify_critical_links()` now returns an object of class
  `"critical_links"` with `print()`, `summary()` and `plot()` methods, so the
  result reports itself instead of having to be picked apart by hand.
  `print()` gives the populations compared, the global-test result and the
  most significant critical edges; `summary()` adds the call, the settings the
  analysis ran with, and the node-level ranking from `get_critical_nodes()`;
  `plot()` is a method interface to `plot_critical_edges()` and takes the
  original populations as its second argument.
* The returned object gained three components alongside the existing three:
  `p_value` (the permutation p-value of the global test on the unmodified
  data, which was previously computed and discarded), `n_edges` and
  `settings`, plus the matched `call`. `critical_edges`, `edges_removed` and
  `modified_populations` keep their names, positions and contents, so existing
  code that extracts them is unaffected.
* Corrected the DOI given for Fraiman and Fraiman (2018) throughout the
  package. The previous value, `10.1038/s41598-018-21688-0`, resolves to an
  unrelated article; the correct one is `10.1038/s41598-018-23152-5`.
* Count arguments are validated. `generate_category_graphs(0.7)` used to fail
  with `attempt to select less than one element in integerOneIndex`,
  `generate_random_graph(0.7)` and `generate_community_graph(n_nodes = 0.7)`
  returned a degenerate 0 x 0 matrix, and a fractional `n_graphs` was silently
  truncated. These now raise an error that names the offending argument. The
  same applies to `batch_size`, `n_permutations`, `alpha` and `a` in
  `identify_critical_links()`.
* `identify_critical_links()` now checks that `populations` is a list of
  non-empty lists of equally sized square binary matrices, as its
  documentation has always required. Weighted matrices were previously passed
  to `fisher.test()`, which rounded the resulting counts and returned a
  meaningless ranking.
* `compute_edge_pvalues()` no longer subscripts out of bounds for single-node
  networks.
* `compute_test_statistic()` returns an unnamed scalar; it used to carry the
  name of the first population. `rank_edges()` resets the row names of its
  output, so the ranking reads 1, 2, 3, ... rather than the positions the
  edges occupied before sorting.
* Internal simplifications following review feedback: `lengths()` in place of
  `sapply(x, length)`, and `Reduce("+", x) / length(x)` in place of an
  accumulation loop in `compute_central_graph()`.
* `.Rbuildignore` now excludes the manuscript sources, so they are no longer
  included in the source tarball.

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
