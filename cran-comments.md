## Test environments

* Local: macOS 26.6.2, aarch64-apple-darwin, R 4.6.1
* win-builder (release and devel)
* R-hub: ubuntu-latest, windows-latest, macos-latest

## R CMD check results

0 errors | 0 warnings | 1 note

* checking HTML version of manual ... NOTE
  Skipping checking math rendering: package 'V8' unavailable

The NOTE is environmental: the check machine has no V8 installation, so the
math-rendering check is skipped. It is unrelated to the package.

## Changes in this version

This release responds to reviewer feedback on a manuscript describing the
package, and corrects a reference error.

* New exported function `global_test()`: the permutation test for a
  difference between populations on its own, returning the statistic, the
  p-value and the null distribution as a classed object with a `print()`
  method. It was previously available only as the first step of
  `identify_critical_links()`.
* `identify_critical_links()` returns an object of class `"critical_links"`
  with `print()`, `summary()` and `plot()` methods. The three components it
  had before keep their names and positions.
* The DOI cited for Fraiman and Fraiman (2018) was wrong in DESCRIPTION, the
  README, NEWS and both vignettes. It resolved to an unrelated article. The
  correct DOI is `10.1038/s41598-018-23152-5`.
* Count arguments are now validated, so mistakes such as
  `generate_category_graphs(0.7)` report which argument is wrong instead of
  failing inside the numerical code or returning a degenerate 0 x 0 matrix.
  Fractional counts and non-binary adjacency matrices are rejected rather than
  silently truncated or rounded.
* `compute_edge_pvalues()` no longer subscripts out of bounds for single-node
  networks.
* `compute_test_statistic()` returns an unnamed scalar and `rank_edges()`
  resets the row names of its result. Neither attribute was documented.
* `.Rbuildignore` now excludes the manuscript sources; the 0.2.1 tarball
  included them by accident, which is why it was substantially larger.

## Downstream dependencies

There are no downstream dependencies.
