## Test environments

* Local: Windows 11 x64, R 4.4.3
* win-builder (release and devel)
* R-hub: ubuntu-latest, windows-latest, macos-latest

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking for future file timestamps ... NOTE
  unable to verify current time

* checking CRAN incoming feasibility ... NOTE
  Days since last update: 2

Both NOTEs are environmental and not related to the package itself.

This is a patch release that fixes the issue reported by the CRAN team on the
noLD (no long double) check.

## Changes in this version

* Fixed the noLD check failure in `compute_edge_pvalues()`. On builds without
  extended (long double) precision, `fisher.test()` could return a p-value
  fractionally greater than 1 due to floating-point rounding, which tripped a
  `0 <= p <= 1` assertion in the package tests. P-values are now clamped to
  `[0, 1]`.
* Added a regression test (using `testthat::with_mocked_bindings()`) that
  reproduces the out-of-range p-value on any platform, guarding against a
  recurrence.

## Downstream dependencies

There are no downstream dependencies.
