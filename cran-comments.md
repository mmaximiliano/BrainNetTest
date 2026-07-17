## Test environments

* Local: macOS 26.5.2 arm64, R 4.6.1
* GitHub Actions: R-devel, release, and oldrel on Ubuntu; R release on macOS
  and Windows
* R-hub: Linux, macOS, and Windows before submission

## R CMD check results

Local:

* 0 errors
* 0 warnings
* 0 notes

## Version 1.0.0

This release intentionally replaces the original unclassed API with validated
S3 input and result objects.

The statistical changes are substantive:

* permutation p-values now include ties and use the Monte Carlo plus-one
  correction;
* undirected distances count each edge once;
* the global statistic uses two-sided randomization extremeness;
* edge-level tests are separate from the global result and use Holm
  multiplicity adjustment by default; and
* the previous adaptively stopped “critical edge” result is replaced by an
  explicitly descriptive ablation path.

The package documentation includes a complete migration map in `NEWS.md`.

The Fraiman and Fraiman (2018) DOI was corrected to
<doi:10.1038/s41598-018-23152-5>.

## Downstream dependencies

There are no downstream dependencies.
