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

This update removes functions that were not yet in use by any downstream
package, simplifying the API before the package gains users.

## Changes in this version

* Removed two unexported-quality plotting functions
  (`plot_graph_with_communities()`, `plot_graphs_grid()`) that added
  external dependencies (`ggplotify`, `gridExtra`) without providing
  significant value over the primary `plot_critical_edges()` function.
* Removed corresponding Suggests, unused imports, tests, and documentation.

## Downstream dependencies

There are no downstream dependencies.
