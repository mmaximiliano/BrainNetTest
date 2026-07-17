# Point-by-point response to the editorial assessment

We appreciate the permission to resubmit. The revision is a breaking package
and manuscript redesign. Paths below are relative to the supplied repository
and paper bundle.

## Journal fit and software emphasis

**Comment:** The submission appeared primarily methodological, did not focus
enough on an open-source implementation for data analysis in general, did not
make suitable use of R classes and methods, and sketched the application
context narrowly.

**Response:** We rebalanced the manuscript around the brain-network problem,
the statistical method, its R implementation, and a real clinical application.
This follows the structure of accepted JSS papers that combine a domain
method, package interface, application, simulation, and computational
benchmark. Sections 2--5 present the clinical data contract, related
neuroimaging software, S3 classes, standard methods, statistical engine, and
executable workflows.

The package now has:

- a validated `brainnet_data` input class;
- a `brainnet_result` output class;
- `print()`, `summary()`, `plot()`, and `as.data.frame()` methods;
- `selected_edges()`, `selected_nodes()`, and `ablation_path()` extractors;
- complete validation of graph values, symmetry, diagonals, dimensions, group
  structure, sample sizes, and node ordering; and
- a frozen 10-function public API (`tests/testthat/test-api.R`).

Section 5 now exercises these methods explicitly: the two-group workflow prints
and summarizes a `brainnet_data` object and its `brainnet_result`, and the
three-group workflow shows the `plot()` method output in a figure. This makes
the class-and-method design visible in the manuscript rather than only in the
package sources.

The application-contract section explains how subject-level networks enter
the package and which upstream steps remain the analyst's responsibility.

To make the application context concrete, Section 6 is now a clinical
brain-network case study using openly shared, preprocessed ABIDE I
resting-state fMRI data. It builds one AAL network per participant and compares
autism and control populations. The observed global statistic is extreme under
label randomization in the full NYU cohort and remains so in a deterministic
sensitivity sample balanced on sex, age, and motion. No individual edge
survives strict multiplicity control, so the paper reports this
preprocessing-conditional population association without making a single-edge
biomarker, causal, or diagnostic claim. The application also reports the
observed runtime on 110 networks, 114 regions, 6,441 edge hypotheses, and 999
assignments.

The `data-application.R` script downloads the public derivatives, records all
preprocessing choices, and builds both cohorts. The manuscript reads only the
small derived object, so it compiles without network access.

## Inferential interpretation

Although not raised as a line-item code defect, our audit found that the
previous strict-tail Monte Carlo estimate and adaptively stopped edge set
could not support the manuscript's inferential language.

Version 1.0 therefore:

- counts each undirected edge once;
- uses inclusive two-sided randomization extremeness;
- applies the Monte Carlo plus-one correction;
- reports finite p-value resolution and ties;
- separates global whole-graph exchangeability inference from edge-marginal
  inference;
- tests the complete edge family with Holm adjustment by default; and
- retains prefix-sum ablation only as a descriptive diagnostic, with no
  first-non-rejecting set or equivalence/causal/minimality claim.

The computational contribution is retained rather than hidden. A new
benchmark verifies exact numerical identity and measures the prefix-sum update
against repeated tail summation, alongside the existing chunked-global-engine
benchmark. Thus the manuscript presents the optimization as a tested software
contribution while keeping its inferential limits explicit.

The exact contract is in `statistical-contract.md` and is enforced by
exact-enumeration and algebraic identity tests in
`tests/testthat/test-brainnet-inference.R`.

## Replication organization and runnable paths

**Comment:** JSS asks for one replication file, or a README when several files
are supplied. The submitted scripts failed because they sourced nonexistent
package paths and used `devtools::load_all("R-code")`.

**Response:** `paper_source/README.md` now provides the complete run order and
output map. `code-full.R` is the authoritative standalone replication script
and regenerates every reported simulation and benchmark result, table, figure,
and numeric LaTeX macro.
`code.R` is a separate standalone reduced study, and `code.html` is generated
with `knitr::spin()`. Neither script sources package internals from relative
`R/` paths or uses `devtools::load_all()`; both use the installed package.

The historical implementation is identified by immutable tag and commit in
`reference/README.md`, rather than by a fragile `source("R/...")` call.

## Manuscript code versus replication code

**Comment:** Code shown in the manuscript differed from `article.R`, and code
did not use JSS `R>` markup.

**Response:** `article.Rnw` is now the authoritative manuscript source and
uses `knitr::render_sweave()` with `R>` and `+` prompts. Displayed examples
call the same version 1.0 API exercised by the replication scripts. All inline
simulation claims are generated from `generated/results.tex` or hidden
CSV-reading chunks in `article.Rnw`; the manuscript build fails when the
required generated outputs are absent.

## Figure mismatch

**Comment:** Figure 1 from the replication materials differed from the
manuscript and lacked the inset.

**Response:** All old manually edited and duplicate figures were removed from
the submission workflow. Simulation and benchmark figures are generated
directly by `code-full.R` under `generated/figures/`. The workflow and ABIDE
clinical figures are produced by the manuscript's reproducible `knitr` chunks.
The ABIDE figure combines an anatomical edge-frequency difference map with
edge-wise evidence. No figure involves an undocumented editing or inset step,
so the manuscript and replication figures cannot diverge.

## Long code listings

**Comment:** Long unexplained listings can be omitted because they belong in
the replication material.

**Response:** The code-listing appendix was removed. The manuscript contains
only short excerpts needed to show the two-group, multi-group, and clinical
workflows. Complete code is in the standalone scripts, package tests, and
vignettes.

## R implementation details

**Comment:** Prefer `lengths(x)` to `sapply(x, length)`, consider `Reduce()`,
and improve input checks such as `generate_category_graphs(0.7)`.

**Response:** The implementation uses `lengths()` for group sizes and
`Reduce()` for central graphs. Shared scalar validators now reject fractional
counts, including `generate_category_graphs(0.7)`, with an argument-specific
error. Generator tests cover fractional node, graph, and community counts;
probability bounds; community sizes; reproducibility; and caller RNG-state
preservation.

## References

**Comment:** The manuscript contained hard-coded references without
bibliography entries.

**Response:** All prose citations now use BibTeX keys. The bibliography was
rebuilt from used references only. The Fraiman and Fraiman DOI was corrected
throughout the package and paper to
`10.1038/s41598-018-23152-5`.

## Verification

The package passes:

- formatting with `styler`;
- package lint with `lintr`;
- the complete `testthat` suite;
- vignette rebuilds;
- PDF manual and HTML manual validation; and
- local `R CMD check --as-cran` with 0 errors, 0 warnings, and 0 notes.

Automatic R-release checks are configured for Linux, macOS, and Windows.
The statistical study has prespecified complete- and partial-null gates in
`simulation-design.md`; failed gates stop `code-full.R` before publication
artifacts are accepted.
