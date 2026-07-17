# BrainNetTest paper replication

This directory contains the manuscript source and complete replication
materials for the brain-network-focused BrainNetTest 1.0.0 JSS resubmission.

## Install the software

From the repository root:

```sh
Rscript -e 'install.packages(c("igraph", "knitr", "rmarkdown", "testthat"))'
R CMD build .
R CMD INSTALL BrainNetTest_1.0.0.tar.gz
```

The replication scripts intentionally use the installed package. They do not
use `devtools::load_all()` or source files from `R/`.
`package-source.txt` records the SHA-256 digest of the exact source archive
that must accompany the JSS upload.
Build the upload from the explicit allowlist in `SUBMISSION-FILES.txt`; do not
zip this working directory indiscriminately.

External package versions are recorded in `renv.lock`. This is an optional
environment record rather than a hidden activation requirement. From
`paper_source/`, it can be restored explicitly with:

```sh
Rscript -e 'install.packages("renv"); renv::restore(project = ".", lockfile = "renv.lock", prompt = FALSE)'
```

## Full published results

Run from this directory:

```sh
Rscript code-full.R
```

`code-full.R` is the authoritative standalone script. It regenerates every
published simulation result, table, figure, inline LaTeX macro, benchmark, and
the session information under `generated/`.

The full study uses:

- 2,000 replicates per complete-null calibration cell;
- 1,000 replicates per partial-null/power cell;
- 2,000 weak-null stress replicates;
- 500 common-data sensitivity replicates; and
- 10 repetitions per network size for both computational benchmarks.

The frozen inferential design, pass/fail criteria, and separately disclosed
computational benchmarks are in `simulation-design.md`.
On the reference macOS system recorded in `generated/sessionInfo.txt`, the
entire script completes well within one hour.

## Reduced reviewer-time replication

```sh
Rscript code.R
Rscript -e 'knitr::spin("code.R")'
```

`code.R` is also standalone and produces `code.html` plus a small set of
outputs under `generated-quick/`. It verifies the package and simulation
pipeline within a short runtime. It does not reproduce the exact Monte Carlo
numbers in the manuscript; only `code-full.R` does that.

## Clinical brain-network application (ABIDE)

The manuscript applies the workflow to a real autism-control study of
resting-state functional connectomes.
The derived, ready-to-use binary networks are already provided under
`application/`, so the manuscript builds without any download. To regenerate
them from the public source (needs internet, about 35 MB on first run):

```sh
Rscript data-application.R
```

`data-application.R` downloads the openly shared, preprocessed ABIDE I
derivatives from the Preprocessed Connectomes Project S3 bucket, builds one
binary connectome per participant over the AAL atlas with a fixed 10% edge
density, validates them with `brainnet_data()`, and creates both the full NYU
cohort and a deterministic sex-, age-, and motion-balanced sensitivity sample.
It writes `application/abide_connectomes.rds` plus provenance. Raw downloads
are cached in `abide_cache/` (not part of the upload). The manuscript reads
only the derived `.rds`.

Data attribution and license: ABIDE I (Di Martino et al. 2014); preprocessed
derivatives from the Neuro Bureau / Preprocessed Connectomes Project (Craddock
et al. 2013); AAL atlas (Tzourio-Mazoyer et al. 2002). ABIDE derivatives are
shared for non-commercial research use under a Creative Commons
Attribution-NonCommercial-ShareAlike license, which the derived connectomes in
`application/` inherit.

## Manuscript

The authoritative manuscript source is `article.Rnw`. After full results have
been generated:

```sh
Rscript -e 'knitr::knit("article.Rnw")'
Rscript -e 'tinytex::latexmk("article.tex", clean = FALSE)'
Rscript verify.R
```

The manuscript reads generated numeric macros from `generated/results.tex`.
Missing generated results cause the build to fail.

Knitting also produces `figure/plot-1.pdf` from the live `plot()` example and
`figure/abide-clinical-1.pdf` from the clinical application. Both are
referenced by `article.tex`, are part of the manuscript upload, and are
regenerated on every knit.

A TeX distribution with `latexmk` is required. One user-scoped option is
TinyTeX:

```sh
Rscript -e 'install.packages("tinytex"); tinytex::install_tinytex()'
```

## Output map

- `generated/tables/complete_null_summary.csv`: global and complete-null edge
  calibration.
- `generated/tables/partial_null_summary.csv`: unconditional localization
  metrics.
- `generated/tables/partial_null_conditional_summary.csv`: separately labeled
  metrics conditional on global rejection.
- `generated/tables/partial_null_gate.csv`: prespecified partial-null error
  criteria.
- `generated/tables/weak_null_summary.csv`: equal-marginal/different-dependence
  scope stress test.
- `generated/tables/sensitivity_summary.csv`: permutation-count sensitivity.
- `generated/tables/benchmark_summary.csv`: measured runtime and input-object
  size.
- `generated/tables/ablation_benchmark_summary.csv`: measured repeated-tail
  and prefix-sum path-update runtime.
- `generated/tables/*_replicates.csv`: raw replicate-level results for complete
  nulls, partial nulls, weak-null stress, sensitivity, and benchmarks.
- `generated/figures/*.pdf`: manuscript figures.
- `generated/results.tex`: inline manuscript macros.
- `generated/sessionInfo.txt`: full execution environment.
- `generated/replication_metadata.csv`: runtime, platform, package version,
  seeds, and source/design hashes.
- `generated/manifest.csv`: immutable sizes and MD5 hashes verified by
  `verify.R`.
- `application/abide_connectomes.rds`: full and balanced derived binary ABIDE
  connectomes read by the manuscript's clinical application.
- `application/abide_provenance.csv`: dataset, atlas, threshold, sample sizes,
  and the connectomes MD5 verified by `verify.R`.
- `application/abide_subjects.csv`: anonymized participant identifiers, group,
  age, sex, motion, and balanced-sample membership.

## Determinism and platform notes

Both simulation scripts initialize L'Ecuyer-CMRG and use deterministic
replicate seeds.
Floating-point timings and PDF bytes can vary across platforms; numeric
statistical outputs are compared with documented tolerances rather than byte
identity. No raw images or direct identifiers are distributed; the application
contains anonymized ABIDE identifiers and derived binary connectomes under the
source data license.

Cross-platform verification uses these rules:

- integer counts, logical decisions, scenario labels, and table dimensions
  must match exactly;
- non-timing statistical columns must agree within absolute tolerance
  `1e-10`;
- timing and input-size columns must be finite, positive, and structurally
  complete but need not match the reference values; and
- figures must be non-empty and contain the expected panels, but PDF hashes
  are platform-specific.

## Historical reference

`reference/README.md` identifies the immutable version 0.2.1 tag used only for
historical comparison. Version 0.2.1 is not part of the 1.0 inferential API.
