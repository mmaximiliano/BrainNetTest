# Preparing a source release for review

This process prepares an R source package and its validation evidence. It does
not create a Git tag, publish a GitHub release, or submit anything to CRAN.
Those are separate decisions after the artifacts have been reviewed.

Choose the exact approved source commit, including all changes intended for the
release. Use its full 40-character lowercase SHA, not a branch name or a tag.
The preparation tool archives that commit, so unrelated local edits are excluded.
Keep the tool and the chosen source commit available in the local repository.

## Manual preparation

Install R, Git, Pandoc, `rcmdcheck`, and all dependencies listed in the chosen
commit's DESCRIPTION, including Suggests. A SHA-256 utility (`sha256sum` or
`shasum`) must also be available. From a checkout containing the preparation tool:

```sh
commit=$(git rev-parse HEAD)
bash .github/scripts/prepare-release.sh "$commit" /tmp/BrainNetTest-release
```

The output directory must be empty. To prepare an earlier approved commit,
replace the value of `commit` with its full SHA. Inspect that commit's DESCRIPTION
when installing dependencies. The tool builds in a temporary extracted snapshot
and checks the resulting tarball using `R CMD check --as-cran --no-manual` with
`NOT_CRAN=true` by default (matching the GitHub Actions workflow so tests gated
by `skip_on_cran()` still run; set `NOT_CRAN=false` explicitly to emulate CRAN's
skipped-test mode, and check `build-metadata.txt` for the recorded value).
Vignettes are built normally. PDF manual generation is excluded; perform that
additional check with a LaTeX installation before a CRAN submission.

An incomplete check, nonzero exit status, errors, or warnings fail preparation
with a diagnostic message pointing to the relevant log. Notes are printed to
stderr, recorded in `check-summary.txt` (and `check-notes.txt` when present),
and must be reviewed.
If a step fails, inspect `build.log`, `environment.log`, `check-console.log`, and
the `check/` directory as available. Use a new empty output directory to retry.

## GitHub Actions

Run **Prepare source release artifacts** manually from the Actions tab and enter
the approved source SHA. The workflow uses release R on Linux and records the
resolved environment. The preparation tool is taken from the workflow revision;
the package is taken from the requested source commit, which may predate the
tool. Both revisions are recorded. The workflow retains artifacts for 30 days,
including available logs when a build or check fails.

## Evidence to review and retain

- The `BrainNetTest_<version>.tar.gz` source archive and `SHA256SUMS`.
- `source-commit.txt`, `source-DESCRIPTION`, and `build-metadata.txt` identify
  the source, preparation tool, commands, and build time.
- `prepare-release.sh` is an exact copy of the tool used for this run.
- `session-info.txt`, `check-session-info.txt`, and `installed-packages.csv`
  record the R/platform information and installed dependency versions.
- `build.log`, `check-console.log`, `check-summary.txt`, and `check/` contain
  the build and test evidence.

Verify the downloaded archive from its output directory:

```sh
sha256sum -c SHA256SUMS       # Linux
shasum -a 256 -c SHA256SUMS # macOS
```

Confirm that the source SHA and version are the intended ones, the check reports
zero errors and warnings, and every note has been reviewed. Download and retain
the evidence before workflow artifacts expire.

These records identify the inputs used in a particular build; they are not a
dependency lockfile. Rebuilding with a different R version, platform, dependency
version, or build time may produce different bytes. The recorded SHA-256 checks
the integrity of the supplied archive, not a promise of identical future builds.
