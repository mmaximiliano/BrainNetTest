#!/usr/bin/env bash
# Prepare reviewable source artifacts from a committed snapshot. No publishing.
set -euo pipefail

if [[ $# -ne 2 || ! "$1" =~ ^[0-9a-f]{40}$ ]]; then
  echo "Usage: $0 <full 40-character lowercase commit SHA> <empty output directory>" >&2
  exit 2
fi
commit=$1
output=$2
script_path=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
resolved=$(git rev-parse --verify "${commit}^{commit}")
if [[ "$resolved" != "$commit" ]]; then
  echo "The requested revision is not an exact commit SHA." >&2
  exit 2
fi
mkdir -p "$output"
output=$(cd "$output" && pwd)
if [[ -n "$(find "$output" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "Output directory must be empty: $output" >&2
  exit 2
fi
if command -v sha256sum >/dev/null 2>&1; then
  hash_command=(sha256sum)
elif command -v shasum >/dev/null 2>&1; then
  hash_command=(shasum -a 256)
else
  echo "Install sha256sum or shasum to record the source archive checksum." >&2
  exit 2
fi

scratch=$(mktemp -d "${TMPDIR:-/tmp}/brainnettest-release.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
mkdir "$scratch/source"
git archive "$commit" | tar -x -C "$scratch/source"
cp "$scratch/source/DESCRIPTION" "$output/source-DESCRIPTION"
printf '%s\n' "$commit" > "$output/source-commit.txt"
cp "$script_path" "$output/prepare-release.sh"
{
  printf 'Source commit: %s\n' "$commit"
  printf 'Preparation tool commit: %s\n' "${RELEASE_TOOL_COMMIT:-not recorded (local invocation)}"
  printf 'Preparation tool checksum: '
  "${hash_command[@]}" "$script_path"
  printf 'Built at (UTC): %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'Build: R CMD build --no-manual <committed source>\n'
  printf 'Check: R CMD check --as-cran --no-manual <source archive>\n'
  printf 'Gate: errors and warnings fail; notes require review.\n'
  printf 'PDF manual generation is excluded from this preparation check.\n'
} > "$output/build-metadata.txt"

Rscript --vanilla - "$output" <<'RSCRIPT' > "$output/environment.log" 2>&1
out <- commandArgs(trailingOnly = TRUE)[1L]
if (!requireNamespace("rcmdcheck", quietly = TRUE)) {
  stop("Install rcmdcheck and the source package's Depends/Imports/Suggests first.")
}
capture.output(sessionInfo(), file = file.path(out, "session-info.txt"))
installed <- installed.packages()
fields <- c("Package", "Version", "Built", "Priority")
write.csv(installed[order(installed[, "Package"]), fields, drop = FALSE],
          file.path(out, "installed-packages.csv"), row.names = FALSE)
RSCRIPT

cd "$output"
R CMD build --no-manual "$scratch/source" > build.log 2>&1
archive=$(Rscript --vanilla -e '
  d <- read.dcf("source-DESCRIPTION")
  cat(paste0(d[1, "Package"], "_", d[1, "Version"], ".tar.gz"))')
test -f "$archive"
"${hash_command[@]}" "$archive" > SHA256SUMS

# Check the built tarball, not the working tree. Keep the complete check folder
# even if the check reports problems, so an unsuccessful run remains reviewable.
Rscript --vanilla - "$archive" <<'RSCRIPT' > check-console.log 2>&1
archive <- commandArgs(trailingOnly = TRUE)[1L]
result <- rcmdcheck::rcmdcheck(archive,
  args = c("--as-cran", "--no-manual"), check_dir = "check", error_on = "never")
capture.output(sessionInfo(), file = "check-session-info.txt")
writeLines(c(paste("Process exit status:", result$status),
             paste("Errors:", length(result$errors)),
             paste("Warnings:", length(result$warnings)),
             paste("Notes:", length(result$notes))), "check-summary.txt")
logs <- list.files("check", pattern = "^00check[.]log$", recursive = TRUE,
                   full.names = TRUE)
completed <- length(logs) == 1L &&
  any(grepl("^\\* DONE$", readLines(logs, warn = FALSE)))
if (!identical(result$status, 0L) || !completed ||
    length(result$errors) || length(result$warnings)) {
  stop("Source check failed; inspect check-summary.txt and the check directory.")
}
RSCRIPT
printf 'Prepared source artifacts in %s\n' "$output"
