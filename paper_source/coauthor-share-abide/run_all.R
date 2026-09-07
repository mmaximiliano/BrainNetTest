script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- sub("^--file=", "", script_arg)
  setwd(dirname(normalizePath(script_path)))
}

steps <- c(
  "00_install_brainnettest.R",
  "01_prepare_abide.R",
  "02_check_balance.R",
  "03_run_analysis.R",
  "99_verify_bundle.R"
)

for (step in steps) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat("RUNNING: ", step, "\n", sep = "")
  cat(strrep("=", 72), "\n", sep = "")
  sys.source(step, envir = globalenv())
}

cat("\nProceso completo. Revisar application/ y results/.\n")
