options(width = 90)

required <- c(
  "article.Rnw",
  "article.tex",
  "article.pdf",
  "figure/plot-1.pdf",
  "figure/abide-clinical-1.pdf",
  "refs.bib",
  "jss.cls",
  "jss.bst",
  "jsslogo.jpg",
  "Sweave.sty",
  "letter_editor.tex",
  "letter_editor.pdf",
  "code-full.R",
  "code.R",
  "code.html",
  "data-application.R",
  "README.md",
  "LICENSE.md",
  "SUBMISSION-FILES.txt",
  "package-source.txt",
  "renv.lock",
  "simulation-design.md",
  "statistical-contract.md",
  "response-to-editor.md",
  "reference/README.md",
  "application/abide_connectomes.rds",
  "application/abide_provenance.csv",
  "application/abide_subjects.csv",
  "generated/results.tex",
  "generated/manifest.csv",
  "generated/replication_metadata.csv",
  "generated/sessionInfo.txt",
  "generated/figures/complete_null_calibration.pdf",
  "generated/figures/localization_performance.pdf",
  "generated/figures/performance.pdf",
  "generated/figures/ablation_performance.pdf",
  "generated/tables/complete_null_summary.csv",
  "generated/tables/partial_null_summary.csv",
  "generated/tables/partial_null_conditional_summary.csv",
  "generated/tables/partial_null_gate.csv",
  "generated/tables/weak_null_summary.csv",
  "generated/tables/sensitivity_summary.csv",
  "generated/tables/benchmark_summary.csv",
  "generated/tables/ablation_benchmark_summary.csv"
)

missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop("Missing required replication files: ", paste(missing, collapse = ", "))
}

if (as.character(utils::packageVersion("BrainNetTest")) != "1.0.0") {
  stop("Verification requires the installed BrainNetTest 1.0.0 package.")
}

metadata <- read.csv(
  "generated/replication_metadata.csv",
  stringsAsFactors = FALSE
)
expected_hashes <- c(
  code_md5 = unname(tools::md5sum("code-full.R")),
  design_md5 = unname(tools::md5sum("simulation-design.md")),
  contract_md5 = unname(tools::md5sum("statistical-contract.md"))
)
recorded_hashes <- unlist(metadata[1L, names(expected_hashes)])
if (!identical(unname(recorded_hashes), unname(expected_hashes))) {
  stop("Replication source hashes differ from the generated metadata.")
}
if (!is.finite(metadata$elapsed_seconds) || metadata$elapsed_seconds <= 0) {
  stop("Replication runtime metadata is invalid.")
}

application_provenance <- read.csv(
  "application/abide_provenance.csv",
  stringsAsFactors = FALSE
)
recorded_connectomes_md5 <- application_provenance$connectomes_md5[1L]
current_connectomes_md5 <- unname(
  tools::md5sum("application/abide_connectomes.rds")
)
if (!identical(recorded_connectomes_md5, current_connectomes_md5)) {
  stop("Derived ABIDE connectomes differ from the recorded provenance hash.")
}
connectomes <- readRDS("application/abide_connectomes.rds")
if (!all(c("populations", "balanced_populations") %in% names(connectomes))) {
  stop("The ABIDE object must contain full and balanced populations.")
}
balanced_sizes <- lengths(connectomes$balanced_populations)
if (length(unique(balanced_sizes)) != 1L || any(balanced_sizes < 2L)) {
  stop("The ABIDE balanced populations must have equal valid sample sizes.")
}
application_subjects <- read.csv(
  "application/abide_subjects.csv",
  stringsAsFactors = FALSE
)
required_subject_columns <- c(
  "FILE_ID", "group", "AGE_AT_SCAN", "SEX", "func_mean_fd", "balanced"
)
if (!all(required_subject_columns %in% names(application_subjects))) {
  stop("The ABIDE subject table is missing required balance fields.")
}
balanced_subjects <- application_subjects[application_subjects$balanced, ]
balanced_group_sizes <- table(balanced_subjects$group)
balanced_sex <- table(balanced_subjects$group, balanced_subjects$SEX)
if (
  !identical(
    unname(as.integer(balanced_group_sizes[names(balanced_sizes)])),
    unname(as.integer(balanced_sizes))
  ) ||
    nrow(balanced_sex) != 2L ||
    !identical(
      unname(as.integer(balanced_sex[1L, ])),
      unname(as.integer(balanced_sex[2L, ]))
    )
) {
  stop("The ABIDE subject table does not match the balanced populations.")
}

source_lines <- readLines("article.Rnw", warn = FALSE)
graphic_lines <- grep(
  "\\{generated/figures/[^}]+\\}",
  source_lines,
  value = TRUE
)
graphic_paths <- sub(
  ".*\\{(generated/figures/[^}]+)\\}.*",
  "\\1",
  graphic_lines
)
missing_graphics <- graphic_paths[!file.exists(graphic_paths)]
if (length(missing_graphics) > 0L) {
  stop(
    "Manuscript references missing graphics: ",
    paste(missing_graphics, collapse = ", ")
  )
}

generated_files <- list.files(
  "generated",
  recursive = TRUE,
  full.names = TRUE
)
generated_files <- generated_files[file.info(generated_files)$isdir == FALSE]
generated_files <- setdiff(generated_files, "generated/manifest.csv")
current_manifest <- data.frame(
  file = sub("^generated/", "", generated_files),
  bytes = file.info(generated_files)$size,
  md5 = unname(tools::md5sum(generated_files)),
  stringsAsFactors = FALSE
)
current_manifest <- current_manifest[order(current_manifest$file), ]
rownames(current_manifest) <- NULL
supplied_manifest <- read.csv(
  "generated/manifest.csv",
  stringsAsFactors = FALSE
)
supplied_manifest <- supplied_manifest[order(supplied_manifest$file), ]
rownames(supplied_manifest) <- NULL
manifest_matches <- identical(
  current_manifest$file,
  supplied_manifest$file
) &&
  identical(
    as.numeric(current_manifest$bytes),
    as.numeric(supplied_manifest$bytes)
  ) &&
  identical(current_manifest$md5, supplied_manifest$md5)
if (!manifest_matches) {
  stop("Generated artifacts do not match `generated/manifest.csv`.")
}

if (file.info("article.pdf")$size <= 0L) {
  stop("The manuscript PDF is empty.")
}

message(
  "Verified ", length(required), " required files, ",
  length(graphic_paths), " manuscript figures, ",
  nrow(current_manifest), " generated artifacts, and the ABIDE ",
  "connectomes hash (", application_provenance$n_control[1L], " control, ",
  application_provenance$n_autism[1L], " autism, ",
  application_provenance$n_nodes[1L], " nodes; balanced sensitivity: ",
  balanced_sizes[[1L]], " per group)."
)
