script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- sub("^--file=", "", script_arg)
  setwd(dirname(normalizePath(script_path)))
}

application_dir <- Sys.getenv("ABIDE_APPLICATION_DIR", "application")
subjects_file <- file.path(application_dir, "abide_subjects.csv")
if (!file.exists(subjects_file)) {
  stop("No se encontró ", subjects_file, ". Correr 01_prepare_abide.R.")
}

subjects <- utils::read.csv(subjects_file, stringsAsFactors = FALSE)
required_columns <- c(
  "FILE_ID", "group", "AGE_AT_SCAN", "SEX", "func_mean_fd",
  "balanced", "balance_id", "balance_stratum"
)
if (!all(required_columns %in% names(subjects))) {
  stop("El subjects file no tiene todas las columnas esperadas.")
}
if (anyDuplicated(subjects$FILE_ID)) {
  stop("Hay FILE_ID duplicados.")
}
if (!setequal(unique(subjects$group), c("autism", "control"))) {
  stop("Se esperaban solamente los grupos autism y control.")
}

dir.create("results", showWarnings = FALSE)

summarize_sample <- function(data, sample_name) {
  rows <- lapply(c("autism", "control"), function(group_name) {
    group <- data[data$group == group_name, ]
    data.frame(
      sample = sample_name,
      group = group_name,
      n = nrow(group),
      male = sum(group$SEX == 1L),
      female = sum(group$SEX == 2L),
      male_percent = 100 * mean(group$SEX == 1L),
      age_mean = mean(group$AGE_AT_SCAN),
      age_sd = stats::sd(group$AGE_AT_SCAN),
      age_min = min(group$AGE_AT_SCAN),
      age_max = max(group$AGE_AT_SCAN),
      mean_fd = mean(group$func_mean_fd),
      mean_fd_sd = stats::sd(group$func_mean_fd),
      mean_fd_min = min(group$func_mean_fd),
      mean_fd_max = max(group$func_mean_fd)
    )
  })
  do.call(rbind, rows)
}

standardized_mean_difference <- function(value, group) {
  means <- tapply(value, group, mean)
  sds <- tapply(value, group, stats::sd)
  unname(
    (means[["autism"]] - means[["control"]]) /
      sqrt(mean(sds^2))
  )
}

balanced_subjects <- subjects[subjects$balanced, ]
summary_table <- rbind(
  summarize_sample(subjects, "full"),
  summarize_sample(balanced_subjects, "balanced")
)

smd_table <- data.frame(
  sample = c("full", "balanced"),
  age_smd_autism_minus_control = c(
    standardized_mean_difference(subjects$AGE_AT_SCAN, subjects$group),
    standardized_mean_difference(
      balanced_subjects$AGE_AT_SCAN,
      balanced_subjects$group
    )
  ),
  mean_fd_smd_autism_minus_control = c(
    standardized_mean_difference(
      subjects$func_mean_fd,
      subjects$group
    ),
    standardized_mean_difference(
      balanced_subjects$func_mean_fd,
      balanced_subjects$group
    )
  )
)

pair_ids <- sort(unique(balanced_subjects$balance_id))
pair_rows <- lapply(pair_ids, function(id) {
  pair <- balanced_subjects[balanced_subjects$balance_id == id, ]
  if (nrow(pair) != 2L ||
    !setequal(pair$group, c("autism", "control")) ||
    length(unique(pair$balance_stratum)) != 1L) {
    stop("Balance pair inválido: ", id)
  }
  control <- pair[pair$group == "control", ]
  autism <- pair[pair$group == "autism", ]
  data.frame(
    balance_id = id,
    balance_stratum = pair$balance_stratum[[1L]],
    control_FILE_ID = control$FILE_ID,
    autism_FILE_ID = autism$FILE_ID,
    same_sex = control$SEX == autism$SEX,
    same_age_bin = floor(control$AGE_AT_SCAN / 5) ==
      floor(autism$AGE_AT_SCAN / 5),
    same_fd_bin = floor(control$func_mean_fd / 0.05) ==
      floor(autism$func_mean_fd / 0.05)
  )
})
pair_table <- do.call(rbind, pair_rows)

if (!all(pair_table$same_sex) ||
  !all(pair_table$same_age_bin) ||
  !all(pair_table$same_fd_bin)) {
  stop("Algún par no respeta el coarsened exact matching.")
}
if (sum(balanced_subjects$group == "autism") !=
  sum(balanced_subjects$group == "control")) {
  stop("La sensitivity sample no quedó 1:1.")
}

expected <- utils::read.csv(
  "EXPECTED_BALANCE.csv",
  stringsAsFactors = FALSE
)
observed_for_check <- summary_table[, names(expected)]
numeric_columns <- setdiff(names(expected), c("sample", "group"))
same_keys <- identical(
  as.character(observed_for_check$sample),
  as.character(expected$sample)
) && identical(
  as.character(observed_for_check$group),
  as.character(expected$group)
)
same_values <- all(vapply(numeric_columns, function(column) {
  isTRUE(all.equal(
    observed_for_check[[column]],
    expected[[column]],
    tolerance = 1e-5
  ))
}, logical(1L)))
if (!same_keys || !same_values) {
  stop("El balance observado no coincide con EXPECTED_BALANCE.csv.")
}

utils::write.csv(
  summary_table,
  "results/balance_summary.csv",
  row.names = FALSE
)
utils::write.csv(
  smd_table,
  "results/balance_smd.csv",
  row.names = FALSE
)
utils::write.csv(
  pair_table,
  "results/balance_pairs.csv",
  row.names = FALSE
)

retention <- data.frame(
  group = c("autism", "control", "total"),
  original_n = c(
    sum(subjects$group == "autism"),
    sum(subjects$group == "control"),
    nrow(subjects)
  ),
  balanced_n = c(
    sum(balanced_subjects$group == "autism"),
    sum(balanced_subjects$group == "control"),
    nrow(balanced_subjects)
  )
)
retention$retained_percent <- 100 *
  retention$balanced_n / retention$original_n
utils::write.csv(
  retention,
  "results/balance_retention.csv",
  row.names = FALSE
)

diagnostics <- capture.output({
  cat("ABIDE balance diagnostics\n")
  cat("=========================\n\n")
  print(summary_table, row.names = FALSE)
  cat("\nStandardized mean differences (autism - control)\n")
  print(smd_table, row.names = FALSE)
  cat("\nRetention\n")
  print(retention, row.names = FALSE)
  cat("\nPairs verified:", nrow(pair_table), "\n")
  cat("Every pair has one subject per group: TRUE\n")
  cat("Every pair has exact sex stratum: TRUE\n")
  cat("Every pair shares 5-year age bin: TRUE\n")
  cat("Every pair shares 0.05 mean-FD bin: TRUE\n")
})
writeLines(diagnostics, "results/balance_diagnostics.txt")
cat(paste(diagnostics, collapse = "\n"), "\n")
