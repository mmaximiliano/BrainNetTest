script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- sub("^--file=", "", script_arg)
  setwd(dirname(normalizePath(script_path)))
}
if (!requireNamespace("BrainNetTest", quietly = TRUE) ||
  as.character(utils::packageVersion("BrainNetTest")) != "1.0.0") {
  stop("La verificación requiere BrainNetTest 1.0.0 instalado.")
}
suppressPackageStartupMessages(library("BrainNetTest"))

manifest_file <- "MANIFEST.csv"
if (!file.exists(manifest_file)) {
  stop("Falta ", manifest_file, ".")
}
manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
expected_columns <- c("file", "bytes", "md5")
if (!identical(names(manifest), expected_columns) ||
  anyDuplicated(manifest$file)) {
  stop("MANIFEST.csv tiene un schema inválido o paths duplicados.")
}

missing <- manifest$file[!file.exists(manifest$file)]
if (length(missing) > 0L) {
  stop("Faltan archivos del manifest: ", paste(missing, collapse = ", "))
}
observed_bytes <- file.info(manifest$file)$size
observed_md5 <- unname(tools::md5sum(manifest$file))
bad_size <- manifest$file[observed_bytes != manifest$bytes]
bad_md5 <- manifest$file[observed_md5 != manifest$md5]
if (length(bad_size) > 0L || length(bad_md5) > 0L) {
  stop(
    "Integrity check falló. Size mismatch: ",
    paste(bad_size, collapse = ", "),
    "; MD5 mismatch: ",
    paste(bad_md5, collapse = ", ")
  )
}

package_info <- readLines("PACKAGE_ARCHIVE.txt", warn = FALSE)
expected_package_md5 <- sub(
  "^MD5: ",
  "",
  grep("^MD5: ", package_info, value = TRUE)
)
observed_package_md5 <- unname(
  tools::md5sum("BrainNetTest_1.0.0.tar.gz")
)
if (!identical(observed_package_md5, expected_package_md5)) {
  stop("El package archive no coincide con PACKAGE_ARCHIVE.txt.")
}

balance <- utils::read.csv("results/balance_summary.csv")
smd <- utils::read.csv("results/balance_smd.csv")
pairs <- utils::read.csv("results/balance_pairs.csv")
retention <- utils::read.csv("results/balance_retention.csv")
global <- utils::read.csv("results/brainnet_global_results.csv")
edges <- utils::read.csv("results/brainnet_balanced_edges.csv")
top_edges <- utils::read.csv("results/brainnet_top_edges.csv")
holm <- utils::read.csv("results/brainnet_holm_selected_edges.csv")

if (nrow(balance) != 4L ||
  nrow(smd) != 2L ||
  nrow(pairs) != 55L ||
  nrow(retention) != 3L ||
  nrow(global) != 2L ||
  nrow(edges) != 6441L ||
  nrow(top_edges) != 10L ||
  nrow(holm) != 0L) {
  stop("Las dimensiones de los resultados no coinciden con lo esperado.")
}
if (!all(pairs$same_sex & pairs$same_age_bin & pairs$same_fd_bin) ||
  !setequal(global$sample, c("full", "balanced")) ||
  !all(global$n_edge_hypotheses == 6441L) ||
  !all(global$edges_per_network == 644L) ||
  !all(global$holm_selected == 0L) ||
  !all(global$by_selected == 0L)) {
  stop("Falló una verificación sustantiva de balance o inferencia.")
}

expected_t <- c(-44.68865, -27.39285)
expected_p <- c(0.001, 0.010)
if (!isTRUE(all.equal(global$statistic_T, expected_t, tolerance = 1e-5)) ||
  !isTRUE(all.equal(global$global_p_value, expected_p, tolerance = 1e-12))) {
  stop("Los resultados globales no coinciden con el paper.")
}

result_object <- readRDS("results/brainnet_results.rds")
if (!inherits(result_object$full, "brainnet_result") ||
  !inherits(result_object$balanced, "brainnet_result")) {
  stop("Los objetos brainnet_result del RDS son inválidos.")
}

for (sample_name in c("full", "balanced")) {
  result <- result_object[[sample_name]]
  row <- global[global$sample == sample_name, ]
  edge_table <- as.data.frame(result)
  by_selected <- sum(stats::p.adjust(edge_table$p_value, "BY") <= 0.05)
  if (!isTRUE(all.equal(
    unname(result$global$statistic),
    row$statistic_T,
    tolerance = 1e-8
  )) ||
    !isTRUE(all.equal(
      result$global$p_value,
      row$global_p_value,
      tolerance = 1e-12
    )) ||
    unname(result$group_sizes[["control"]]) != row$n_control ||
    unname(result$group_sizes[["autism"]]) != row$n_autism ||
    result$n_nodes != row$n_nodes ||
    result$edge_inference$family_size != row$n_edge_hypotheses ||
    result$global$n_assignments != row$random_assignments ||
    sum(edge_table$p_value < 0.001) != row$raw_p_below_0_001 ||
    sum(edge_table$selected) != row$holm_selected ||
    by_selected != row$by_selected) {
    stop("El RDS no coincide con el CSV global para ", sample_name, ".")
  }
}

rds_edges <- as.data.frame(result_object$balanced)
comparison_columns <- c(
  "node1", "node2", "risk_difference", "effect_range",
  "p_value", "p_adjusted", "selected", "rank"
)
if (!isTRUE(all.equal(
  rds_edges[, comparison_columns],
  edges[, comparison_columns],
  tolerance = 1e-12,
  check.attributes = FALSE
))) {
  stop("El edge table del RDS no coincide con el CSV determinístico.")
}

pdf_header <- rawToChar(readBin(
  "results/edge_effects.pdf",
  what = "raw",
  n = 4L
))
if (file.size("results/edge_effects.pdf") <= 0L ||
  !identical(pdf_header, "%PDF")) {
  stop("Los objetos RDS/PDF de resultados son inválidos.")
}

message(
  "Bundle verificado: ", nrow(manifest),
  " archivos con size+MD5; 55 pairs; 6.441 edge tests; 0 Holm/BY selections."
)
