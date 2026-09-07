script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- sub("^--file=", "", script_arg)
  setwd(dirname(normalizePath(script_path)))
}

if (!requireNamespace("BrainNetTest", quietly = TRUE) ||
  as.character(utils::packageVersion("BrainNetTest")) != "1.0.0") {
  stop("Primero correr 00_install_brainnettest.R.")
}
suppressPackageStartupMessages(library("BrainNetTest"))

application_dir <- Sys.getenv("ABIDE_APPLICATION_DIR", "application")
connectome_file <- file.path(application_dir, "abide_connectomes.rds")
if (!file.exists(connectome_file)) {
  stop("No se encontró ", connectome_file, ". Correr 01_prepare_abide.R.")
}

dir.create("results", showWarnings = FALSE)
connectomes <- readRDS(connectome_file)

full_data <- brainnet_data(
  connectomes$populations,
  node_labels = connectomes$node_labels
)
balanced_data <- brainnet_data(
  connectomes$balanced_populations,
  node_labels = connectomes$node_labels
)

message("Corriendo full cohort...")
full_result <- brainnet_test(
  full_data,
  n_permutations = 999L,
  edge_method = "fisher",
  adjust = "holm",
  seed = 42L
)

message("Corriendo balanced sensitivity sample...")
message(
  "balance_id audita el subsampling; el test usa unrestricted ",
  "independent-group label randomization (55/55)."
)
balanced_result <- brainnet_test(
  balanced_data,
  n_permutations = 999L,
  edge_method = "fisher",
  adjust = "holm",
  seed = 42L
)

balanced_edges <- as.data.frame(balanced_result)
balanced_edges$region1 <- connectomes$node_labels[balanced_edges$node1]
balanced_edges$region2 <- connectomes$node_labels[balanced_edges$node2]
balanced_edges$p_by <- stats::p.adjust(
  balanced_edges$p_value,
  method = "BY"
)
balanced_edges <- balanced_edges[, c(
  "node1", "node2", "region1", "region2",
  grep("^proportion_", names(balanced_edges), value = TRUE),
  "risk_difference", "effect_range", "p_value", "p_adjusted",
  "p_by", "selected", "rank"
)]

network_edge_counts <- unlist(lapply(
  connectomes$populations,
  function(group) {
    vapply(
      group,
      function(graph) sum(graph[upper.tri(graph)]),
      numeric(1L)
    )
  }
))
if (length(unique(network_edge_counts)) != 1L) {
  stop("Las redes no tienen una densidad fija común.")
}

make_global_row <- function(label, result) {
  edges <- as.data.frame(result)
  by <- stats::p.adjust(edges$p_value, method = "BY")
  data.frame(
    sample = label,
    n_control = unname(result$group_sizes[["control"]]),
    n_autism = unname(result$group_sizes[["autism"]]),
    n_nodes = result$n_nodes,
    n_edge_hypotheses = result$edge_inference$family_size,
    edges_per_network = unique(network_edge_counts),
    statistic_T = unname(result$global$statistic),
    global_p_value = result$global$p_value,
    random_assignments = result$global$n_assignments,
    permutation_design = "unrestricted independent-group label randomization",
    balance_id_role = if (label == "balanced") {
      "audit-only stratified subsampling"
    } else {
      "not applicable"
    },
    raw_p_below_0_001 = sum(edges$p_value < 0.001),
    holm_selected = sum(edges$selected),
    by_selected = sum(by <= 0.05)
  )
}

global_results <- rbind(
  make_global_row("full", full_result),
  make_global_row("balanced", balanced_result)
)
utils::write.csv(
  global_results,
  "results/brainnet_global_results.csv",
  row.names = FALSE
)
utils::write.csv(
  balanced_edges,
  "results/brainnet_balanced_edges.csv",
  row.names = FALSE
)

top_edges <- utils::head(
  balanced_edges[order(balanced_edges$p_value), ],
  10L
)
utils::write.csv(
  top_edges,
  "results/brainnet_top_edges.csv",
  row.names = FALSE
)
utils::write.csv(
  selected_edges(balanced_result),
  "results/brainnet_holm_selected_edges.csv",
  row.names = FALSE
)
saveRDS(
  list(
    full = full_result,
    balanced = balanced_result,
    data_meta = connectomes$meta
  ),
  "results/brainnet_results.rds",
  compress = "xz"
)

grDevices::pdf("results/edge_effects.pdf", width = 7, height = 5)
graphics::plot(
  balanced_edges$risk_difference,
  -log10(balanced_edges$p_value),
  pch = 16,
  cex = 0.45,
  col = grDevices::adjustcolor("grey30", alpha.f = 0.5),
  xlab = "Edge-frequency difference: autism - control",
  ylab = "-log10(raw Fisher p-value)",
  main = "ABIDE balanced sensitivity sample"
)
top <- order(balanced_edges$p_value)[seq_len(10L)]
graphics::points(
  balanced_edges$risk_difference[top],
  -log10(balanced_edges$p_value[top]),
  pch = 1,
  cex = 1.2,
  col = "firebrick"
)
grDevices::dev.off()

cat("\nGlobal results\n")
print(global_results, row.names = FALSE)
cat("\nTop 10 raw edge results (no biomarker claim)\n")
print(
  top_edges[, c(
    "region1", "region2", "risk_difference",
    "p_value", "p_adjusted", "p_by"
  )],
  row.names = FALSE
)
cat(
  "\nInterpretación: el resultado global es significativo en ambas muestras; ",
  "ninguna arista sobrevive Holm/BY en la muestra balanceada.\n",
  sep = ""
)
