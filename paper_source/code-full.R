###########################################################################
## Full replication for:
## BrainNetTest: Global and Edge-Wise Inference for Binary Network
## Populations
##
## Run from this directory after installing BrainNetTest 1.0.0.
## This script regenerates every simulation result, table, figure, inline
## LaTeX macro, and the session information used by the manuscript.
###########################################################################

options(width = 90, useFancyQuotes = FALSE)
RNGkind("L'Ecuyer-CMRG")

suppressPackageStartupMessages(library("BrainNetTest"))

if (as.character(utils::packageVersion("BrainNetTest")) != "1.0.0") {
  stop("Full replication requires the installed BrainNetTest 1.0.0 package.")
}

replication_started_at <- Sys.time()
replication_started_elapsed <- proc.time()[["elapsed"]]
alpha <- 0.05
n_permutations <- 199L
master_seed <- 20260714L
n_complete <- 2000L
n_partial <- 1000L
n_weak <- 2000L
n_sensitivity <- 500L

final_generated_dir <- "generated"
generated_dir <- tempfile(
  pattern = "brainnettest-generated-",
  tmpdir = "."
)
on.exit({
  if (dir.exists(generated_dir)) {
    unlink(generated_dir, recursive = TRUE, force = TRUE)
  }
}, add = TRUE)
figures_dir <- file.path(generated_dir, "figures")
tables_dir <- file.path(generated_dir, "tables")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

edge_index <- function(n_nodes) {
  which(upper.tri(matrix(FALSE, n_nodes, n_nodes)), arr.ind = TRUE)
}

graph_from_edges <- function(n_nodes, edge_values) {
  graph <- matrix(0L, nrow = n_nodes, ncol = n_nodes)
  graph[upper.tri(graph)] <- edge_values
  graph + t(graph)
}

draw_independent_graphs <- function(n_graphs, n_nodes, probabilities) {
  n_edges <- n_nodes * (n_nodes - 1L) / 2L
  if (length(probabilities) == 1L) {
    probabilities <- rep(probabilities, n_edges)
  }
  replicate(
    n_graphs,
    graph_from_edges(
      n_nodes,
      stats::rbinom(n_edges, 1L, probabilities)
    ),
    simplify = FALSE
  )
}

draw_latent_graphs <- function(n_graphs, n_nodes, mean_density = 0.3,
                               concentration = 20,
                               override_edges = integer(),
                               override_probability = NULL) {
  n_edges <- n_nodes * (n_nodes - 1L) / 2L
  shape_1 <- mean_density * concentration
  shape_2 <- (1 - mean_density) * concentration

  lapply(seq_len(n_graphs), function(index) {
    probability <- rep(
      stats::rbeta(1L, shape_1, shape_2),
      n_edges
    )
    if (length(override_edges) > 0L) {
      probability[override_edges] <- override_probability
    }
    graph_from_edges(
      n_nodes,
      stats::rbinom(n_edges, 1L, probability)
    )
  })
}

draw_fixed_density_graphs <- function(n_graphs, n_nodes, n_present) {
  n_edges <- n_nodes * (n_nodes - 1L) / 2L
  replicate(
    n_graphs,
    {
      values <- integer(n_edges)
      values[sample.int(n_edges, n_present)] <- 1L
      graph_from_edges(n_nodes, values)
    },
    simplify = FALSE
  )
}

fit_groups <- function(populations, seed, permutations = n_permutations) {
  data <- brainnet_data(populations)
  brainnet_test(
    data,
    alpha = alpha,
    n_permutations = permutations,
    exact = "never",
    edge_method = "fisher",
    adjust = "holm",
    seed = seed
  )
}

wilson_upper <- function(successes, trials, level = 0.95) {
  z <- stats::qnorm(level)
  estimate <- successes / trials
  denominator <- 1 + z^2 / trials
  center <- estimate + z^2 / (2 * trials)
  radius <- z * sqrt(
    estimate * (1 - estimate) / trials +
      z^2 / (4 * trials^2)
  )
  (center + radius) / denominator
}

bootstrap_mean_upper <- function(values, seed, n_boot = 2000L,
                                 level = 0.95) {
  set.seed(seed)
  means <- replicate(
    n_boot,
    mean(sample(values, replace = TRUE))
  )
  mean(values) + stats::qnorm(level) * stats::sd(means)
}

selection_metrics <- function(selected, truth) {
  true_positive <- sum(selected & truth)
  false_positive <- sum(selected & !truth)
  n_selected <- sum(selected)
  union_size <- sum(selected | truth)

  c(
    tpr = true_positive / sum(truth),
    precision = true_positive / max(n_selected, 1L),
    fdp = false_positive / max(n_selected, 1L),
    fwer = as.numeric(false_positive > 0L),
    set_size_error = n_selected - sum(truth),
    jaccard = if (union_size == 0L) 1 else true_positive / union_size,
    n_selected = n_selected
  )
}

summarize_binary <- function(values) {
  c(
    estimate = mean(values),
    lower = stats::qbinom(0.025, length(values), mean(values)) /
      length(values),
    upper = stats::qbinom(0.975, length(values), mean(values)) /
      length(values)
  )
}

scenario_seed <- function(scenario_index, replicate_index) {
  master_seed + scenario_index * 100000L + replicate_index
}

###########################################################################
## Complete-null calibration
###########################################################################

complete_scenarios <- list(
  independent_balanced = function() {
    list(
      GroupA = draw_independent_graphs(10L, 6L, 0.30),
      GroupB = draw_independent_graphs(10L, 6L, 0.30)
    )
  },
  independent_unbalanced = function() {
    list(
      GroupA = draw_independent_graphs(7L, 6L, 0.30),
      GroupB = draw_independent_graphs(13L, 6L, 0.30)
    )
  },
  latent_density = function() {
    list(
      GroupA = draw_latent_graphs(10L, 6L),
      GroupB = draw_latent_graphs(10L, 6L)
    )
  },
  fixed_density = function() {
    list(
      GroupA = draw_fixed_density_graphs(10L, 6L, 5L),
      GroupB = draw_fixed_density_graphs(10L, 6L, 5L)
    )
  },
  sparse_three_group = function() {
    list(
      GroupA = draw_independent_graphs(8L, 6L, 0.05),
      GroupB = draw_independent_graphs(8L, 6L, 0.05),
      GroupC = draw_independent_graphs(8L, 6L, 0.05)
    )
  },
  dense_three_group = function() {
    list(
      GroupA = draw_independent_graphs(8L, 6L, 0.95),
      GroupB = draw_independent_graphs(8L, 6L, 0.95),
      GroupC = draw_independent_graphs(8L, 6L, 0.95)
    )
  }
)

complete_rows <- vector(
  "list",
  length(complete_scenarios) * n_complete
)
row_index <- 0L

for (scenario_index in seq_along(complete_scenarios)) {
  scenario_name <- names(complete_scenarios)[[scenario_index]]
  message("Complete null: ", scenario_name)

  for (replicate_index in seq_len(n_complete)) {
    seed <- scenario_seed(scenario_index, replicate_index)
    set.seed(seed)
    result <- fit_groups(
      complete_scenarios[[scenario_index]](),
      seed = seed + 50000L
    )
    by_selected <- stats::p.adjust(
      result$edges$p_value,
      method = "BY"
    ) <= alpha

    row_index <- row_index + 1L
    complete_rows[[row_index]] <- data.frame(
      scenario = scenario_name,
      replicate = replicate_index,
      global_reject = result$global$p_value <= alpha,
      holm_fwer = any(result$edges$selected),
      by_fdp = as.numeric(any(by_selected))
    )
  }
}

complete_results <- do.call(rbind, complete_rows)
write.csv(
  complete_results,
  file.path(tables_dir, "complete_null_replicates.csv"),
  row.names = FALSE
)

complete_summary <- do.call(rbind, lapply(
  split(complete_results, complete_results$scenario),
  function(cell) {
    data.frame(
      scenario = cell$scenario[[1L]],
      n_replications = nrow(cell),
      global_rejection = mean(cell$global_reject),
      global_upper_95 = wilson_upper(
        sum(cell$global_reject),
        nrow(cell)
      ),
      holm_fwer = mean(cell$holm_fwer),
      holm_upper_95 = wilson_upper(sum(cell$holm_fwer), nrow(cell)),
      by_fdr = mean(cell$by_fdp)
    )
  }
))
rownames(complete_summary) <- NULL
write.csv(
  complete_summary,
  file.path(tables_dir, "complete_null_summary.csv"),
  row.names = FALSE
)

if (any(complete_summary$global_upper_95 > 0.065 |
  complete_summary$holm_upper_95 > 0.065)) {
  stop("A prespecified complete-null calibration criterion failed.")
}

###########################################################################
## Partial-null power and localization
###########################################################################

partial_scenarios <- list(
  localized_positive = function() {
    n_nodes <- 8L
    n_edges <- n_nodes * (n_nodes - 1L) / 2L
    truth <- seq_len(6L)
    probability_a <- rep(0.20, n_edges)
    probability_b <- probability_a
    probability_b[truth] <- 0.60
    list(
      populations = list(
        GroupA = draw_independent_graphs(20L, n_nodes, probability_a),
        GroupB = draw_independent_graphs(20L, n_nodes, probability_b)
      ),
      truth = seq_len(n_edges) %in% truth
    )
  },
  scattered_mixed = function() {
    n_nodes <- 8L
    n_edges <- n_nodes * (n_nodes - 1L) / 2L
    truth <- c(1L, 7L, 13L, 19L, 24L, 28L)
    probability_a <- rep(0.20, n_edges)
    probability_b <- probability_a
    probability_b[truth[seq_len(3L)]] <- 0.60
    probability_a[truth[4:6]] <- 0.60
    probability_b[truth[4:6]] <- 0.20
    list(
      populations = list(
        GroupA = draw_independent_graphs(20L, n_nodes, probability_a),
        GroupB = draw_independent_graphs(20L, n_nodes, probability_b)
      ),
      truth = seq_len(n_edges) %in% truth
    )
  },
  latent_dependent = function() {
    n_nodes <- 8L
    n_edges <- n_nodes * (n_nodes - 1L) / 2L
    truth <- seq_len(6L)
    list(
      populations = list(
        GroupA = draw_latent_graphs(
          20L,
          n_nodes,
          mean_density = 0.20,
          override_edges = truth,
          override_probability = 0.20
        ),
        GroupB = draw_latent_graphs(
          20L,
          n_nodes,
          mean_density = 0.20,
          override_edges = truth,
          override_probability = 0.60
        )
      ),
      truth = seq_len(n_edges) %in% truth
    )
  },
  three_group = function() {
    n_nodes <- 8L
    n_edges <- n_nodes * (n_nodes - 1L) / 2L
    truth <- seq_len(6L)
    probabilities <- lapply(c(0.20, 0.40, 0.65), function(value) {
      output <- rep(0.20, n_edges)
      output[truth] <- value
      output
    })
    list(
      populations = list(
        GroupA = draw_independent_graphs(15L, n_nodes, probabilities[[1L]]),
        GroupB = draw_independent_graphs(15L, n_nodes, probabilities[[2L]]),
        GroupC = draw_independent_graphs(15L, n_nodes, probabilities[[3L]])
      ),
      truth = seq_len(n_edges) %in% truth
    )
  }
)

partial_rows <- vector(
  "list",
  length(partial_scenarios) * n_partial * 3L
)
row_index <- 0L

for (scenario_index in seq_along(partial_scenarios)) {
  scenario_name <- names(partial_scenarios)[[scenario_index]]
  message("Partial null: ", scenario_name)

  for (replicate_index in seq_len(n_partial)) {
    seed <- scenario_seed(20L + scenario_index, replicate_index)
    set.seed(seed)
    generated <- partial_scenarios[[scenario_index]]()
    result <- fit_groups(
      generated$populations,
      seed = seed + 50000L
    )
    selections <- list(
      unadjusted = result$edges$p_value <= alpha,
      Holm = result$edges$p_adjusted <= alpha,
      BY = stats::p.adjust(result$edges$p_value, "BY") <= alpha
    )

    for (procedure in names(selections)) {
      metrics <- selection_metrics(
        selections[[procedure]],
        generated$truth
      )
      row_index <- row_index + 1L
      partial_rows[[row_index]] <- data.frame(
        scenario = scenario_name,
        replicate = replicate_index,
        procedure = procedure,
        global_reject = result$global$p_value <= alpha,
        t(metrics)
      )
    }
  }
}

partial_results <- do.call(rbind, partial_rows)
write.csv(
  partial_results,
  file.path(tables_dir, "partial_null_replicates.csv"),
  row.names = FALSE
)

summarize_metrics <- function(data) {
  data.frame(
    n_replications = nrow(data),
    global_power = mean(data$global_reject),
    tpr = mean(data$tpr),
    precision = mean(data$precision),
    fdp = mean(data$fdp),
    fwer = mean(data$fwer),
    mean_set_size_error = mean(data$set_size_error),
    jaccard = mean(data$jaccard),
    mean_selected = mean(data$n_selected)
  )
}

partial_summary <- do.call(rbind, lapply(
  split(
    partial_results,
    interaction(
      partial_results$scenario,
      partial_results$procedure,
      drop = TRUE
    )
  ),
  function(cell) {
    cbind(
      data.frame(
        scenario = cell$scenario[[1L]],
        procedure = cell$procedure[[1L]]
      ),
      summarize_metrics(cell)
    )
  }
))
rownames(partial_summary) <- NULL

partial_conditional <- do.call(rbind, lapply(
  split(
    partial_results[partial_results$global_reject, ],
    interaction(
      partial_results$scenario[partial_results$global_reject],
      partial_results$procedure[partial_results$global_reject],
      drop = TRUE
    )
  ),
  function(cell) {
    cbind(
      data.frame(
        scenario = cell$scenario[[1L]],
        procedure = cell$procedure[[1L]]
      ),
      summarize_metrics(cell)
    )
  }
))
rownames(partial_conditional) <- NULL

write.csv(
  partial_summary,
  file.path(tables_dir, "partial_null_summary.csv"),
  row.names = FALSE
)
write.csv(
  partial_conditional,
  file.path(tables_dir, "partial_null_conditional_summary.csv"),
  row.names = FALSE
)

holm_cells <- split(
  partial_results[partial_results$procedure == "Holm", ],
  partial_results$scenario[partial_results$procedure == "Holm"]
)
holm_gate <- vapply(holm_cells, function(cell) {
  wilson_upper(sum(cell$fwer), nrow(cell))
}, numeric(1L))

by_cells <- split(
  partial_results[partial_results$procedure == "BY", ],
  partial_results$scenario[partial_results$procedure == "BY"]
)
by_gate <- vapply(seq_along(by_cells), function(index) {
  bootstrap_mean_upper(
    by_cells[[index]]$fdp,
    seed = master_seed + 900000L + index
  )
}, numeric(1L))
names(by_gate) <- names(by_cells)

partial_gate <- data.frame(
  scenario = names(holm_gate),
  holm_fwer_upper_95 = unname(holm_gate),
  by_fdr_upper_95 = unname(by_gate[names(holm_gate)])
)
write.csv(
  partial_gate,
  file.path(tables_dir, "partial_null_gate.csv"),
  row.names = FALSE
)

if (any(partial_gate$holm_fwer_upper_95 > 0.065 |
  partial_gate$by_fdr_upper_95 > 0.065)) {
  stop("A prespecified partial-null edge criterion failed.")
}

###########################################################################
## Weak-null scope stress test
###########################################################################

weak_rows <- vector("list", n_weak)
message("Weak-null dependence stress test")
for (replicate_index in seq_len(n_weak)) {
  seed <- scenario_seed(40L, replicate_index)
  set.seed(seed)
  n_nodes <- 6L
  n_edges <- n_nodes * (n_nodes - 1L) / 2L
  group_a <- lapply(seq_len(10L), function(index) {
    graph_from_edges(n_nodes, rep(stats::rbinom(1L, 1L, 0.5), n_edges))
  })
  group_b <- draw_independent_graphs(10L, n_nodes, 0.5)
  result <- fit_groups(
    list(GroupA = group_a, GroupB = group_b),
    seed = seed + 50000L
  )
  weak_rows[[replicate_index]] <- data.frame(
    replicate = replicate_index,
    global_reject = result$global$p_value <= alpha,
    holm_fwer = any(result$edges$selected)
  )
}

weak_results <- do.call(rbind, weak_rows)
weak_summary <- data.frame(
  n_replications = nrow(weak_results),
  global_rejection = mean(weak_results$global_reject),
  holm_fwer = mean(weak_results$holm_fwer),
  holm_upper_95 = wilson_upper(
    sum(weak_results$holm_fwer),
    nrow(weak_results)
  )
)
write.csv(
  weak_results,
  file.path(tables_dir, "weak_null_replicates.csv"),
  row.names = FALSE
)
write.csv(
  weak_summary,
  file.path(tables_dir, "weak_null_summary.csv"),
  row.names = FALSE
)

if (weak_summary$holm_upper_95 > 0.065) {
  stop("Holm edge inference failed the weak-null stress criterion.")
}

###########################################################################
## Permutation-count sensitivity
###########################################################################

sensitivity_counts <- c(99L, 199L, 999L, 4999L)
sensitivity_rows <- vector(
  "list",
  n_sensitivity * length(sensitivity_counts)
)
row_index <- 0L
message("Permutation-count sensitivity")

for (replicate_index in seq_len(n_sensitivity)) {
  seed <- scenario_seed(50L, replicate_index)
  set.seed(seed)
  generated <- partial_scenarios$localized_positive()

  for (permutations in sensitivity_counts) {
    result <- fit_groups(
      generated$populations,
      seed = seed +
        match(permutations, sensitivity_counts) * 10000000L,
      permutations = permutations
    )
    row_index <- row_index + 1L
    sensitivity_rows[[row_index]] <- data.frame(
      replicate = replicate_index,
      n_permutations = permutations,
      global_p_value = result$global$p_value,
      global_reject = result$global$p_value <= alpha,
      n_selected = sum(result$edges$selected),
      elapsed_seconds = result$elapsed_seconds
    )
  }
}

sensitivity_results <- do.call(rbind, sensitivity_rows)
sensitivity_summary <- do.call(rbind, lapply(
  split(sensitivity_results, sensitivity_results$n_permutations),
  function(cell) {
    data.frame(
      n_permutations = cell$n_permutations[[1L]],
      global_power = mean(cell$global_reject),
      median_p_value = stats::median(cell$global_p_value),
      mean_selected = mean(cell$n_selected),
      median_seconds = stats::median(cell$elapsed_seconds)
    )
  }
))
rownames(sensitivity_summary) <- NULL
write.csv(
  sensitivity_results,
  file.path(tables_dir, "sensitivity_replicates.csv"),
  row.names = FALSE
)
write.csv(
  sensitivity_summary,
  file.path(tables_dir, "sensitivity_summary.csv"),
  row.names = FALSE
)

###########################################################################
## Performance benchmark
###########################################################################

benchmark_nodes <- c(16L, 32L, 64L, 128L, 256L)
benchmark_repetitions <- 10L
benchmark_rows <- list()
row_index <- 0L
message("Performance benchmark")

measure_per_call <- function(code, minimum_elapsed = 0.25,
                             maximum_repetitions = 500L) {
  code()
  pilot <- unname(system.time(code())[["elapsed"]])
  repetitions <- min(
    maximum_repetitions,
    max(1L, ceiling(minimum_elapsed / max(pilot, 0.001)))
  )
  elapsed <- unname(system.time({
    for (index in seq_len(repetitions)) {
      code()
    }
  })[["elapsed"]])
  list(
    seconds = elapsed / repetitions,
    repetitions = repetitions
  )
}

for (n_nodes in benchmark_nodes) {
  for (benchmark_index in seq_len(benchmark_repetitions)) {
    seed <- scenario_seed(60L + n_nodes, benchmark_index)
    set.seed(seed)
    populations <- list(
      GroupA = draw_independent_graphs(20L, n_nodes, 0.30),
      GroupB = draw_independent_graphs(20L, n_nodes, 0.45)
    )
    data <- brainnet_data(populations)
    encoded <- BrainNetTest:::.brainnet_encode(data)
    assignments <- BrainNetTest:::.with_preserved_seed(
      seed + 50000L,
      function() {
        BrainNetTest:::.sample_group_assignments(
          data$group_sizes,
          499L
        )
      }
    )

    optimized_code <- function() {
      BrainNetTest:::.brainnet_scores_for_assignments(
        encoded$edge_matrix,
        assignments,
        data$group_sizes
      )
    }
    direct_code <- function() {
      apply(assignments, 1L, function(labels) {
        counts <- BrainNetTest:::.brainnet_group_counts(
          encoded$edge_matrix,
          labels,
          n_groups = length(data$group_sizes)
        )
        BrainNetTest:::.brainnet_score(counts, data$group_sizes)
      })
    }

    optimized <- optimized_code()
    direct <- direct_code()
    if (!isTRUE(all.equal(optimized, direct, tolerance = 1e-10))) {
      stop("Optimized and direct assignment scores differ.")
    }

    method_order <- sample(c("optimized", "direct"))
    timings <- list()
    for (method in method_order) {
      timings[[method]] <- measure_per_call(
        if (method == "optimized") optimized_code else direct_code
      )
    }

    row_index <- row_index + 1L
    benchmark_rows[[row_index]] <- data.frame(
      n_nodes = n_nodes,
      n_edges = n_nodes * (n_nodes - 1L) / 2L,
      n_graphs = sum(data$group_sizes),
      n_permutations = 499L,
      optimized_seconds = timings$optimized$seconds,
      direct_seconds = timings$direct$seconds,
      speedup = timings$direct$seconds / timings$optimized$seconds,
      optimized_batch_repetitions = timings$optimized$repetitions,
      direct_batch_repetitions = timings$direct$repetitions,
      input_megabytes = as.numeric(
        object.size(encoded$edge_matrix) + object.size(assignments)
      ) / 1024^2
    )
  }
}

benchmark_results <- do.call(rbind, benchmark_rows)
benchmark_summary <- do.call(rbind, lapply(
  split(benchmark_results, benchmark_results$n_nodes),
  function(cell) {
    data.frame(
      n_nodes = cell$n_nodes[[1L]],
      n_edges = cell$n_edges[[1L]],
      median_optimized_seconds = stats::median(cell$optimized_seconds),
      optimized_iqr = stats::IQR(cell$optimized_seconds),
      median_direct_seconds = stats::median(cell$direct_seconds),
      median_speedup = stats::median(cell$speedup),
      input_megabytes = stats::median(cell$input_megabytes)
    )
  }
))
rownames(benchmark_summary) <- NULL
write.csv(
  benchmark_results,
  file.path(tables_dir, "benchmark_replicates.csv"),
  row.names = FALSE
)
write.csv(
  benchmark_summary,
  file.path(tables_dir, "benchmark_summary.csv"),
  row.names = FALSE
)

###########################################################################
## Prefix-sum ablation benchmark
###########################################################################

ablation_nodes <- c(16L, 24L, 32L, 40L, 48L, 56L, 64L)
ablation_assignments <- 99L
ablation_repetitions <- 10L
ablation_benchmark_rows <- list()
row_index <- 0L
message("Prefix-sum ablation benchmark")

for (n_nodes in ablation_nodes) {
  for (benchmark_index in seq_len(ablation_repetitions)) {
    seed <- scenario_seed(400L + n_nodes, benchmark_index)
    set.seed(seed)
    populations <- list(
      GroupA = draw_independent_graphs(20L, n_nodes, 0.30),
      GroupB = draw_independent_graphs(20L, n_nodes, 0.45)
    )
    data <- brainnet_data(populations)
    encoded <- BrainNetTest:::.brainnet_encode(data)
    assignments <- BrainNetTest:::.with_preserved_seed(
      seed + 50000L,
      function() {
        BrainNetTest:::.sample_group_assignments(
          data$group_sizes,
          ablation_assignments
        )
      }
    )
    n_edges <- ncol(encoded$edge_matrix)
    observed_counts <- BrainNetTest:::.brainnet_group_counts(
      encoded$edge_matrix,
      encoded$group,
      n_groups = length(data$group_sizes)
    )
    observed_contributions <- BrainNetTest:::.brainnet_edge_contributions(
      observed_counts,
      data$group_sizes
    )
    edge_order <- order(
      -abs(observed_contributions),
      seq_len(n_edges),
      method = "radix"
    )
    assignment_contributions <- t(vapply(
      seq_len(nrow(assignments)),
      function(assignment_index) {
        counts <- BrainNetTest:::.brainnet_group_counts(
          encoded$edge_matrix,
          assignments[assignment_index, ],
          n_groups = length(data$group_sizes)
        )
        BrainNetTest:::.brainnet_edge_contributions(
          counts,
          data$group_sizes
        )
      },
      numeric(n_edges)
    ))

    prefix_code <- function() {
      scores <- rowSums(assignment_contributions)
      t(vapply(
        seq_len(nrow(assignment_contributions)),
        function(assignment_index) {
          c(
            scores[[assignment_index]],
            scores[[assignment_index]] - cumsum(
              assignment_contributions[assignment_index, edge_order]
            )
          )
        },
        numeric(n_edges + 1L)
      ))
    }
    repeated_tail_code <- function() {
      t(vapply(
        seq_len(nrow(assignment_contributions)),
        function(assignment_index) {
          ordered <- assignment_contributions[
            assignment_index,
            edge_order
          ]
          c(
            sum(ordered),
            vapply(
              seq_len(n_edges),
              function(n_removed) {
                if (n_removed == n_edges) {
                  return(0)
                }
                sum(ordered[-seq_len(n_removed)])
              },
              numeric(1L)
            )
          )
        },
        numeric(n_edges + 1L)
      ))
    }

    prefix <- prefix_code()
    repeated <- repeated_tail_code()
    if (!isTRUE(all.equal(prefix, repeated, tolerance = 1e-10))) {
      stop("Prefix and repeated-tail ablation paths differ.")
    }

    method_order <- sample(c("prefix", "repeated"))
    timings <- list()
    for (method in method_order) {
      timings[[method]] <- measure_per_call(
        if (method == "prefix") prefix_code else repeated_tail_code
      )
    }

    row_index <- row_index + 1L
    ablation_benchmark_rows[[row_index]] <- data.frame(
      n_nodes = n_nodes,
      n_edges = n_edges,
      n_assignments = ablation_assignments,
      prefix_seconds = timings$prefix$seconds,
      repeated_tail_seconds = timings$repeated$seconds,
      speedup = timings$repeated$seconds / timings$prefix$seconds,
      prefix_batch_repetitions = timings$prefix$repetitions,
      repeated_batch_repetitions = timings$repeated$repetitions
    )
  }
}

ablation_benchmark_results <- do.call(rbind, ablation_benchmark_rows)
ablation_benchmark_summary <- do.call(rbind, lapply(
  split(ablation_benchmark_results, ablation_benchmark_results$n_nodes),
  function(cell) {
    data.frame(
      n_nodes = cell$n_nodes[[1L]],
      n_edges = cell$n_edges[[1L]],
      median_prefix_seconds = stats::median(cell$prefix_seconds),
      prefix_iqr = stats::IQR(cell$prefix_seconds),
      median_repeated_tail_seconds = stats::median(
        cell$repeated_tail_seconds
      ),
      median_speedup = stats::median(cell$speedup)
    )
  }
))
rownames(ablation_benchmark_summary) <- NULL
write.csv(
  ablation_benchmark_results,
  file.path(tables_dir, "ablation_benchmark_replicates.csv"),
  row.names = FALSE
)
write.csv(
  ablation_benchmark_summary,
  file.path(tables_dir, "ablation_benchmark_summary.csv"),
  row.names = FALSE
)

###########################################################################
## Figures
###########################################################################

pdf(file.path(figures_dir, "complete_null_calibration.pdf"), 9, 5)
op <- par(mar = c(8, 4, 2, 1))
positions <- seq_len(nrow(complete_summary))
plot(
  positions,
  complete_summary$global_rejection,
  ylim = c(0, max(0.08, complete_summary$global_upper_95)),
  xaxt = "n",
  xlab = "",
  ylab = "Rejection / family-wise error rate",
  pch = 19,
  col = "steelblue"
)
segments(
  positions,
  0,
  positions,
  complete_summary$global_upper_95,
  col = "steelblue"
)
points(
  positions + 0.12,
  complete_summary$holm_fwer,
  pch = 17,
  col = "firebrick"
)
segments(
  positions + 0.12,
  0,
  positions + 0.12,
  complete_summary$holm_upper_95,
  col = "firebrick"
)
abline(h = alpha, lty = 2)
axis(
  1,
  at = positions,
  labels = complete_summary$scenario,
  las = 2,
  cex.axis = 0.8
)
legend(
  "topright",
  legend = c("Global", "Holm FWER", "Nominal alpha"),
  pch = c(19, 17, NA),
  lty = c(NA, NA, 2),
  col = c("steelblue", "firebrick", "black"),
  bty = "n"
)
par(op)
dev.off()

pdf(file.path(figures_dir, "localization_performance.pdf"), 9, 5)
holm_summary <- partial_summary[partial_summary$procedure == "Holm", ]
by_summary <- partial_summary[partial_summary$procedure == "BY", ]
matrix_values <- rbind(
  Holm_TPR = holm_summary$tpr,
  Holm_precision = holm_summary$precision,
  BY_TPR = by_summary$tpr,
  BY_precision = by_summary$precision
)
barplot(
  matrix_values,
  beside = TRUE,
  names.arg = holm_summary$scenario,
  las = 2,
  ylim = c(0, 1),
  ylab = "Mean rate",
  legend.text = rownames(matrix_values),
  args.legend = list(x = "topright", bty = "n", cex = 0.8),
  col = c("steelblue", "lightblue", "firebrick", "mistyrose")
)
dev.off()

pdf(file.path(figures_dir, "performance.pdf"), 7, 5)
plot(
  benchmark_summary$n_edges,
  benchmark_summary$median_direct_seconds,
  type = "b",
  log = "xy",
  pch = 17,
  col = "firebrick",
  ylim = range(c(
    benchmark_summary$median_direct_seconds,
    benchmark_summary$median_optimized_seconds
  )),
  xlab = "Number of possible edges",
  ylab = "Median elapsed seconds"
)
lines(
  benchmark_summary$n_edges,
  benchmark_summary$median_optimized_seconds,
  type = "b",
  pch = 19,
  col = "steelblue"
)
legend(
  "topleft",
  legend = c("Direct assignment loop", "Chunked matrix engine"),
  pch = c(17, 19),
  lty = 1,
  col = c("firebrick", "steelblue"),
  bty = "n"
)
dev.off()

pdf(file.path(figures_dir, "ablation_performance.pdf"), 7, 5)
plot(
  ablation_benchmark_summary$n_edges,
  ablation_benchmark_summary$median_repeated_tail_seconds,
  type = "b",
  log = "xy",
  pch = 17,
  col = "firebrick",
  ylim = range(c(
    ablation_benchmark_summary$median_repeated_tail_seconds,
    ablation_benchmark_summary$median_prefix_seconds
  )),
  xlab = "Number of possible edges",
  ylab = "Median path-update seconds"
)
lines(
  ablation_benchmark_summary$n_edges,
  ablation_benchmark_summary$median_prefix_seconds,
  type = "b",
  pch = 19,
  col = "steelblue"
)
legend(
  "topleft",
  legend = c("Repeated tail sums", "Prefix-sum update"),
  pch = c(17, 19),
  lty = 1,
  col = c("firebrick", "steelblue"),
  bty = "n"
)
dev.off()

###########################################################################
## Inline manuscript macros and session information
###########################################################################

best_null_global <- max(complete_summary$global_rejection)
best_null_holm <- max(complete_summary$holm_fwer)
mean_holm_tpr <- mean(holm_summary$tpr)
mean_holm_precision <- mean(holm_summary$precision)
largest_speedup <- max(benchmark_summary$median_speedup)
largest_ablation_speedup <- max(
  ablation_benchmark_summary$median_speedup
)

macro_lines <- c(
  sprintf("\\newcommand{\\CompleteNullReps}{%d}", n_complete),
  sprintf("\\newcommand{\\PartialNullReps}{%d}", n_partial),
  sprintf("\\newcommand{\\MaxGlobalNullRate}{%.3f}", best_null_global),
  sprintf("\\newcommand{\\MaxHolmNullFWER}{%.3f}", best_null_holm),
  sprintf("\\newcommand{\\MeanHolmTPR}{%.3f}", mean_holm_tpr),
  sprintf("\\newcommand{\\MeanHolmPrecision}{%.3f}", mean_holm_precision),
  sprintf("\\newcommand{\\LargestMeasuredSpeedup}{%.1f}", largest_speedup),
  sprintf(
    "\\newcommand{\\LargestAblationSpeedup}{%.0f}",
    largest_ablation_speedup
  )
)
writeLines(macro_lines, file.path(generated_dir, "results.tex"))

writeLines(
  capture.output(sessionInfo()),
  file.path(generated_dir, "sessionInfo.txt")
)

replication_finished_at <- Sys.time()
replication_elapsed_seconds <- unname(
  proc.time()[["elapsed"]] - replication_started_elapsed
)
pdflatex <- unname(Sys.which("pdflatex"))
tex_engine <- if (nzchar(pdflatex)) {
  system2(pdflatex, "--version", stdout = TRUE)[[1L]]
} else {
  NA_character_
}
metadata <- data.frame(
  started_at_utc = format(
    replication_started_at,
    tz = "UTC",
    usetz = TRUE
  ),
  finished_at_utc = format(
    replication_finished_at,
    tz = "UTC",
    usetz = TRUE
  ),
  elapsed_seconds = replication_elapsed_seconds,
  package_version = as.character(
    utils::packageVersion("BrainNetTest")
  ),
  r_version = R.version.string,
  platform = R.version$platform,
  machine = unname(Sys.info()[["machine"]]),
  logical_cores = parallel::detectCores(logical = TRUE),
  master_seed = master_seed,
  code_md5 = unname(tools::md5sum("code-full.R")),
  design_md5 = unname(tools::md5sum("simulation-design.md")),
  contract_md5 = unname(tools::md5sum("statistical-contract.md")),
  tex_engine = tex_engine,
  stringsAsFactors = FALSE
)
write.csv(
  metadata,
  file.path(generated_dir, "replication_metadata.csv"),
  row.names = FALSE
)

generated_files <- list.files(
  generated_dir,
  recursive = TRUE,
  full.names = TRUE
)
generated_files <- generated_files[file.info(generated_files)$isdir == FALSE]
generated_files <- setdiff(
  generated_files,
  file.path(generated_dir, "manifest.csv")
)
manifest <- data.frame(
  file = substring(generated_files, nchar(generated_dir) + 2L),
  bytes = file.info(generated_files)$size,
  md5 = unname(tools::md5sum(generated_files)),
  stringsAsFactors = FALSE
)
write.csv(
  manifest,
  file.path(generated_dir, "manifest.csv"),
  row.names = FALSE
)

if (dir.exists(final_generated_dir)) {
  unlink(final_generated_dir, recursive = TRUE, force = TRUE)
}
if (!file.rename(generated_dir, final_generated_dir)) {
  stop("Could not promote validated outputs into `generated/`.")
}

message(
  "Full replication completed successfully in ",
  round(replication_elapsed_seconds, 1),
  " seconds."
)
