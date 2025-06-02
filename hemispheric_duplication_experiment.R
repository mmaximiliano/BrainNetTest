# Hemispheric Duplication Experiment
# Stress-test identify_critical_links() with synthetic brain networks
# ──────────────────────────────────────────────────────────────────────────────

library(parallel)
source("R/identify_critical_links.R")

# 1. Graph Generation Model ───────────────────────────────────────────────────

#' Compute ring distance within hemisphere
#' @param i,j Node indices within same hemisphere
#' @param hemi_size Size of hemisphere (N/2)
ring_distance <- function(i, j, hemi_size) {
  diff <- abs(i - j)
  min(diff, hemi_size - diff)
}

#' Solve for lambda_N to achieve target density ≈ 0.15
#' @param N Total number of nodes
#' @param target_density Target edge density (default 0.15)
solve_lambda <- function(N, target_density = 0.15) {
  hemi_size <- N / 2

  # Objective function: expected density - target
  obj_fn <- function(lambda) {
    total_prob <- 0
    count <- 0
    for (i in 1:(hemi_size-1)) {
      for (j in (i+1):hemi_size) {
        d_ij <- ring_distance(i, j, hemi_size)
        total_prob <- total_prob + exp(-lambda * d_ij)
        count <- count + 1
      }
    }
    expected_density <- total_prob / count
    expected_density - target_density
  }

  # Simple bisection search
  lambda_low <- 0.01
  lambda_high <- 5.0

  for (iter in 1:20) {
    lambda_mid <- (lambda_low + lambda_high) / 2
    if (abs(obj_fn(lambda_mid)) < 1e-6) break

    if (obj_fn(lambda_mid) > 0) {
      lambda_low <- lambda_mid
    } else {
      lambda_high <- lambda_mid
    }
  }

  lambda_mid
}

# Pre-compute lambda values for all N
LAMBDA_CACHE <- list()
for (N in c(10, 100, 1000, 10000)) {
  LAMBDA_CACHE[[as.character(N)]] <- solve_lambda(N)
}

#' Generate hemisphere adjacency matrix using ring topology
#' @param hemi_size Size of hemisphere (N/2)
#' @param lambda Decay parameter
generate_hemisphere <- function(hemi_size, lambda) {
  adj <- matrix(0, hemi_size, hemi_size)

  for (i in 1:(hemi_size-1)) {
    for (j in (i+1):hemi_size) {
      d_ij <- ring_distance(i, j, hemi_size)
      p_ij <- exp(-lambda * d_ij)
      if (runif(1) < p_ij) {
        adj[i, j] <- adj[j, i] <- 1
      }
    }
  }

  adj
}

#' Generate Population A (hemispheric duplication)
#' @param N Total nodes
#' @param n_graphs Number of graphs per population
generate_population_A <- function(N, n_graphs) {
  hemi_size <- N / 2
  lambda <- LAMBDA_CACHE[[as.character(N)]]

  graphs <- vector("list", n_graphs)

  for (g in 1:n_graphs) {
    # Generate left hemisphere
    left_hemi <- generate_hemisphere(hemi_size, lambda)

    # Full graph: copy left to right, zero inter-hemisphere
    full_adj <- matrix(0, N, N)
    full_adj[1:hemi_size, 1:hemi_size] <- left_hemi
    full_adj[(hemi_size+1):N, (hemi_size+1):N] <- left_hemi  # DUPLICATE

    graphs[[g]] <- full_adj
  }

  graphs
}

#' Generate Population B (independent hemispheres)
#' @param N Total nodes
#' @param n_graphs Number of graphs per population
generate_population_B <- function(N, n_graphs) {
  hemi_size <- N / 2
  lambda <- LAMBDA_CACHE[[as.character(N)]]

  graphs <- vector("list", n_graphs)

  for (g in 1:n_graphs) {
    # Generate both hemispheres independently
    left_hemi <- generate_hemisphere(hemi_size, lambda)
    right_hemi <- generate_hemisphere(hemi_size, lambda)

    # Full graph: independent hemispheres, zero inter-hemisphere
    full_adj <- matrix(0, N, N)
    full_adj[1:hemi_size, 1:hemi_size] <- left_hemi
    full_adj[(hemi_size+1):N, (hemi_size+1):N] <- right_hemi  # INDEPENDENT

    graphs[[g]] <- full_adj
  }

  graphs
}

# 2. Ground Truth and Evaluation ──────────────────────────────────────────────

#' Get ground truth critical edges (right hemisphere upper triangle)
#' @param N Total number of nodes
get_ground_truth <- function(N) {
  hemi_size <- N / 2
  right_start <- hemi_size + 1

  edges <- list()
  for (i in right_start:(N-1)) {
    for (j in (i+1):N) {
      edges[[length(edges) + 1]] <- c(i, j)
    }
  }

  edges
}

#' Compute confusion matrix metrics
#' @param predicted_edges List of predicted critical edges
#' @param true_edges List of ground truth critical edges
#' @param N Total number of nodes
compute_metrics <- function(predicted_edges, true_edges, N) {
  # Convert to sets of edge keys
  make_edge_key <- function(edge_list) {
    if (is.null(edge_list) || length(edge_list) == 0) return(character(0))
    sapply(edge_list, function(e) paste(sort(e), collapse="-"))
  }

  pred_keys <- make_edge_key(predicted_edges)
  true_keys <- make_edge_key(true_edges)

  # All possible edges (upper triangle)
  all_edges <- list()
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      all_edges[[length(all_edges) + 1]] <- c(i, j)
    }
  }
  all_keys <- make_edge_key(all_edges)

  # Confusion matrix
  TP <- length(intersect(pred_keys, true_keys))
  FP <- length(setdiff(pred_keys, true_keys))
  FN <- length(setdiff(true_keys, pred_keys))
  TN <- length(all_keys) - TP - FP - FN

  # Metrics
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  precision <- if (TP + FP == 0) 0 else TP / (TP + FP)
  recall <- if (TP + FN == 0) 0 else TP / (TP + FN)
  f1 <- if (precision + recall == 0) 0 else 2 * precision * recall / (precision + recall)

  # Matthews Correlation Coefficient
  mcc_denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- if (mcc_denom == 0) 0 else (TP * TN - FP * FN) / mcc_denom

  list(
    TP = TP, FP = FP, TN = TN, FN = FN,
    accuracy = accuracy,
    precision = precision,
    recall = recall,
    f1 = f1,
    mcc = mcc
  )
}

# 3. Main Experiment Runner ───────────────────────────────────────────────────

#' Run single experiment for given N
#' @param N Number of nodes
#' @param run_id Unique run identifier for seeding
#' @param n_bootstrap Number of bootstrap samples
run_single_experiment <- function(N, run_id, n_bootstrap = 1000) {
  set.seed(run_id)

  # Determine sample sizes
  n_graphs <- if (N <= 1000) 30 else 15

  # Generate populations
  start_time <- Sys.time()

  pop_A <- generate_population_A(N, n_graphs)
  pop_B <- generate_population_B(N, n_graphs)
  populations <- list(A = pop_A, B = pop_B)

  # Run algorithm
  tryCatch({
    result <- identify_critical_links(
      populations = populations,
      alpha = 0.05,
      method = "fisher",
      adjust_method = "none",
      batch_size = 1,
      n_bootstrap = n_bootstrap,
      seed = run_id + 1000
    )

    end_time <- Sys.time()
    runtime <- as.numeric(end_time - start_time, units = "secs")

    # Extract results
    predicted_edges <- if (!is.null(result$critical_edges)) {
      lapply(1:nrow(result$critical_edges), function(i) {
        c(result$critical_edges$node1[i], result$critical_edges$node2[i])
      })
    } else {
      list()
    }

    # Get ground truth and compute metrics
    true_edges <- get_ground_truth(N)
    metrics <- compute_metrics(predicted_edges, true_edges, N)

    # Check for residual signal
    residual_result <- NULL
    residual_p <- NA
    global_significant <- !is.null(result$critical_edges)

    if (global_significant) {
      tryCatch({
        residual_result <- identify_critical_links(
          populations = result$modified_populations,
          alpha = 0.05,
          method = "fisher",
          adjust_method = "none",
          batch_size = 1,
          n_bootstrap = n_bootstrap,
          seed = run_id + 2000
        )
        residual_p <- ifelse(is.null(residual_result$critical_edges), 1, 0)
      }, error = function(e) {
        residual_p <<- NA
      })
    }

    return(list(
      N = N,
      run_id = run_id,
      runtime = runtime,
      global_significant = global_significant,
      TP = metrics$TP,
      FP = metrics$FP,
      TN = metrics$TN,
      FN = metrics$FN,
      accuracy = metrics$accuracy,
      precision = metrics$precision,
      recall = metrics$recall,
      f1 = metrics$f1,
      mcc = metrics$mcc,
      n_predicted = length(predicted_edges),
      n_true = length(true_edges),
      residual_p = residual_p,
      success = TRUE,
      error = NULL
    ))

  }, error = function(e) {
    end_time <- Sys.time()
    runtime <- as.numeric(end_time - start_time, units = "secs")

    return(list(
      N = N,
      run_id = run_id,
      runtime = runtime,
      global_significant = FALSE,
      TP = 0, FP = 0, TN = 0, FN = 0,
      accuracy = 0, precision = 0, recall = 0, f1 = 0, mcc = 0,
      n_predicted = 0, n_true = length(get_ground_truth(N)),
      residual_p = NA,
      success = FALSE,
      error = as.character(e)
    ))
  })
}

#' Run full hemispheric duplication experiment
#' @param N_values Vector of node counts to test
#' @param n_runs Number of Monte Carlo runs per N value
#' @param max_runtime_minutes Maximum runtime for N=10000
#' @param verbose Whether to print progress
#' @return Data frame with all results
run_hemispheric_experiment <- function(N_values = c(10, 100, 1000, 10000),
                                     n_runs = 200,
                                     max_runtime_minutes = 30,
                                     verbose = TRUE) {

  results <- list()
  result_idx <- 1

  if (verbose) {
    cat("Starting Hemispheric Duplication Experiment\n")
    cat("============================================\n")
  }

  for (N in N_values) {
    # Determine sample sizes
    n_graphs <- if (N <= 1000) 30 else 15

    if (verbose) {
      cat(sprintf("\nTesting N = %d with %d graphs per population\n", N, n_graphs))
      cat(sprintf("Lambda = %.4f\n", LAMBDA_CACHE[[as.character(N)]]))
    }

    # Adjust bootstrap for large N if needed
    n_bootstrap <- if (N == 10000) 500 else 1000

    # Track statistics for this N
    successful_runs <- 0
    significant_runs <- 0
    total_tp <- 0
    total_fp <- 0
    total_fn <- 0
    total_tn <- 0
    runtimes <- numeric(0)

    for (run_id in 1:n_runs) {
      if (verbose && run_id %% 50 == 0) {
        cat(sprintf("  Progress: %d/%d\n", run_id, n_runs))
      }

      result <- run_single_experiment(N, run_id, n_bootstrap)
      results[[result_idx]] <- result
      result_idx <- result_idx + 1

      # Track statistics
      if (result$success) {
        successful_runs <- successful_runs + 1
        if (result$global_significant) {
          significant_runs <- significant_runs + 1
        }
        total_tp <- total_tp + result$TP
        total_fp <- total_fp + result$FP
        total_fn <- total_fn + result$FN
        total_tn <- total_tn + result$TN
        runtimes <- c(runtimes, result$runtime)
      }
    }

    # Print summary for this N
    if (verbose && successful_runs > 0) {
      success_rate <- successful_runs / n_runs
      power <- significant_runs / successful_runs
      overall_recall <- if (total_tp + total_fn > 0) total_tp / (total_tp + total_fn) else 0
      overall_precision <- if (total_tp + total_fp > 0) total_tp / (total_tp + total_fp) else 1
      median_runtime <- median(runtimes)

      cat(sprintf("\n  Results for N = %d:\n", N))
      cat(sprintf("    Success rate: %.1f%% (%d/%d)\n",
                  success_rate * 100, successful_runs, n_runs))
      cat(sprintf("    Power (global significance): %.3f\n", power))
      cat(sprintf("    Overall Recall: %.3f\n", overall_recall))
      cat(sprintf("    Overall Precision: %.3f\n", overall_precision))
      cat(sprintf("    Median Runtime: %.2f seconds\n", median_runtime))
      cat(sprintf("    Ground truth edges: %d, Total predicted: %d\n",
                  result$n_true, total_tp + total_fp))
    }
  }

  # Convert to data frame
  df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      N = x$N,
      run_id = x$run_id,
      runtime = x$runtime,
      global_significant = x$global_significant,
      TP = x$TP,
      FP = x$FP,
      TN = x$TN,
      FN = x$FN,
      accuracy = x$accuracy,
      precision = x$precision,
      recall = x$recall,
      f1 = x$f1,
      mcc = x$mcc,
      n_predicted = x$n_predicted,
      n_true = x$n_true,
      residual_p = x$residual_p,
      success = x$success,
      error = ifelse(is.null(x$error), "", x$error),
      stringsAsFactors = FALSE
    )
  }))

  if (verbose) {
    cat("\nExperiment completed!\n")
  }

  return(df)
}

# 4. Metrics and Acceptance Criteria ──────────────────────────────────────────

#' Compute performance metrics from results
#' @param results Data frame from run_hemispheric_experiment
#' @return Data frame with aggregated metrics per N
compute_hemispheric_metrics <- function(results) {

  # Add derived metrics
  results$power <- results$global_significant
  results$fdr <- 1 - results$precision
  results$no_residual <- ifelse(is.na(results$residual_p), 0, results$residual_p > 0.05)

  # Aggregate by N
  library(dplyr)

  metrics <- results %>%
    group_by(N) %>%
    summarise(
      n_experiments = n(),
      success_rate = mean(success, na.rm = TRUE),
      power = mean(power, na.rm = TRUE),
      recall_median = median(recall, na.rm = TRUE),
      recall_q25 = quantile(recall, 0.25, na.rm = TRUE),
      recall_q75 = quantile(recall, 0.75, na.rm = TRUE),
      precision_median = median(precision, na.rm = TRUE),
      precision_q25 = quantile(precision, 0.25, na.rm = TRUE),
      precision_q75 = quantile(precision, 0.75, na.rm = TRUE),
      fdr_median = median(fdr, na.rm = TRUE),
      fdr_q25 = quantile(fdr, 0.25, na.rm = TRUE),
      fdr_q75 = quantile(fdr, 0.75, na.rm = TRUE),
      accuracy_median = median(accuracy, na.rm = TRUE),
      f1_median = median(f1, na.rm = TRUE),
      mcc_median = median(mcc, na.rm = TRUE),
      mcc_q25 = quantile(mcc, 0.25, na.rm = TRUE),
      mcc_q75 = quantile(mcc, 0.75, na.rm = TRUE),
      no_residual_rate = mean(no_residual, na.rm = TRUE),
      median_runtime = median(runtime, na.rm = TRUE),
      q95_runtime = quantile(runtime, 0.95, na.rm = TRUE),
      .groups = "drop"
    )

  return(metrics)
}

#' Check acceptance criteria for hemispheric experiment
#' @param metrics Data frame from compute_hemispheric_metrics
#' @return List with pass/fail status and details
check_hemispheric_acceptance <- function(metrics) {

  # Define runtime limits (in seconds)
  runtime_limits <- c("10" = 5, "100" = 30, "1000" = 180, "10000" = 1800)

  # Tolerance for floating point comparisons
  tolerance <- 1e-6

  results <- list()
  overall_pass <- TRUE

  for (i in 1:nrow(metrics)) {
    N <- metrics$N[i]

    # Check acceptance criteria
    power_pass <- !is.na(metrics$power[i]) && (metrics$power[i] >= (0.95 - tolerance))
    recall_pass <- !is.na(metrics$recall_median[i]) && (metrics$recall_median[i] >= (0.80 - tolerance))
    fdr_pass <- !is.na(metrics$fdr_median[i]) && (metrics$fdr_median[i] <= (0.10 + tolerance))
    residual_pass <- !is.na(metrics$no_residual_rate[i]) && (metrics$no_residual_rate[i] >= (0.95 - tolerance))
    runtime_pass <- !is.na(metrics$q95_runtime[i]) && metrics$q95_runtime[i] <= runtime_limits[as.character(N)]

    condition_pass <- power_pass & recall_pass & fdr_pass & residual_pass & runtime_pass

    results[[as.character(N)]] <- list(
      N = N,
      power_pass = power_pass,
      recall_pass = recall_pass,
      fdr_pass = fdr_pass,
      residual_pass = residual_pass,
      runtime_pass = runtime_pass,
      overall_pass = condition_pass,
      power = ifelse(is.na(metrics$power[i]), 0, metrics$power[i]),
      recall_median = ifelse(is.na(metrics$recall_median[i]), 0, metrics$recall_median[i]),
      fdr_median = ifelse(is.na(metrics$fdr_median[i]), 1, metrics$fdr_median[i]),
      mcc_median = ifelse(is.na(metrics$mcc_median[i]), 0, metrics$mcc_median[i]),
      no_residual_rate = ifelse(is.na(metrics$no_residual_rate[i]), 0, metrics$no_residual_rate[i]),
      runtime_95th = ifelse(is.na(metrics$q95_runtime[i]), Inf, metrics$q95_runtime[i]),
      runtime_limit = runtime_limits[as.character(N)]
    )

    overall_pass <- overall_pass & condition_pass
  }

  return(list(
    overall_pass = overall_pass,
    conditions = results
  ))
}

#' Print hemispheric validation results
#' @param acceptance_results Results from check_hemispheric_acceptance
print_hemispheric_results <- function(acceptance_results) {

  cat("=== HEMISPHERIC DUPLICATION EXPERIMENT RESULTS ===\n\n")

  if (acceptance_results$overall_pass) {
    cat("✓ OVERALL RESULT: PASS\n\n")
  } else {
    cat("✗ OVERALL RESULT: FAIL\n\n")
  }

  cat("Detailed results by N:\n")
  cat("────────────────────────────────────────────────────────\n")

  for (condition in acceptance_results$conditions) {
    status <- if (condition$overall_pass) "✓ PASS" else "✗ FAIL"

    cat(sprintf("N=%d: %s\n", condition$N, status))

    cat(sprintf("  Power (global p ≤ 0.05): %.3f (≥0.95: %s)\n",
                condition$power,
                if (condition$power_pass) "✓" else "✗"))

    cat(sprintf("  Median Recall: %.3f (≥0.80: %s)\n",
                condition$recall_median,
                if (condition$recall_pass) "✓" else "✗"))

    cat(sprintf("  Median FDR: %.3f (≤0.10: %s)\n",
                condition$fdr_median,
                if (condition$fdr_pass) "✓" else "✗"))

    cat(sprintf("  Median MCC: %.3f (reported)\n", condition$mcc_median))

    cat(sprintf("  No residual: %.3f (≥0.95: %s)\n",
                condition$no_residual_rate,
                if (condition$residual_pass) "✓" else "✗"))

    cat(sprintf("  Runtime 95th: %.2fs (≤%.0fs: %s)\n",
                condition$runtime_95th, condition$runtime_limit,
                if (condition$runtime_pass) "✓" else "✗"))

    cat("\n")
  }
}

#' Save hemispheric experiment results to CSV files
#' @param results Raw results data frame
#' @param metrics Aggregated metrics data frame
#' @param acceptance Acceptance criteria results
#' @param output_dir Directory to save files
#' @param master_seed Master seed used
#' @param quick_test Whether this was a quick test
save_hemispheric_results <- function(results, metrics, acceptance, output_dir = ".",
                                   master_seed = 42, quick_test = FALSE) {

  # Create timestamp for unique filenames
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  test_type <- if (quick_test) "quick" else "full"

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save detailed results
  results_file <- file.path(output_dir,
                           sprintf("hemispheric_experiment_results_%s_%s_seed%d.csv",
                                   test_type, timestamp, master_seed))

  # Add derived metrics
  results_with_metrics <- results
  results_with_metrics$power <- results$global_significant
  results_with_metrics$fdr <- 1 - results$precision
  results_with_metrics$no_residual <- ifelse(is.na(results$residual_p), 0, results$residual_p > 0.05)

  write.csv(results_with_metrics, results_file, row.names = FALSE)
  cat(sprintf("  Detailed results saved to: %s\n", results_file))

  # Save aggregated metrics
  metrics_file <- file.path(output_dir,
                           sprintf("hemispheric_experiment_metrics_%s_%s_seed%d.csv",
                                   test_type, timestamp, master_seed))
  write.csv(metrics, metrics_file, row.names = FALSE)
  cat(sprintf("  Aggregated metrics saved to: %s\n", metrics_file))

  # Save acceptance criteria summary
  acceptance_file <- file.path(output_dir,
                              sprintf("hemispheric_experiment_acceptance_%s_%s_seed%d.csv",
                                      test_type, timestamp, master_seed))

  # Convert acceptance results to data frame
  acceptance_df <- do.call(rbind, lapply(acceptance$conditions, function(x) {
    data.frame(
      N = x$N,
      overall_pass = x$overall_pass,
      power = x$power,
      power_pass = x$power_pass,
      recall_median = x$recall_median,
      recall_pass = x$recall_pass,
      fdr_median = x$fdr_median,
      fdr_pass = x$fdr_pass,
      mcc_median = x$mcc_median,
      no_residual_rate = x$no_residual_rate,
      residual_pass = x$residual_pass,
      runtime_95th = x$runtime_95th,
      runtime_limit = x$runtime_limit,
      runtime_pass = x$runtime_pass,
      stringsAsFactors = FALSE
    )
  }))

  # Add overall summary
  acceptance_df$experiment_timestamp <- timestamp
  acceptance_df$master_seed <- master_seed
  acceptance_df$test_type <- test_type
  acceptance_df$overall_experiment_pass <- acceptance$overall_pass

  write.csv(acceptance_df, acceptance_file, row.names = FALSE)
  cat(sprintf("  Acceptance criteria saved to: %s\n", acceptance_file))

  # Create a summary report
  summary_file <- file.path(output_dir,
                           sprintf("hemispheric_experiment_summary_%s_%s_seed%d.txt",
                                   test_type, timestamp, master_seed))

  sink(summary_file)
  cat("=== HEMISPHERIC DUPLICATION EXPERIMENT SUMMARY ===\n")
  cat(sprintf("Timestamp: %s\n", Sys.time()))
  cat(sprintf("Master Seed: %d\n", master_seed))
  cat(sprintf("Test Type: %s\n", test_type))
  cat(sprintf("Overall Pass: %s\n", if (acceptance$overall_pass) "YES" else "NO"))
  cat("\n")
  print_hemispheric_results(acceptance)
  sink()

  cat(sprintf("  Summary report saved to: %s\n", summary_file))
}

# 5. Main Validation Function ─────────────────────────────────────────────────

#' Main function to run the complete hemispheric duplication validation
#' @param quick_test If TRUE, run a reduced test for development
#' @param nodes Maximum N value to test (10, 100, 1000, or 10000)
#' @param master_seed Master seed for reproducibility
#' @param output_csv If TRUE, save results to CSV files
#' @param output_dir Directory to save CSV files
#' @param ... Additional arguments passed to run_hemispheric_experiment
main_hemispheric_validation <- function(quick_test = FALSE, nodes = 10000, master_seed = 42,
                                       output_csv = TRUE, output_dir = ".", ...) {

  # Set master seed for reproducibility
  set.seed(master_seed)
  cat(sprintf("Using master seed: %d\n", master_seed))

  # Determine N values to test
  if (quick_test) {
    if (nodes == 10) {
      N_values <- c(10)
      n_runs <- 10
    } else if (nodes == 100) {
      N_values <- c(10, 100)
      n_runs <- 20
    } else if (nodes == 1000) {
      N_values <- c(10, 100, 1000)
      n_runs <- 50
    } else {
      N_values <- c(10, 100, 1000, 10000)
      n_runs <- 10
    }
    cat("Running quick test...\n")
  } else {
    if (nodes == 10) {
      N_values <- c(10)
    } else if (nodes == 100) {
      N_values <- c(10, 100)
    } else if (nodes == 1000) {
      N_values <- c(10, 100, 1000)
    } else {
      N_values <- c(10, 100, 1000, 10000)
    }
    n_runs <- 200
    cat("Running full validation...\n")
  }

  # Run the experiment
  results <- run_hemispheric_experiment(
    N_values = N_values,
    n_runs = n_runs,
    verbose = TRUE,
    ...
  )

  cat("\nComputing metrics...\n")
  metrics <- compute_hemispheric_metrics(results)

  # Print detailed metrics summary
  cat("\n=== DETAILED METRICS SUMMARY ===\n")
  for (i in 1:nrow(metrics)) {
    N <- metrics$N[i]
    cat(sprintf("\nN = %d:\n", N))
    cat(sprintf("  Experiments: %d (%.1f%% success)\n",
                metrics$n_experiments[i], metrics$success_rate[i] * 100))
    cat(sprintf("  Power: %.3f\n", metrics$power[i]))
    cat(sprintf("  Recall: %.3f [%.3f, %.3f]\n",
                metrics$recall_median[i], metrics$recall_q25[i], metrics$recall_q75[i]))
    cat(sprintf("  Precision: %.3f [%.3f, %.3f]\n",
                metrics$precision_median[i], metrics$precision_q25[i], metrics$precision_q75[i]))
    cat(sprintf("  FDR: %.3f [%.3f, %.3f]\n",
                metrics$fdr_median[i], metrics$fdr_q25[i], metrics$fdr_q75[i]))
    cat(sprintf("  MCC: %.3f [%.3f, %.3f]\n",
                metrics$mcc_median[i], metrics$mcc_q25[i], metrics$mcc_q75[i]))
    cat(sprintf("  No residual rate: %.3f\n", metrics$no_residual_rate[i]))
    cat(sprintf("  Runtime: %.2fs (95th: %.2fs)\n",
                metrics$median_runtime[i], metrics$q95_runtime[i]))
  }

  cat("\nChecking acceptance criteria...\n")
  acceptance <- check_hemispheric_acceptance(metrics)

  print_hemispheric_results(acceptance)

  # Save results to CSV if requested
  if (output_csv) {
    cat("\nSaving results to CSV files...\n")
    save_hemispheric_results(results, metrics, acceptance, output_dir, master_seed, quick_test)
  }

  return(list(
    results = results,
    metrics = metrics,
    acceptance = acceptance,
    master_seed = master_seed
  ))
}

# 6. Example Usage ────────────────────────────────────────────────────────────

if (TRUE) {  # Set to TRUE to run
  # Quick test for development
  quick_results <- main_hemispheric_validation(quick_test = TRUE, nodes = 100, master_seed = 42)

  # Full validation (uncomment to run)
  # full_results <- main_hemispheric_validation(master_seed = 42)
}

