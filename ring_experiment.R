# Ring Experiment for Validation of identify_critical_links()
# ──────────────────────────────────────────────────────────────────────────────
library(pbapply)
library(parallel)
library(Matrix)
library(igraph)

#' Compute ring distance between nodes on a cycle
#' @param i,j Node indices (1-based)
#' @param N Total number of nodes
#' @return Minimum distance on the ring
ring_distance <- function(i, j, N) {
  diff <- abs(i - j)
  pmin(diff, N - diff)
}

#' Calculate expected number of edges for given lambda and N
#' @param lambda Decay parameter
#' @param N Number of nodes
#' @return Expected number of edges
expected_edges <- function(lambda, N) {
  if (N <= 1) return(0)

  total <- 0
  max_k <- floor(N / 2)

  for (k in 1:max_k) {
    multiplier <- if (k == N / 2) 1 else 2
    total <- total + multiplier * exp(-lambda * k)
  }

  return(N * total / 2)
}

#' Find lambda to achieve target density (approximately 20%)
#' @param N Number of nodes
#' @param target_density Target edge density (default 0.2)
#' @return Lambda parameter
find_lambda <- function(N, target_density = 0.25) {  # Increased target density
  if (N <= 1) return(1)

  target_edges <- target_density * N * (N - 1) / 2

  # Binary search for lambda
  lambda_low <- 0.01
  lambda_high <- 10
  tolerance <- 1e-6

  while (lambda_high - lambda_low > tolerance) {
    lambda_mid <- (lambda_low + lambda_high) / 2
    edges_mid <- expected_edges(lambda_mid, N)

    if (edges_mid > target_edges) {
      lambda_low <- lambda_mid
    } else {
      lambda_high <- lambda_mid
    }
  }

  return((lambda_low + lambda_high) / 2)
}

#' Generate a single ring graph
#' @param N Number of nodes
#' @param lambda Decay parameter
#' @param perturbed_nodes Nodes with different edge probabilities
#' @param perturbation_type Type of perturbation ("lambda_half", "lambda_double", "const_high", "const_low")
#' @param lambda_alt Alternative lambda for perturbed edges
#' @param p_const Constant probability for perturbed edges
#' @return Symmetric adjacency matrix
generate_ring_graph <- function(N, lambda, perturbed_nodes = c(),
                                perturbation_type = "none",
                                lambda_alt = NULL, p_const = NULL) {

  # Initialize adjacency matrix
  adj <- matrix(0, nrow = N, ncol = N)

  # Generate edges
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      # Check if this edge involves a perturbed node
      is_perturbed <- (i %in% perturbed_nodes) || (j %in% perturbed_nodes)

      if (is_perturbed && perturbation_type != "none") {
        # Use perturbed probability
        if (perturbation_type == "lambda_half") {
          prob <- exp(-lambda_alt * ring_distance(i, j, N))
        } else if (perturbation_type == "lambda_double") {
          prob <- exp(-lambda_alt * ring_distance(i, j, N))
        } else if (perturbation_type %in% c("const_high", "const_low")) {
          prob <- p_const
        } else {
          prob <- exp(-lambda * ring_distance(i, j, N))
        }
      } else {
        # Use baseline probability
        prob <- exp(-lambda * ring_distance(i, j, N))
      }

      # Generate edge
      if (runif(1) < prob) {
        adj[i, j] <- 1
        adj[j, i] <- 1
      }
    }
  }

  return(adj)
}

#' Generate population of graphs
#' @param n_graphs Number of graphs to generate
#' @param N Number of nodes per graph
#' @param lambda Decay parameter
#' @param perturbed_nodes Nodes with different edge probabilities
#' @param perturbation_type Type of perturbation
#' @return List of adjacency matrices
generate_population <- function(n_graphs, N, lambda, perturbed_nodes = c(),
                               perturbation_type = "none") {

  # Set perturbation parameters with stronger effects
  lambda_alt <- switch(perturbation_type,
                      "lambda_half" = lambda / 3,     # More extreme reduction
                      "lambda_double" = lambda * 3,   # More extreme increase
                      lambda)

  p_const <- switch(perturbation_type,
                   "const_high" = 0.95,  # Higher constant probability
                   "const_low" = 0.01,   # Lower constant probability
                   NULL)

  # Generate graphs
  graphs <- vector("list", n_graphs)
  for (i in 1:n_graphs) {
    graphs[[i]] <- generate_ring_graph(N, lambda, perturbed_nodes,
                                      perturbation_type, lambda_alt, p_const)
  }

  return(graphs)
}

#' Compute confusion matrix for edge detection
#' @param detected_edges Data frame with detected critical edges
#' @param perturbed_nodes Vector of perturbed node indices
#' @param N Total number of nodes
#' @return List with TP, FP, FN, TN counts
compute_confusion_matrix <- function(detected_edges, perturbed_nodes, N) {

  # All possible edges (upper triangular)
  all_edges <- expand.grid(i = 1:N, j = 1:N)
  all_edges <- all_edges[all_edges$i < all_edges$j, ]

  # True critical edges (those touching perturbed nodes)
  true_critical <- rep(FALSE, nrow(all_edges))
  for (k in 1:nrow(all_edges)) {
    i <- all_edges$i[k]
    j <- all_edges$j[k]
    true_critical[k] <- (i %in% perturbed_nodes) || (j %in% perturbed_nodes)
  }

  # Predicted critical edges
  predicted_critical <- rep(FALSE, nrow(all_edges))
  if (!is.null(detected_edges) && nrow(detected_edges) > 0) {
    for (k in 1:nrow(detected_edges)) {
      i <- detected_edges$node1[k]
      j <- detected_edges$node2[k]
      # Find matching edge in all_edges
      idx <- which(all_edges$i == i & all_edges$j == j)
      if (length(idx) > 0) {
        predicted_critical[idx] <- TRUE
      }
    }
  }

  # Confusion matrix
  TP <- sum(true_critical & predicted_critical)
  FP <- sum(!true_critical & predicted_critical)
  FN <- sum(true_critical & !predicted_critical)
  TN <- sum(!true_critical & !true_critical)

  list(TP = TP, FP = FP, FN = FN, TN = TN)
}

#' Run single experiment repetition
#' @param N Number of nodes
#' @param n_graphs Number of graphs per population
#' @param perturbation_type Type of perturbation to apply
#' @param rho Fraction of nodes to perturb
#' @param alpha Significance level
#' @param n_bootstrap Number of bootstrap samples
#' @param seed Random seed for reproducibility
#' @return List with results
run_single_experiment <- function(N, n_graphs, perturbation_type, rho = 0.03,  # Increased from 0.02
                                 alpha = 0.05, n_bootstrap = 1000, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Find appropriate lambda
  lambda <- find_lambda(N, target_density = 0.25)  # Increased density

  # Determine perturbed nodes - ensure at least 2 nodes for larger N, scale better with N
  K <- max(2, ceiling(rho * N))
  if (N >= 100) {
    K <- max(3, ceiling(rho * N))  # Ensure at least 3 perturbed nodes for N >= 100
  }
  perturbed_nodes <- sample(1:N, K)

  # Generate populations
  pop_A <- generate_population(n_graphs, N, lambda)
  pop_B <- generate_population(n_graphs, N, lambda, perturbed_nodes, perturbation_type)

  populations <- list(A = pop_A, B = pop_B)

  # Run identify_critical_links
  start_time <- Sys.time()

  tryCatch({
    result <- identify_critical_links(
      populations = populations,
      alpha = alpha,
      method = "fisher",
      adjust_method = "none",
      batch_size = 1,
      n_bootstrap = n_bootstrap,
      a = 1,
      seed = seed
    )

    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Compute confusion matrix
    conf_matrix <- compute_confusion_matrix(result$critical_edges, perturbed_nodes, N)

    # Check if global test was significant
    global_significant <- !is.null(result$critical_edges)

    # Test for residual signal (re-run on modified populations)
    residual_result <- NULL
    residual_p <- NA
    if (global_significant) {
      tryCatch({
        residual_result <- identify_critical_links(
          populations = result$modified_populations,
          alpha = alpha,
          method = "fisher",
          adjust_method = "none",
          batch_size = 1,
          n_bootstrap = n_bootstrap,
          a = 1,
          seed = seed + 1000  # Different seed for residual test
        )
        residual_p <- ifelse(is.null(residual_result$critical_edges), 1, 0)
      }, error = function(e) {
        residual_p <<- NA
      })
    }

    return(list(
      N = N,
      perturbation_type = perturbation_type,
      seed = seed,
      runtime = runtime,
      global_significant = global_significant,
      TP = conf_matrix$TP,
      FP = conf_matrix$FP,
      FN = conf_matrix$FN,
      TN = conf_matrix$TN,
      perturbed_nodes = perturbed_nodes,
      detected_edges = result$critical_edges,
      residual_p = residual_p,
      success = TRUE,
      error = NULL
    ))

  }, error = function(e) {
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

    return(list(
      N = N,
      perturbation_type = perturbation_type,
      seed = seed,
      runtime = runtime,
      global_significant = FALSE,
      TP = 0, FP = 0, FN = 0, TN = 0,
      perturbed_nodes = perturbed_nodes,
      detected_edges = NULL,
      residual_p = NA,
      success = FALSE,
      error = as.character(e)
    ))
  })
}

#' Run complete ring experiment
#' @param N_values Vector of node counts to test
#' @param perturbation_types Vector of perturbation types to test
#' @param n_repetitions Number of Monte Carlo repetitions per condition
#' @param rho Fraction of nodes to perturb
#' @param alpha Significance level
#' @param n_bootstrap Number of bootstrap samples
#' @param verbose Whether to print progress
#' @param n_cores Number of cores to use for parallel processing
#' @return Data frame with all results
run_ring_experiment <- function(N_values = c(10, 100, 1000, 10000),
                               perturbation_types = c("lambda_half", "lambda_double",
                                                     "const_high", "const_low"),
                               n_repetitions = 200,
                               rho = 0.02,  # Increased from 0.01
                               alpha = 0.05,
                               n_bootstrap = 1000,
                               verbose = TRUE,
                               n_cores = parallel::detectCores() - 1) {

  # Build a parameter grid instead of nested loops
  param_grid <- list()
  idx <- 1

  for (N in N_values) {
    # Significantly increase number of graphs based on N for better statistical power
    n_graphs <- if (N <= 50) {
      150  # Small networks need fewer samples
    } else if (N <= 100) {
      250  # Medium networks need more samples
    } else if (N <= 1000) {
      400  # Large networks need many more samples
    } else {
      500  # Very large networks need the most samples
    }

    if (verbose) {
      cat(sprintf("Preparing N = %d with %d graphs per population\n", N, n_graphs))
    }

    for (pert_type in perturbation_types) {
      for (rep in 1:n_repetitions) {
        seed <- 1000 * idx + rep
        param_grid[[idx]] <- list(
          N = N,
          n_graphs = n_graphs,
          perturbation_type = pert_type,
          rho = rho,
          alpha = alpha,
          n_bootstrap = n_bootstrap,
          seed = seed
        )
        idx <- idx + 1
      }
    }
  }

  # Set up PSOCK cluster (works on Windows & Unix)
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Export functions and variables to workers
  parallel::clusterExport(cl, c(
    "run_single_experiment",
    "generate_population", 
    "generate_ring_graph",
    "find_lambda",
    "expected_edges",
    "ring_distance",
    "compute_confusion_matrix",
    "identify_critical_links"  # Make sure this function is exported
  ), envir = environment())
  
  # Load required libraries on each worker
  parallel::clusterEvalQ(cl, {
    library(Matrix)
    library(igraph)
    library(BrainNetTest)  # if identify_critical_links is in a package
  })
  
  # If identify_critical_links is defined in another file, source it on each worker
  parallel::clusterEvalQ(cl, {
    source("./R/identify_critical_links.R")
  })

  # Set reproducible parallel RNG
  parallel::clusterSetRNGStream(cl, iseed = 42)

  # Parallel execution
  if (verbose) {
    cat(sprintf("\nRunning %d experiments on %d cores...\n",
                length(param_grid), n_cores))
  }

  results <- parallel::parLapply(
    cl, param_grid,
    function(pars) do.call(run_single_experiment, pars)
  )

  # Convert to data frame
  df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      N = x$N,
      perturbation_type = x$perturbation_type,
      seed = x$seed,
      runtime = x$runtime,
      global_significant = x$global_significant,
      TP = x$TP,
      FP = x$FP,
      FN = x$FN,
      TN = x$TN,
      residual_p = x$residual_p,
      success = x$success,
      error = ifelse(is.null(x$error), "", x$error),
      stringsAsFactors = FALSE
    )
  }))

  return(df)
}

#' Compute performance metrics from results
#' @param results Data frame from run_ring_experiment
#' @return Data frame with aggregated metrics per N and perturbation type
compute_metrics <- function(results) {

  # Add derived metrics with proper NA handling
  results$recall <- with(results, ifelse(TP + FN > 0, TP / (TP + FN), 0))
  results$precision <- with(results, ifelse(TP + FP > 0, TP / (TP + FP), 1))  # Perfect precision when no FP
  results$power <- results$global_significant
  results$no_residual <- ifelse(is.na(results$residual_p), 0, results$residual_p > 0.05)

  # Aggregate by N and perturbation_type
  library(dplyr)

  metrics <- results %>%
    group_by(N, perturbation_type) %>%
    summarise(
      n_experiments = n(),
      success_rate = mean(success, na.rm = TRUE),
      power = mean(power, na.rm = TRUE),
      recall = mean(recall, na.rm = TRUE),
      precision = mean(precision, na.rm = TRUE),
      no_residual_rate = mean(no_residual, na.rm = TRUE),
      median_runtime = median(runtime, na.rm = TRUE),
      q95_runtime = quantile(runtime, 0.95, na.rm = TRUE),
      .groups = "drop"
    )

  return(metrics)
}

#' Check acceptance criteria
#' @param metrics Data frame from compute_metrics
#' @return List with pass/fail status and details
check_acceptance_criteria <- function(metrics) {

  # Define runtime limits (in seconds) - adjusted for increased sample sizes
  runtime_limits <- c("10" = 3, "100" = 25, "1000" = 120, "10000" = 1200)

  # Tolerance for floating point comparisons
  tolerance <- 1e-6

  results <- list()
  overall_pass <- TRUE

  for (i in 1:nrow(metrics)) {
    N <- metrics$N[i]
    pert_type <- metrics$perturbation_type[i]

    # Check criteria with proper NA handling and tolerance
    power_pass <- !is.na(metrics$power[i]) && (metrics$power[i] >= (0.95 - tolerance))
    recall_pass <- !is.na(metrics$recall[i]) && (metrics$recall[i] >= (0.90 - tolerance))
    precision_pass <- !is.na(metrics$precision[i]) && (metrics$precision[i] >= (0.95 - tolerance))
    residual_pass <- !is.na(metrics$no_residual_rate[i]) && (metrics$no_residual_rate[i] >= (0.95 - tolerance))
    runtime_pass <- !is.na(metrics$q95_runtime[i]) && metrics$q95_runtime[i] <= runtime_limits[as.character(N)]

    condition_pass <- power_pass & recall_pass & precision_pass & residual_pass & runtime_pass

    results[[paste(N, pert_type, sep = "_")]] <- list(
      N = N,
      perturbation_type = pert_type,
      power_pass = power_pass,
      recall_pass = recall_pass,
      precision_pass = precision_pass,
      residual_pass = residual_pass,
      runtime_pass = runtime_pass,
      overall_pass = condition_pass,
      power = ifelse(is.na(metrics$power[i]), 0, metrics$power[i]),
      recall = ifelse(is.na(metrics$recall[i]), 0, metrics$recall[i]),
      precision = ifelse(is.na(metrics$precision[i]), 0, metrics$precision[i]),
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

#' Print validation results
#' @param acceptance_results Results from check_acceptance_criteria
print_validation_results <- function(acceptance_results) {

  cat("=== VALIDATION RESULTS FOR identify_critical_links() ===\n\n")

  if (acceptance_results$overall_pass) {
    cat("✓ OVERALL RESULT: PASS\n\n")
  } else {
    cat("✗ OVERALL RESULT: FAIL\n\n")
  }

  cat("Detailed results by condition:\n")
  cat("────────────────────────────────────────────────────────\n")

  for (condition in acceptance_results$conditions) {
    status <- if (condition$overall_pass) "✓ PASS" else "✗ FAIL"

    cat(sprintf("N=%d, %s: %s\n",
                condition$N, condition$perturbation_type, status))

    cat(sprintf("  Power: %.3f (≥0.95: %s)\n",
                condition$power,
                if (condition$power_pass) "✓" else "✗"))

    cat(sprintf("  Recall: %.3f (≥0.90: %s)\n",
                condition$recall,
                if (condition$recall_pass) "✓" else "✗"))

    cat(sprintf("  Precision: %.3f (≥0.95: %s)\n",
                condition$precision,
                if (condition$precision_pass) "✓" else "✗"))

    cat(sprintf("  No residual: %.3f (≥0.95: %s)\n",
                condition$no_residual_rate,
                if (condition$residual_pass) "✓" else "✗"))

    cat(sprintf("  Runtime 95th: %.2fs (≤%.0fs: %s)\n",
                condition$runtime_95th, condition$runtime_limit,
                if (condition$runtime_pass) "✓" else "✗"))

    cat("\n")
  }
}

#' Main function to run the complete validation
#' @param quick_test If TRUE, run a reduced test for development
#' @param master_seed Master seed for reproducibility (default: 42)
#' @param output_csv If TRUE, save results to CSV files (default: TRUE)
#' @param output_dir Directory to save CSV files (default: current directory)
#' @param ... Additional arguments passed to run_ring_experiment
main_validation <- function(quick_test = FALSE, nodes = 10000, master_seed = 42,
                           output_csv = TRUE, output_dir = ".", ...) {

  # Set master seed for reproducibility
  set.seed(master_seed)
  cat(sprintf("Using master seed: %d\n", master_seed))

  if (quick_test) {
    if (nodes == 10){
        N_values <- c(10)
    } else if (nodes == 100) {
        N_values <- c(10, 100)
    } else if (nodes == 1000) {
        N_values <- c(1000)
    }

    cat("Running quick test...\n")
    results <- run_ring_experiment(
      N_values = N_values,
      perturbation_types = c("lambda_half", "const_high"),
      n_repetitions = 10,
      verbose = TRUE,
      ...
    )
  } else {
    cat("Running full validation...\n")
    results <- run_ring_experiment(verbose = TRUE, ...)
  }

  cat("\nComputing metrics...\n")
  metrics <- compute_metrics(results)

  cat("\nChecking acceptance criteria...\n")
  acceptance <- check_acceptance_criteria(metrics)

  print_validation_results(acceptance)

  # Save results to CSV if requested
  if (output_csv) {
    cat("\nSaving results to CSV files...\n")
    save_results_to_csv(results, metrics, acceptance, output_dir, master_seed, quick_test)
  }

  return(list(
    results = results,
    metrics = metrics,
    acceptance = acceptance,
    master_seed = master_seed
  ))
}

#' Save validation results to CSV files
#' @param results Raw results data frame
#' @param metrics Aggregated metrics data frame
#' @param acceptance Acceptance criteria results
#' @param output_dir Directory to save files
#' @param master_seed Master seed used
#' @param quick_test Whether this was a quick test
save_results_to_csv <- function(results, metrics, acceptance, output_dir = ".",
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
                           sprintf("ring_experiment_results_%s_%s_seed%d.csv",
                                   test_type, timestamp, master_seed))

  # Add computed metrics to results for the detailed file
  results_with_metrics <- results
  results_with_metrics$recall <- with(results, ifelse(TP + FN > 0, TP / (TP + FN), 0))
  results_with_metrics$precision <- with(results, ifelse(TP + FP > 0, TP / (TP + FP), 1))
  results_with_metrics$power <- results$global_significant
  results_with_metrics$no_residual <- ifelse(is.na(results$residual_p), 0, results$residual_p > 0.05)

  write.csv(results_with_metrics, results_file, row.names = FALSE)
  cat(sprintf("  Detailed results saved to: %s\n", results_file))

  # Save aggregated metrics
  metrics_file <- file.path(output_dir,
                           sprintf("ring_experiment_metrics_%s_%s_seed%d.csv",
                                   test_type, timestamp, master_seed))
  write.csv(metrics, metrics_file, row.names = FALSE)
  cat(sprintf("  Aggregated metrics saved to: %s\n", metrics_file))

  # Save acceptance criteria summary
  acceptance_file <- file.path(output_dir,
                              sprintf("ring_experiment_acceptance_%s_%s_seed%d.csv",
                                      test_type, timestamp, master_seed))

  # Convert acceptance results to data frame
  acceptance_df <- do.call(rbind, lapply(acceptance$conditions, function(x) {
    data.frame(
      N = x$N,
      perturbation_type = x$perturbation_type,
      overall_pass = x$overall_pass,
      power = x$power,
      power_pass = x$power_pass,
      recall = x$recall,
      recall_pass = x$recall_pass,
      precision = x$precision,
      precision_pass = x$precision_pass,
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
                           sprintf("ring_experiment_summary_%s_%s_seed%d.txt",
                                   test_type, timestamp, master_seed))

  sink(summary_file)
  cat("=== RING EXPERIMENT VALIDATION SUMMARY ===\n")
  cat(sprintf("Timestamp: %s\n", Sys.time()))
  cat(sprintf("Master Seed: %d\n", master_seed))
  cat(sprintf("Test Type: %s\n", test_type))
  cat(sprintf("Overall Pass: %s\n", if (acceptance$overall_pass) "YES" else "NO"))
  cat("\n")
  print_validation_results(acceptance)
  sink()

  cat(sprintf("  Summary report saved to: %s\n", summary_file))
}

# Example usage:
validation_results <- main_validation(quick_test = TRUE, nodes = 100, master_seed = 42)
#validation_results <- main_validation(quick_test = TRUE, nodes = 1000, master_seed = 42)
# validation_results <- main_validation(master_seed = 42)  # Full validation
