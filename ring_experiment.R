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

  # Compute distance matrix using vectorized operation
  D <- outer(1:N, 1:N, ring_distance, N)

  # Compute baseline probability matrix
  P <- exp(-lambda * D)

  # Handle perturbations if any
  if (length(perturbed_nodes) > 0 && perturbation_type != "none") {
    # Create perturbation mask - edges involving perturbed nodes
    perturb_mask <- outer(1:N, 1:N, function(i, j) {
      (i %in% perturbed_nodes) | (j %in% perturbed_nodes)
    })

    if (perturbation_type == "lambda_half") {
      P_perturbed <- exp(-lambda_alt * D)
      P[perturb_mask] <- P_perturbed[perturb_mask]
    } else if (perturbation_type == "lambda_double") {
      P_perturbed <- exp(-lambda_alt * D)
      P[perturb_mask] <- P_perturbed[perturb_mask]
    } else if (perturbation_type %in% c("const_high", "const_low")) {
      P[perturb_mask] <- p_const
    }
  }

  # Get upper triangular indices (avoid diagonal and lower triangle)
  upper_tri_idx <- which(upper.tri(P, diag = FALSE))

  # Generate edges only for upper triangle
  n_edges <- length(upper_tri_idx)
  edge_probs <- P[upper_tri_idx]
  edges <- rbinom(n_edges, 1, edge_probs)

  # Initialize adjacency matrix
  adj <- matrix(0, nrow = N, ncol = N)

  # Set edges in upper triangle
  adj[upper_tri_idx] <- edges

  # Make symmetric (copy upper triangle to lower triangle)
  adj <- adj + t(adj)

  return(adj)
}

#' Generate population of graphs
#' @param n_graphs Number of graphs to generate
#' @param N Number of nodes per graph
#' @param lambda Decay parameter
#' @param perturbed_nodes Nodes with different edge probabilities
#' @param perturbation_type Type of perturbation
#' @param lambda_mult Specific lambda multiplier (for lambda_half, lambda_double)
#' @param p_const_value Specific constant probability value (for const_high, const_low)
#' @return List of adjacency matrices
generate_population <- function(n_graphs, N, lambda, perturbed_nodes = c(),
                               perturbation_type = "none",
                               lambda_mult = NULL, p_const_value = NULL) {

  # Define available parameter values (used if not explicitly provided)
  lambda_multipliers <- switch(perturbation_type,
                      "lambda_half" = c(0.5, 0.33, 0.25),  # Various reduction levels
                      "lambda_double" = c(2, 3, 4),        # Various increase levels
                      1)  # Default case

  p_const_values <- switch(perturbation_type,
                   "const_high" = c(0.85, 0.90, 0.95),  # Different high probabilities
                   "const_low" = c(0.05, 0.03, 0.01),   # Different low probabilities
                   NULL)

  # Use provided parameter value or default to the first/last
  if (perturbation_type %in% c("lambda_half", "lambda_double")) {
    if (is.null(lambda_mult)) {
      lambda_mult <- lambda_multipliers[1]  # Default to first multiplier
    }
    lambda_alt <- lambda * lambda_mult
  } else {
    lambda_alt <- lambda  # Default for other perturbation types
  }

  if (perturbation_type %in% c("const_high", "const_low")) {
    if (is.null(p_const_value)) {
      p_const_value <- p_const_values[length(p_const_values)]  # Default to last value
    }
    p_const <- p_const_value
  } else {
    p_const <- NULL  # Default for other perturbation types
  }

  # Generate graphs
  graphs <- vector("list", n_graphs)

  # Use a fixed parameter value for all graphs in this population
  # (lambda_alt and p_const are now passed directly to generate_population)
  for (i in 1:n_graphs) {
    graphs[[i]] <- generate_ring_graph(N, lambda, perturbed_nodes,
                                      perturbation_type, lambda_alt, p_const)
  }

  # Add parameter information to the returned list
  attr(graphs, "lambda_alt") <- lambda_alt
  attr(graphs, "p_const") <- p_const

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
#' @param lambda_mult Specific lambda multiplier (for lambda_half, lambda_double)
#' @param p_const_value Specific constant probability value (for const_high, const_low)
#' @return List with comprehensive experiment results
run_single_experiment <- function(N, n_graphs, perturbation_type, rho = 0.03,
                                 alpha = 0.05, n_bootstrap = 1000, seed = NULL,
                                 lambda_mult = NULL, p_const_value = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # ---- Parameter setup and recording ----
  experiment_params <- list(
    N = N,
    n_graphs = n_graphs,
    perturbation_type = perturbation_type,
    rho = rho,
    alpha = alpha,
    n_bootstrap = n_bootstrap,
    seed = seed,
    lambda_mult = lambda_mult,
    p_const_value = p_const_value
  )

  # Find appropriate lambda
  lambda <- find_lambda(N, target_density = 0.25)
  experiment_params$lambda_base <- lambda

  # Determine perturbed nodes
  K <- max(2, ceiling(rho * N))
  if (N >= 100) {
    K <- max(3, ceiling(rho * N))
  }
  perturbed_nodes <- sample(1:N, K)
  experiment_params$perturbed_nodes <- perturbed_nodes
  experiment_params$n_perturbed_nodes <- length(perturbed_nodes)

  # ---- Generate ground truth data ----
  # Get all edges involving perturbed nodes (true critical links)
  true_critical_links <- list()
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (i %in% perturbed_nodes || j %in% perturbed_nodes) {
        true_critical_links[[length(true_critical_links) + 1]] <- c(i, j)
      }
    }
  }
  experiment_params$true_critical_links <- true_critical_links
  experiment_params$n_true_critical_links <- length(true_critical_links)

  # ---- Generate populations ----
  pop_A <- generate_population(n_graphs, N, lambda)

  # Store actual lambda value used for perturbation
  if (perturbation_type %in% c("lambda_half", "lambda_double")) {
    if (is.null(lambda_mult)) {
      lambda_multipliers <- switch(perturbation_type,
                                 "lambda_half" = c(0.5, 0.33, 0.25),
                                 "lambda_double" = c(2, 3, 4),
                                 1)
      lambda_mult <- lambda_multipliers[1]
    }
    experiment_params$lambda_alt <- lambda * lambda_mult
  }

  # Store actual constant probability used
  if (perturbation_type %in% c("const_high", "const_low")) {
    if (is.null(p_const_value)) {
      p_const_values <- switch(perturbation_type,
                             "const_high" = c(0.85, 0.90, 0.95),
                             "const_low" = c(0.05, 0.03, 0.01))
      p_const_value <- p_const_values[length(p_const_values)]
    }
    experiment_params$p_const_value <- p_const_value
  }

  pop_B <- generate_population(n_graphs, N, lambda, perturbed_nodes, perturbation_type,
                              lambda_mult, p_const_value)
  populations <- list(A = pop_A, B = pop_B)

  # ---- Run algorithm and collect detailed results ----
  start_time <- Sys.time()

  # Initialize result container
  result_data <- list(
    parameters = experiment_params,
    timing = list(start_time = start_time),
    algorithm_execution = list(),
    confusion_matrix = list(),
    evaluation = list(),
    success = FALSE
  )

  tryCatch({
    # Set a hook to capture internal algorithm metrics if available
    options(BrainNetTest.debug = TRUE)

    # Run the algorithm
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

    # Capture timing information
    end_time <- Sys.time()
    result_data$timing$end_time <- end_time
    result_data$timing$runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Store algorithm execution details
    result_data$algorithm_execution$global_test_p <- ifelse(is.null(result$debug_info$global_p),
                                                           NA,
                                                           result$debug_info$global_p)
    result_data$algorithm_execution$global_significant <- !is.null(result$critical_edges)
    result_data$algorithm_execution$initial_T <- ifelse(is.null(result$debug_info$initial_T),
                                                       NA,
                                                       result$debug_info$initial_T)

    # Store T value trajectory if available
    if (!is.null(result$debug_info$T_trajectory)) {
      result_data$algorithm_execution$T_trajectory <- result$debug_info$T_trajectory
    }

    # Store number of iterations performed
    result_data$algorithm_execution$n_iterations <- ifelse(is.null(result$debug_info$n_iterations),
                                                          0,
                                                          result$debug_info$n_iterations)

    # Store all edges identified and removed
    predicted_edges <- if (!is.null(result$critical_edges)) {
      lapply(1:nrow(result$critical_edges), function(i) {
        c(result$critical_edges$node1[i], result$critical_edges$node2[i])
      })
    } else {
      list()
    }
    result_data$algorithm_execution$predicted_edges <- predicted_edges
    result_data$algorithm_execution$n_predicted_edges <- length(predicted_edges)

    # Store p-value trajectory if available
    if (!is.null(result$debug_info$p_trajectory)) {
      result_data$algorithm_execution$p_trajectory <- result$debug_info$p_trajectory
    }

    # Compute confusion matrix
    conf_matrix <- compute_confusion_matrix(result$critical_edges, perturbed_nodes, N)
    result_data$confusion_matrix <- conf_matrix

    # Calculate additional evaluation metrics
    result_data$evaluation$recall <- if (conf_matrix$TP + conf_matrix$FN > 0)
                                     conf_matrix$TP / (conf_matrix$TP + conf_matrix$FN) else 0
    result_data$evaluation$precision <- if (conf_matrix$TP + conf_matrix$FP > 0)
                                       conf_matrix$TP / (conf_matrix$TP + conf_matrix$FP) else 1
    result_data$evaluation$f1 <- if (result_data$evaluation$recall + result_data$evaluation$precision > 0)
                                2 * (result_data$evaluation$precision * result_data$evaluation$recall) /
                                (result_data$evaluation$precision + result_data$evaluation$recall) else 0
    result_data$evaluation$accuracy <- (conf_matrix$TP + conf_matrix$TN) /
                                      (conf_matrix$TP + conf_matrix$TN + conf_matrix$FP + conf_matrix$FN)

    # Compute Matthews Correlation Coefficient
    denom <- sqrt((conf_matrix$TP + conf_matrix$FP) * (conf_matrix$TP + conf_matrix$FN) *
                  (conf_matrix$TN + conf_matrix$FP) * (conf_matrix$TN + conf_matrix$FN))
    result_data$evaluation$mcc <- if (denom > 0) (conf_matrix$TP * conf_matrix$TN - conf_matrix$FP * conf_matrix$FN) / denom else 0

    # Store links correctly and incorrectly identified
    true_keys <- sapply(true_critical_links, function(e) paste(sort(e), collapse="-"))
    pred_keys <- sapply(predicted_edges, function(e) paste(sort(e), collapse="-"))

    result_data$evaluation$correctly_identified <- intersect(true_keys, pred_keys)
    result_data$evaluation$incorrectly_identified <- setdiff(pred_keys, true_keys)
    result_data$evaluation$missed_critical <- setdiff(true_keys, pred_keys)

    # Test for residual signal
    residual_result <- NULL
    result_data$evaluation$residual_test_performed <- FALSE

    if (result_data$algorithm_execution$global_significant) {
      tryCatch({
        residual_result <- identify_critical_links(
          populations = result$modified_populations,
          alpha = alpha,
          method = "fisher",
          adjust_method = "none",
          batch_size = 1,
          n_bootstrap = n_bootstrap,
          a = 1,
          seed = seed + 1000
        )
        result_data$evaluation$residual_test_performed <- TRUE
        result_data$evaluation$residual_p <- ifelse(is.null(residual_result$debug_info$global_p),
                                                   1,
                                                   residual_result$debug_info$global_p)
        result_data$evaluation$no_residual <- is.null(residual_result$critical_edges)
      }, error = function(e) {
        result_data$evaluation$residual_test_error <- as.character(e)
      })
    }

    result_data$success <- TRUE

    # Print real-time summary
    cat(sprintf("[N=%d, %s] TP=%d FP=%d FN=%d | Recall=%.2f Prec=%.2f | Runtime=%.2fs\n",
                N, perturbation_type, conf_matrix$TP, conf_matrix$FP, conf_matrix$FN,
                result_data$evaluation$recall, result_data$evaluation$precision,
                result_data$timing$runtime))

    return(result_data)

  }, error = function(e) {
    # Capture error information
    end_time <- Sys.time()
    result_data$timing$end_time <- end_time
    result_data$timing$runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
    result_data$error <- as.character(e)

    # Print real-time error summary
    cat(sprintf("[N=%d, %s] ERROR: %s | Runtime=%.2fs\n",
                N, perturbation_type, substr(result_data$error, 1, 50),
                result_data$timing$runtime))

    return(result_data)
  })
}

#' Create organized directory structure for experiment outputs
#' @param base_dir Base directory for all experiments
#' @param experiment_id Unique experiment identifier
#' @return List with paths to various output directories
create_experiment_structure <- function(base_dir, experiment_id) {
  # Create main experiment directory
  exp_dir <- file.path(base_dir, experiment_id)

  # Create subdirectories
  dirs <- list(
    root = exp_dir,
    metadata = file.path(exp_dir, "metadata"),
    by_N = file.path(exp_dir, "by_network_size"),
    by_perturbation = file.path(exp_dir, "by_perturbation"),
    individual_runs = file.path(exp_dir, "individual_runs"),
    summaries = file.path(exp_dir, "summaries"),
    reports = file.path(exp_dir, "reports"),
    plots = file.path(exp_dir, "plots")
  )

  # Create all directories
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }

  return(dirs)
}

#' Generate descriptive filename for individual run
#' @param N Network size
#' @param pert_type Perturbation type
#' @param rep Repetition number
#' @param param_value Parameter value (lambda_mult or p_const)
#' @param run_id Global run ID
#' @return Formatted filename
generate_run_filename <- function(N, pert_type, rep, param_value = NULL, run_id) {
  # Base filename components
  base <- sprintf("N%05d_%s_rep%03d", N, pert_type, rep)

  # Add parameter value if applicable
  if (!is.null(param_value)) {
    if (pert_type %in% c("lambda_half", "lambda_double")) {
      param_str <- sprintf("_lam%.3f", param_value)
    } else if (pert_type %in% c("const_high", "const_low")) {
      param_str <- sprintf("_p%.3f", param_value)
    } else {
      param_str <- ""
    }
    base <- paste0(base, param_str)
  }

  # Add run ID
  filename <- sprintf("%s_run%05d.rds", base, run_id)
  return(filename)
}

#' Create experiment manifest file
#' @param dirs Directory structure from create_experiment_structure
#' @param experiment_metadata Metadata about the experiment
create_experiment_manifest <- function(dirs, experiment_metadata) {
  manifest <- list(
    created = Sys.time(),
    experiment_metadata = experiment_metadata,
    directory_structure = dirs,
    file_patterns = list(
      individual_runs = "N{size}_{perturbation}_rep{rep}_run{id}.rds",
      summaries = "{level}_summary.csv",
      reports = "{type}_report.{ext}"
    )
  )

  manifest_file <- file.path(dirs$metadata, "experiment_manifest.json")
  jsonlite::write_json(manifest, manifest_file, pretty = TRUE, auto_unbox = TRUE)

  # Also save as RDS for easy R access
  saveRDS(manifest, file.path(dirs$metadata, "experiment_manifest.rds"))
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
#' @param output_dir Directory to save intermediate results
#' @param lambda_multipliers Custom lambda multipliers for lambda-based perturbations
#' @param p_const_values Custom probability values for constant perturbations
#' @return Data frame with all results
run_ring_experiment <- function(N_values = c(10, 100, 1000, 10000),
                               perturbation_types = c("lambda_half", "lambda_double",
                                                     "const_high", "const_low"),
                               n_repetitions = 200,
                               rho = 0.02,
                               alpha = 0.05,
                               n_bootstrap = 1000,
                               verbose = TRUE,
                               n_cores = parallel::detectCores() - 1,
                               output_dir = ".",
                               lambda_multipliers = NULL,
                               p_const_values = NULL) {

  # Define default parameter values if not provided
  if (is.null(lambda_multipliers)) {
    lambda_multipliers <- list(
      "lambda_half" = c(0.5, 0.33, 0.25),
      "lambda_double" = c(2, 3, 4)
    )
  }

  if (is.null(p_const_values)) {
    p_const_values <- list(
      "const_high" = c(0.85, 0.90, 0.95),
      "const_low" = c(0.05, 0.03, 0.01)
    )
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Create timestamp for this experiment run
  experiment_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  experiment_id <- paste0("ring_exp_", experiment_timestamp)

  # Create organized directory structure
  dirs <- create_experiment_structure(output_dir, experiment_id)

  # Create experiment metadata
  experiment_metadata <- list(
    experiment_id = experiment_id,
    timestamp = experiment_timestamp,
    N_values = N_values,
    perturbation_types = perturbation_types,
    n_repetitions = n_repetitions,
    rho = rho,
    alpha = alpha,
    n_bootstrap = n_bootstrap,
    n_cores = n_cores,
    lambda_multipliers = lambda_multipliers,
    p_const_values = p_const_values
  )

  # Save experiment metadata
  metadata_file <- file.path(dirs$metadata, "experiment_config.rds")
  saveRDS(experiment_metadata, metadata_file)

  # Also save as JSON for easy inspection
  jsonlite::write_json(experiment_metadata,
                      file.path(dirs$metadata, "experiment_config.json"),
                      pretty = TRUE, auto_unbox = TRUE)

  # Create experiment manifest
  create_experiment_manifest(dirs, experiment_metadata)

  # Create results tracking data frame
  results_summary <- data.frame()

  # Initialize progress tracking
  experiment_start_time <- Sys.time()
  last_progress_update <- experiment_start_time
  progress_update_interval <- 600  # 10 minutes in seconds

  # Print experiment header
  if (verbose) {
    cat("\n========================================================\n")
    cat("RING EXPERIMENT FOR identify_critical_links() VALIDATION\n")
    cat("========================================================\n")
    cat(sprintf("Experiment ID: %s\n", experiment_id))
    cat(sprintf("Output Directory: %s\n", dirs$root))
    cat(sprintf("Started: %s\n", format(experiment_start_time, "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("Parameters: N=%s, Perturbations=%s, Repetitions=%d\n",
                paste(N_values, collapse=","),
                paste(perturbation_types, collapse=","),
                n_repetitions))
    cat("--------------------------------------------------------\n\n")
  }

  # Build a parameter grid with all combinations
  param_grid <- list()
  idx <- 1

  # Count total parameter combinations for progress tracking
  total_param_combinations <- 0
  for (N in N_values) {
    for (pert_type in perturbation_types) {
      if (pert_type %in% c("lambda_half", "lambda_double")) {
        total_param_combinations <- total_param_combinations + length(lambda_multipliers[[pert_type]])
      } else if (pert_type %in% c("const_high", "const_low")) {
        total_param_combinations <- total_param_combinations + length(p_const_values[[pert_type]])
      } else {
        total_param_combinations <- total_param_combinations + 1
      }
    }
  }
  total_experiments <- total_param_combinations * n_repetitions

  if (verbose) {
    cat(sprintf("Total experiments to run: %d\n", total_experiments))
    cat(sprintf("Estimated experiments per minute (per core): ~%.1f\n",
                ifelse(total_experiments <= 100, 2, 0.5)))
    estimated_minutes <- total_experiments / (n_cores * ifelse(total_experiments <= 100, 2, 0.5))
    cat(sprintf("Estimated runtime: %.1f minutes (%.1f hours)\n\n",
                estimated_minutes, estimated_minutes / 60))
  }

  for (N in N_values) {
    # Create subdirectory for this N value
    N_dir <- file.path(dirs$by_N, sprintf("N_%05d", N))
    if (!dir.exists(N_dir)) {
      dir.create(N_dir, recursive = TRUE)
    }

    # Determine number of graphs based on N for better statistical power
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
      # Create subdirectory for this perturbation type
      pert_dir <- file.path(dirs$by_perturbation, pert_type)
      if (!dir.exists(pert_dir)) {
        dir.create(pert_dir, recursive = TRUE)
      }

      # Get applicable parameter values for this perturbation type
      if (pert_type %in% c("lambda_half", "lambda_double")) {
        param_values <- lambda_multipliers[[pert_type]]
        param_type <- "lambda_mult"
      } else if (pert_type %in% c("const_high", "const_low")) {
        param_values <- p_const_values[[pert_type]]  # Fixed: was perturbation_type
        param_type <- "p_const_value"
      } else {
        param_values <- list(NULL)  # Just one run with default parameters
        param_type <- "none"
      }

      # For each parameter value, do n_repetitions
      for (param_value in param_values) {
        # Print parameter value header
        if (verbose) {
          if (param_type == "lambda_mult") {
            cat(sprintf("\n  Testing %s with λ multiplier = %.2f\n", pert_type, param_value))
          } else if (param_type == "p_const_value") {
            cat(sprintf("\n  Testing %s with p = %.2f\n", pert_type, param_value))
          } else {
            cat(sprintf("\n  Testing %s\n", pert_type))
          }
        }

        for (rep in 1:n_repetitions) {
          seed <- 1000 * idx + rep

          # Create parameter list with appropriate named parameter
          params <- list(
            N = N,
            n_graphs = n_graphs,
            perturbation_type = pert_type,
            rho = rho,
            alpha = alpha,
            n_bootstrap = n_bootstrap,
            seed = seed
          )

          # Add the specific parameter being tested
          if (param_type == "lambda_mult") {
            params$lambda_mult <- param_value
          } else if (param_type == "p_const_value") {
            params$p_const_value <- param_value
          }

          params$experiment_id <- experiment_id
          params$run_id <- idx
          params$rep_id <- rep
          params$dirs <- dirs  # Pass directory structure

          param_grid[[idx]] <- params
          idx <- idx + 1
        }
      }
    }
  }

  # Set up parallel processing
  if (n_cores > 1) {
    if (verbose) {
      cat(sprintf("\nSetting up parallel processing with %d cores\n", n_cores))
    }

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
      "identify_critical_links",
      "generate_run_filename"
    ), envir = environment())

    # Load required libraries on each worker
    parallel::clusterEvalQ(cl, {
      library(Matrix)
      library(igraph)
      suppressPackageStartupMessages(library(BrainNetTest))
    })

    # If identify_critical_links is defined in another file, source it on each worker
    parallel::clusterEvalQ(cl, {
      if (!requireNamespace("BrainNetTest", quietly = TRUE)) {
        source("./R/identify_critical_links.R")
      }
    })

    # Set reproducible parallel RNG
    parallel::clusterSetRNGStream(cl, iseed = 42)

    # Parallel execution
    if (verbose) {
      cat(sprintf("\nRunning %d experiments across %d cores...\n\n",
                  length(param_grid), n_cores))
    }

    # Process parameters in batches for real-time results
    batch_size <- min(20, length(param_grid))  # Changed from 100 to 20
    n_batches <- ceiling(length(param_grid) / batch_size)

    all_results <- list()
    experiments_completed <- 0

    for (batch in 1:n_batches) {
      batch_start_time <- Sys.time()
      start_idx <- (batch - 1) * batch_size + 1
      end_idx <- min(batch * batch_size, length(param_grid))
      batch_experiments <- end_idx - start_idx + 1

      if (verbose && n_batches > 1) {
        cat(sprintf("\n[BATCH %d/%d] Starting experiments %d-%d (batch size: %d)\n",
                    batch, n_batches, start_idx, end_idx, batch_experiments))
      }

      # Extract current batch of parameters
      current_batch <- param_grid[start_idx:end_idx]

      # Store tracking parameters separately but don't pass them to run_single_experiment
      tracking_info <- lapply(current_batch, function(params) {
        list(
          experiment_id = params$experiment_id,
          run_id = params$run_id,
          rep_id = params$rep_id,
          dirs = params$dirs,
          N = params$N,
          perturbation_type = params$perturbation_type,
          lambda_mult = params$lambda_mult,
          p_const_value = params$p_const_value
        )
      })

      # Remove tracking parameters before passing to run_single_experiment
      execution_params <- lapply(current_batch, function(params) {
        params$experiment_id <- NULL
        params$run_id <- NULL
        params$rep_id <- NULL
        params$dirs <- NULL
        params
      })

      # Run current batch in parallel with clean parameters
      batch_results <- parallel::parLapply(
        cl, execution_params,
        function(pars) do.call(run_single_experiment, pars)
      )

      # Process and save results from this batch
      for (i in seq_along(batch_results)) {
        result <- batch_results[[i]]
        # Get tracking info from our saved list instead of from parameters
        info <- tracking_info[[i]]
        run_id <- info$run_id

        # Generate descriptive filename
        filename <- generate_run_filename(
          N = info$N,
          pert_type = info$perturbation_type,
          rep = info$rep_id,
          param_value = ifelse(is.null(info$lambda_mult), info$p_const_value, info$lambda_mult),
          run_id = run_id
        )

        # Save in multiple locations for easy access
        # 1. In individual_runs directory with full filename
        result_file <- file.path(info$dirs$individual_runs, filename)
        saveRDS(result, result_file)

        # 2. In by_N subdirectory
        N_subdir <- file.path(info$dirs$by_N, sprintf("N_%05d", info$N))
        result_file_N <- file.path(N_subdir, filename)
        saveRDS(result, result_file_N)

        # 3. In by_perturbation subdirectory
        pert_subdir <- file.path(info$dirs$by_perturbation, info$perturbation_type)
        result_file_pert <- file.path(pert_subdir, filename)
        saveRDS(result, result_file_pert)

        # Extract summary for the tracking data frame
        if (result$success) {
          summary_row <- data.frame(
            experiment_id = experiment_id,
            run_id = run_id,
            N = result$parameters$N,
            perturbation_type = result$parameters$perturbation_type,
            lambda_base = result$parameters$lambda_base,
            lambda_mult = ifelse(is.null(result$parameters$lambda_mult),
                                NA, result$parameters$lambda_mult),
            p_const = ifelse(is.null(result$parameters$p_const_value),
                            NA, result$parameters$p_const_value),
            global_significant = result$algorithm_execution$global_significant,
            TP = result$confusion_matrix$TP,
            FP = result$confusion_matrix$FP,
            FN = result$confusion_matrix$FN,
            TN = result$confusion_matrix$TN,
            recall = result$evaluation$recall,
            precision = result$evaluation$precision,
            f1 = result$evaluation$f1,
            mcc = result$evaluation$mcc,
            n_predicted = result$algorithm_execution$n_predicted_edges,
            n_true = result$parameters$n_true_critical_links,
            runtime = result$timing$runtime,
            no_residual = ifelse(is.null(result$evaluation$no_residual),
                                NA, result$evaluation$no_residual),
            error = NA
          )
        } else {
          # For failed runs, create a basic summary row
          summary_row <- data.frame(
            experiment_id = experiment_id,
            run_id = run_id,
            N = result$parameters$N,
            perturbation_type = result$parameters$perturbation_type,
            lambda_base = ifelse(is.null(result$parameters$lambda_base),
                                NA, result$parameters$lambda_base),
            lambda_mult = ifelse(is.null(result$parameters$lambda_mult),
                                NA, result$parameters$lambda_mult),
            p_const = ifelse(is.null(result$parameters$p_const_value),
                            NA, result$parameters$p_const_value),
            global_significant = FALSE,
            TP = NA, FP = NA, FN = NA, TN = NA,
            recall = NA, precision = NA, f1 = NA, mcc = NA,
            n_predicted = 0,
            n_true = ifelse(is.null(result$parameters$n_true_critical_links),
                            NA, result$parameters$n_true_critical_links),
            runtime = result$timing$runtime,
            no_residual = NA,
            error = ifelse(is.null(result$error), NA, result$error)
          )
        }

        # Append to the results summary data frame
        results_summary <- rbind(results_summary, summary_row)
      }

      # Update experiment count
      experiments_completed <- experiments_completed + batch_experiments
      batch_runtime <- as.numeric(difftime(Sys.time(), batch_start_time, units = "secs"))

      # Save intermediate summary results in summaries directory
      summary_file <- file.path(dirs$summaries, "running_summary.csv")
      write.csv(results_summary, summary_file, row.names = FALSE)

      # Also save summaries by N and perturbation type
      for (N in unique(results_summary$N)) {
        N_summary <- subset(results_summary, N == N)
        N_summary_file <- file.path(dirs$summaries, sprintf("summary_N_%05d.csv", N))
        write.csv(N_summary, N_summary_file, row.names = FALSE)
      }

      for (pert in unique(results_summary$perturbation_type)) {
        pert_summary <- subset(results_summary, perturbation_type == pert)
        pert_summary_file <- file.path(dirs$summaries, sprintf("summary_%s.csv", pert))
        write.csv(pert_summary, pert_summary_file, row.names = FALSE)
      }

      # Print batch completion summary
      successful <- sum(!is.na(results_summary$global_significant))
      significant <- sum(results_summary$global_significant, na.rm = TRUE)
      avg_recall <- mean(results_summary$recall, na.rm = TRUE)
      avg_precision <- mean(results_summary$precision, na.rm = TRUE)

      if (verbose && n_batches > 1) {
        cat(sprintf("\n[BATCH %d/%d COMPLETE] Runtime: %.1f seconds (%.2f sec/experiment)\n",
                    batch, n_batches, batch_runtime, batch_runtime / batch_experiments))
        cat(sprintf("  Progress: %d/%d experiments (%.1f%%) completed\n",
                    experiments_completed, total_experiments,
                    100 * experiments_completed / total_experiments))
        cat(sprintf("  Cumulative: %.1f%% success, %.1f%% power, Avg recall=%.3f, Avg precision=%.3f\n",
                    100 * successful / nrow(results_summary),
                    100 * significant / successful,
                    avg_recall, avg_precision))

        # Estimate time remaining
        total_elapsed <- as.numeric(difftime(Sys.time(), experiment_start_time, units = "mins"))
        experiments_per_minute <- experiments_completed / total_elapsed
        minutes_remaining <- (total_experiments - experiments_completed) / experiments_per_minute

        cat(sprintf("  Time elapsed: %.1f minutes | Estimated remaining: %.1f minutes (%.1f hours)\n",
                    total_elapsed, minutes_remaining, minutes_remaining / 60))
        cat(sprintf("  Estimated completion: %s\n",
                    format(Sys.time() + minutes_remaining * 60, "%Y-%m-%d %H:%M:%S")))
      }

      # Check if it's time for a periodic update
      current_time <- Sys.time()
      time_since_last_update <- as.numeric(difftime(current_time, last_progress_update, units = "secs"))

      if (verbose && time_since_last_update >= progress_update_interval) {
        cat("\n--------------------------------------------------------\n")
        cat(sprintf("[PERIODIC UPDATE] %s\n", format(current_time, "%Y-%m-%d %H:%M:%S")))
        cat(sprintf("  Total progress: %d/%d experiments (%.1f%%)\n",
                    experiments_completed, total_experiments,
                    100 * experiments_completed / total_experiments))

        # Performance metrics
        total_elapsed_mins <- as.numeric(difftime(current_time, experiment_start_time, units = "mins"))
        cat(sprintf("  Runtime: %.1f minutes (%.1f hours)\n",
                    total_elapsed_mins, total_elapsed_mins / 60))
        cat(sprintf("  Average speed: %.2f experiments/minute\n",
                    experiments_completed / total_elapsed_mins))

        # Current status
        if (experiments_completed < total_experiments) {
          current_N <- param_grid[[min(experiments_completed + 1, length(param_grid))]]$N
          current_pert <- param_grid[[min(experiments_completed + 1, length(param_grid))]]$perturbation_type
          cat(sprintf("  Currently processing: N=%d, %s\n", current_N, current_pert))
        }

        cat("--------------------------------------------------------\n\n")
        last_progress_update <- current_time
      }

      # Store batch results
      all_results <- c(all_results, batch_results)
    }
  } else {
    # Sequential execution (for debugging or single-core machines)
    if (verbose) {
      cat("\nRunning experiments sequentially...\n\n")
    }

    all_results <- list()
    experiments_completed <- 0

    for (i in seq_along(param_grid)) {
      experiment_start <- Sys.time()

      # Store tracking info separately
      tracking_info <- list(
        experiment_id = param_grid[[i]]$experiment_id,
        run_id = param_grid[[i]]$run_id,
        rep_id = param_grid[[i]]$rep_id,
        dirs = param_grid[[i]]$dirs,
        N = param_grid[[i]]$N,
        perturbation_type = param_grid[[i]]$perturbation_type,
        lambda_mult = param_grid[[i]]$lambda_mult,
        p_const_value = param_grid[[i]]$p_const_value
      )

      # Remove tracking params
      execution_params <- param_grid[[i]]
      execution_params$experiment_id <- NULL
      execution_params$run_id <- NULL
      execution_params$rep_id <- NULL
      execution_params$dirs <- NULL

      # Execute with clean parameters
      result <- do.call(run_single_experiment, execution_params)
      all_results[[i]] <- result

      # Generate descriptive filename
      run_id <- tracking_info$run_id
      filename <- generate_run_filename(
        N = tracking_info$N,
        pert_type = tracking_info$perturbation_type,
        rep = tracking_info$rep_id,
        param_value = ifelse(is.null(tracking_info$lambda_mult),
                           tracking_info$p_const_value,
                           tracking_info$lambda_mult),
        run_id = run_id
      )

      # Save in multiple locations
      # 1. In individual_runs directory
      result_file <- file.path(tracking_info$dirs$individual_runs, filename)
      saveRDS(result, result_file)

      # 2. In by_N subdirectory
      N_subdir <- file.path(tracking_info$dirs$by_N, sprintf("N_%05d", tracking_info$N))
      if (!dir.exists(N_subdir)) dir.create(N_subdir, recursive = TRUE)
      result_file_N <- file.path(N_subdir, filename)
      saveRDS(result, result_file_N)

      # 3. In by_perturbation subdirectory
      pert_subdir <- file.path(tracking_info$dirs$by_perturbation, tracking_info$perturbation_type)
      if (!dir.exists(pert_subdir)) dir.create(pert_subdir, recursive = TRUE)
      result_file_pert <- file.path(pert_subdir, filename)
      saveRDS(result, result_file_pert)

      # Extract and save summary info
      if (result$success) {
        summary_row <- data.frame(
          experiment_id = experiment_id,
          run_id = run_id,
          N = result$parameters$N,
          perturbation_type = result$parameters$perturbation_type,
          lambda_base = result$parameters$lambda_base,
          lambda_mult = ifelse(is.null(result$parameters$lambda_mult),
                              NA, result$parameters$lambda_mult),
          p_const = ifelse(is.null(result$parameters$p_const_value),
                          NA, result$parameters$p_const_value),
          global_significant = result$algorithm_execution$global_significant,
          TP = result$confusion_matrix$TP,
          FP = result$confusion_matrix$FP,
          FN = result$confusion_matrix$FN,
          TN = result$confusion_matrix$TN,
          recall = result$evaluation$recall,
          precision = result$evaluation$precision,
          f1 = result$evaluation$f1,
          mcc = result$evaluation$mcc,
          n_predicted = result$algorithm_execution$n_predicted_edges,
          n_true = result$parameters$n_true_critical_links,
          runtime = result$timing$runtime,
          no_residual = ifelse(is.null(result$evaluation$no_residual),
                              NA, result$evaluation$no_residual),
          error = NA
        )
      } else {
        summary_row <- data.frame(
          experiment_id = experiment_id,
          run_id = run_id,
          N = result$parameters$N,
          perturbation_type = result$parameters$perturbation_type,
          lambda_base = ifelse(is.null(result$parameters$lambda_base),
                              NA, result$parameters$lambda_base),
          lambda_mult = ifelse(is.null(result$parameters$lambda_mult),
                              NA, result$parameters$lambda_mult),
          p_const = ifelse(is.null(result$parameters$p_const_value),
                          NA, result$parameters$p_const_value),
          global_significant = FALSE,
          TP = NA, FP = NA, FN = NA, TN = NA,
          recall = NA, precision = NA, f1 = NA, mcc = NA,
          n_predicted = 0,
          n_true = ifelse(is.null(result$parameters$n_true_critical_links),
                          NA, result$parameters$n_true_critical_links),
          runtime = result$timing$runtime,
          no_residual = NA,
          error = ifelse(is.null(result$error), NA, result$error)
        )
      }

      # Append to results summary
      results_summary <- rbind(results_summary, summary_row)
      experiments_completed <- experiments_completed + 1

      # Save intermediate summary results
      summary_file <- file.path(dirs$summaries, "running_summary.csv")
      write.csv(results_summary, summary_file, row.names = FALSE)

      # Check for periodic updates in sequential mode
      current_time <- Sys.time()
      time_since_last_update <- as.numeric(difftime(current_time, last_progress_update, units = "secs"))

      if (verbose && (experiments_completed %% 10 == 0 || time_since_last_update >= progress_update_interval)) {
        cat(sprintf("\n[PROGRESS] %d/%d experiments completed (%.1f%%)\n",
                    experiments_completed, total_experiments,
                    100 * experiments_completed / total_experiments))

        total_elapsed <- as.numeric(difftime(current_time, experiment_start_time, units = "mins"))
        experiments_per_minute <- experiments_completed / total_elapsed
        minutes_remaining <- (total_experiments - experiments_completed) / experiments_per_minute

        cat(sprintf("  Time elapsed: %.1f minutes | Estimated remaining: %.1f minutes\n",
                    total_elapsed, minutes_remaining))

        if (time_since_last_update >= progress_update_interval) {
          last_progress_update <- current_time
        }
      }
    }
  }

  # Create final consolidated results
  if (verbose) {
    cat("\n========================================================\n")
    cat("EXPERIMENT COMPLETED!\n")
    cat("========================================================\n")
    total_runtime <- as.numeric(difftime(Sys.time(), experiment_start_time, units = "mins"))
    cat(sprintf("Total runtime: %.1f minutes (%.2f hours)\n", total_runtime, total_runtime / 60))
    cat(sprintf("Average speed: %.2f experiments/minute\n", total_experiments / total_runtime))
    cat(sprintf("Output directory: %s\n", dirs$root))
    cat("\nFinalizing results...\n")
  }

  # Save complete summary
  final_summary_file <- file.path(dirs$summaries, "complete_summary.csv")
  write.csv(results_summary, final_summary_file, row.names = FALSE)

  # Create experiment report
  report_file <- file.path(dirs$reports, "experiment_report.txt")
  sink(report_file)

  cat("=================================================================\n")
  cat("RING EXPERIMENT SUMMARY REPORT\n")
  cat("=================================================================\n")
  cat(sprintf("Experiment ID: %s\n", experiment_id))
  cat(sprintf("Started: %s\n", format(as.POSIXct(experiment_timestamp, format="%Y%m%d_%H%M%S"), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("Completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat("=================================================================\n\n")

  # Summarize by N and perturbation type
  cat("PERFORMANCE BY CONDITION:\n")
  cat("-------------------------------------------------------------------\n")

  for (N in N_values) {
    cat(sprintf("\nResults for N = %d:\n", N))
    cat("-------------------------------------------------------------------\n")

    for (pert_type in perturbation_types) {
      condition_results <- subset(results_summary,
                                  N == N & perturbation_type == pert_type)

      if (nrow(condition_results) > 0) {
        successful <- sum(!is.na(condition_results$global_significant))
        significant <- sum(condition_results$global_significant, na.rm = TRUE)
        power <- if (successful > 0) significant / successful else 0

        avg_recall <- mean(condition_results$recall, na.rm = TRUE)
        avg_precision <- mean(condition_results$precision, na.rm = TRUE)
        avg_f1 <- mean(condition_results$f1, na.rm = TRUE)

        median_recall <- median(condition_results$recall, na.rm = TRUE)
        median_precision <- median(condition_results$precision, na.rm = TRUE)
        median_f1 <- median(condition_results$f1, na.rm = TRUE)

        avg_runtime <- mean(condition_results$runtime, na.rm = TRUE)

        cat(sprintf("  %s (n=%d):\n", pert_type, nrow(condition_results)))
        cat(sprintf("    Success rate: %.1f%% (%d/%d)\n",
                    100 * successful / nrow(condition_results),
                    successful, nrow(condition_results)))
        cat(sprintf("    Power: %.3f (%d/%d)\n",
                    power, significant, successful))
        cat(sprintf("    Recall: %.3f (median: %.3f)\n",
                    avg_recall, median_recall))
        cat(sprintf("    Precision: %.3f (median: %.3f)\n",
                    avg_precision, median_precision))
        cat(sprintf("    F1 Score: %.3f (median: %.3f)\n",
                    avg_f1, median_f1))
        cat(sprintf("    Average Runtime: %.2f seconds\n", avg_runtime))
        cat("\n")
      }
    }
  }

  # Overall summary
  cat("=================================================================\n")
  cat("OVERALL EXPERIMENT SUMMARY:\n")
  cat("=================================================================\n")

  total_runs <- nrow(results_summary)
  successful_runs <- sum(!is.na(results_summary$global_significant))
  failed_runs <- total_runs - successful_runs
  significant_runs <- sum(results_summary$global_significant, na.rm = TRUE)
  overall_power <- if (successful_runs > 0) significant_runs / successful_runs else 0

  cat(sprintf("Total experiments: %d\n", total_runs))
  cat(sprintf("Successful runs: %d (%.1f%%)\n",
              successful_runs, 100 * successful_runs / total_runs))
  cat(sprintf("Failed runs: %d (%.1f%%)\n",
              failed_runs, 100 * failed_runs / total_runs))
  cat(sprintf("Overall power: %.3f (%d/%d)\n",
              overall_power, significant_runs, successful_runs))

  # Performance metrics
  avg_recall <- mean(results_summary$recall, na.rm = TRUE)
  avg_precision <- mean(results_summary$precision, na.rm = TRUE)
  avg_f1 <- mean(results_summary$f1, na.rm = TRUE)

  cat(sprintf("Average Recall: %.3f\n", avg_recall))
  cat(sprintf("Average Precision: %.3f\n", avg_precision))
  cat(sprintf("Average F1 Score: %.3f\n", avg_f1))

  # Runtime information
  avg_runtime <- mean(results_summary$runtime, na.rm = TRUE)
  med_runtime <- median(results_summary$runtime, na.rm = TRUE)
  max_runtime <- max(results_summary$runtime, na.rm = TRUE)

  cat(sprintf("Runtime - Average: %.2f seconds, Median: %.2f seconds, Max: %.2f seconds\n",
              avg_runtime, med_runtime, max_runtime))

  sink()

  # Save performance metrics as JSON for easy parsing
  performance_metrics <- list()
  for (N in N_values) {
    performance_metrics[[paste0("N_", N)]] <- list()
    for (pert_type in perturbation_types) {
      condition_results <- subset(results_summary,
                                  N == N & perturbation_type == pert_type)
      if (nrow(condition_results) > 0) {
        successful <- sum(!is.na(condition_results$global_significant))
        significant <- sum(condition_results$global_significant, na.rm = TRUE)

        performance_metrics[[paste0("N_", N)]][[pert_type]] <- list(
          n_runs = nrow(condition_results),
          success_rate = successful / nrow(condition_results),
          power = if (successful > 0) significant / successful else 0,
          recall = list(
            mean = mean(condition_results$recall, na.rm = TRUE),
            median = median(condition_results$recall, na.rm = TRUE),
            sd = sd(condition_results$recall, na.rm = TRUE)
          ),
          precision = list(
            mean = mean(condition_results$precision, na.rm = TRUE),
            median = median(condition_results$precision, na.rm = TRUE),
            sd = sd(condition_results$precision, na.rm = TRUE)
          ),
          f1 = list(
            mean = mean(condition_results$f1, na.rm = TRUE),
            median = median(condition_results$f1, na.rm = TRUE),
            sd = sd(condition_results$f1, na.rm = TRUE)
          ),
          runtime = list(
            mean = mean(condition_results$runtime, na.rm = TRUE),
            median = median(condition_results$runtime, na.rm = TRUE),
            max = max(condition_results$runtime, na.rm = TRUE)
          )
        )
      }
    }
  }

  jsonlite::write_json(performance_metrics,
                      file.path(dirs$reports, "performance_metrics.json"),
                      pretty = TRUE, auto_unbox = TRUE)

  # Create an index file for easy navigation
  index_file <- file.path(dirs$root, "index.txt")
  cat("EXPERIMENT INDEX\n", file = index_file)
  cat("================\n\n", file = index_file, append = TRUE)
  cat(sprintf("Experiment ID: %s\n", experiment_id), file = index_file, append = TRUE)
  cat(sprintf("Created: %s\n\n", format(experiment_start_time, "%Y-%m-%d %H:%M:%S")),
      file = index_file, append = TRUE)
  cat("Directory Structure:\n", file = index_file, append = TRUE)
  cat("-------------------\n", file = index_file, append = TRUE)
  cat("metadata/          - Experiment configuration and manifest\n", file = index_file, append = TRUE)
  cat("by_network_size/   - Results organized by N value\n", file = index_file, append = TRUE)
  cat("by_perturbation/   - Results organized by perturbation type\n", file = index_file, append = TRUE)
  cat("individual_runs/   - All individual run results\n", file = index_file, append = TRUE)
  cat("summaries/         - CSV summaries at various levels\n", file = index_file, append = TRUE)
  cat("reports/           - Text and JSON reports\n", file = index_file, append = TRUE)
  cat("plots/             - Directory for generated plots\n", file = index_file, append = TRUE)
  cat("\nKey Files:\n", file = index_file, append = TRUE)
  cat("----------\n", file = index_file, append = TRUE)
  cat("metadata/experiment_manifest.json  - Complete experiment documentation\n", file = index_file, append = TRUE)
  cat("summaries/complete_summary.csv     - All results in one CSV\n", file = index_file, append = TRUE)
  cat("reports/experiment_report.txt      - Human-readable summary\n", file = index_file, append = TRUE)
  cat("reports/performance_metrics.json   - Machine-readable metrics\n", file = index_file, append = TRUE)

  if (verbose) {
    cat(sprintf("\nDetailed report saved to: %s\n", report_file))
    cat(sprintf("Summary data saved to: %s\n", final_summary_file))
    cat(sprintf("Performance metrics saved to: %s\n",
                file.path(dirs$reports, "performance_metrics.json")))
    cat(sprintf("Experiment root directory: %s\n", dirs$root))
    cat("\nExperiment completed successfully!\n")
  }

  # Add directory paths to results for downstream use
  attr(results_summary, "experiment_dirs") <- dirs
  attr(results_summary, "experiment_id") <- experiment_id

  return(results_summary)
}

#' Main function to run the complete validation
#' @param quick_test If TRUE, run a reduced test for development
#' @param nodes Size of network to test (10, 100, 1000, or 10000)
#' @param master_seed Master seed for reproducibility (default: 42)
#' @param output_csv If TRUE, save results to CSV files (default: TRUE)
#' @param output_dir Directory to save CSV files (default: current directory)
#' @param perturbation_types Perturbation types to test
#' @param lambda_multipliers Custom lambda multipliers for lambda-based perturbations
#' @param p_const_values Custom probability values for constant perturbations
#' @param ... Additional arguments passed to run_ring_experiment
#' @return List with experiment results
main_validation <- function(quick_test = FALSE, nodes = 10000, master_seed = 42,
                           output_csv = TRUE, output_dir = "ring_experiments",
                           perturbation_types = NULL,
                           lambda_multipliers = NULL,
                           p_const_values = NULL,
                           ...) {

  # Set master seed for reproducibility
  set.seed(master_seed)
  cat(sprintf("Using master seed: %d\n", master_seed))

  # Create main output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Configure N_values based on the nodes parameter
  if (quick_test) {
    if (nodes == 10) {
      N_values <- c(10)
    } else if (nodes == 100) {
      N_values <- c(10, 100)
    } else if (nodes == 1000) {
      N_values <- c(1000)
    } else {
      N_values <- c(10000)
    }

    # Default perturbation types for quick test if not specified
    if (is.null(perturbation_types)) {
      perturbation_types <- c("lambda_half", "const_high")
    }

    n_repetitions <- 10

    cat("Running quick test...\n")
  } else {
    # Full validation
    N_values <- c(10, 100, 1000, 10000)

    # Default perturbation types for full validation if not specified
    if (is.null(perturbation_types)) {
      perturbation_types <- c("lambda_half", "lambda_double", "const_high", "const_low")
    }

    n_repetitions <- 200

    cat("Running full validation...\n")
  }

  # Run the experiment
  results <- run_ring_experiment(
    N_values = N_values,
    perturbation_types = perturbation_types,
    n_repetitions = n_repetitions,
    lambda_multipliers = lambda_multipliers,
    p_const_values = p_const_values,
    verbose = TRUE,
    output_dir = output_dir,
    ...
  )

  cat("\nExperiment completed!\n")

  return(list(
    results = results,
    master_seed = master_seed,
    N_values = N_values,
    perturbation_types = perturbation_types,
    quick_test = quick_test
  ))
}

# Example usage:
# Basic quick test with default parameters (nodes=100)
#validation_results <- main_validation(quick_test = TRUE, nodes = 100, master_seed = 42,
#                                    output_dir = "ring_experiments")

# Quick test with custom lambda multipliers for lambda_half perturbation
#custom_lambda_test <- main_validation(
#  quick_test = TRUE,
#  nodes = 100,
#  master_seed = 43,
#  perturbation_types = "lambda_half",
#  lambda_multipliers = list("lambda_half" = c(0.6, 0.4, 0.2))  # Custom multipliers
#)

# Quick test with custom constant probabilities for const_high perturbation
#custom_prob_test <- main_validation(
#  quick_test = TRUE,
#  nodes = 10,
#  master_seed = 44,
#  perturbation_types = "const_high",
#  p_const_values = list("const_high" = c(0.80, 0.90, 0.99))  # Custom probabilities
#)

# Full validation with all perturbation types and default parameters
# Uncomment to run the full validation (takes a long time)
 full_validation <- main_validation(
   master_seed = 45,
   output_dir = "results/full_validation"
 )

# Full validation with custom parameters for all perturbation types
# Uncomment to run a comprehensive test
# comprehensive_test <- main_validation(
#   master_seed = 46,
#   N_values = c(50, 500),  # Custom network sizes
#   n_repetitions = 20,     # Fewer repetitions for faster results
#   lambda_multipliers = list(
#     "lambda_half" = c(0.5, 0.25),
#     "lambda_double" = c(2, 4)
#   ),
#   p_const_values = list(
#     "const_high" = c(0.85, 0.95),
#     "const_low" = c(0.05, 0.01)
#   ),
#   output_dir = "results/custom_params"
# )
