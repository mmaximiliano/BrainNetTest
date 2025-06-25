#' Analyze ring experiment results from organized directory structure
#' @param experiment_dir Path to the experiment directory
#' @param output_plots Whether to generate and save plots
#' @return List with analysis results and generated plots
analyze_ring_experiment <- function(experiment_dir, output_plots = TRUE) {
  
  # Verify directory structure
  if (!dir.exists(experiment_dir)) {
    stop("Experiment directory does not exist: ", experiment_dir)
  }
  
  # Define expected subdirectories
  dirs <- list(
    metadata = file.path(experiment_dir, "metadata"),
    summary = file.path(experiment_dir, "summary"),
    raw_results = file.path(experiment_dir, "raw_results"),
    by_n = file.path(experiment_dir, "results_by_n"),
    by_perturbation = file.path(experiment_dir, "results_by_perturbation"),
    reports = file.path(experiment_dir, "reports"),
    plots = file.path(experiment_dir, "plots")
  )
  
  # Load metadata
  metadata_file <- file.path(dirs$metadata, "experiment_config.rds")
  if (!file.exists(metadata_file)) {
    stop("Metadata file not found")
  }
  metadata <- readRDS(metadata_file)
  
  # Load summary data
  summary_file <- file.path(dirs$summary, "all_results.csv")
  if (!file.exists(summary_file)) {
    stop("Summary file not found")
  }
  results <- read.csv(summary_file, stringsAsFactors = FALSE)
  
  # Create analysis results list
  analysis <- list(
    metadata = metadata,
    summary = results,
    by_n = list(),
    by_perturbation = list(),
    overall = list()
  )
  
  # Analyze by N
  for (N in unique(results$N)) {
    n_results <- results[results$N == N, ]
    
    analysis$by_n[[as.character(N)]] <- list(
      n_experiments = nrow(n_results),
      success_rate = mean(!is.na(n_results$global_significant)),
      power = mean(n_results$global_significant, na.rm = TRUE),
      recall = list(
        mean = mean(n_results$recall, na.rm = TRUE),
        median = median(n_results$recall, na.rm = TRUE),
        sd = sd(n_results$recall, na.rm = TRUE)
      ),
      precision = list(
        mean = mean(n_results$precision, na.rm = TRUE),
        median = median(n_results$precision, na.rm = TRUE),
        sd = sd(n_results$precision, na.rm = TRUE)
      ),
      runtime = list(
        mean = mean(n_results$runtime, na.rm = TRUE),
        median = median(n_results$runtime, na.rm = TRUE),
        max = max(n_results$runtime, na.rm = TRUE)
      )
    )
  }
  
  # Analyze by perturbation type
  for (pert in unique(results$perturbation_type)) {
    pert_results <- results[results$perturbation_type == pert, ]
    
    analysis$by_perturbation[[pert]] <- list(
      n_experiments = nrow(pert_results),
      success_rate = mean(!is.na(pert_results$global_significant)),
      power = mean(pert_results$global_significant, na.rm = TRUE),
      recall = list(
        mean = mean(pert_results$recall, na.rm = TRUE),
        median = median(pert_results$recall, na.rm = TRUE),
        sd = sd(pert_results$recall, na.rm = TRUE)
      ),
      precision = list(
        mean = mean(pert_results$precision, na.rm = TRUE),
        median = median(pert_results$precision, na.rm = TRUE),
        sd = sd(pert_results$precision, na.rm = TRUE)
      )
    )
  }
  
  # Generate plots if requested
  if (output_plots) {
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    
    # Power by N and perturbation type
    power_plot <- results %>%
      group_by(N, perturbation_type) %>%
      summarise(power = mean(global_significant, na.rm = TRUE),
                .groups = "drop") %>%
      ggplot(aes(x = factor(N), y = power, fill = perturbation_type)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = "Statistical Power by Network Size and Perturbation Type",
           x = "Network Size (N)",
           y = "Power",
           fill = "Perturbation") +
      theme_minimal() +
      scale_y_continuous(limits = c(0, 1))
    
    ggsave(file.path(dirs$plots, "power_by_condition.png"), 
           power_plot, width = 10, height = 6)
    
    # Recall vs Precision scatter
    recall_precision_plot <- results %>%
      filter(!is.na(recall) & !is.na(precision)) %>%
      ggplot(aes(x = recall, y = precision, color = perturbation_type)) +
      geom_point(alpha = 0.6) +
      facet_wrap(~N, scales = "free") +
      labs(title = "Recall vs Precision by Network Size",
           x = "Recall",
           y = "Precision") +
      theme_minimal()
    
    ggsave(file.path(dirs$plots, "recall_vs_precision.png"), 
           recall_precision_plot, width = 12, height = 8)
    
    # Runtime distribution
    runtime_plot <- results %>%
      ggplot(aes(x = factor(N), y = runtime)) +
      geom_boxplot(aes(fill = perturbation_type)) +
      scale_y_log10() +
      labs(title = "Runtime Distribution by Network Size",
           x = "Network Size (N)",
           y = "Runtime (seconds, log scale)") +
      theme_minimal()
    
    ggsave(file.path(dirs$plots, "runtime_distribution.png"), 
           runtime_plot, width = 10, height = 6)
    
    # F1 score heatmap
    f1_summary <- results %>%
      group_by(N, perturbation_type) %>%
      summarise(mean_f1 = mean(f1, na.rm = TRUE),
                .groups = "drop")
    
    f1_heatmap <- ggplot(f1_summary, aes(x = factor(N), y = perturbation_type)) +
      geom_tile(aes(fill = mean_f1)) +
      geom_text(aes(label = sprintf("%.3f", mean_f1))) +
      scale_fill_gradient2(low = "red", mid = "yellow", high = "green", 
                          midpoint = 0.5, limits = c(0, 1)) +
      labs(title = "Mean F1 Score Heatmap",
           x = "Network Size (N)",
           y = "Perturbation Type",
           fill = "Mean F1") +
      theme_minimal()
    
    ggsave(file.path(dirs$plots, "f1_heatmap.png"), 
           f1_heatmap, width = 8, height = 6)
  }
  
  # Save analysis results
  analysis_file <- file.path(dirs$reports, "analysis_results.rds")
  saveRDS(analysis, analysis_file)
  
  # Generate analysis report
  report_file <- file.path(dirs$reports, "analysis_report.txt")
  sink(report_file)
  
  cat("RING EXPERIMENT ANALYSIS REPORT\n")
  cat("================================\n\n")
  cat(sprintf("Experiment ID: %s\n", metadata$experiment_id))
  cat(sprintf("Analysis Date: %s\n\n", Sys.Date()))
  
  cat("OVERALL RESULTS:\n")
  cat("----------------\n")
  cat(sprintf("Total experiments: %d\n", nrow(results)))
  cat(sprintf("Overall success rate: %.1f%%\n", 
              100 * mean(!is.na(results$global_significant))))
  cat(sprintf("Overall power: %.3f\n", 
              mean(results$global_significant, na.rm = TRUE)))
  
  cat("\nRESULTS BY NETWORK SIZE:\n")
  cat("------------------------\n")
  for (N in names(analysis$by_n)) {
    n_analysis <- analysis$by_n[[N]]
    cat(sprintf("\nN = %s:\n", N))
    cat(sprintf("  Experiments: %d\n", n_analysis$n_experiments))
    cat(sprintf("  Power: %.3f\n", n_analysis$power))
    cat(sprintf("  Recall: %.3f ± %.3f\n", 
                n_analysis$recall$mean, n_analysis$recall$sd))
    cat(sprintf("  Precision: %.3f ± %.3f\n", 
                n_analysis$precision$mean, n_analysis$precision$sd))
    cat(sprintf("  Runtime: %.2fs (median), %.2fs (max)\n", 
                n_analysis$runtime$median, n_analysis$runtime$max))
  }
  
  cat("\nRESULTS BY PERTURBATION TYPE:\n")
  cat("-----------------------------\n")
  for (pert in names(analysis$by_perturbation)) {
    pert_analysis <- analysis$by_perturbation[[pert]]
    cat(sprintf("\n%s:\n", pert))
    cat(sprintf("  Experiments: %d\n", pert_analysis$n_experiments))
    cat(sprintf("  Power: %.3f\n", pert_analysis$power))
    cat(sprintf("  Recall: %.3f ± %.3f\n", 
                pert_analysis$recall$mean, pert_analysis$recall$sd))
    cat(sprintf("  Precision: %.3f ± %.3f\n", 
                pert_analysis$precision$mean, pert_analysis$precision$sd))
  }
  
  sink()
  
  cat(sprintf("\nAnalysis complete. Results saved to:\n"))
  cat(sprintf("  - Analysis data: %s\n", analysis_file))
  cat(sprintf("  - Analysis report: %s\n", report_file))
  if (output_plots) {
    cat(sprintf("  - Plots directory: %s\n", dirs$plots))
  }
  
  return(analysis)
}

#' Load specific experiment runs for detailed analysis
#' @param experiment_dir Path to the experiment directory
#' @param N Network size to analyze
#' @param perturbation_type Perturbation type to analyze
#' @param max_runs Maximum number of runs to load
#' @return List of detailed run results
load_experiment_runs <- function(experiment_dir, N = NULL, 
                                perturbation_type = NULL, 
                                max_runs = 100) {
  
  # Determine which directory to load from
  if (!is.null(N) && is.null(perturbation_type)) {
    # Load from N-specific directory
    run_dir <- file.path(experiment_dir, "results_by_n", sprintf("N_%05d", N))
  } else if (is.null(N) && !is.null(perturbation_type)) {
    # Load from perturbation-specific directory
    run_dir <- file.path(experiment_dir, "results_by_perturbation", perturbation_type)
  } else {
    # Load from raw results
    run_dir <- file.path(experiment_dir, "raw_results")
  }
  
  if (!dir.exists(run_dir)) {
    stop("Results directory not found: ", run_dir)
  }
  
  # Get list of result files
  result_files <- list.files(run_dir, pattern = "\\.rds$", full.names = TRUE)
  
  # Limit number of files if specified
  if (length(result_files) > max_runs) {
    result_files <- result_files[1:max_runs]
  }
  
  # Load results
  results <- lapply(result_files, readRDS)
  
  # Filter if both N and perturbation_type are specified
  if (!is.null(N) && !is.null(perturbation_type)) {
    results <- Filter(function(r) {
      r$parameters$N == N && r$parameters$perturbation_type == perturbation_type
    }, results)
  }
  
  return(results)
}
