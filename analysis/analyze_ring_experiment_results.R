# analyze_ring_experiment_results.R
# Comprehensive analysis of ring experiment results for identify_critical_links

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)
library(viridis)

# Function to read and combine all summary files
read_experiment_results <- function(results_dir) {
  # Find all summary CSV files
  summary_files <- list.files(
    path = file.path(results_dir, "summaries"),
    pattern = "^summary_.*\\.csv$",
    full.names = TRUE
  )
  
  # Read and combine all files
  all_results <- map_df(summary_files, read_csv, show_col_types = FALSE)
  
  # Calculate actual lambda values
  all_results <- all_results %>%
    mutate(
      lambda_actual = lambda_base * lambda_mult,
      # Create a unique identifier for each parameter combination
      param_id = paste(N, perturbation_type, lambda_actual, sep = "_")
    )
  
  return(all_results)
}

# Function to aggregate results across multiple runs
aggregate_results <- function(results_df) {
  results_df %>%
    group_by(N, perturbation_type, lambda_actual) %>%
    summarise(
      # Average metrics
      TP_mean = mean(TP, na.rm = TRUE),
      FP_mean = mean(FP, na.rm = TRUE),
      FN_mean = mean(FN, na.rm = TRUE),
      TN_mean = mean(TN, na.rm = TRUE),
      recall_mean = mean(recall, na.rm = TRUE),
      precision_mean = mean(precision, na.rm = TRUE),
      f1_mean = mean(f1, na.rm = TRUE),
      mcc_mean = mean(mcc, na.rm = TRUE),
      n_predicted_mean = mean(n_predicted, na.rm = TRUE),
      n_true_mean = mean(n_true, na.rm = TRUE),
      runtime_mean = mean(runtime, na.rm = TRUE),
      
      # Standard deviations
      TP_sd = sd(TP, na.rm = TRUE),
      FP_sd = sd(FP, na.rm = TRUE),
      FN_sd = sd(FN, na.rm = TRUE),
      TN_sd = sd(TN, na.rm = TRUE),
      recall_sd = sd(recall, na.rm = TRUE),
      precision_sd = sd(precision, na.rm = TRUE),
      f1_sd = sd(f1, na.rm = TRUE),
      mcc_sd = sd(mcc, na.rm = TRUE),
      n_predicted_sd = sd(n_predicted, na.rm = TRUE),
      runtime_sd = sd(runtime, na.rm = TRUE),
      
      # Count of runs
      n_runs = n(),
      .groups = "drop"
    ) %>%
    mutate(
      # Calculate derived metrics
      predicted_ratio = n_predicted_mean / n_true_mean,
      tp_fn_tradeoff = TP_mean - FN_mean
    )
}

# Function to create plots for a single N value
create_plots_for_n <- function(data_n, n_value, output_dir) {
  # Filter for lambda_half and lambda_double only
  data_n <- data_n %>%
    filter(perturbation_type %in% c("lambda_half", "lambda_double"))
  
  # Define color palette
  color_palette <- scale_color_viridis_d(end = 0.8)
  
  # 1. Lambda vs TP
  p1 <- ggplot(data_n, aes(x = lambda_actual, y = TP_mean, color = as.factor(lambda_actual))) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = TP_mean - TP_sd, ymax = TP_mean + TP_sd), width = 0.02) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    labs(
      title = sprintf("True Positives vs Lambda (N = %d)", n_value),
      subtitle = "Average number of correctly identified critical links",
      x = "Lambda Value (log scale)",
      y = "True Positives",
      color = "Lambda"
    ) +
    theme_minimal() +
    color_palette +
    theme(legend.position = "right")
  
  # 2. Lambda vs FP
  p2 <- ggplot(data_n, aes(x = lambda_actual, y = FP_mean, color = as.factor(lambda_actual))) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = FP_mean - FP_sd, ymax = FP_mean + FP_sd), width = 0.02) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    labs(
      title = sprintf("False Positives vs Lambda (N = %d)", n_value),
      subtitle = "Average number of incorrectly identified critical links",
      x = "Lambda Value (log scale)",
      y = "False Positives",
      color = "Lambda"
    ) +
    theme_minimal() +
    color_palette
  
  # 3. Lambda vs FN
  p3 <- ggplot(data_n, aes(x = lambda_actual, y = FN_mean, color = as.factor(lambda_actual))) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = FN_mean - FN_sd, ymax = FN_mean + FN_sd), width = 0.02) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    labs(
      title = sprintf("False Negatives vs Lambda (N = %d)", n_value),
      subtitle = "Average number of missed critical links",
      x = "Lambda Value (log scale)",
      y = "False Negatives",
      color = "Lambda"
    ) +
    theme_minimal() +
    color_palette
  
  # 4. Lambda vs TN
  p4 <- ggplot(data_n, aes(x = lambda_actual, y = TN_mean, color = as.factor(lambda_actual))) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = TN_mean - TN_sd, ymax = TN_mean + TN_sd), width = 0.02) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    labs(
      title = sprintf("True Negatives vs Lambda (N = %d)", n_value),
      subtitle = "Average number of correctly identified non-critical links",
      x = "Lambda Value (log scale)",
      y = "True Negatives",
      color = "Lambda"
    ) +
    theme_minimal() +
    color_palette
  
  # 5. Lambda vs Recall
  p5 <- ggplot(data_n, aes(x = lambda_actual, y = recall_mean, color = as.factor(lambda_actual))) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = recall_mean - recall_sd, ymax = recall_mean + recall_sd), width = 0.02) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(
      title = sprintf("Recall vs Lambda (N = %d)", n_value),
      subtitle = "Proportion of actual critical links correctly identified",
      x = "Lambda Value (log scale)",
      y = "Recall",
      color = "Lambda"
    ) +
    theme_minimal() +
    color_palette
  
  # 6. Lambda vs Precision
  p6 <- ggplot(data_n, aes(x = lambda_actual, y = precision_mean, color = as.factor(lambda_actual))) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = precision_mean - precision_sd, ymax = precision_mean + precision_sd), width = 0.02) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(
      title = sprintf("Precision vs Lambda (N = %d)", n_value),
      subtitle = "Proportion of predicted critical links that are correct",
      x = "Lambda Value (log scale)",
      y = "Precision",
      color = "Lambda"
    ) +
    theme_minimal() +
    color_palette
  
  # 7. Lambda vs F1
  p7 <- ggplot(data_n, aes(x = lambda_actual, y = f1_mean, color = as.factor(lambda_actual))) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = f1_mean - f1_sd, ymax = f1_mean + f1_sd), width = 0.02) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(
      title = sprintf("F1 Score vs Lambda (N = %d)", n_value),
      subtitle = "Harmonic mean of precision and recall",
      x = "Lambda Value (log scale)",
      y = "F1 Score",
      color = "Lambda"
    ) +
    theme_minimal() +
    color_palette
  
  # 8. Lambda vs MCC
  p8 <- ggplot(data_n, aes(x = lambda_actual, y = mcc_mean, color = as.factor(lambda_actual))) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mcc_mean - mcc_sd, ymax = mcc_mean + mcc_sd), width = 0.02) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    scale_y_continuous(limits = c(-1, 1)) +
    labs(
      title = sprintf("Matthews Correlation Coefficient vs Lambda (N = %d)", n_value),
      subtitle = "Balanced measure of binary classification quality",
      x = "Lambda Value (log scale)",
      y = "MCC",
      color = "Lambda"
    ) +
    theme_minimal() +
    color_palette
  
  # 9. Lambda vs Predicted Ratio
  p9 <- ggplot(data_n, aes(x = lambda_actual, y = predicted_ratio, color = as.factor(lambda_actual))) +
    geom_point(size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    labs(
      title = sprintf("Predicted/True Ratio vs Lambda (N = %d)", n_value),
      subtitle = "Ratio of predicted to actual number of critical links",
      x = "Lambda Value (log scale)",
      y = "Predicted/True Ratio",
      color = "Lambda"
    ) +
    theme_minimal() +
    color_palette
  
  # 10. Lambda vs TP-FN Tradeoff
  p10 <- ggplot(data_n, aes(x = lambda_actual, y = tp_fn_tradeoff, color = as.factor(lambda_actual))) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    labs(
      title = sprintf("TP-FN Tradeoff vs Lambda (N = %d)", n_value),
      subtitle = "Difference between true positives and false negatives",
      x = "Lambda Value (log scale)",
      y = "TP - FN",
      color = "Lambda"
    ) +
    theme_minimal() +
    color_palette
  
  # 11. Outlier Analysis
  p11 <- ggplot(data_n, aes(x = as.factor(lambda_actual), y = runtime_mean)) +
    geom_boxplot(aes(fill = as.factor(lambda_actual))) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
    scale_y_log10() +
    labs(
      title = sprintf("Runtime Distribution vs Lambda (N = %d)", n_value),
      subtitle = "Computational time analysis with outlier detection",
      x = "Lambda Value",
      y = "Runtime (seconds, log scale)",
      fill = "Lambda"
    ) +
    theme_minimal() +
    scale_fill_viridis_d(end = 0.8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save individual plots
  plot_dir <- file.path(output_dir, sprintf("N_%05d", n_value))
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  ggsave(file.path(plot_dir, "01_lambda_vs_tp.png"), p1, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "02_lambda_vs_fp.png"), p2, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "03_lambda_vs_fn.png"), p3, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "04_lambda_vs_tn.png"), p4, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "05_lambda_vs_recall.png"), p5, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "06_lambda_vs_precision.png"), p6, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "07_lambda_vs_f1.png"), p7, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "08_lambda_vs_mcc.png"), p8, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "09_lambda_vs_predicted_ratio.png"), p9, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "10_lambda_vs_tp_fn_tradeoff.png"), p10, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "11_outlier_analysis.png"), p11, width = 10, height = 6, dpi = 300)
  
  # Create combined plot
  combined <- (p1 + p2) / (p3 + p4) / (p5 + p6) / (p7 + p8) / (p9 + p10) +
    plot_annotation(
      title = sprintf("Ring Experiment Results Summary (N = %d)", n_value),
      subtitle = sprintf("Based on %d runs per parameter combination", unique(data_n$n_runs)),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  ggsave(
    file.path(plot_dir, "00_combined_summary.png"), 
    combined, 
    width = 20, 
    height = 30, 
    dpi = 300
  )
}

# Function to create combined plots across all N values
create_combined_plots <- function(aggregated_data, output_dir) {
  # Filter for lambda_half and lambda_double only
  aggregated_data <- aggregated_data %>%
    filter(perturbation_type %in% c("lambda_half", "lambda_double"))
  
  combined_dir <- file.path(output_dir, "combined_all_N")
  dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. F1 Score across all N
  p_f1_all <- ggplot(aggregated_data, aes(x = lambda_actual, y = f1_mean, color = as.factor(N))) +
    geom_point(size = 3) +
    geom_line(aes(group = N), alpha = 0.5) +
    geom_errorbar(aes(ymin = f1_mean - f1_sd, ymax = f1_mean + f1_sd), width = 0.02, alpha = 0.5) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(
      title = "F1 Score vs Lambda for Different Network Sizes",
      subtitle = "Performance comparison across network scales",
      x = "Lambda Value (log scale)",
      y = "F1 Score",
      color = "Network Size (N)"
    ) +
    theme_minimal() +
    scale_color_viridis_d(end = 0.8)
  
  # 2. MCC across all N
  p_mcc_all <- ggplot(aggregated_data, aes(x = lambda_actual, y = mcc_mean, color = as.factor(N))) +
    geom_point(size = 3) +
    geom_line(aes(group = N), alpha = 0.5) +
    geom_errorbar(aes(ymin = mcc_mean - mcc_sd, ymax = mcc_mean + mcc_sd), width = 0.02, alpha = 0.5) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    scale_y_continuous(limits = c(-1, 1)) +
    labs(
      title = "Matthews Correlation Coefficient vs Lambda for Different Network Sizes",
      subtitle = "Balanced classification performance across network scales",
      x = "Lambda Value (log scale)",
      y = "MCC",
      color = "Network Size (N)"
    ) +
    theme_minimal() +
    scale_color_viridis_d(end = 0.8)
  
  # 3. Recall across all N
  p_recall_all <- ggplot(aggregated_data, aes(x = lambda_actual, y = recall_mean, color = as.factor(N))) +
    geom_point(size = 3) +
    geom_line(aes(group = N), alpha = 0.5) +
    geom_errorbar(aes(ymin = recall_mean - recall_sd, ymax = recall_mean + recall_sd), width = 0.02, alpha = 0.5) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(
      title = "Recall vs Lambda for Different Network Sizes",
      subtitle = "Sensitivity across network scales",
      x = "Lambda Value (log scale)",
      y = "Recall",
      color = "Network Size (N)"
    ) +
    theme_minimal() +
    scale_color_viridis_d(end = 0.8)
  
  # 4. Precision across all N
  p_precision_all <- ggplot(aggregated_data, aes(x = lambda_actual, y = precision_mean, color = as.factor(N))) +
    geom_point(size = 3) +
    geom_line(aes(group = N), alpha = 0.5) +
    geom_errorbar(aes(ymin = precision_mean - precision_sd, ymax = precision_mean + precision_sd), width = 0.02, alpha = 0.5) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(
      title = "Precision vs Lambda for Different Network Sizes",
      subtitle = "Positive predictive value across network scales",
      x = "Lambda Value (log scale)",
      y = "Precision",
      color = "Network Size (N)"
    ) +
    theme_minimal() +
    scale_color_viridis_d(end = 0.8)
  
  # 5. Runtime scaling
  p_runtime_all <- ggplot(aggregated_data, aes(x = N, y = runtime_mean, color = as.factor(lambda_actual))) +
    geom_point(size = 3) +
    geom_line(aes(group = lambda_actual), alpha = 0.5) +
    geom_errorbar(aes(ymin = runtime_mean - runtime_sd, ymax = runtime_mean + runtime_sd), width = 0.02, alpha = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "Runtime Scaling with Network Size",
      subtitle = "Computational complexity analysis",
      x = "Network Size (N, log scale)",
      y = "Runtime (seconds, log scale)",
      color = "Lambda Value"
    ) +
    theme_minimal() +
    scale_color_viridis_d(end = 0.8)
  
  # 6. Predicted ratio across all N
  p_ratio_all <- ggplot(aggregated_data, aes(x = lambda_actual, y = predicted_ratio, color = as.factor(N))) +
    geom_point(size = 3) +
    geom_line(aes(group = N), alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
    labs(
      title = "Predicted/True Ratio vs Lambda for Different Network Sizes",
      subtitle = "Over/under-prediction across network scales",
      x = "Lambda Value (log scale)",
      y = "Predicted/True Ratio",
      color = "Network Size (N)"
    ) +
    theme_minimal() +
    scale_color_viridis_d(end = 0.8)
  
  # Save plots
  ggsave(file.path(combined_dir, "f1_score_all_n.png"), p_f1_all, width = 12, height = 8, dpi = 300)
  ggsave(file.path(combined_dir, "mcc_all_n.png"), p_mcc_all, width = 12, height = 8, dpi = 300)
  ggsave(file.path(combined_dir, "recall_all_n.png"), p_recall_all, width = 12, height = 8, dpi = 300)
  ggsave(file.path(combined_dir, "precision_all_n.png"), p_precision_all, width = 12, height = 8, dpi = 300)
  ggsave(file.path(combined_dir, "runtime_scaling.png"), p_runtime_all, width = 12, height = 8, dpi = 300)
  ggsave(file.path(combined_dir, "predicted_ratio_all_n.png"), p_ratio_all, width = 12, height = 8, dpi = 300)
  
  # Create summary dashboard
  dashboard <- (p_f1_all + p_mcc_all) / (p_recall_all + p_precision_all) / (p_ratio_all + p_runtime_all) +
    plot_annotation(
      title = "Ring Experiment Results: Cross-Network Comparison",
      subtitle = "Performance metrics across different network sizes and lambda values",
      theme = theme(plot.title = element_text(size = 18, face = "bold"))
    )
  
  ggsave(
    file.path(combined_dir, "dashboard_all_metrics.png"), 
    dashboard, 
    width = 24, 
    height = 18, 
    dpi = 300
  )
}

# Main analysis function
analyze_ring_experiment <- function(experiment_dir) {
  # Read results
  cat("Reading experiment results...\n")
  results <- read_experiment_results(experiment_dir)
  
  # Aggregate results
  cat("Aggregating results across runs...\n")
  aggregated <- aggregate_results(results)
  
  # Create output directory
  output_dir <- file.path(experiment_dir, "analysis_plots")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create plots for each N value
  n_values <- unique(aggregated$N)
  cat(sprintf("Creating plots for %d different network sizes...\n", length(n_values)))
  
  for (n_val in n_values) {
    cat(sprintf("  Processing N = %d...\n", n_val))
    data_n <- aggregated %>% filter(N == n_val)
    create_plots_for_n(data_n, n_val, output_dir)
  }
  
  # Create combined plots
  cat("Creating combined plots across all N values...\n")
  create_combined_plots(aggregated, output_dir)
  
  # Save aggregated data for further analysis
  write_csv(aggregated, file.path(output_dir, "aggregated_results.csv"))
  
  # Create summary report
  create_summary_report(aggregated, output_dir)
  
  cat(sprintf("\nAnalysis complete! Results saved to: %s\n", output_dir))
}

# Function to create a summary report
create_summary_report <- function(aggregated_data, output_dir) {
  report_lines <- c(
    "# Ring Experiment Analysis Summary",
    "",
    sprintf("Generated: %s", Sys.time()),
    "",
    "## Experiment Overview",
    "",
    sprintf("- Network sizes tested: %s", paste(sort(unique(aggregated_data$N)), collapse = ", ")),
    sprintf("- Lambda values tested: %s", paste(sort(unique(aggregated_data$lambda_actual)), collapse = ", ")),
    sprintf("- Total parameter combinations: %d", nrow(aggregated_data)),
    "",
    "## Best Performing Parameters",
    ""
  )
  
  # Find best parameters for each N
  best_by_n <- aggregated_data %>%
    group_by(N) %>%
    slice_max(f1_mean, n = 1) %>%
    select(N, lambda_actual, f1_mean, mcc_mean, recall_mean, precision_mean)
  
  report_lines <- c(report_lines,
    "### By F1 Score:",
    "",
    knitr::kable(best_by_n, format = "markdown", digits = 3)
  )
  
  # Write report
  writeLines(report_lines, file.path(output_dir, "analysis_summary.md"))
}

# Example usage
if (interactive()) {
  # Replace with your actual experiment directory
  experiment_dir <- "results/full_validation/ring_exp_20250627_120647"
  analyze_ring_experiment(experiment_dir)
}
