# plot_experiment_results.R
# Comprehensive analysis and visualization of identify_critical_links experiment results
# Author: [Your Name]
# Date: [Sys.Date()]

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)
library(stringr)
library(reshape2)
library(ggvenn)
library(gridExtra)
library(viridis)

# Set paths - check multiple possible locations
# First, determine the base path
if (file.exists("../results/full_validation/")) {
  results_base <- "../results/full_validation/"
} else if (file.exists("./results/full_validation/")) {
  results_base <- "./results/full_validation/"
} else if (file.exists("~/Desktop/repos/BrainNetTest/results/full_validation/")) {
  results_base <- "~/Desktop/repos/BrainNetTest/results/full_validation/"
} else {
  # Try to find the most recent experiment directory
  possible_paths <- c(
    "../results/",
    "./results/",
    "~/Desktop/repos/BrainNetTest/results/"
  )
  
  for (path in possible_paths) {
    if (dir.exists(path)) {
      # Look for directories starting with "ring_exp_"
      exp_dirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
      ring_dirs <- exp_dirs[grepl("ring_exp_", basename(exp_dirs))]
      
      if (length(ring_dirs) > 0) {
        # Use the most recent one (by name, which includes timestamp)
        results_base <- ring_dirs[order(basename(ring_dirs), decreasing = TRUE)[1]]
        cat(sprintf("Found experiment directory: %s\n", results_base))
        break
      }
    }
  }
  
  if (!exists("results_base")) {
    stop("Could not find results directory. Please specify the correct path.")
  }
}

# Set the results directory
results_dir <- results_base
plots_dir <- './plots/full_validation/'
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Helper: Recursively find all summary CSV files
find_csv_files <- function(path, pattern = '\\.csv$') {
  # First, check if the path exists
  if (!dir.exists(path)) {
    stop(paste("Directory does not exist:", path))
  }
  
  # List all files recursively
  all_files <- list.files(path, pattern = pattern, recursive = TRUE, full.names = TRUE)
  
  # Filter for files in summaries subdirectories
  summary_files <- all_files[grepl('/summaries/', all_files)]
  
  # If no files found in summaries subdirectories, look for any CSV files
  if (length(summary_files) == 0) {
    cat("No CSV files found in summaries subdirectories. Looking for any CSV files...\n")
    summary_files <- all_files
  }
  
  return(summary_files)
}

# 1. Load all experiment summary CSV files
csv_files <- find_csv_files(results_dir)

# Debug: Print found files
if (length(csv_files) > 0) {
  cat(sprintf("Found %d CSV files:\n", length(csv_files)))
  for (f in head(csv_files, 5)) {
    cat(sprintf("  - %s\n", f))
  }
  if (length(csv_files) > 5) {
    cat(sprintf("  ... and %d more files\n", length(csv_files) - 5))
  }
} else {
  # Try to find the correct path by looking for the experiment directory
  exp_dirs <- list.dirs(results_dir, recursive = TRUE)
  summaries_dirs <- exp_dirs[grepl('/summaries$', exp_dirs)]
  
  if (length(summaries_dirs) > 0) {
    cat("Found summaries directories:\n")
    for (d in summaries_dirs) {
      cat(sprintf("  - %s\n", d))
      csv_in_dir <- list.files(d, pattern = '\\.csv$', full.names = TRUE)
      if (length(csv_in_dir) > 0) {
        cat(sprintf("    Contains %d CSV files\n", length(csv_in_dir)))
      }
    }
    stop('CSV files found but not loaded. Please check the directory structure.')
  } else {
    stop('No summary CSV files found in results_dir. Please check the path.')
  }
}

# Load all CSVs into a single data frame
results_df <- bind_rows(lapply(csv_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  # Add source file info
  df$source_file <- basename(f)
  df
}))

# Remove rows with missing N or perturbation_type
results_df <- results_df %>% filter(!is.na(N), !is.na(perturbation_type))

# Create a combined parameter column for lambda-based and const-based perturbations
results_df <- results_df %>%
  mutate(
    param_value = case_when(
      perturbation_type %in% c("lambda_half", "lambda_double") ~ lambda_mult,
      perturbation_type %in% c("const_high", "const_low") ~ p_const,
      TRUE ~ NA_real_
    ),
    param_type = case_when(
      perturbation_type %in% c("lambda_half", "lambda_double") ~ "lambda_multiplier",
      perturbation_type %in% c("const_high", "const_low") ~ "probability",
      TRUE ~ "none"
    )
  )

# Calculate summary statistics by group
summary_stats <- results_df %>%
  group_by(perturbation_type, N, param_value) %>%
  summarise(
    n_runs = n(),
    mean_precision = mean(precision, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE),
    mean_f1 = mean(f1, na.rm = TRUE),
    sd_precision = sd(precision, na.rm = TRUE),
    sd_recall = sd(recall, na.rm = TRUE),
    sd_f1 = sd(f1, na.rm = TRUE),
    mean_runtime = mean(runtime, na.rm = TRUE),
    success_rate = mean(global_significant, na.rm = TRUE),
    .groups = 'drop'
  )

# 2. Performance metrics vs parameter values (for perturbations with parameters)
param_results <- results_df %>%
  filter(!is.na(param_value))

if (nrow(param_results) > 0) {
  g_param_f1 <- ggplot(param_results, aes(x = param_value, y = f1, color = factor(N))) +
    geom_boxplot(aes(group = interaction(param_value, N)), alpha = 0.7) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.3) +
    facet_wrap(~ perturbation_type, scales = "free_x") +
    labs(title = 'F1 Score vs Parameter Value',
         x = 'Parameter Value', y = 'F1 Score', color = 'Network Size (N)') +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(file.path(plots_dir, 'f1_vs_param_value.png'), g_param_f1, width = 10, height = 6)
}

# 3. Precision-Recall curves
pr_summary <- results_df %>%
  group_by(perturbation_type, N) %>%
  arrange(recall) %>%
  mutate(
    recall_bin = cut(recall, breaks = seq(0, 1, 0.1), include.lowest = TRUE),
    precision_bin = cut(precision, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
  )

g_pr <- ggplot(results_df, aes(x = recall, y = precision, color = factor(N))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  facet_wrap(~ perturbation_type) +
  labs(title = 'Precision-Recall Curves', x = 'Recall', y = 'Precision', 
       color = 'Network Size (N)') +
  theme_bw() +
  xlim(0, 1) + ylim(0, 1) +
  geom_abline(intercept = 1, slope = -1, linetype = "dashed", alpha = 0.5)
ggsave(file.path(plots_dir, 'precision_recall_curves.png'), g_pr, width = 10, height = 8)

# 4. F1 Score heatmap (N vs perturbation type)
f1_heatmap_data <- results_df %>%
  group_by(N, perturbation_type) %>%
  summarise(mean_f1 = mean(f1, na.rm = TRUE), .groups = 'drop')

g_heatmap <- ggplot(f1_heatmap_data, aes(x = factor(N), y = perturbation_type, fill = mean_f1)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", mean_f1)), color = "white", size = 4) +
  scale_fill_viridis(name = "Mean F1\nScore", limits = c(0, 1)) +
  labs(title = 'Mean F1 Score by Network Size and Perturbation Type',
       x = 'Network Size (N)', y = 'Perturbation Type') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plots_dir, 'f1_heatmap.png'), g_heatmap, width = 8, height = 6)

# 5. Runtime scaling analysis
g_runtime <- ggplot(results_df, aes(x = N, y = runtime, color = perturbation_type)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = TRUE) +
  scale_y_log10() +
  scale_x_log10() +
  labs(title = 'Runtime Scaling with Network Size',
       x = 'Network Size (N, log scale)', 
       y = 'Runtime (seconds, log scale)',
       color = 'Perturbation Type') +
  theme_bw()
ggsave(file.path(plots_dir, 'runtime_scaling.png'), g_runtime, width = 10, height = 6)

# 6. Confusion matrix components over N
conf_matrix_long <- results_df %>%
  select(N, perturbation_type, TP, FP, FN, TN) %>%
  gather(key = "metric", value = "count", TP, FP, FN, TN) %>%
  group_by(N, perturbation_type, metric) %>%
  summarise(mean_count = mean(count, na.rm = TRUE), .groups = 'drop')

g_conf <- ggplot(conf_matrix_long, aes(x = factor(N), y = mean_count, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ perturbation_type, scales = "free_y") +
  labs(title = 'Confusion Matrix Components by Network Size',
       x = 'Network Size (N)', y = 'Mean Count', fill = 'Metric') +
  theme_bw() +
  scale_fill_manual(values = c("TP" = "darkgreen", "TN" = "lightgreen", 
                              "FP" = "darkred", "FN" = "orange"))
ggsave(file.path(plots_dir, 'confusion_matrix_components.png'), g_conf, width = 12, height = 8)

# 7. Success rate and power analysis
power_stats <- results_df %>%
  group_by(perturbation_type, N) %>%
  summarise(
    success_rate = mean(!is.na(global_significant)),
    power = mean(global_significant, na.rm = TRUE),
    n_experiments = n(),
    .groups = 'drop'
  )

g_power <- ggplot(power_stats, aes(x = factor(N), y = power, fill = perturbation_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.2f", power)), 
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = 'Statistical Power by Network Size and Perturbation Type',
       x = 'Network Size (N)', y = 'Power (Proportion of Significant Tests)',
       fill = 'Perturbation Type') +
  theme_bw() +
  ylim(0, 1.1)
ggsave(file.path(plots_dir, 'power_analysis.png'), g_power, width = 10, height = 6)

# 8. Distribution plots for each perturbation type
for (pert_type in unique(results_df$perturbation_type)) {
  pert_data <- results_df %>% filter(perturbation_type == pert_type)
  
  # Create multi-panel plot
  p1 <- ggplot(pert_data, aes(x = precision, fill = factor(N))) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    labs(title = paste(pert_type, "- Precision Distribution"), x = "Precision", y = "Count") +
    theme_bw() + theme(legend.position = "none")
  
  p2 <- ggplot(pert_data, aes(x = recall, fill = factor(N))) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    labs(title = paste(pert_type, "- Recall Distribution"), x = "Recall", y = "Count") +
    theme_bw() + theme(legend.position = "none")
  
  p3 <- ggplot(pert_data, aes(x = f1, fill = factor(N))) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    labs(title = paste(pert_type, "- F1 Distribution"), x = "F1 Score", y = "Count") +
    theme_bw() + theme(legend.position = "bottom")
  
  p4 <- ggplot(pert_data, aes(x = mcc, fill = factor(N))) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    labs(title = paste(pert_type, "- MCC Distribution"), x = "MCC", y = "Count", 
         fill = "Network Size") +
    theme_bw() + theme(legend.position = "bottom")
  
  combined <- grid.arrange(p1, p2, p3, p4, ncol = 2)
  ggsave(file.path(plots_dir, paste0('distributions_', pert_type, '.png')), 
         combined, width = 12, height = 10)
}

# 9. Edge detection accuracy (true positive rate vs false positive rate)
tpr_fpr_data <- results_df %>%
  mutate(
    TPR = TP / (TP + FN),  # True Positive Rate (Sensitivity/Recall)
    FPR = FP / (FP + TN)   # False Positive Rate
  ) %>%
  filter(!is.na(TPR) & !is.na(FPR))

g_roc <- ggplot(tpr_fpr_data, aes(x = FPR, y = TPR, color = factor(N))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ perturbation_type) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  labs(title = 'ROC-like Curves (TPR vs FPR)',
       x = 'False Positive Rate', y = 'True Positive Rate',
       color = 'Network Size (N)') +
  theme_bw() +
  xlim(0, 1) + ylim(0, 1)
ggsave(file.path(plots_dir, 'roc_curves.png'), g_roc, width = 10, height = 8)

# 10. Parameter sensitivity analysis (if multiple parameter values tested)
param_sensitivity <- results_df %>%
  filter(!is.na(param_value)) %>%
  group_by(perturbation_type, N, param_value) %>%
  summarise(
    mean_f1 = mean(f1, na.rm = TRUE),
    se_f1 = sd(f1, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

if (nrow(param_sensitivity) > 0) {
  g_sensitivity <- ggplot(param_sensitivity, 
                         aes(x = param_value, y = mean_f1, color = factor(N))) +
    geom_point(size = 3) +
    geom_line() +
    geom_errorbar(aes(ymin = mean_f1 - se_f1, ymax = mean_f1 + se_f1), width = 0.05) +
    facet_wrap(~ perturbation_type, scales = "free_x") +
    labs(title = 'Parameter Sensitivity Analysis',
         x = 'Parameter Value', y = 'Mean F1 Score (± SE)',
         color = 'Network Size (N)') +
    theme_bw()
  ggsave(file.path(plots_dir, 'parameter_sensitivity.png'), g_sensitivity, width = 10, height = 6)
}

# 11. Lambda value effects on performance metrics
# Calculate actual lambda values for lambda-based perturbations
lambda_results <- results_df %>%
  filter(perturbation_type %in% c("lambda_half", "lambda_double")) %>%
  mutate(
    actual_lambda = lambda_base * lambda_mult,
    # Create a label for lambda multiplier
    lambda_label = case_when(
      lambda_mult == 0.5 ~ "λ × 0.5",
      lambda_mult == 2 ~ "λ × 2",
      lambda_mult == 3 ~ "λ × 3",
      lambda_mult == 4 ~ "λ × 4",
      lambda_mult == 0.33 ~ "λ × 0.33",
      lambda_mult == 0.25 ~ "λ × 0.25",
      TRUE ~ paste0("λ × ", lambda_mult)
    )
  )

# Calculate mean performance metrics for each lambda value and N
lambda_summary <- lambda_results %>%
  group_by(N, actual_lambda, lambda_mult, lambda_label) %>%
  summarise(
    n_runs = n(),
    mean_precision = mean(precision, na.rm = TRUE),
    se_precision = sd(precision, na.rm = TRUE) / sqrt(n()),
    mean_recall = mean(recall, na.rm = TRUE),
    se_recall = sd(recall, na.rm = TRUE) / sqrt(n()),
    mean_f1 = mean(f1, na.rm = TRUE),
    se_f1 = sd(f1, na.rm = TRUE) / sqrt(n()),
    mean_mcc = mean(mcc, na.rm = TRUE),
    se_mcc = sd(mcc, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Convert to long format for faceting
lambda_long <- lambda_summary %>%
  select(N, actual_lambda, lambda_mult, lambda_label, starts_with("mean_")) %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "metric",
    values_to = "value",
    names_prefix = "mean_"
  ) %>%
  mutate(metric = factor(metric, 
                        levels = c("precision", "recall", "f1", "mcc"),
                        labels = c("Precision", "Recall", "F1 Score", "MCC")))

# Add standard errors
lambda_se_long <- lambda_summary %>%
  select(N, actual_lambda, lambda_mult, lambda_label, starts_with("se_")) %>%
  pivot_longer(
    cols = starts_with("se_"),
    names_to = "metric",
    values_to = "se",
    names_prefix = "se_"
  ) %>%
  mutate(metric = factor(metric, 
                        levels = c("precision", "recall", "f1", "mcc"),
                        labels = c("Precision", "Recall", "F1 Score", "MCC")))

# Combine mean values and standard errors
lambda_plot_data <- lambda_long %>%
  left_join(lambda_se_long, by = c("N", "actual_lambda", "lambda_mult", "lambda_label", "metric"))

# Create the plot with actual lambda values on x-axis
g_lambda_performance <- ggplot(lambda_plot_data, 
                              aes(x = actual_lambda, y = value, color = factor(N))) +
  geom_point(size = 3) +
  geom_line(alpha = 0.7) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.02, alpha = 0.5) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2) +
  scale_x_continuous(trans = "log10") +
  labs(title = 'Performance Metrics vs Lambda Values',
       x = 'Lambda Value (log scale)', 
       y = 'Score',
       color = 'Network Size (N)') +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plots_dir, 'lambda_performance_metrics.png'), 
       g_lambda_performance, width = 12, height = 10)

# Alternative plot with lambda multipliers on x-axis (categorical)
g_lambda_mult_performance <- ggplot(lambda_plot_data, 
                                   aes(x = lambda_label, y = value, fill = factor(N))) +
  geom_boxplot(data = lambda_results %>%
                 mutate(lambda_label = case_when(
                   lambda_mult == 0.5 ~ "λ × 0.5",
                   lambda_mult == 2 ~ "λ × 2",
                   lambda_mult == 3 ~ "λ × 3",
                   lambda_mult == 4 ~ "λ × 4",
                   lambda_mult == 0.33 ~ "λ × 0.33",
                   lambda_mult == 0.25 ~ "λ × 0.25",
                   TRUE ~ paste0("λ × ", lambda_mult)
                 )) %>%
                 select(N, lambda_label, precision, recall, f1, mcc) %>%
                 pivot_longer(cols = c(precision, recall, f1, mcc),
                            names_to = "metric",
                            values_to = "value") %>%
                 mutate(metric = factor(metric, 
                                      levels = c("precision", "recall", "f1", "mcc"),
                                      labels = c("Precision", "Recall", "F1 Score", "MCC"))),
               aes(x = lambda_label, y = value, fill = factor(N)),
               alpha = 0.6, outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.75), size = 2, shape = 21, color = "black") +
  facet_wrap(~ metric, scales = "free_y", nrow = 2) +
  labs(title = 'Performance Metrics vs Lambda Multipliers',
       x = 'Lambda Multiplier', 
       y = 'Score',
       fill = 'Network Size (N)') +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plots_dir, 'lambda_multiplier_performance.png'), 
       g_lambda_mult_performance, width = 12, height = 10)

# Create a detailed table of lambda effects
lambda_table <- lambda_summary %>%
  arrange(N, lambda_mult) %>%
  mutate(
    precision = sprintf("%.3f ± %.3f", mean_precision, se_precision),
    recall = sprintf("%.3f ± %.3f", mean_recall, se_recall),
    f1 = sprintf("%.3f ± %.3f", mean_f1, se_f1),
    mcc = sprintf("%.3f ± %.3f", mean_mcc, se_mcc)
  ) %>%
  select(N, lambda_mult, actual_lambda, n_runs, precision, recall, f1, mcc)

write.csv(lambda_table, file.path(plots_dir, 'lambda_effects_table.csv'), row.names = FALSE)

# Print summary to console
cat("\nLambda Effects Summary:\n")
cat("======================\n")
print(lambda_table)

# 8. Save summary statistics with additional metrics
overall_stats <- results_df %>%
  group_by(perturbation_type, N, param_value) %>%
  summarise(
    n_runs = n(),
    mean_precision = mean(precision, na.rm=TRUE),
    mean_recall = mean(recall, na.rm=TRUE),
    mean_f1 = mean(f1, na.rm=TRUE),
    mean_mcc = mean(mcc, na.rm=TRUE),
    mean_runtime = mean(runtime, na.rm=TRUE),
    power = mean(global_significant, na.rm=TRUE),
    mean_tp = mean(TP, na.rm=TRUE),
    mean_fp = mean(FP, na.rm=TRUE),
    mean_fn = mean(FN, na.rm=TRUE),
    mean_tn = mean(TN, na.rm=TRUE),
    .groups = 'drop')
write.csv(overall_stats, file.path(plots_dir, 'summary_stats.csv'), row.names=FALSE)

# Create a summary report
report_file <- file.path(plots_dir, 'plot_summary_report.txt')
sink(report_file)
cat("RING EXPERIMENT VISUALIZATION SUMMARY\n")
cat("=====================================\n\n")
cat(sprintf("Total experiments analyzed: %d\n", nrow(results_df)))
cat(sprintf("Network sizes tested: %s\n", paste(sort(unique(results_df$N)), collapse=", ")))
cat(sprintf("Perturbation types: %s\n", paste(unique(results_df$perturbation_type), collapse=", ")))
cat("\nBest performing conditions (by mean F1 score):\n")
cat("---------------------------------------------\n")
top_conditions <- overall_stats %>%
  arrange(desc(mean_f1)) %>%
  head(10)
print(top_conditions)
sink()

cat('All plots and summary statistics saved to', plots_dir, '\n')
cat('Summary report saved to', report_file, '\n')

# ---
# EXPLANATION OF PLOTS (for user reference):
#
# 1. Precision/Recall/F1 vs Lambda: Shows how detection performance varies with network sparsity/density.
# 2. Precision/Recall/F1 vs Network Size: Shows scalability and robustness to network size.
# 3. Confusion Matrix Heatmap: Visualizes the types of errors made by the algorithm.
# 4. F1 by Perturbation Type: Compares which perturbations are easier/harder to detect.
# 5. Runtime vs N: Shows computational cost scaling.
# 6. p-value Distribution: Checks calibration of edge tests (if available).
# 7. Summary stats: For further custom analysis.
#
# You can extend this script to add more plots (e.g., Venn diagrams for edge overlap) as needed.
