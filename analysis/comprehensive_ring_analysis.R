# comprehensive_ring_analysis.R
# Comprehensive multi-perspective analysis of ring experiment results
# Generates extensive visualizations to understand algorithm performance

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)
library(stringr)
library(reshape2)
library(gridExtra)
library(viridis)
library(ggridges)
library(corrplot)
library(pheatmap)
library(scales)
library(ggforce)
library(ggrepel)
library(tibble)  # Add this for column_to_rownames

# Set up paths and directories
source("analysis/plot_experiment_results.R")  # Reuse path finding logic

# Create subdirectory for comprehensive plots
comprehensive_plots_dir <- file.path(plots_dir, 'comprehensive_analysis')
if (!dir.exists(comprehensive_plots_dir)) dir.create(comprehensive_plots_dir, recursive = TRUE)

# Helper function to save plots with consistent formatting
save_plot <- function(plot, filename, width = 10, height = 8, ...) {
  ggsave(file.path(comprehensive_plots_dir, filename), plot, width = width, height = height, ...)
}

cat("Starting comprehensive ring experiment analysis...\n")
cat(sprintf("Total experiments: %d\n", nrow(results_df)))
cat(sprintf("Output directory: %s\n", comprehensive_plots_dir))

# ============================================================================
# 1. PERFORMANCE LANDSCAPE ANALYSIS
# ============================================================================
cat("\n1. Creating performance landscape visualizations...\n")

# 1.1 3D-like performance surface plot (F1 score as function of N and parameter)
performance_surface <- results_df %>%
  filter(!is.na(param_value)) %>%
  group_by(N, perturbation_type, param_value) %>%
  summarise(
    mean_f1 = mean(f1, na.rm = TRUE),
    .groups = 'drop'
  )

g_surface <- ggplot(performance_surface, aes(x = log10(N), y = param_value)) +
  geom_tile(aes(fill = mean_f1)) +
  geom_text(aes(label = sprintf("%.2f", mean_f1)), size = 3) +
  facet_wrap(~ perturbation_type, scales = "free_y") +
  scale_fill_viridis(name = "Mean F1", limits = c(0, 1)) +
  labs(
    title = "Performance Landscape: F1 Score across Network Size and Parameters",
    subtitle = "Shows how algorithm performance varies with both network size and perturbation strength",
    x = "Network Size (log10 scale)",
    y = "Parameter Value"
  ) +
  theme_minimal()
save_plot(g_surface, "01_performance_landscape.png", width = 12, height = 8)

# 1.2 Ridge plot showing F1 distribution evolution
g_ridge <- ggplot(results_df, aes(x = f1, y = factor(N), fill = factor(N))) +
  geom_density_ridges(alpha = 0.7, scale = 2) +
  facet_wrap(~ perturbation_type) +
  scale_fill_viridis_d() +
  labs(
    title = "F1 Score Distribution Evolution across Network Sizes",
    subtitle = "Ridge plots show how performance distributions change with network size",
    x = "F1 Score",
    y = "Network Size",
    fill = "Network Size"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
save_plot(g_ridge, "02_f1_distribution_evolution.png", width = 12, height = 10)

# ============================================================================
# 2. ERROR ANALYSIS AND PATTERNS
# ============================================================================
cat("\n2. Analyzing error patterns...\n")

# 2.1 Error composition analysis
error_composition <- results_df %>%
  mutate(
    total_errors = FP + FN,
    error_rate = total_errors / (TP + TN + FP + FN),
    fp_ratio = FP / (FP + FN + 0.001),  # Avoid division by zero
    fn_ratio = FN / (FP + FN + 0.001)
  ) %>%
  filter(!is.na(error_rate))

g_error_comp <- ggplot(error_composition, aes(x = factor(N), y = error_rate)) +
  geom_violin(aes(fill = perturbation_type), alpha = 0.6) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~ perturbation_type) +
  labs(
    title = "Error Rate Distribution by Network Size",
    subtitle = "Shows total error rate (FP + FN) as proportion of all predictions",
    x = "Network Size",
    y = "Error Rate",
    fill = "Perturbation Type"
  ) +
  theme_bw()
save_plot(g_error_comp, "03_error_rate_distribution.png", width = 12, height = 8)

# 2.2 False positive vs false negative trade-off
g_fp_fn_tradeoff <- ggplot(error_composition, aes(x = FP, y = FN, color = factor(N))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ perturbation_type, scales = "free") +
  labs(
    title = "False Positive vs False Negative Trade-off",
    subtitle = "Examines the relationship between different types of errors",
    x = "False Positives",
    y = "False Negatives",
    color = "Network Size"
  ) +
  theme_bw()
save_plot(g_fp_fn_tradeoff, "04_fp_fn_tradeoff.png", width = 12, height = 8)

# ============================================================================
# 3. PARAMETER OPTIMIZATION ANALYSIS
# ============================================================================
cat("\n3. Analyzing parameter optimization...\n")

# 3.1 Optimal parameter identification
param_optimization <- results_df %>%
  filter(!is.na(param_value)) %>%
  group_by(N, perturbation_type, param_value) %>%
  summarise(
    mean_f1 = mean(f1, na.rm = TRUE),
    mean_mcc = mean(mcc, na.rm = TRUE),
    mean_precision = mean(precision, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE),
    n_runs = n(),
    .groups = 'drop'
  ) %>%
  group_by(N, perturbation_type) %>%
  mutate(
    is_optimal_f1 = mean_f1 == max(mean_f1),
    is_optimal_mcc = mean_mcc == max(mean_mcc, na.rm = TRUE)
  )

g_param_opt <- ggplot(param_optimization, aes(x = param_value, y = mean_f1)) +
  geom_line(aes(color = factor(N)), size = 1.2) +
  geom_point(aes(color = factor(N), size = is_optimal_f1)) +
  scale_size_manual(values = c(3, 6), guide = "none") +
  facet_wrap(~ perturbation_type, scales = "free_x") +
  labs(
    title = "Parameter Optimization Curves",
    subtitle = "Larger points indicate optimal parameter values for each network size",
    x = "Parameter Value",
    y = "Mean F1 Score",
    color = "Network Size"
  ) +
  theme_bw()
save_plot(g_param_opt, "05_parameter_optimization.png", width = 12, height = 8)

# 3.2 Parameter sensitivity heatmap
param_sensitivity_matrix <- param_optimization %>%
  select(N, perturbation_type, param_value, mean_f1) %>%
  unite("condition", perturbation_type, param_value) %>%
  spread(condition, mean_f1) %>%
  tibble::column_to_rownames("N")  # Use explicit namespace

# Create heatmap
pheatmap(
  param_sensitivity_matrix,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  main = "Parameter Sensitivity Heatmap\nF1 scores across all parameter combinations",
  filename = file.path(comprehensive_plots_dir, "06_parameter_sensitivity_heatmap.png"),
  width = 14,
  height = 8,
  cellwidth = 30,
  cellheight = 30
)

# ============================================================================
# 4. STATISTICAL POWER AND SIGNIFICANCE
# ============================================================================
cat("\n4. Analyzing statistical power...\n")

# 4.1 Power analysis with confidence intervals
power_analysis <- results_df %>%
  group_by(N, perturbation_type, param_value) %>%
  summarise(
    n_total = n(),
    n_significant = sum(global_significant, na.rm = TRUE),
    power = n_significant / n_total,
    power_se = sqrt(power * (1 - power) / n_total),
    power_ci_lower = power - 1.96 * power_se,
    power_ci_upper = power + 1.96 * power_se,
    .groups = 'drop'
  )

g_power_ci <- ggplot(power_analysis, aes(x = factor(N), y = power, color = perturbation_type)) +
  geom_point(position = position_dodge(0.3), size = 3) +
  geom_errorbar(
    aes(ymin = power_ci_lower, ymax = power_ci_upper),
    position = position_dodge(0.3),
    width = 0.2
  ) +
  geom_hline(yintercept = 0.8, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Statistical Power Analysis with 95% Confidence Intervals",
    subtitle = "Dashed line indicates 80% power threshold",
    x = "Network Size",
    y = "Statistical Power",
    color = "Perturbation Type"
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1))
save_plot(g_power_ci, "07_power_analysis_ci.png", width = 12, height = 8)

# 4.2 Sample size requirements analysis
# Calculate minimum N needed for 80% power
min_n_for_power <- power_analysis %>%
  filter(power >= 0.8) %>%
  group_by(perturbation_type) %>%
  summarise(
    min_N = min(N),
    .groups = 'drop'
  )

g_sample_size <- ggplot(power_analysis, aes(x = N, y = power, color = perturbation_type)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.8, linetype = "dashed", alpha = 0.5) +
  geom_vline(
    data = min_n_for_power,
    aes(xintercept = min_N, color = perturbation_type),
    linetype = "dotted",
    size = 1
  ) +
  scale_x_log10() +
  labs(
    title = "Sample Size Requirements for Adequate Power",
    subtitle = "Vertical lines show minimum N for 80% power",
    x = "Network Size (log scale)",
    y = "Statistical Power",
    color = "Perturbation Type"
  ) +
  theme_bw()
save_plot(g_sample_size, "08_sample_size_requirements.png", width = 10, height = 6)

# ============================================================================
# 5. RUNTIME AND COMPUTATIONAL EFFICIENCY
# ============================================================================
cat("\n5. Analyzing computational efficiency...\n")

# 5.1 Runtime scaling analysis with fitted curves
runtime_models <- results_df %>%
  group_by(perturbation_type) %>%
  do(
    linear = lm(log(runtime) ~ log(N), data = .),
    quadratic = lm(log(runtime) ~ log(N) + I(log(N)^2), data = .)
  )

# Create predictions for smooth curves
N_seq <- seq(min(results_df$N), max(results_df$N), length.out = 100)
runtime_predictions <- runtime_models %>%
  rowwise() %>%
  do(
    data.frame(
      N = N_seq,
      perturbation_type = .$perturbation_type,
      linear_pred = exp(predict(.$linear, newdata = data.frame(N = N_seq))),
      quadratic_pred = exp(predict(.$quadratic, newdata = data.frame(N = N_seq)))
    )
  )

g_runtime_scaling <- ggplot(results_df, aes(x = N, y = runtime)) +
  geom_point(aes(color = perturbation_type), alpha = 0.3) +
  geom_line(
    data = runtime_predictions,
    aes(x = N, y = quadratic_pred, color = perturbation_type),
    size = 1.2
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Runtime Scaling Analysis",
    subtitle = "Fitted curves show O(N²) scaling behavior",
    x = "Network Size (log scale)",
    y = "Runtime in seconds (log scale)",
    color = "Perturbation Type"
  ) +
  theme_bw()
save_plot(g_runtime_scaling, "09_runtime_scaling_fitted.png", width = 10, height = 8)

# 5.2 Efficiency metric: F1 score per second
efficiency_analysis <- results_df %>%
  mutate(
    efficiency = f1 / runtime,
    edges_per_second = (N * (N - 1) / 2) / runtime
  ) %>%
  filter(!is.na(efficiency), efficiency < Inf)

g_efficiency <- ggplot(efficiency_analysis, aes(x = factor(N), y = efficiency)) +
  geom_violin(aes(fill = perturbation_type), alpha = 0.6) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~ perturbation_type) +
  labs(
    title = "Algorithm Efficiency: F1 Score per Second",
    subtitle = "Higher values indicate better performance-to-runtime ratio",
    x = "Network Size",
    y = "F1 Score / Runtime (sec⁻¹)",
    fill = "Perturbation Type"
  ) +
  theme_bw()
save_plot(g_efficiency, "10_algorithm_efficiency.png", width = 12, height = 8)

# ============================================================================
# 6. ROBUSTNESS AND STABILITY ANALYSIS
# ============================================================================
cat("\n6. Analyzing robustness and stability...\n")

# 6.1 Coefficient of variation for performance metrics
stability_analysis <- results_df %>%
  group_by(N, perturbation_type) %>%
  summarise(
    cv_precision = sd(precision, na.rm = TRUE) / mean(precision, na.rm = TRUE),
    cv_recall = sd(recall, na.rm = TRUE) / mean(recall, na.rm = TRUE),
    cv_f1 = sd(f1, na.rm = TRUE) / mean(f1, na.rm = TRUE),
    cv_runtime = sd(runtime, na.rm = TRUE) / mean(runtime, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  gather(metric, cv, -N, -perturbation_type) %>%
  mutate(
    metric = str_replace(metric, "cv_", ""),
    metric = factor(metric, levels = c("precision", "recall", "f1", "runtime"))
  )

g_stability <- ggplot(stability_analysis, aes(x = factor(N), y = cv, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ perturbation_type) +
  labs(
    title = "Algorithm Stability: Coefficient of Variation",
    subtitle = "Lower values indicate more consistent performance",
    x = "Network Size",
    y = "Coefficient of Variation",
    fill = "Metric"
  ) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2")
save_plot(g_stability, "11_algorithm_stability.png", width = 12, height = 8)

# 6.2 Outlier analysis
outlier_analysis <- results_df %>%
  group_by(N, perturbation_type) %>%
  mutate(
    f1_zscore = abs((f1 - mean(f1, na.rm = TRUE)) / sd(f1, na.rm = TRUE)),
    is_outlier = f1_zscore > 2
  ) %>%
  summarise(
    outlier_rate = mean(is_outlier, na.rm = TRUE),
    n_outliers = sum(is_outlier, na.rm = TRUE),
    n_total = n(),
    .groups = 'drop'
  )

g_outliers <- ggplot(outlier_analysis, aes(x = factor(N), y = outlier_rate, fill = perturbation_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    aes(label = sprintf("%d/%d", n_outliers, n_total)),
    position = position_dodge(0.9),
    vjust = -0.5,
    size = 3
  ) +
  labs(
    title = "Outlier Analysis: F1 Score Outliers (|z-score| > 2)",
    subtitle = "Shows proportion of runs with unusual performance",
    x = "Network Size",
    y = "Outlier Rate",
    fill = "Perturbation Type"
  ) +
  theme_bw() +
  scale_y_continuous(labels = percent_format())
save_plot(g_outliers, "12_outlier_analysis.png", width = 10, height = 6)

# ============================================================================
# 7. COMPARATIVE ANALYSIS
# ============================================================================
cat("\n7. Creating comparative visualizations...\n")

# 7.1 Head-to-head comparison matrix
comparison_matrix <- results_df %>%
  group_by(N, perturbation_type) %>%
  summarise(mean_f1 = mean(f1, na.rm = TRUE), .groups = 'drop') %>%
  spread(perturbation_type, mean_f1) %>%
  tibble::column_to_rownames("N")  # Use explicit namespace

# Calculate pairwise differences
n_pert <- ncol(comparison_matrix)
diff_matrix <- matrix(0, n_pert, n_pert)
rownames(diff_matrix) <- colnames(comparison_matrix)
colnames(diff_matrix) <- colnames(comparison_matrix)

for (i in 1:n_pert) {
  for (j in 1:n_pert) {
    diff_matrix[i, j] <- mean(comparison_matrix[, i] - comparison_matrix[, j], na.rm = TRUE)
  }
}

# Check for any NA/NaN/Inf values and replace with 0
if (any(!is.finite(diff_matrix))) {
  warning("Found non-finite values in difference matrix, replacing with 0")
  diff_matrix[!is.finite(diff_matrix)] <- 0
}

# Plot difference matrix
png(file.path(comprehensive_plots_dir, "13_perturbation_comparison_matrix.png"), 
    width = 800, height = 800)

# Use try-catch in case corrplot still has issues
tryCatch({
  corrplot(
    diff_matrix,
    method = "color",
    type = "full",
    order = "original",  # Changed from "hclust" to avoid clustering issues
    title = "Mean F1 Score Differences Between Perturbation Types",
    mar = c(0, 0, 2, 0),
    addCoef.col = "black",
    tl.col = "black",
    tl.srt = 45
  )
}, error = function(e) {
  # Fallback: create a simple heatmap if corrplot fails
  cat("corrplot failed, using simple heatmap instead\n")
  # Reset the plot
  dev.off()
  png(file.path(comprehensive_plots_dir, "13_perturbation_comparison_matrix.png"), 
      width = 800, height = 800)
  
  # Create a simple heatmap
  image(1:n_pert, 1:n_pert, diff_matrix, 
        xlab = "", ylab = "", axes = FALSE,
        main = "Mean F1 Score Differences Between Perturbation Types")
  axis(1, at = 1:n_pert, labels = colnames(diff_matrix), las = 2)
  axis(2, at = 1:n_pert, labels = rownames(diff_matrix), las = 2)
  
  # Add values
  for (i in 1:n_pert) {
    for (j in 1:n_pert) {
      text(i, j, sprintf("%.3f", diff_matrix[i, j]))
    }
  }
})
dev.off()

# 7.2 Radar chart for multi-metric comparison
radar_data <- results_df %>%
  group_by(perturbation_type) %>%
  summarise(
    Precision = mean(precision, na.rm = TRUE),
    Recall = mean(recall, na.rm = TRUE),
    `F1 Score` = mean(f1, na.rm = TRUE),
    MCC = mean(mcc, na.rm = TRUE),
    Power = mean(global_significant, na.rm = TRUE),
    `Speed (1/runtime)` = 1 / mean(runtime, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate_if(is.numeric, ~ (. - min(.)) / (max(.) - min(.)))  # Normalize to 0-1

# Convert to long format for ggplot
radar_long <- radar_data %>%
  gather(metric, value, -perturbation_type)

# Create a polar bar chart as an alternative to radar chart
g_radar <- ggplot(radar_long, aes(x = metric, y = value, fill = perturbation_type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  coord_polar(theta = "x") +
  labs(
    title = "Multi-Metric Performance Comparison",
    subtitle = "Polar chart shows normalized performance across multiple dimensions",
    fill = "Perturbation Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80")
  ) +
  ylim(0, 1.1)  # Add some space at the top

save_plot(g_radar, "14_radar_comparison.png", width = 10, height = 10)

# Alternative: Create a regular bar chart if polar doesn't work well
g_metrics_comparison <- ggplot(radar_long, aes(x = metric, y = value, fill = perturbation_type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  labs(
    title = "Multi-Metric Performance Comparison",
    subtitle = "Normalized performance scores across multiple dimensions",
    x = "Metric",
    y = "Normalized Score (0-1)",
    fill = "Perturbation Type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  ylim(0, 1.1)

save_plot(g_metrics_comparison, "14b_metrics_comparison_bars.png", width = 10, height = 6)

# ============================================================================
# 8. EDGE DETECTION PATTERNS
# ============================================================================
cat("\n8. Analyzing edge detection patterns...\n")

# 8.1 True positive rate by edge position
edge_pattern_analysis <- results_df %>%
  mutate(
    tpr = TP / (TP + FN),
    edge_density = (TP + FP) / (N * (N - 1) / 2),
    detection_ratio = n_predicted / n_true
  ) %>%
  filter(!is.na(tpr))

g_detection_ratio <- ggplot(edge_pattern_analysis, aes(x = detection_ratio, fill = perturbation_type)) +
  geom_histogram(alpha = 0.7, bins = 30) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  facet_grid(perturbation_type ~ N, scales = "free_y") +
  labs(
    title = "Edge Detection Ratio Distribution",
    subtitle = "Ratio of predicted to true critical edges (red line = perfect detection)",
    x = "Detection Ratio (Predicted / True)",
    y = "Count",
    fill = "Perturbation Type"
  ) +
  theme_bw() +
  theme(legend.position = "none")
save_plot(g_detection_ratio, "15_edge_detection_ratio.png", width = 14, height = 10)

# 8.2 Precision-Recall trade-off by parameter value
pr_by_param <- results_df %>%
  filter(!is.na(param_value)) %>%
  group_by(perturbation_type, param_value) %>%
  summarise(
    mean_precision = mean(precision, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE),
    se_precision = sd(precision, na.rm = TRUE) / sqrt(n()),
    se_recall = sd(recall, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

g_pr_param <- ggplot(pr_by_param, aes(x = mean_recall, y = mean_precision)) +
  geom_point(aes(color = param_value), size = 4) +
  geom_errorbar(
    aes(ymin = mean_precision - se_precision, ymax = mean_precision + se_precision),
    width = 0.01, alpha = 0.5
  ) +
  geom_errorbarh(
    aes(xmin = mean_recall - se_recall, xmax = mean_recall + se_recall),
    height = 0.01, alpha = 0.5
  ) +
  geom_text_repel(aes(label = round(param_value, 2)), size = 3) +
  facet_wrap(~ perturbation_type) +
  scale_color_viridis() +
  labs(
    title = "Precision-Recall Trade-off by Parameter Value",
    subtitle = "Points labeled with parameter values, error bars show ±1 SE",
    x = "Mean Recall",
    y = "Mean Precision",
    color = "Parameter\nValue"
  ) +
  theme_bw() +
  coord_fixed() +
  xlim(0, 1) + ylim(0, 1)
save_plot(g_pr_param, "16_pr_tradeoff_by_parameter.png", width = 12, height = 8)

# ============================================================================
# 9. SUMMARY DASHBOARD
# ============================================================================
cat("\n9. Creating summary dashboard...\n")

# Create a comprehensive summary dashboard
# Select key metrics for each perturbation type and network size
dashboard_data <- results_df %>%
  group_by(perturbation_type, N) %>%
  summarise(
    `Runs` = n(),
    `Power` = mean(global_significant, na.rm = TRUE),
    `F1` = mean(f1, na.rm = TRUE),
    `MCC` = mean(mcc, na.rm = TRUE),
    `Runtime` = mean(runtime, na.rm = TRUE),
    .groups = 'drop'
  )

# Create individual plots for dashboard
p1 <- ggplot(dashboard_data, aes(x = factor(N), y = Power, fill = perturbation_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Statistical Power", x = "N", y = "Power") +
  theme_minimal() + theme(legend.position = "none")

p2 <- ggplot(dashboard_data, aes(x = factor(N), y = F1, fill = perturbation_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mean F1 Score", x = "N", y = "F1") +
  theme_minimal() + theme(legend.position = "none")

p3 <- ggplot(dashboard_data, aes(x = factor(N), y = MCC, fill = perturbation_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mean MCC", x = "N", y = "MCC") +
  theme_minimal() + theme(legend.position = "bottom")

p4 <- ggplot(dashboard_data, aes(x = factor(N), y = Runtime, fill = perturbation_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_log10() +
  labs(title = "Mean Runtime (log scale)", x = "N", y = "Runtime (s)") +
  theme_minimal() + theme(legend.position = "bottom")

dashboard <- grid.arrange(p1, p2, p3, p4, ncol = 2,
                         top = "Ring Experiment Summary Dashboard")
save_plot(dashboard, "17_summary_dashboard.png", width = 12, height = 10)

# ============================================================================
# 10. DETAILED SUMMARY REPORT
# ============================================================================
cat("\n10. Generating detailed summary report...\n")

report_file <- file.path(comprehensive_plots_dir, "comprehensive_analysis_report.txt")
sink(report_file)

cat("COMPREHENSIVE RING EXPERIMENT ANALYSIS REPORT\n")
cat("============================================\n\n")
cat(sprintf("Analysis Date: %s\n", Sys.Date()))
cat(sprintf("Total Experiments Analyzed: %d\n", nrow(results_df)))
cat(sprintf("Network Sizes: %s\n", paste(sort(unique(results_df$N)), collapse = ", ")))
cat(sprintf("Perturbation Types: %s\n\n", paste(unique(results_df$perturbation_type), collapse = ", ")))

cat("KEY FINDINGS:\n")
cat("-------------\n\n")

# Best performing conditions
cat("1. BEST PERFORMING CONDITIONS (by F1 score):\n")
best_conditions <- results_df %>%
  group_by(perturbation_type, N, param_value) %>%
  summarise(
    mean_f1 = mean(f1, na.rm = TRUE),
    mean_precision = mean(precision, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE),
    n_runs = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_f1)) %>%
  head(10)
print(best_conditions)

cat("\n2. OPTIMAL PARAMETERS BY NETWORK SIZE:\n")
optimal_params <- param_optimization %>%
  filter(is_optimal_f1) %>%
  select(N, perturbation_type, param_value, mean_f1) %>%
  arrange(N, perturbation_type)
print(optimal_params)

cat("\n3. COMPUTATIONAL EFFICIENCY SUMMARY:\n")
efficiency_summary <- results_df %>%
  group_by(N) %>%
  summarise(
    mean_runtime = mean(runtime, na.rm = TRUE),
    edges_per_second = mean((N * (N - 1) / 2) / runtime, na.rm = TRUE),
    .groups = 'drop'
  )
print(efficiency_summary)

cat("\n4. STATISTICAL POWER BY CONDITION:\n")
power_summary <- power_analysis %>%
  select(N, perturbation_type, power, power_ci_lower, power_ci_upper) %>%
  arrange(N, perturbation_type)
print(power_summary)

sink()

# Create plot descriptions file
descriptions_file <- file.path(comprehensive_plots_dir, "plot_descriptions.txt")
sink(descriptions_file)

cat("COMPREHENSIVE ANALYSIS PLOT DESCRIPTIONS\n")
cat("=======================================\n\n")

cat("1. performance_landscape.png\n")
cat("   - Shows F1 score as a heatmap across network size and parameter values\n")
cat("   - Helps identify optimal parameter regions for different network sizes\n\n")

cat("2. f1_distribution_evolution.png\n")
cat("   - Ridge plots showing how F1 score distributions change with network size\n")
cat("   - Reveals whether performance becomes more or less variable as N increases\n\n")

cat("3. error_rate_distribution.png\n")
cat("   - Violin plots of total error rates by condition\n")
cat("   - Shows which conditions produce more reliable results\n\n")

cat("4. fp_fn_tradeoff.png\n")
cat("   - Scatter plots examining the relationship between false positives and false negatives\n")
cat("   - Reveals whether the algorithm tends to over- or under-detect edges\n\n")

cat("5. parameter_optimization.png\n")
cat("   - Line plots showing how mean F1 varies with parameter values\n")
cat("   - Optimal values are highlighted with larger points\n\n")

cat("6. parameter_sensitivity_heatmap.png\n")
cat("   - Heatmap showing F1 scores for all parameter combinations\n")
cat("   - Useful for identifying robust parameter regions\n\n")

cat("7. power_analysis_ci.png\n")
cat("   - Statistical power with 95% confidence intervals\n")
cat("   - Shows reliability of the algorithm for detecting true differences\n\n")

cat("8. sample_size_requirements.png\n")
cat("   - Shows minimum network size needed for adequate statistical power\n")
cat("   - Helps determine appropriate network sizes for future studies\n\n")

cat("9. runtime_scaling_fitted.png\n")
cat("   - Log-log plot with fitted curves showing computational complexity\n")
cat("   - Confirms expected O(N²) scaling behavior\n\n")

cat("10. algorithm_efficiency.png\n")
cat("    - Shows F1 score achieved per second of computation\n")
cat("    - Useful for comparing cost-effectiveness across conditions\n\n")

cat("11. algorithm_stability.png\n")
cat("    - Coefficient of variation for different metrics\n")
cat("    - Lower values indicate more consistent/stable performance\n\n")

cat("12. outlier_analysis.png\n")
cat("    - Shows frequency of unusual results (|z-score| > 2)\n")
cat("    - Helps identify conditions prone to erratic behavior\n\n")

cat("13. perturbation_comparison_matrix.png\n")
cat("    - Pairwise differences in mean F1 between perturbation types\n")
cat("    - Quickly shows which perturbations are easier/harder to detect\n\n")

cat("14. radar_comparison.png\n")
cat("    - Multi-dimensional comparison across six performance metrics\n")
cat("    - Provides holistic view of algorithm performance\n\n")

cat("15. edge_detection_ratio.png\n")
cat("    - Distribution of predicted/true edge ratios\n")
cat("    - Shows tendency to over- or under-detect across conditions\n\n")

cat("16. pr_tradeoff_by_parameter.png\n")
cat("    - Precision-recall scatter plot with parameter values labeled\n")
cat("    - Reveals how parameter tuning affects the precision-recall trade-off\n\n")

cat("17. summary_dashboard.png\n")
cat("    - Four-panel dashboard summarizing key metrics\n")
cat("    - Provides quick overview of main results\n\n")

sink()

cat("\n============================================\n")
cat("Comprehensive analysis complete!\n")
cat(sprintf("Generated %d plots in: %s\n", 17, comprehensive_plots_dir))
cat(sprintf("Detailed report saved to: %s\n", report_file))
cat(sprintf("Plot descriptions saved to: %s\n", descriptions_file))
cat("============================================\n")
