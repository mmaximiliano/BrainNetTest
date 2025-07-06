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
    pattern = "^summary_N_00010.*\\.csv$",
    full.names = TRUE
  )
  
  # Read and combine all files
  all_results <- map_df(summary_files, read_csv, show_col_types = FALSE)
  
  # Calculate actual lambda values
  all_results <- all_results %>%
    mutate(
      # Handle const_* types where lambda_mult is NA
      lambda_actual = case_when(
        perturbation_type %in% c("const_high", "const_low") ~ p_const, # the edge probability is constant
        TRUE ~ exp(-lambda_base * lambda_mult * 1) # probability depends on distance and lambda
      ),
      # Round lambda_actual to 3 decimals for display
      lambda_rounded = round(lambda_actual, 3),
      # Create a unique identifier for each parameter combination
      param_id = paste(N, perturbation_type, p_const, lambda_actual, sep = "_")
    )
  
  return(all_results)
}

# Function to aggregate results across multiple runs
aggregate_results <- function(results_df) {
  results_df %>%
    group_by(N, perturbation_type, lambda_base, lambda_mult, p_const, lambda_actual, lambda_rounded) %>%
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
  # Separate constant and non-constant perturbations
  data_const <- data_n %>% filter(perturbation_type %in% c("const_high", "const_low"))
  data_non_const <- data_n %>% filter(!perturbation_type %in% c("const_high", "const_low"))
  
  # Define color palette
  color_palette <- scale_color_viridis_d(end = 0.8)
  
  # Define common theme with white background
  theme_white <- theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA)
    )
  
  # Create plots for both subsets
  plot_subsets <- list(
    list(data = data_const, suffix = "const", title_suffix = " (Constant Perturbations)"),
    list(data = data_non_const, suffix = "non_const", title_suffix = " (Non-Constant Perturbations)")
  )
  
  # Create plot directory
  plot_dir <- file.path(output_dir, sprintf("N_%05d", n_value))
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (subset_info in plot_subsets) {
    subset_data <- subset_info$data
    suffix <- subset_info$suffix
    title_suffix <- subset_info$title_suffix
    
    # Skip if no data for this subset
    if (nrow(subset_data) == 0) next
    
    # Calculate additional rate metrics
    subset_data <- subset_data %>%
      mutate(
        # True Positive Rate (TPR / Recall) = TP / (TP + FN)
        TPR = TP_mean / (TP_mean + FN_mean),
        # False Negative Rate (FNR) = FN / (TP + FN)
        FNR = FN_mean / (TP_mean + FN_mean),
        # False Positive Rate (FPR) = FP / (FP + TN)
        FPR = FP_mean / (FP_mean + TN_mean),
        # True Negative Rate (TNR / Specificity) = TN / (FP + TN)
        TNR = TN_mean / (FP_mean + TN_mean)
      )
    
    # 1. Lambda vs TPR (instead of TP)
    p1 <- ggplot(subset_data, aes(x = lambda_actual, y = TPR, color = as.factor(lambda_rounded))) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = pmax(0, TPR - recall_sd), ymax = pmin(1, TPR + recall_sd)), width = 0.02, linewidth = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      labs(
        title = sprintf("True Positive Rate vs P(edge) (N = %d)%s", n_value, title_suffix),
        subtitle = "Proportion of actual critical links correctly identified (Recall)",
        x = "P(edge) in log scale",
        y = "True Positive Rate",
        color = "P(edge)"
      ) +
      theme_white +
      color_palette +
      theme(legend.position = "right")
    
    # 2. Lambda vs FPR (instead of FP)
    p2 <- ggplot(subset_data, aes(x = lambda_actual, y = FPR, color = as.factor(lambda_rounded))) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = pmax(0, FPR - FPR), ymax = pmin(1, FPR + FPR)), width = 0.02, linewidth = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      labs(
        title = sprintf("False Positive Rate vs P(edge) (N = %d)%s", n_value, title_suffix),
        subtitle = "Proportion of non-critical links incorrectly identified as critical",
        x = "P(edge) in log scale",
        y = "False Positive Rate",
        color = "P(edge)"
      ) +
      theme_white +
      color_palette
    
    # 3. Lambda vs FNR (instead of FN)
    p3 <- ggplot(subset_data, aes(x = lambda_actual, y = FNR, color = as.factor(lambda_rounded))) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = pmax(0, FNR - FNR), ymax = pmin(1, FNR + FNR)), width = 0.02, linewidth = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      labs(
        title = sprintf("False Negative Rate vs P(edge) (N = %d)%s", n_value, title_suffix),
        subtitle = "Proportion of actual critical links missed",
        x = "P(edge) in log scale",
        y = "False Negative Rate",
        color = "P(edge)"
      ) +
      theme_white +
      color_palette
    
    # 4. Lambda vs TNR (instead of TN)
    p4 <- ggplot(subset_data, aes(x = lambda_actual, y = TNR, color = as.factor(lambda_rounded))) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = pmax(0, TNR - TNR), ymax = pmin(1, TNR + TNR)), width = 0.02, linewidth = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      labs(
        title = sprintf("True Negative Rate vs P(edge) (N = %d)%s", n_value, title_suffix),
        subtitle = "Proportion of non-critical links correctly identified (Specificity)",
        x = "P(edge) in log scale",
        y = "True Negative Rate",
        color = "P(edge)"
      ) +
      theme_white +
      color_palette
    
    # 5. Lambda vs Recall - already using rates, so keep this plot
    p5 <- ggplot(subset_data, aes(x = lambda_actual, y = recall_mean, color = as.factor(lambda_rounded))) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = pmax(0, recall_mean - recall_sd), ymax = pmin(1, recall_mean + recall_sd)), width = 0.02, linewidth = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      labs(
        title = sprintf("Recall vs P(edge) (N = %d)%s", n_value, title_suffix),
        subtitle = "Proportion of actual critical links correctly identified",
        x = "P(edge) in log scale",
        y = "Recall",
        color = "P(edge)"
      ) +
      theme_white +
      color_palette
    
    # 6. Lambda vs Precision - already using rates, so keep this plot
    p6 <- ggplot(subset_data, aes(x = lambda_actual, y = precision_mean, color = as.factor(lambda_rounded))) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = pmax(0, precision_mean - precision_sd), ymax = pmin(1, precision_mean + precision_sd)), width = 0.02, linewidth = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      labs(
        title = sprintf("Precision vs P(edge) (N = %d)%s", n_value, title_suffix),
        subtitle = "Proportion of predicted critical links that are correct",
        x = "P(edge) in log scale",
        y = "Precision",
        color = "P(edge)"
      ) +
      theme_white +
      color_palette
    
    # 7. Lambda vs F1
    p7 <- ggplot(subset_data, aes(x = lambda_actual, y = f1_mean, color = as.factor(lambda_rounded))) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = pmax(0, f1_mean - f1_sd), ymax = pmin(1, f1_mean + f1_sd)), width = 0.02, linewidth = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      labs(
        title = sprintf("F1 Score vs P(edge) (N = %d)%s", n_value, title_suffix),
        subtitle = "Harmonic mean of precision and recall",
        x = "P(edge) in log scale",
        y = "F1 Score",
        color = "P(edge)"
      ) +
      theme_white +
      color_palette
    
    # 8. Lambda vs MCC
    p8 <- ggplot(subset_data, aes(x = lambda_actual, y = mcc_mean, color = as.factor(lambda_rounded))) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = pmax(-1, mcc_mean - mcc_sd), ymax = pmin(1, mcc_mean + mcc_sd)), width = 0.02, linewidth = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(limits = c(-1, 1)) +
      labs(
        title = sprintf("Matthews Correlation Coefficient vs P(edge) (N = %d)%s", n_value, title_suffix),
        subtitle = "Balanced measure of binary classification quality",
        x = "P(edge) in log scale",
        y = "MCC",
        color = "P(edge)"
      ) +
      theme_white +
      color_palette
    
    # 9. Lambda vs Predicted Ratio
    p9 <- ggplot(subset_data, aes(x = lambda_actual, y = predicted_ratio, color = as.factor(lambda_rounded))) +
      geom_point(size = 3) +
      geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      labs(
        title = sprintf("Predicted/True Ratio vs P(edge) (N = %d)%s", n_value, title_suffix),
        subtitle = "Ratio of predicted to actual number of critical links",
        x = "P(edge) in log scale",
        y = "Predicted/True Ratio",
        color = "P(edge)"
      ) +
      theme_white +
      color_palette
    
    # 10. Lambda vs TP-FN Tradeoff
    p10 <- ggplot(subset_data, aes(x = lambda_actual, y = tp_fn_tradeoff, color = as.factor(lambda_rounded))) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      labs(
        title = sprintf("TP-FN Tradeoff vs P(edge) (N = %d)%s", n_value, title_suffix),
        subtitle = "Difference between true positives and false negatives",
        x = "P(edge) in log scale",
        y = "TP - FN",
        color = "P(edge)"
      ) +
      theme_white +
      color_palette
    
    # 11. Outlier Analysis - Modified to handle zeros
    p11 <- ggplot(subset_data, aes(x = as.factor(lambda_rounded), y = runtime_mean)) +
      geom_boxplot(aes(fill = as.factor(lambda_rounded))) +
      geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
      # Use pseudo_log to handle zeros or very small values properly
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), 
                        labels = scales::label_number()) +
      labs(
        title = sprintf("Runtime Distribution vs P(edge) (N = %d)%s", n_value, title_suffix),
        subtitle = "Computational time analysis with outlier detection",
        x = "P(edge)",
        y = "Runtime (seconds, pseudo-log scale)",
        fill = "P(edge)"
      ) +
      theme_white +
      scale_fill_viridis_d(end = 0.8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Function to create a confusion matrix plot using rates
    create_confusion_matrix_plot <- function(subset_data, n_value, title_suffix) {
      # Calculate confusion matrix rates (normalize by row so each row sums to 1)
      # True Positive Rate (Sensitivity/Recall) and False Negative Rate
      TPR <- subset_data$TP_mean / (subset_data$TP_mean + subset_data$FN_mean)
      FNR <- subset_data$FN_mean / (subset_data$TP_mean + subset_data$FN_mean)
      
      # False Positive Rate and True Negative Rate (Specificity)
      FPR <- subset_data$FP_mean / (subset_data$FP_mean + subset_data$TN_mean)
      TNR <- subset_data$TN_mean / (subset_data$FP_mean + subset_data$TN_mean)
      
      # Create a data frame for plotting
      cm_data <- data.frame(
        lambda_actual = rep(subset_data$lambda_actual, each = 4),
        lambda_rounded = rep(subset_data$lambda_rounded, each = 4),
        true_class = rep(c("Positive", "Negative"), each = 2, times = nrow(subset_data)),
        predicted_class = rep(c("Positive", "Negative"), times = 2 * nrow(subset_data)),
        rate = c(rbind(TPR, FNR, FPR, TNR))  # Interleave rates to match structure
      )
      
      # Define color palette with better contrast
      heatmap_colors <- scale_fill_gradient2(
        low = "white", 
        mid = "#ffcc99", 
        high = "#ff6600", 
        midpoint = 0.5,
        limits = c(0, 1),
        name = "Rate"
      )
      
      # Create a custom labeller function that explicitly formats the labels
      p_edge_labeller <- labeller(lambda_rounded = function(x) paste("P(edge):", x))
      
      # Create a faceted plot for each lambda value
      confusion_plot <- ggplot(cm_data, aes(x = predicted_class, y = true_class, fill = rate)) +
        geom_tile(color = "black", linewidth = 0.5) +
        facet_wrap(~ lambda_rounded, labeller = p_edge_labeller, ncol = 4) +
        geom_text(aes(label = sprintf("%.2f", rate)), color = "black", size = 3) +
        heatmap_colors +
        scale_x_discrete(position = "top") +
        labs(
          title = sprintf("Confusion Matrix Rates vs P(edge) (N = %d)%s", n_value, title_suffix),
          subtitle = "Rows are normalized to sum to 1 (TPR+FNR=1, FPR+TNR=1)",
          x = "Predicted Class",
          y = "True Class"
        ) +
        theme_minimal() +
        theme(
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          axis.text = element_text(size = 10),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "#f0f0f0", color = "black")
        )
      
      return(confusion_plot)
    }
    
    # Create confusion matrix plot
    p_confusion <- create_confusion_matrix_plot(subset_data, n_value, title_suffix)
    
    # Save confusion matrix plot
    ggsave(file.path(plot_dir, sprintf("12_confusion_matrix_%s.png", suffix)), p_confusion, 
           width = 12, height = 8, dpi = 300)
    
    # Save individual plots
    ggsave(file.path(plot_dir, sprintf("01_lambda_vs_tpr_%s.png", suffix)), p1, width = 10, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, sprintf("02_lambda_vs_fpr_%s.png", suffix)), p2, width = 10, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, sprintf("03_lambda_vs_fnr_%s.png", suffix)), p3, width = 10, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, sprintf("04_lambda_vs_tnr_%s.png", suffix)), p4, width = 10, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, sprintf("05_lambda_vs_recall_%s.png", suffix)), p5, width = 10, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, sprintf("06_lambda_vs_precision_%s.png", suffix)), p6, width = 10, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, sprintf("07_lambda_vs_f1_%s.png", suffix)), p7, width = 10, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, sprintf("08_lambda_vs_mcc_%s.png", suffix)), p8, width = 10, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, sprintf("09_lambda_vs_predicted_ratio_%s.png", suffix)), p9, width = 10, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, sprintf("10_lambda_vs_tp_fn_tradeoff_%s.png", suffix)), p10, width = 10, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, sprintf("11_outlier_analysis_%s.png", suffix)), p11, width = 10, height = 6, dpi = 300)
    
    # Create combined plot
    combined <- (p1 + p2) / (p3 + p4) / (p5 + p6) / (p7 + p8) +
      plot_annotation(
        title = sprintf("Ring Experiment Results Summary (N = %d)%s", n_value, title_suffix),
        subtitle = sprintf("Based on %d runs per parameter combination", unique(subset_data$n_runs)),
        theme = theme(plot.title = element_text(size = 16, face = "bold"))
      )
    
    ggsave(
      file.path(plot_dir, sprintf("00_combined_summary_%s.png", suffix)), 
      combined, 
      width = 20, 
      height = 30, 
      dpi = 300
    )
    
    # Create a separate combined plot that includes the confusion matrix
    combined_with_cm <- (p1 + p2) / (p3 + p4) / (p_confusion) +
      plot_annotation(
        title = sprintf("Ring Experiment Confusion Matrix Analysis (N = %d)%s", n_value, title_suffix),
        subtitle = sprintf("Based on %d runs per parameter combination", unique(subset_data$n_runs)),
        theme = theme(plot.title = element_text(size = 16, face = "bold"))
      )
    
    ggsave(
      file.path(plot_dir, sprintf("00_confusion_analysis_%s.png", suffix)), 
      combined_with_cm, 
      width = 20, 
      height = 20, 
      dpi = 300
    )
  }
}

# Function to create combined plots across all N values
create_combined_plots <- function(aggregated_data, output_dir) {
  # Separate constant and non-constant perturbations
  data_const <- aggregated_data %>% filter(perturbation_type %in% c("const_high", "const_low"))
  data_non_const <- aggregated_data %>% filter(!perturbation_type %in% c("const_high", "const_low"))
  
  # Calculate additional rate metrics for both datasets
  if(nrow(data_const) > 0) {
    data_const <- data_const %>%
      mutate(
        TPR = TP_mean / (TP_mean + FN_mean),
        FNR = FN_mean / (TP_mean + FN_mean),
        FPR = FP_mean / (FP_mean + TN_mean),
        TNR = TN_mean / (FP_mean + TN_mean)
      )
  }
  
  if(nrow(data_non_const) > 0) {
    data_non_const <- data_non_const %>%
      mutate(
        TPR = TP_mean / (TP_mean + FN_mean),
        FNR = FN_mean / (TP_mean + FN_mean),
        FPR = FP_mean / (FP_mean + TN_mean),
        TNR = TN_mean / (FP_mean + TN_mean)
      )
  }
  
  combined_dir <- file.path(output_dir, "combined_all_N")
  dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define common theme with white background
  theme_white <- theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA)
    )
  
  # Create plots for both subsets
  plot_subsets <- list(
    list(data = data_const, suffix = "const", title_suffix = " (Constant Perturbations)"),
    list(data = data_non_const, suffix = "non_const", title_suffix = " (Non-Constant Perturbations)")
  )
  
  for (subset_info in plot_subsets) {
    subset_data <- subset_info$data
    suffix <- subset_info$suffix
    title_suffix <- subset_info$title_suffix
    
    # Skip if no data for this subset
    if (nrow(subset_data) == 0) next
    
    # Create cross-network plots for rates instead of counts
    # TPR across all N
    p_tpr_all <- ggplot(subset_data, aes(x = lambda_actual, y = TPR, color = as.factor(N))) +
      geom_point(size = 3) +
      geom_line(aes(group = as.factor(N)), alpha = 0.5, linewidth = 0.7) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      labs(
        title = paste0("True Positive Rate vs P(edge) for Different Network Sizes", title_suffix),
        subtitle = "Proportion of actual critical links correctly identified across network scales",
        x = "P(edge) in log scale",
        y = "True Positive Rate",
        color = "Network Size (N)"
      ) +
      theme_white +
      scale_color_viridis_d(end = 0.8)
    
    # FPR across all N
    p_fpr_all <- ggplot(subset_data, aes(x = lambda_actual, y = FPR, color = as.factor(N))) +
      geom_point(size = 3) +
      geom_line(aes(group = as.factor(N)), alpha = 0.5, linewidth = 0.7) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      labs(
        title = paste0("False Positive Rate vs P(edge) for Different Network Sizes", title_suffix),
        subtitle = "Proportion of non-critical links incorrectly identified as critical across network scales",
        x = "P(edge) in log scale",
        y = "False Positive Rate",
        color = "Network Size (N)"
      ) +
      theme_white +
      scale_color_viridis_d(end = 0.8)
    
    # 1. F1 Score across all N - Add this missing plot
    p_f1_all <- ggplot(subset_data, aes(x = lambda_actual, y = f1_mean, color = as.factor(N))) +
      geom_point(size = 3) +
      geom_line(aes(group = as.factor(N)), alpha = 0.5, linewidth = 0.7) +
      geom_errorbar(aes(ymin = f1_mean - f1_sd, ymax = f1_mean + f1_sd), width = 0.02, alpha = 0.5, linewidth = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      labs(
        title = paste0("F1 Score vs P(edge) for Different Network Sizes", title_suffix),
        subtitle = "Performance comparison across network scales",
        x = "P(edge) in log scale",
        y = "F1 Score",
        color = "Network Size (N)"
      ) +
      theme_white +
      scale_color_viridis_d(end = 0.8)
    
    # 2. MCC across all N - Add this missing plot
    p_mcc_all <- ggplot(subset_data, aes(x = lambda_actual, y = mcc_mean, color = as.factor(N))) +
      geom_point(size = 3) +
      geom_line(aes(group = as.factor(N)), alpha = 0.5, linewidth = 0.7) +
      geom_errorbar(aes(ymin = mcc_mean - mcc_sd, ymax = mcc_mean + mcc_sd), width = 0.02, alpha = 0.5, linewidth = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(limits = c(-1, 1)) +
      labs(
        title = paste0("Matthews Correlation Coefficient vs P(edge) for Different Network Sizes", title_suffix),
        subtitle = "Balanced classification performance across network scales",
        x = "P(edge) in log scale",
        y = "MCC",
        color = "Network Size (N)"
      ) +
      theme_white +
      scale_color_viridis_d(end = 0.8)
    
    # Continue with existing plots like recall, precision, etc.
    # 3. Recall across all N
    p_recall_all <- ggplot(subset_data, aes(x = lambda_actual, y = recall_mean, color = as.factor(N))) +
      geom_point(size = 3) +
      geom_line(aes(group = as.factor(N)), alpha = 0.5, linewidth = 0.7) +
      geom_errorbar(aes(ymin = recall_mean - recall_sd, ymax = recall_mean + recall_sd), width = 0.02, alpha = 0.5, linewidth = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      labs(
        title = paste0("Recall vs P(edge) for Different Network Sizes", title_suffix),
        subtitle = "Sensitivity across network scales",
        x = "P(edge) in log scale",
        y = "Recall",
        color = "Network Size (N)"
      ) +
      theme_white +
      scale_color_viridis_d(end = 0.8)
    
    # 4. Precision across all N
    p_precision_all <- ggplot(subset_data, aes(x = lambda_actual, y = precision_mean, color = as.factor(N))) +
      geom_point(size = 3) +
      geom_line(aes(group = as.factor(N)), alpha = 0.5, linewidth = 0.7) +
      geom_errorbar(aes(ymin = precision_mean - precision_sd, ymax = precision_mean + precision_sd), width = 0.02, alpha = 0.5, linewidth = 0.5) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      labs(
        title = paste0("Precision vs P(edge) for Different Network Sizes", title_suffix),
        subtitle = "Positive predictive value across network scales",
        x = "P(edge) in log scale",
        y = "Precision",
        color = "Network Size (N)"
      ) +
      theme_white +
      scale_color_viridis_d(end = 0.8)
    
    # 5. Runtime scaling - Modified to handle zeros or very small values
    p_runtime_all <- ggplot(subset_data, aes(x = N, y = runtime_mean, color = as.factor(lambda_rounded))) +
      geom_point(size = 3) +
      geom_line(aes(group = as.factor(lambda_rounded)), alpha = 0.5, linewidth = 0.7) +
      geom_errorbar(aes(ymin = pmax(0.001, runtime_mean - runtime_sd), 
                        ymax = runtime_mean + runtime_sd), width = 0.02, alpha = 0.5, linewidth = 0.5) +
      scale_x_log10() +
      # Use pseudo_log to handle zeros or very small values properly
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                        labels = scales::label_number()) +
      labs(
        title = paste0("Runtime Scaling with Network Size", title_suffix),
        subtitle = "Computational complexity analysis",
        x = "Network Size (N, log scale)",
        y = "Runtime (seconds, pseudo-log scale)",
        color = "P(edge)"
      ) +
      theme_white +
      scale_color_viridis_d(end = 0.8)
    
    # 6. Predicted ratio across all N
    p_ratio_all <- ggplot(subset_data, aes(x = lambda_actual, y = predicted_ratio, color = as.factor(N))) +
      geom_point(size = 3) +
      geom_line(aes(group = as.factor(N)), alpha = 0.5, linewidth = 0.7) +
      geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5, linewidth = 0.7) +
      scale_x_log10(labels = scales::number_format(accuracy = 0.001)) +
      labs(
        title = paste0("Predicted/True Ratio vs Lambda for Different Network Sizes", title_suffix),
        subtitle = "Over/under-prediction across network scales",
        x = "Lambda Value (log scale)",
        y = "Predicted/True Ratio",
        color = "Network Size (N)"
      ) +
      theme_white +
      scale_color_viridis_d(end = 0.8)
    
    # Save plots
    ggsave(file.path(combined_dir, sprintf("tpr_all_n_%s.png", suffix)), p_tpr_all, width = 12, height = 8, dpi = 300)
    ggsave(file.path(combined_dir, sprintf("fpr_all_n_%s.png", suffix)), p_fpr_all, width = 12, height = 8, dpi = 300)
    ggsave(file.path(combined_dir, sprintf("f1_score_all_n_%s.png", suffix)), p_f1_all, width = 12, height = 8, dpi = 300)
    ggsave(file.path(combined_dir, sprintf("mcc_all_n_%s.png", suffix)), p_mcc_all, width = 12, height = 8, dpi = 300)
    ggsave(file.path(combined_dir, sprintf("recall_all_n_%s.png", suffix)), p_recall_all, width = 12, height = 8, dpi = 300)
    ggsave(file.path(combined_dir, sprintf("precision_all_n_%s.png", suffix)), p_precision_all, width = 12, height = 8, dpi = 300)
    ggsave(file.path(combined_dir, sprintf("runtime_scaling_%s.png", suffix)), p_runtime_all, width = 12, height = 8, dpi = 300)
    ggsave(file.path(combined_dir, sprintf("predicted_ratio_all_n_%s.png", suffix)), p_ratio_all, width = 12, height = 8, dpi = 300)
    
    # Create summary dashboard
    dashboard <- (p_tpr_all + p_fpr_all) / (p_recall_all + p_precision_all) / (p_f1_all + p_mcc_all) +
      plot_annotation(
        title = paste0("Ring Experiment Results: Cross-Network Comparison", title_suffix),
        subtitle = "Performance metrics across different network sizes and edge probabilities",
        theme = theme(plot.title = element_text(size = 18, face = "bold"))
      )
    
    ggsave(
      file.path(combined_dir, sprintf("dashboard_all_metrics_%s.png", suffix)), 
      dashboard, 
      width = 24, 
      height = 18, 
      dpi = 300
    )
  }
}

# Main analysis function
analyze_ring_experiment <- function(experiment_dir) {
  # Read results
  cat("Reading experiment results...\n")
  results <- read_experiment_results(experiment_dir)
  
  # Aggregate results
  cat("Aggregating results across runs...\n")
  aggregated <- aggregate_results(results)
  
  # Create output directory in the same location as this script
  script_dir <- dirname(sys.frame(1)$ofile)
  if (is.null(script_dir)) {
    # If running interactively, use the current working directory
    script_dir <- getwd()
  }
  output_dir <- file.path(script_dir, "ring_experiment_analysis_plots")
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
    sprintf("- Lambda values tested: %s", paste(sort(unique(round(aggregated_data$lambda_actual, 3))), collapse = ", ")),
    sprintf("- Perturbation types: %s", paste(sort(unique(aggregated_data$perturbation_type)), collapse = ", ")),
    sprintf("- Total parameter combinations: %d", nrow(aggregated_data)),
    "",
    "## Best Performing Parameters",
    ""
  )
  
  # Find best parameters for each N
  best_by_n <- aggregated_data %>%
    group_by(N) %>%
    slice_max(f1_mean, n = 1) %>%
    select(N, perturbation_type, lambda_rounded, f1_mean, mcc_mean, recall_mean, precision_mean)
  
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
