#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Set working directory and define paths
results_dir <- "/faststorage/project/EcoGenetics/people/jilong/local_adaptation/swpf_threshold/results/EntNic_lr_distribution"
output_dir <- "/faststorage/project/EcoGenetics/people/jilong/local_adaptation/local_adaptation_arthropods_denmark/swpf_threshold"

# Read the LR summary statistics
lr_data <- read.table(file.path(results_dir, "LR_summary_statistics.tsv"), 
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Clean population names for better visualization
lr_data$pop_clean <- gsub("EntNic_", "", lr_data$population)
lr_data$pop_clean <- gsub("-", "_", lr_data$pop_clean)

# Create output directory for plots
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. Basic distribution statistics comparison
# ============================================================================

cat("Creating basic distribution comparison plots...\n")

# Plot 1: Mean vs Median LR across populations
p1 <- ggplot(lr_data, aes(x = reorder(pop_clean, mean_LR))) +
  geom_col(aes(y = mean_LR), fill = "steelblue", alpha = 0.7, width = 0.8) +
  geom_point(aes(y = median_LR), color = "red", size = 1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mean LR (bars) vs Median LR (red points) by Population",
       x = "Population", y = "Likelihood Ratio",
       subtitle = "EntNic Species - SweepFinder2 Results") +
  scale_y_continuous(limits = c(0, max(lr_data$mean_LR) * 1.05))

ggsave(file.path(output_dir, "plots", "mean_vs_median_LR.png"), p1, 
       width = 14, height = 8, dpi = 300)

# Plot 2: Distribution of high LR sites (>1, >5, >10)
high_lr_data <- lr_data %>%
  select(pop_clean, pct_LR_gt_1, pct_LR_gt_5, pct_LR_gt_10) %>%
  pivot_longer(cols = starts_with("pct_LR_gt"), 
               names_to = "threshold", values_to = "percentage") %>%
  mutate(threshold = case_when(
    threshold == "pct_LR_gt_1" ~ "LR > 1",
    threshold == "pct_LR_gt_5" ~ "LR > 5",
    threshold == "pct_LR_gt_10" ~ "LR > 10"
  ))

p2 <- ggplot(high_lr_data, aes(x = reorder(pop_clean, percentage), 
                               y = percentage, fill = threshold)) +
  geom_col(position = "dodge", alpha = 0.8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  labs(title = "Percentage of Sites with High Likelihood Ratios",
       x = "Population", y = "Percentage of Sites (%)",
       fill = "Threshold") +
  scale_fill_manual(values = c("lightblue", "orange", "darkred"))

ggsave(file.path(output_dir, "plots", "high_LR_sites_percentage.png"), p2, 
       width = 16, height = 8, dpi = 300)

# ============================================================================
# 2. Percentile-based distribution visualization
# ============================================================================

cat("Creating percentile distribution plots...\n")

# Extract percentile columns (q0.1_LR to q99.9_LR)
percentile_cols <- grep("^q[0-9]+\\.[0-9]_LR$", names(lr_data), value = TRUE)
percentile_data <- lr_data %>%
  select(pop_clean, all_of(percentile_cols)) %>%
  pivot_longer(cols = all_of(percentile_cols), 
               names_to = "percentile", values_to = "LR_value") %>%
  mutate(percentile_num = as.numeric(gsub("q([0-9]+\\.[0-9])_LR", "\\1", percentile)))

# Plot 3: Percentile curves for all populations
p3 <- ggplot(percentile_data, aes(x = percentile_num, y = LR_value, color = pop_clean)) +
  geom_line(alpha = 0.7, size = 0.5) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "LR Distribution Percentiles Across All Populations",
       x = "Percentile", y = "Likelihood Ratio",
       subtitle = "Each line represents one population") +
  scale_y_log10() +
  annotation_logticks(sides = "l")

ggsave(file.path(output_dir, "plots", "percentile_curves_all_populations.png"), p3, 
       width = 12, height = 8, dpi = 300)

# Plot 4: Focus on high percentiles (95th-99.9th)
high_percentile_data <- percentile_data %>%
  filter(percentile_num >= 95)

p4 <- ggplot(high_percentile_data, aes(x = percentile_num, y = LR_value, color = pop_clean)) +
  geom_line(alpha = 0.8, size = 0.8) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "High Percentile LR Distribution (95th-99.9th percentiles)",
       x = "Percentile", y = "Likelihood Ratio",
       subtitle = "Focus on tail of distribution where selection signals are strongest") +
  scale_y_log10() +
  annotation_logticks(sides = "l")

ggsave(file.path(output_dir, "plots", "high_percentiles_focus.png"), p4, 
       width = 12, height = 8, dpi = 300)

# ============================================================================
# 3. Population comparison heatmaps
# ============================================================================

cat("Creating heatmap visualizations...\n")

# Plot 5: Heatmap of key percentiles across populations
key_percentiles <- c("q50.0_LR", "q90.0_LR", "q95.0_LR", "q99.0_LR", "q99.9_LR")
heatmap_data <- lr_data %>%
  select(pop_clean, all_of(key_percentiles)) %>%
  pivot_longer(cols = all_of(key_percentiles), 
               names_to = "percentile", values_to = "LR_value") %>%
  mutate(percentile = case_when(
    percentile == "q50.0_LR" ~ "50th",
    percentile == "q90.0_LR" ~ "90th", 
    percentile == "q95.0_LR" ~ "95th",
    percentile == "q99.0_LR" ~ "99th",
    percentile == "q99.9_LR" ~ "99.9th"
  ))

p5 <- ggplot(heatmap_data, aes(x = percentile, y = reorder(pop_clean, LR_value), 
                               fill = log10(LR_value + 1))) +
  geom_tile() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Heatmap of Key LR Percentiles by Population",
       x = "Percentile", y = "Population",
       fill = "log10(LR + 1)") +
  scale_fill_gradient(low = "white", high = "darkblue")

ggsave(file.path(output_dir, "plots", "percentile_heatmap.png"), p5, 
       width = 10, height = 16, dpi = 300)

# ============================================================================
# 4. Summary statistics
# ============================================================================

cat("Generating summary statistics...\n")

# Create summary table
summary_stats <- lr_data %>%
  select(pop_clean, total_sites, mean_LR, median_LR, max_LR, 
         pct_LR_gt_1, pct_LR_gt_5, pct_LR_gt_10) %>%
  arrange(desc(mean_LR))

write.table(summary_stats, 
            file.path(output_dir, "plots", "population_LR_summary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print top 10 populations by mean LR
cat("\nTop 10 populations by mean LR:\n")
print(head(summary_stats, 10))

# Print population with highest max LR
max_lr_pop <- lr_data[which.max(lr_data$max_LR), ]
cat("\nPopulation with highest maximum LR:\n")
cat(sprintf("%s: Max LR = %.2f\n", max_lr_pop$pop_clean, max_lr_pop$max_LR))

# ============================================================================
# 5. Boxplot comparison of key statistics
# ============================================================================

# Plot 6: Boxplot of basic statistics
basic_stats <- lr_data %>%
  select(pop_clean, mean_LR, median_LR, max_LR) %>%
  pivot_longer(cols = c(mean_LR, median_LR, max_LR),
               names_to = "statistic", values_to = "value") %>%
  mutate(statistic = case_when(
    statistic == "mean_LR" ~ "Mean",
    statistic == "median_LR" ~ "Median", 
    statistic == "max_LR" ~ "Maximum"
  ))

p6 <- ggplot(basic_stats, aes(x = statistic, y = value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Distribution of LR Statistics Across Populations",
       x = "Statistic", y = "Likelihood Ratio") +
  scale_y_log10() +
  annotation_logticks(sides = "l")

ggsave(file.path(output_dir, "plots", "LR_statistics_boxplot.png"), p6, 
       width = 10, height = 8, dpi = 300)

cat("\nVisualization complete!\n")
cat("All plots saved to:", file.path(output_dir, "plots"), "\n")
cat("Summary table saved to:", file.path(output_dir, "plots", "population_LR_summary.txt"), "\n")