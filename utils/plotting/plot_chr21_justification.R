#!/usr/bin/env Rscript
# Diagnostic plots justifying the symmetric MIN_CLUSTER_SUPPORT filter
# and characterizing the chr21 co-occurrence graph.
#
# Usage: Rscript plot_chr21_justification.R [input_dir] [output_dir]
# Expects the four CSVs produced by prepare_plot_data.sh in input_dir.

suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "."
output_dir <- if (length(args) >= 2) args[2] else "."
label      <- if (length(args) >= 3) args[3] else "chr21"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

N_INDIVIDUALS <- 2548
MIN_SUPPORT   <- 25
MAX_SUPPORT   <- N_INDIVIDUALS - MIN_SUPPORT  # symmetric ceiling, derived not chosen

theme_set(theme_minimal(base_size = 13))

## ---- Plot 1: node support distribution, with the symmetric filter's
## kept region shaded, on a log-y scale (heavy-tailed distribution) ----

support_dist <- read.csv(file.path(input_dir, "support_distribution.csv"))

p1 <- ggplot(support_dist, aes(x = support, y = node_count)) +
    annotate("rect", xmin = MIN_SUPPORT, xmax = MAX_SUPPORT,
             ymin = 0, ymax = Inf, fill = "steelblue", alpha = 0.08) +
    geom_col(width = 1, fill = "grey40") +
    geom_vline(xintercept = MIN_SUPPORT, linetype = "dashed", color = "firebrick") +
    geom_vline(xintercept = MAX_SUPPORT, linetype = "dashed", color = "firebrick") +
    scale_y_log10(labels = scales::comma) +
    labs(
        title = paste(label, "haploblock cluster support distribution"),
        subtitle = paste0(
            "Shaded region kept by symmetric filter: min(support, N-support) >= ", MIN_SUPPORT,
            "  (", MIN_SUPPORT, " to ", MAX_SUPPORT, " of ", N_INDIVIDUALS, " individuals)"
        ),
        x = "Support (individuals carrying the cluster)",
        y = "Number of nodes (log scale)"
    )

ggsave(file.path(output_dir, "01_support_distribution.png"), p1, width = 9, height = 5.5, dpi = 150)

## ---- Plot 2: diminishing-returns curve -- avg_k vs MIN threshold ----
## This is the plot that most directly justifies MIN_CLUSTER_SUPPORT=25:
## shows how little average k grows as the threshold is relaxed further.

curve_data <- read.csv(file.path(input_dir, "support_threshold_curve.csv"))
curve_data <- curve_data[order(curve_data$min_threshold), ]

max_possible_k <- max(curve_data$avg_k)
k_at_25 <- curve_data$avg_k[curve_data$min_threshold == MIN_SUPPORT]

p2 <- ggplot(curve_data, aes(x = min_threshold, y = avg_k)) +
    geom_line(color = "grey30") +
    geom_point(size = 0.6, alpha = 0.4) +
    geom_vline(xintercept = MIN_SUPPORT, linetype = "dashed", color = "firebrick") +
    geom_hline(yintercept = k_at_25, linetype = "dotted", color = "firebrick", alpha = 0.6) +
    scale_x_log10() +
    labs(
        title = "Diminishing returns: average per-individual k vs. MIN_CLUSTER_SUPPORT",
        subtitle = sprintf(
            "At threshold=%d: avg_k=%.1f (%.1f%% of the max achievable %.1f at threshold=1)",
            MIN_SUPPORT, k_at_25, 100 * k_at_25 / max_possible_k, max_possible_k
        ),
        x = "MIN_CLUSTER_SUPPORT threshold (log scale)",
        y = "Average surviving nodes per individual (k)"
    )

ggsave(file.path(output_dir, "02_threshold_diminishing_returns.png"), p2, width = 9, height = 5.5, dpi = 150)

## ---- Plot 3: edge weight distribution (log-log) ----

edge_weights <- read.csv(file.path(input_dir, "edge_weight_distribution.csv"))

p3 <- ggplot(edge_weights, aes(x = weight, y = edge_count)) +
    geom_point(size = 0.8, alpha = 0.5, color = "darkgreen") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(
        title = paste(label, "edge weight distribution"),
        subtitle = "Long-tailed: many weak co-occurrences, few very strong ones",
        x = "Edge weight (individuals sharing both clusters, log scale)",
        y = "Number of edges (log scale)"
    )

ggsave(file.path(output_dir, "03_edge_weight_distribution.png"), p3, width = 9, height = 5.5, dpi = 150)

## ---- Plot 4: node degree distribution ----
## Watch for a spike at the maximum possible degree (N_nodes - 1) --
## that's the mega-hub signature discussed alongside this analysis.

node_degree <- read.csv(file.path(input_dir, "node_degree_distribution.csv"))

p4 <- ggplot(node_degree, aes(x = degree)) +
    geom_histogram(bins = 60, fill = "purple4", alpha = 0.7) +
    labs(
        title = paste(label, "node degree distribution"),
        subtitle = "A spike near the maximum possible degree suggests population-stratification-driven hub nodes, not LD",
        x = "Degree (number of distinct co-occurrence partners)",
        y = "Number of nodes"
    )

ggsave(file.path(output_dir, "04_node_degree_distribution.png"), p4, width = 9, height = 5.5, dpi = 150)

cat("Saved 4 plots to:", output_dir, "\n")
