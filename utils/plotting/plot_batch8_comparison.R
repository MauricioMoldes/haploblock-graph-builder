#!/usr/bin/env Rscript
# Cross-chromosome comparison plots for the 8-chromosome batch.
#
# Usage: Rscript plot_batch8_comparison.R [summary_csv] [output_dir]
# Expects chromosome_summary.csv from split_and_summarize_batch8.sh.

suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
summary_csv <- if (length(args) >= 1) args[1] else "chromosome_summary.csv"
output_dir  <- if (length(args) >= 2) args[2] else "."
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 13))

summary_data <- read.csv(summary_csv)

# n_nodes from chromosome_summary.csv is the UNFILTERED count (from
# nodes.csv, before MIN_CLUSTER_SUPPORT). The hub-check plot needs the
# count of nodes actually PRESENT in edges.csv (post-filter) instead --
# otherwise max_degree/n_nodes compares a filtered graph's hub against
# an unrelated, much larger denominator and badly understates hub
# dominance. graph_node_counts.csv supplies the correct figure,
# computed directly from each chromosome's edges.csv.
graph_counts_path <- file.path(dirname(summary_csv), "graph_node_counts.csv")
if (file.exists(graph_counts_path)) {
    graph_counts <- read.csv(graph_counts_path, header = FALSE,
                              col.names = c("chromosome", "n_nodes_in_graph"))
    summary_data <- merge(summary_data, graph_counts, by = "chromosome")
} else {
    warning("graph_node_counts.csv not found -- hub-check plot will use the ",
            "unfiltered nodes.csv count instead, which UNDERSTATES hub dominance. ",
            "Run the graph_node_counts.csv generation step first for a correct plot.")
    summary_data$n_nodes_in_graph <- summary_data$n_nodes
}

# Known chromosome sizes (Mb, GRCh38) and the mean-k-per-individual
# figures reported directly by the pipeline's own workload diagnostic
# in the batch8 run log -- not recomputed here, since the log already
# gives the exact figures the actual algorithm used.
chr_info <- data.frame(
    chromosome = c("chr13", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"),
    size_mb    = c(114,     90,      83,      80,      59,      64,      47,      51),
    mean_k     = c(1888.1,  1800.0,  1674.8,  1674.4,  1357.2,  1459.0,  820.2,   928.4)
)

merged <- merge(summary_data, chr_info, by = "chromosome")

## ---- Plot 1: does mean_k scale with size^2, as assumed throughout? ----

p1 <- ggplot(merged, aes(x = size_mb, y = mean_k, label = chromosome)) +
    geom_point(size = 3, color = "steelblue") +
    geom_text(vjust = -1, size = 3.5) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE,
                color = "grey50", linetype = "dashed", linewidth = 0.6) +
    labs(
        title = "Mean per-individual k vs. chromosome size",
        subtitle = "Dashed line: quadratic fit -- tests the size^2 cost-scaling assumption used throughout this project",
        x = "Chromosome size (Mb)",
        y = "Mean k (nodes/individual, within that chromosome only)"
    )

ggsave(file.path(output_dir, "batch8_01_k_vs_size.png"), p1, width = 8, height = 6, dpi = 150)

## ---- Plot 2: node count vs. edge count per chromosome ----

p2 <- ggplot(merged, aes(x = n_nodes, y = n_edges, label = chromosome, size = size_mb)) +
    geom_point(color = "darkgreen", alpha = 0.7) +
    geom_text(vjust = -1.2, size = 3.5, show.legend = FALSE) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    labs(
        title = "Node count vs. edge count per chromosome",
        subtitle = "Point size = chromosome size (Mb)",
        x = "Number of filtered nodes",
        y = "Number of edges",
        size = "Size (Mb)"
    )

ggsave(file.path(output_dir, "batch8_02_nodes_vs_edges.png"), p2, width = 8, height = 6, dpi = 150)

## ---- Plot 3: mean edge weight and mean degree, per chromosome ----

p3 <- ggplot(merged, aes(x = reorder(chromosome, -size_mb))) +
    geom_col(aes(y = mean_degree), fill = "purple4", alpha = 0.7) +
    labs(
        title = "Mean node degree per chromosome",
        subtitle = "Chromosomes ordered by size (largest first)",
        x = "Chromosome",
        y = "Mean degree"
    )

ggsave(file.path(output_dir, "batch8_03_mean_degree_by_chr.png"), p3, width = 8, height = 6, dpi = 150)

## ---- Plot 4: max degree per chromosome (mega-hub check) ----
## If max_degree approaches n_nodes for that chromosome, that's the
## same mega-hub signature found in chr21 -- worth checking it holds
## (or doesn't) across all 8.

merged$max_degree_fraction <- merged$max_degree / merged$n_nodes_in_graph

p4 <- ggplot(merged, aes(x = reorder(chromosome, -size_mb), y = max_degree_fraction)) +
    geom_col(fill = "firebrick", alpha = 0.7) +
    scale_y_continuous(labels = scales::percent) +
    labs(
        title = "Max node degree as a fraction of total nodes, per chromosome",
        subtitle = "Values near 100% indicate a mega-hub node connected to nearly every other node in that chromosome",
        x = "Chromosome",
        y = "Max degree / n_nodes"
    )

ggsave(file.path(output_dir, "batch8_04_hub_check_by_chr.png"), p4, width = 8, height = 6, dpi = 150)

cat("Saved 4 plots to:", output_dir, "\n")
cat("\nMerged summary table:\n")
print(merged[, c("chromosome", "size_mb", "n_nodes", "n_nodes_in_graph", "n_edges", "mean_k", "mean_weight", "mean_degree", "max_degree_fraction")])
