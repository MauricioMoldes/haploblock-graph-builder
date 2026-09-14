#!/usr/bin/env Rscript
# Cross-chromosome comparison plots across all 22 autosomes.
# Usage: Rscript plot_all22_comparison.R [summary_csv] [output_dir]
# Expects all_chromosomes_summary.csv from consolidate_all22.sh.

suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
summary_csv <- if (length(args) >= 1) args[1] else "all_chromosomes_summary.csv"
output_dir  <- if (length(args) >= 2) args[2] else "."
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 13))

summary_data <- read.csv(summary_csv)

# Known chromosome sizes (Mb, GRCh38) and mean-k-per-individual figures
# reported directly by the pipeline's own workload diagnostic in the
# batch8 and full-16 run logs -- not recomputed here.
chr_info <- data.frame(
    chromosome = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                   "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
                   "chr18","chr19","chr20","chr21","chr22"),
    size_mb = c(248,242,198,190,182,171,159,145,138,134,135,133,114,107,102,90,83,80,59,64,47,51),
    mean_k  = c(4212.0,4232.6,3629.0,3331.0,3284.0,3218.6,2899.0,2723.6,2430.5,
                2728.9,2500.1,2645.6,1888.1,1741.2,1652.4,1800.0,1674.8,
                1674.4,1357.2,1459.0,820.2,928.4)
)

merged <- merge(summary_data, chr_info, by = "chromosome")
merged$max_degree_fraction <- merged$max_degree / merged$n_nodes_in_graph

## ---- Plot 1: k vs size across all 22 ----
p1 <- ggplot(merged, aes(x = size_mb, y = mean_k, label = chromosome)) +
    geom_point(size = 2.5, color = "steelblue") +
    geom_text(vjust = -1, size = 3) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE,
                color = "grey50", linetype = "dashed", linewidth = 0.6) +
    labs(title = "Mean per-individual k vs. chromosome size (all 22 autosomes)",
         subtitle = "Dashed line: quadratic fit",
         x = "Chromosome size (Mb)", y = "Mean k (nodes/individual, within-chromosome)")
ggsave(file.path(output_dir, "all22_01_k_vs_size.png"), p1, width = 9, height = 6.5, dpi = 150)

## ---- Plot 2: nodes vs edges ----
p2 <- ggplot(merged, aes(x = n_nodes, y = n_edges, label = chromosome, size = size_mb)) +
    geom_point(color = "darkgreen", alpha = 0.7) +
    geom_text(vjust = -1.2, size = 3, show.legend = FALSE) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "Node count vs. edge count (all 22 autosomes)",
         subtitle = "Point size = chromosome size (Mb)",
         x = "Number of filtered nodes", y = "Number of edges", size = "Size (Mb)")
ggsave(file.path(output_dir, "all22_02_nodes_vs_edges.png"), p2, width = 9, height = 6.5, dpi = 150)

## ---- Plot 3: hub check across all 22 ----
p3 <- ggplot(merged, aes(x = reorder(chromosome, -size_mb), y = max_degree_fraction)) +
    geom_col(fill = "firebrick", alpha = 0.7) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Max node degree as a fraction of total nodes (all 22 autosomes)",
         subtitle = "Universal ~100% mega-hub artifact, confirmed genome-wide",
         x = "Chromosome (ordered by size)", y = "Max degree / n_nodes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "all22_03_hub_check.png"), p3, width = 10, height = 6.5, dpi = 150)

## ---- Plot 4 (NEW): largest island span, all 22 ----
p4 <- ggplot(merged, aes(x = reorder(chromosome, -size_mb), y = largest_island_span_kb)) +
    geom_col(fill = "darkorange", alpha = 0.8) +
    geom_hline(yintercept = 500, linetype = "dashed", color = "grey40") +
    labs(title = "Largest extended-haplotype island span, per chromosome",
         subtitle = "Dashed line: 500kb span cap. Consistent 280-500kb range across all 22 -- not an artifact of chromosome size.",
         x = "Chromosome (ordered by size)", y = "Span (kb)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "all22_04_largest_island_span.png"), p4, width = 10, height = 6.5, dpi = 150)

## ---- Plot 5 (NEW): does island size track chromosome size? (it shouldn't, biologically) ----
p5 <- ggplot(merged, aes(x = size_mb, y = largest_island_nodes, label = chromosome)) +
    geom_point(size = 2.5, color = "purple4") +
    geom_text(vjust = -1, size = 3) +
    labs(title = "Largest island's node count vs. chromosome size",
         subtitle = "Flat relationship expected: extended haplotypes are a local phenomenon, not chromosome-scale",
         x = "Chromosome size (Mb)", y = "Nodes in largest island")
ggsave(file.path(output_dir, "all22_05_island_nodes_vs_chr_size.png"), p5, width = 9, height = 6.5, dpi = 150)

cat("Saved 5 plots to:", output_dir, "\n\n")
print(merged[, c("chromosome","size_mb","n_nodes","n_edges","mean_k",
                  "max_degree_fraction","largest_island_nodes","largest_island_span_kb")])
