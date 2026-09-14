#!/usr/bin/env Rscript
# Two additional diagnostic plots:
#   - position heatmap: does edge strength concentrate near the
#     genomic diagonal (real LD) or scatter (stratification/hubs)?
#   - node-link graph of the top-N edges BY LIFT (not raw weight) --
#     lift = weight / expected-under-independence, so this surfaces
#     genuine enrichment rather than mega-hub nodes that are just
#     individually common.
#
# Usage: Rscript plot_graph_extras.R [input_dir] [output_dir] [label]
# Expects position_heatmap.csv and top_edges_by_lift.csv in input_dir
# (from bin_position_heatmap.sh and compute_lift_and_top_edges.sh).

suppressMessages({
    library(ggplot2)
    library(igraph)
    library(ggraph)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "."
output_dir <- if (length(args) >= 2) args[2] else "."
label      <- if (length(args) >= 3) args[3] else "chr21"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 13))

## ---- Plot 5: position-binned adjacency heatmap ----

heatmap_data <- read.csv(file.path(input_dir, "position_heatmap.csv"))

# mirror across the diagonal for a symmetric, readable heatmap
mirrored <- heatmap_data
mirrored$tmp <- mirrored$bin_i
mirrored$bin_i <- mirrored$bin_j
mirrored$bin_j <- mirrored$tmp
mirrored$tmp <- NULL
heatmap_full <- rbind(heatmap_data, mirrored[mirrored$bin_i != mirrored$bin_j, ])

p5 <- ggplot(heatmap_full, aes(x = bin_i / 1e6, y = bin_j / 1e6, fill = total_weight)) +
    geom_tile() +
    scale_fill_viridis_c(trans = "log10", labels = scales::comma) +
    coord_fixed() +
    labs(
        title = paste(label, "position-binned edge weight heatmap"),
        subtitle = "Concentration near the diagonal = local LD. Off-diagonal blocks = stratification/hub signature.",
        x = "Position (Mb)",
        y = "Position (Mb)",
        fill = "Total\nweight"
    )

ggsave(file.path(output_dir, "05_position_heatmap.png"), p5, width = 8, height = 7, dpi = 150)

## ---- Plot 6: node-link graph of top edges by lift ----

top_edges <- read.csv(file.path(input_dir, "top_edges_by_lift.csv"))

g <- graph_from_data_frame(top_edges[, c("source", "target", "weight", "lift")], directed = FALSE)

p6 <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(edge_width = lift, edge_alpha = lift), color = "steelblue") +
    geom_node_point(size = 2, color = "grey20") +
    scale_edge_width(range = c(0.3, 2.5)) +
    scale_edge_alpha(range = c(0.3, 0.9)) +
    theme_void() +
    labs(
        title = paste(label, "- top", nrow(top_edges), "edges by lift"),
        subtitle = "Lift = observed / expected-under-independence co-occurrence -- surfaces real enrichment, not just common nodes"
    ) +
    theme(legend.position = "bottom")

ggsave(file.path(output_dir, "06_top_lift_network.png"), p6, width = 10, height = 9, dpi = 150)

cat("Saved 2 plots to:", output_dir, "\n")
