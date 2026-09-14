#!/bin/sh
# Consolidates chromosome_summary.csv + graph_node_counts.csv +
# largest-island data from both graph_output_batch8/ (8 chromosomes)
# and graph_output_full/by_chromosome/ (14 chromosomes) into one
# all_chromosomes_summary.csv covering all 22 autosomes.
#
# Usage: ./consolidate_all22.sh <batch8_dir> <full16_dir> <output_csv>

set -e

BATCH8_DIR="$1"
FULL16_DIR="$2"
OUT="${3:-all_chromosomes_summary.csv}"

echo "chromosome,n_nodes,n_edges,mean_weight,max_weight,mean_degree,max_degree,n_nodes_in_graph,largest_island_nodes,largest_island_span_kb" > "$OUT"

join_row() {
    chr="$1"
    summary_file="$2"
    node_counts_file="$3"
    island_file="$4"

    summary_row=$(awk -F, -v c="$chr" '$1==c{print}' "$summary_file")
    node_count=$(awk -F, -v c="$chr" '$1==c{print $2}' "$node_counts_file")
    island_row=$(awk -F, -v c="$chr" '$1==c{print $2","$3}' "$island_file")

    echo "${summary_row},${node_count},${island_row}"
}

for chr in chr13 chr16 chr17 chr18 chr19 chr20 chr21 chr22; do
    join_row "$chr" "$BATCH8_DIR/chromosome_summary.csv" "$BATCH8_DIR/graph_node_counts.csv" "$BATCH8_DIR/largest_island_corrected.csv" >> "$OUT"
done

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr14 chr15; do
    join_row "$chr" "$FULL16_DIR/chromosome_summary.csv" "$FULL16_DIR/graph_node_counts_full16.csv" "$FULL16_DIR/largest_island_corrected.csv" >> "$OUT"
done

echo "Done."
cat "$OUT"
