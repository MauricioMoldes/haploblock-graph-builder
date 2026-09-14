#!/bin/sh
# Runs compute_lift_and_top_edges.sh, detect_islands.py, and
# bin_position_heatmap.sh across the 16 remaining chromosomes
# (chr1-12, chr14, chr15, chrX, chrY), then aggregates a per-chromosome
# density comparison table -- same logic as run_toolkit_all_chr.sh,
# generalized with the full GRCh38 size table.
#
# Usage: ./run_toolkit_full16.sh <base_dir> [lift_threshold] [min_island_nodes]
# Expects <base_dir>/<chr>/{nodes.csv,edges.csv} for each chromosome
# (from split_and_summarize_batch8.sh run against graph_output_full).

set -e

BASE_DIR="${1:-.}"
LIFT_THRESHOLD="${2:-5}"
MIN_ISLAND_NODES="${3:-3}"

CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr14 chr15 chrX chrY"
SCRIPT_DIR="$(dirname "$0")"

echo "chromosome,size_mb,n_enriched_edges,edges_per_mb,n_islands,largest_island_nodes,largest_island_kb" > "$BASE_DIR/density_comparison_full16.csv"

for chr in $CHROMOSOMES; do
    chr_dir="$BASE_DIR/$chr"
    if [ ! -f "$chr_dir/nodes.csv" ] || [ ! -f "$chr_dir/edges.csv" ]; then
        echo "SKIP $chr: nodes.csv/edges.csv not found in $chr_dir" >&2
        continue
    fi

    echo "=== $chr ==="
    date

    "$SCRIPT_DIR/compute_lift_and_top_edges.sh" "$chr_dir" 200 "$LIFT_THRESHOLD" 2548

    "$SCRIPT_DIR/bin_position_heatmap.sh" "$chr_dir" 500000

    python3 "$SCRIPT_DIR/detect_islands.py" "$chr_dir/edges_lift_above_threshold.csv" "$MIN_ISLAND_NODES" "" 500000 0.3 > "$chr_dir/islands.csv"

    # Full GRCh38 chromosome sizes (Mb)
    case "$chr" in
        chr1)  size_mb=248 ;;
        chr2)  size_mb=242 ;;
        chr3)  size_mb=198 ;;
        chr4)  size_mb=190 ;;
        chr5)  size_mb=182 ;;
        chr6)  size_mb=171 ;;
        chr7)  size_mb=159 ;;
        chr8)  size_mb=145 ;;
        chr9)  size_mb=138 ;;
        chr10) size_mb=134 ;;
        chr11) size_mb=135 ;;
        chr12) size_mb=133 ;;
        chr14) size_mb=107 ;;
        chr15) size_mb=102 ;;
        chrX)  size_mb=156 ;;
        chrY)  size_mb=57 ;;
        *) size_mb=1 ;;
    esac

    n_enriched=$(($(wc -l < "$chr_dir/edges_lift_above_threshold.csv") - 1))
    n_islands=$(($(wc -l < "$chr_dir/islands.csv") - 1))

    if [ "$n_islands" -gt 0 ]; then
        largest=$(awk -F, 'NR==2{print $3","$8}' "$chr_dir/islands.csv")
        largest_nodes=$(echo "$largest" | cut -d, -f1)
        largest_kb=$(echo "$largest" | cut -d, -f2)
    else
        largest_nodes=0
        largest_kb=0
    fi

    edges_per_mb=$(awk -v n="$n_enriched" -v mb="$size_mb" 'BEGIN{printf "%.2f", n/mb}')

    echo "$chr,$size_mb,$n_enriched,$edges_per_mb,$n_islands,$largest_nodes,$largest_kb" >> "$BASE_DIR/density_comparison_full16.csv"
done

echo
echo "=== Density comparison across the 16 remaining chromosomes ==="
cat "$BASE_DIR/density_comparison_full16.csv"
