#!/bin/sh
# Runs compute_lift_and_top_edges.sh, detect_islands.py, and
# bin_position_heatmap.sh across all 8 chromosome subdirectories, then
# aggregates a per-chromosome "real signal density" comparison table
# -- count of genuinely enriched (lift-filtered) edges normalized by
# chromosome size, comparable across all 8.
#
# Usage: ./run_toolkit_all_chr.sh <base_dir> [lift_threshold] [min_island_nodes]
# Expects <base_dir>/<chr>/{nodes.csv,edges.csv} for each of the 8 chromosomes
# (from split_and_summarize_batch8.sh).

set -e

BASE_DIR="${1:-.}"
LIFT_THRESHOLD="${2:-5}"
MIN_ISLAND_NODES="${3:-3}"

CHROMOSOMES="chr13 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
# GRCh38 sizes in Mb, for density normalization
SCRIPT_DIR="$(dirname "$0")"

echo "chromosome,size_mb,n_enriched_edges,edges_per_mb,n_islands,largest_island_nodes,largest_island_kb" > "$BASE_DIR/density_comparison.csv"

for chr in $CHROMOSOMES; do
    chr_dir="$BASE_DIR/$chr"
    if [ ! -f "$chr_dir/nodes.csv" ] || [ ! -f "$chr_dir/edges.csv" ]; then
        echo "SKIP $chr: nodes.csv/edges.csv not found in $chr_dir" >&2
        continue
    fi

    echo "=== $chr ==="

    "$SCRIPT_DIR/compute_lift_and_top_edges.sh" "$chr_dir" 200 "$LIFT_THRESHOLD" 2548

    "$SCRIPT_DIR/bin_position_heatmap.sh" "$chr_dir" 500000

    python3 "$SCRIPT_DIR/detect_islands.py" "$chr_dir/edges_lift_above_threshold.csv" "$MIN_ISLAND_NODES" > "$chr_dir/islands.csv"

    case "$chr" in
        chr13) size_mb=114 ;;
        chr16) size_mb=90 ;;
        chr17) size_mb=83 ;;
        chr18) size_mb=80 ;;
        chr19) size_mb=59 ;;
        chr20) size_mb=64 ;;
        chr21) size_mb=47 ;;
        chr22) size_mb=51 ;;
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

    echo "$chr,$size_mb,$n_enriched,$edges_per_mb,$n_islands,$largest_nodes,$largest_kb" >> "$BASE_DIR/density_comparison.csv"
done

echo
echo "=== Density comparison across all 8 chromosomes ==="
cat "$BASE_DIR/density_comparison.csv"
