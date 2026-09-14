#!/bin/sh
# Bins every edge's source/target by genomic position midpoint, and
# aggregates total weight per bin-pair -- for a Hi-C-style adjacency
# heatmap. Real LD should show strength concentrated near the
# diagonal (nearby positions); off-diagonal concentration is the
# stratification/hub signature discussed alongside this analysis.
#
# Usage: ./bin_position_heatmap.sh <dir> [bin_size_bp]
# Expects edges.csv in <dir>. Writes <dir>/position_heatmap.csv

set -e

DIR="${1:-.}"
BIN_SIZE="${2:-500000}"

cd "$DIR"

AWK=awk
command -v mawk >/dev/null 2>&1 && AWK=mawk

echo "Binning edges by genomic position (bin size: $BIN_SIZE bp)..."
echo "bin_i,bin_j,total_weight,edge_count" > position_heatmap.csv

"$AWK" -F, -v binsize="$BIN_SIZE" '
function midpoint_bin(node_id,    n, parts, coords, mid) {
    n = split(node_id, parts, "_")
    n = split(parts[2], coords, "-")
    mid = (coords[1] + coords[2]) / 2
    return int(mid / binsize) * binsize
}
NR==1 { next }
{
    bi = midpoint_bin($1)
    bj = midpoint_bin($2)
    # canonicalize so (A,B) and (B,A) land in the same cell
    if (bi > bj) { tmp = bi; bi = bj; bj = tmp }
    key = bi "," bj
    weight_sum[key] += $3
    count[key] += 1
}
END {
    for (key in weight_sum) {
        print key "," weight_sum[key] "," count[key]
    }
}
' edges.csv >> position_heatmap.csv

echo "Done."
wc -l position_heatmap.csv
