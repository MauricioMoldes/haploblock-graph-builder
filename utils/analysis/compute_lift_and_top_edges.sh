#!/bin/sh
# Computes lift = weight / (support(A) * support(B) / N) for every edge.
# Writes THREE outputs in one pass over edges.csv:
#   node_support.csv                -- per-node support (individuals carrying it)
#   top_edges_by_lift.csv           -- top-N edges by lift (quick inspection)
#   edges_lift_above_threshold.csv  -- ALL edges with lift >= threshold
#                                       (for island detection / density comparisons,
#                                        more complete than top-N alone)
#
# Usage: ./compute_lift_and_top_edges.sh <dir> [top_n] [lift_threshold] [n_individuals]

set -e

DIR="${1:-.}"
TOP_N="${2:-200}"
LIFT_THRESHOLD="${3:-5}"
N_INDIVIDUALS="${4:-2548}"

cd "$DIR"

AWK=awk
command -v mawk >/dev/null 2>&1 && AWK=mawk

TMP_ALL="./.all_edges_with_lift.tmp"

echo "Computing per-node support from nodes.csv..."
echo "node_id,support" > node_support.csv
"$AWK" -F, 'NR==1{next} {c=0; for(i=3;i<=NF;i++) c+=$i; print $1","c}' nodes.csv >> node_support.csv

echo "Computing lift for every edge (threshold=$LIFT_THRESHOLD, top_n=$TOP_N)..."
echo "source,target,weight,lift" > edges_lift_above_threshold.csv

"$AWK" -F, -v n="$N_INDIVIDUALS" -v thresh="$LIFT_THRESHOLD" -v thresh_out="edges_lift_above_threshold.csv" '
{ gsub(/\r/, "") }
NR==FNR {
    if (FNR == 1) next
    support[$1] = $2
    next
}
FNR==1 { next }
{
    sa = support[$1]; sb = support[$2]
    if (sa == "" || sb == "") next
    expected = (sa * sb) / n
    if (expected <= 0) next
    lift = $3 / expected
    line = $1","$2","$3","lift
    print line
    if (lift >= thresh) print line >> thresh_out
}
' node_support.csv edges.csv > "$TMP_ALL"

echo "source,target,weight,lift" > top_edges_by_lift.csv
LC_ALL=C sort -t, -k4 -rn -T . "$TMP_ALL" | head -n "$TOP_N" >> top_edges_by_lift.csv
rm -f "$TMP_ALL"

echo "Done."
wc -l node_support.csv top_edges_by_lift.csv edges_lift_above_threshold.csv
