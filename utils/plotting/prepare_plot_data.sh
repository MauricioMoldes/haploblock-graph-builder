#!/bin/sh
# Converts the raw exploration outputs already sitting in
# graph_output_chr21/ into tidy CSVs for plotting with ggplot2.
#
# Usage: run from inside graph_output_chr21/ (or pass the dir as $1)
#
# Inputs expected (already generated in earlier exploration steps):
#   chr21_support_histogram.txt   -- from: sort -n | uniq -c  (count, support)
#   edge_weight_histogram.txt     -- from: sort -n | uniq -c  (count, weight)
#   node_degree.txt               -- from: sort | uniq -c | sort -rn (degree, node_id)
#
# Outputs (tidy CSVs, ready for read.csv() in R):
#   support_distribution.csv      support,node_count
#   support_threshold_curve.csv   min_threshold,nodes_kept,avg_k   (the full diminishing-returns curve)
#   edge_weight_distribution.csv  weight,edge_count
#   node_degree_distribution.csv  degree

set -e

DIR="${1:-.}"
cd "$DIR"

N_INDIVIDUALS=2548

AWK=awk
command -v mawk >/dev/null 2>&1 && AWK=mawk

echo "Preparing support_distribution.csv..."
echo "support,node_count" > support_distribution.csv
"$AWK" '{print $2","$1}' chr21_support_histogram.txt >> support_distribution.csv

echo "Preparing support_threshold_curve.csv (full cumulative curve)..."
echo "min_threshold,nodes_kept,avg_k" > support_threshold_curve.csv
sort -k2,2 -rn chr21_support_histogram.txt | "$AWK" -v n="$N_INDIVIDUALS" '
{
    cum_nodes += $1
    cum_incid += $1 * $2
    support = $2
    printf "%d,%d,%.4f\n", support, cum_nodes, cum_incid / n
}
' >> support_threshold_curve.csv

echo "Preparing edge_weight_distribution.csv..."
echo "weight,edge_count" > edge_weight_distribution.csv
"$AWK" '{print $2","$1}' edge_weight_histogram.txt >> edge_weight_distribution.csv

echo "Preparing node_degree_distribution.csv..."
echo "degree" > node_degree_distribution.csv
"$AWK" '{print $1}' node_degree.txt >> node_degree_distribution.csv

echo "Done. Files ready for R:"
ls -la support_distribution.csv support_threshold_curve.csv edge_weight_distribution.csv node_degree_distribution.csv
