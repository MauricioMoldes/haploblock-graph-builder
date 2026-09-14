#!/bin/sh
# Splits the batch8 nodes.csv/edges.csv by chromosome (reusing the same
# logic as split_by_chromosome.sh) AND computes per-chromosome summary
# stats in the same pass, writing one consolidated CSV for cross-
# chromosome comparison plots.
#
# Usage:
#   ./split_and_summarize_batch8.sh <outdir> <nodes.csv> <edges.csv>
#
# Output: <outdir>/<chr>/nodes.csv, <outdir>/<chr>/edges.csv (as before)
#         <outdir>/chromosome_summary.csv  (one row per chromosome)

set -e

OUTDIR="$1"
NODES="$2"
EDGES="$3"

if [ -z "$OUTDIR" ] || [ -z "$NODES" ] || [ -z "$EDGES" ]; then
    echo "Usage: $0 <outdir> <nodes.csv> <edges.csv>" >&2
    exit 1
fi

AWK=awk
command -v mawk >/dev/null 2>&1 && AWK=mawk

mkdir -p "$OUTDIR"

echo "Splitting $NODES by chromosome (node counts)..."
"$AWK" -F, -v outdir="$OUTDIR" '
NR == 1 { header = $0; next }
{
    chr = substr($1, 1, index($1, "_") - 1)
    if (!(chr in seen)) {
        system("mkdir -p \"" outdir "/" chr "\"")
        outfile[chr] = outdir "/" chr "/nodes.csv"
        print header > outfile[chr]
        seen[chr] = 1
    }
    print $0 > outfile[chr]
    node_count[chr]++
}
END {
    for (chr in node_count) print chr "\t" node_count[chr] > "/tmp/.batch8_node_counts.tsv"
}
' "$NODES"

echo "Splitting $EDGES by chromosome, with weight/degree stats per chromosome..."
"$AWK" -F, -v outdir="$OUTDIR" '
NR == 1 { header = $0; next }
{
    gsub(/\r/, "")
    src_chr = substr($1, 1, index($1, "_") - 1)
    tgt_chr = substr($2, 1, index($2, "_") - 1)
    if (src_chr != tgt_chr) {
        print "WARNING: cross-chromosome edge (" src_chr " <-> " tgt_chr ") skipped in stats" > "/dev/stderr"
        next
    }
    chr = src_chr
    if (!(chr in seen)) {
        system("mkdir -p \"" outdir "/" chr "\"")
        outfile[chr] = outdir "/" chr "/edges.csv"
        print header > outfile[chr]
        seen[chr] = 1
    }
    print $0 > outfile[chr]

    w = $3 + 0
    edge_count[chr]++
    weight_sum[chr] += w
    if (!(chr in weight_max) || w > weight_max[chr]) weight_max[chr] = w

    degree[chr","$1]++
    degree[chr","$2]++
}
END {
    for (k in degree) {
        split(k, parts, ",")
        c = parts[1]
        deg_sum[c] += degree[k]
        deg_count[c]++
        if (degree[k] > deg_max[c]) deg_max[c] = degree[k]
    }
    print "chromosome,n_edges,mean_weight,max_weight,mean_degree,max_degree" > outdir "/chromosome_summary_partial.csv"
    for (chr in edge_count) {
        mean_deg = (deg_count[chr] > 0) ? deg_sum[chr] / deg_count[chr] : 0
        printf "%s,%d,%.2f,%d,%.2f,%d\n", chr, edge_count[chr], weight_sum[chr]/edge_count[chr], weight_max[chr], mean_deg, deg_max[chr] >> outdir "/chromosome_summary_partial.csv"
    }
}
' "$EDGES"

echo "Merging node counts into final summary..."
echo "chromosome,n_nodes,n_edges,mean_weight,max_weight,mean_degree,max_degree" > "$OUTDIR/chromosome_summary.csv"
"$AWK" -F'\t' '
NR==FNR { nodes[$1] = $2; next }
FNR==1 { next }
{
    split($0, f, ",")
    chr = f[1]
    print chr "," nodes[chr] "," f[2] "," f[3] "," f[4] "," f[5] "," f[6]
}
' /tmp/.batch8_node_counts.tsv "$OUTDIR/chromosome_summary_partial.csv" >> "$OUTDIR/chromosome_summary.csv"

rm -f /tmp/.batch8_node_counts.tsv "$OUTDIR/chromosome_summary_partial.csv"

echo "Done. Summary:"
cat "$OUTDIR/chromosome_summary.csv"
