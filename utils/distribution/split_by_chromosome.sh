#!/bin/sh
# Splits a combined nodes.csv and/or edges.csv into per-chromosome
# subdirectories, e.g.:
#   outdir/chr21/nodes.csv
#   outdir/chr21/edges.csv
#   outdir/chr22/nodes.csv
#   outdir/chr22/edges.csv
#
# Relies on chromosome being derivable from column 1 (id / source),
# which holds as long as edges were generated with the default
# same-chromosome-only setting (CROSS_CHROMOSOME_EDGES unset/0).
# If that setting was overridden, this script will warn on any edge
# row whose source and target chromosomes disagree, since such a row
# can only be correctly filed under ONE chromosome folder here.
#
# Usage:
#   ./split_by_chromosome.sh <outdir> [nodes.csv] [edges.csv]
#
# Either input file may be omitted (pass "" to skip it).

set -e

OUTDIR="$1"
NODES="$2"
EDGES="$3"

if [ -z "$OUTDIR" ]; then
    echo "Usage: $0 <outdir> [nodes.csv] [edges.csv]" >&2
    exit 1
fi

AWK=awk
command -v mawk >/dev/null 2>&1 && AWK=mawk

mkdir -p "$OUTDIR"

if [ -n "$NODES" ]; then
    echo "Splitting $NODES by chromosome..."
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
        count[chr]++
    }
    END {
        for (chr in count) print "  " chr ": " count[chr] " nodes"
    }
    ' "$NODES"
fi

if [ -n "$EDGES" ]; then
    echo "Splitting $EDGES by chromosome..."
    "$AWK" -F, -v outdir="$OUTDIR" '
    NR == 1 { header = $0; next }
    {
        src_chr = substr($1, 1, index($1, "_") - 1)
        tgt_chr = substr($2, 1, index($2, "_") - 1)
        if (src_chr != tgt_chr) {
            print "WARNING: cross-chromosome edge (" src_chr " <-> " tgt_chr \
                  ") filed under " src_chr " -- CROSS_CHROMOSOME_EDGES was likely enabled" \
                  " for this run" > "/dev/stderr"
        }
        chr = src_chr
        if (!(chr in seen)) {
            system("mkdir -p \"" outdir "/" chr "\"")
            outfile[chr] = outdir "/" chr "/edges.csv"
            print header > outfile[chr]
            seen[chr] = 1
        }
        print $0 > outfile[chr]
        count[chr]++
    }
    END {
        for (chr in count) print "  " chr ": " count[chr] " edges"
    }
    ' "$EDGES"
fi

echo "Done. Output under: $OUTDIR/<chromosome>/"
