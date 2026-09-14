#!/bin/sh
# Stages and compresses nodes.csv/edges.csv (+ toolkit outputs) for all
# 22 autosomes into a distribution-ready directory tree, matching the
# layout described in README.txt.
#
# Usage: ./compress_for_distribution.sh <staging_dir> <batch8_dir> <full16_dir>

set -e

STAGING="$1"
BATCH8_DIR="$2"
FULL16_DIR="$3"

if [ -z "$STAGING" ] || [ -z "$BATCH8_DIR" ] || [ -z "$FULL16_DIR" ]; then
    echo "Usage: $0 <staging_dir> <batch8_dir> <full16_dir>" >&2
    exit 1
fi

ZIP=gzip
command -v pigz >/dev/null 2>&1 && ZIP=pigz
echo "Using compressor: $ZIP"

mkdir -p "$STAGING"

BATCH8_CHRS="chr13 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
FULL16_CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr14 chr15"

compress_chr() {
    chr="$1"
    src_dir="$2"

    out_dir="$STAGING/$chr"
    mkdir -p "$out_dir"

    echo "=== $chr ==="
    date

    for f in nodes.csv edges.csv edges_lift_above_threshold.csv islands.csv position_heatmap.csv top_edges_by_lift.csv; do
        if [ -f "$src_dir/$f" ]; then
            "$ZIP" -c "$src_dir/$f" > "$out_dir/${f}.gz"
            echo "  compressed $f -> ${f}.gz ($(du -h "$out_dir/${f}.gz" | cut -f1))"
        else
            echo "  SKIP: $src_dir/$f not found"
        fi
    done
}

for chr in $BATCH8_CHRS; do
    compress_chr "$chr" "$BATCH8_DIR/$chr"
done

for chr in $FULL16_CHRS; do
    compress_chr "$chr" "$FULL16_DIR/$chr"
done

echo
echo "All done. Staged tree:"
du -sh "$STAGING"/*
