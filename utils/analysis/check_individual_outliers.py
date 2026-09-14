#!/usr/bin/env python3
"""
Checks whether any individuals are consistent outliers (unusually high
or low node counts) across MULTIPLE chromosomes at once -- as opposed
to being an outlier on just one. A consistent cross-chromosome outlier
is more likely a technical/QC signal (sequencing depth, phasing
quality for that individual) than genuine biology, since real biology
wouldn't typically make someone an outlier on every single chromosome
simultaneously.

Uses the checkpoint pickle directly (individual -> filtered node set),
so no need to re-derive anything from nodes.csv/edges.csv.

Usage:
    python3 check_individual_outliers.py <individual_nodes_filtered.pkl> [z_threshold]
"""
import pickle
import statistics
import sys
from collections import defaultdict


def chromosome_of(node_id):
    return node_id.split("_", 1)[0]


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    pkl_path = sys.argv[1]
    z_threshold = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0

    with open(pkl_path, "rb") as f:
        individual_nodes = pickle.load(f)

    # per_chr_k[chr][individual] = k for that individual on that chromosome
    per_chr_k = defaultdict(dict)

    for ind, nodes in individual_nodes.items():
        by_chr = defaultdict(int)
        for node in nodes:
            by_chr[chromosome_of(node)] += 1
        for chrom, k in by_chr.items():
            per_chr_k[chrom][ind] = k

    chromosomes = sorted(per_chr_k.keys())

    # z-score of each individual's k, within each chromosome
    z_scores = defaultdict(dict)  # z_scores[individual][chrom] = z
    for chrom, ind_k in per_chr_k.items():
        values = list(ind_k.values())
        mean_k = statistics.mean(values)
        stdev_k = statistics.stdev(values) if len(values) > 1 else 0
        if stdev_k == 0:
            continue
        for ind, k in ind_k.items():
            z_scores[ind][chrom] = (k - mean_k) / stdev_k

    # average z-score per individual, across chromosomes they appear in
    avg_z = {}
    consistency = {}  # fraction of chromosomes where |z| >= threshold
    for ind, chrom_z in z_scores.items():
        zs = list(chrom_z.values())
        avg_z[ind] = statistics.mean(zs)
        consistency[ind] = sum(1 for z in zs if abs(z) >= z_threshold) / len(zs)

    print(f"Individuals: {len(avg_z)}")
    print(f"Chromosomes: {chromosomes}")
    print(f"z_threshold: {z_threshold}\n")

    # flag individuals who are outliers (|z| >= threshold) on MOST chromosomes
    flagged = [
        (ind, avg_z[ind], consistency[ind])
        for ind in avg_z
        if consistency[ind] >= 0.5
    ]
    flagged.sort(key=lambda x: abs(x[1]), reverse=True)

    print(f"Individuals flagged as outliers on >=50% of chromosomes: {len(flagged)}\n")
    print(f"{'individual':<15}{'avg_z':<10}{'frac_chr_outlier':<18}")
    for ind, az, frac in flagged[:30]:
        print(f"{ind:<15}{az:<10.2f}{frac:<18.2f}")


if __name__ == "__main__":
    main()
