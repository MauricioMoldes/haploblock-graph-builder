#!/usr/bin/env python3
"""
Use case 1: distribution of a phenotype-associated SNP's haploblock
cluster across ancestry groups.

The SNP -> cluster mapping (TARGET_NODE) is a step you supply: find
which haploblock region your SNP of interest falls in (region ranges
are visible in nodes.csv's high_dim_edge column), then determine which
cluster within that region carries the effect allele (requires the
phased sequence data / VCF -- not produced by this pipeline directly).

Usage:
    python3 usecase1_snp_distribution.py \
        --nodes nodes.csv --phenotypes phenotypes.csv \
        --node chr21_1000000-1274441_cluster7 \
        --label "height-assoc SNP rs_TOY (effect allele)"
"""
import argparse
import csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", required=True, help="Path to nodes.csv")
    parser.add_argument("--phenotypes", required=True, help="Path to phenotypes.csv")
    parser.add_argument("--node", required=True, help="Target node id (SNP's cluster)")
    parser.add_argument("--phenotype-field", default="ancestry",
                         help="Which phenotype column value to group by (default: ancestry)")
    parser.add_argument("--label", default="", help="Free-text label for the SNP/trait")
    args = parser.parse_args()

    # --- who carries the target cluster? ---
    with open(args.nodes) as f:
        reader = csv.reader(f)
        header = next(reader)
        individuals = header[2:]

        carriers = set()
        found = False
        for row in reader:
            if row[0] == args.node:
                carriers = {ind for ind, flag in zip(individuals, row[2:]) if flag == "1"}
                found = True
                break

    if not found:
        print(f"Node '{args.node}' not found in {args.nodes}")
        return

    # --- join against phenotype table ---
    groups = {}
    with open(args.phenotypes) as f:
        for row in csv.DictReader(f):
            if row["phenotype"] == args.phenotype_field:
                groups[row["individual_id"]] = row["value"]

    totals = {}
    carrier_counts = {}
    for ind, grp in groups.items():
        totals[grp] = totals.get(grp, 0) + 1
        if ind in carriers:
            carrier_counts[grp] = carrier_counts.get(grp, 0) + 1

    label = f"  ({args.label})" if args.label else ""
    print(f"Node: {args.node}{label}")
    print(f"Total carriers: {len(carriers)}\n")
    print(f"{'Group':<12}{'Carriers':<10}{'Total':<8}{'Frequency':<10}")
    for grp in sorted(totals):
        c = carrier_counts.get(grp, 0)
        t = totals[grp]
        print(f"{grp:<12}{c:<10}{t:<8}{c/t:<10.2f}")


if __name__ == "__main__":
    main()
