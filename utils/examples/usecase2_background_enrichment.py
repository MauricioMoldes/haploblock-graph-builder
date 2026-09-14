#!/usr/bin/env python3
"""
Use case 2: enrichment of background haploblock clusters among carriers
of a rare penetrant SNP, versus non-carriers.

Carrier list comes from an external source (a VCF/callset for the
causal gene) -- rare penetrant variants are below the resolution this
pipeline clusters at, so they are never in nodes.csv/edges.csv
themselves. Supply carriers as a comma-separated list or a file with
one individual_id per line.

Usage:
    python3 usecase2_background_enrichment.py \
        --nodes nodes.csv \
        --carriers HG10001,HG10002 \
        --exclude-region chr21_1000000-1274441
"""
import argparse
import csv
import os


def load_carriers(spec):
    if os.path.isfile(spec):
        with open(spec) as f:
            return {line.strip() for line in f if line.strip()}
    return {x.strip() for x in spec.split(",") if x.strip()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", required=True, help="Path to nodes.csv")
    parser.add_argument("--carriers", required=True,
                         help="Comma-separated individual_ids, or a path to a file with one id per line")
    parser.add_argument("--exclude-region", default=None,
                         help="high_dim_edge region to exclude (the causal locus's own haploblock)")
    parser.add_argument("--top", type=int, default=20, help="How many top-enriched nodes to print")
    args = parser.parse_args()

    carriers = load_carriers(args.carriers)

    with open(args.nodes) as f:
        reader = csv.reader(f)
        header = next(reader)
        individuals = header[2:]

        rows = []
        for row in reader:
            node, region = row[0], row[1]
            if args.exclude_region and region == args.exclude_region:
                continue
            node_inds = {ind for ind, flag in zip(individuals, row[2:]) if flag == "1"}
            rows.append((node, node_inds))

    non_carriers = set(individuals) - carriers

    results = []
    for node, node_inds in rows:
        carrier_freq = len(node_inds & carriers) / len(carriers) if carriers else 0
        noncarrier_freq = len(node_inds & non_carriers) / len(non_carriers) if non_carriers else 0
        enrichment = carrier_freq - noncarrier_freq
        if carrier_freq > 0:
            results.append((node, carrier_freq, noncarrier_freq, enrichment))

    results.sort(key=lambda r: r[3], reverse=True)

    print(f"Carriers ({len(carriers)}): {sorted(carriers)}")
    if args.exclude_region:
        print(f"Excluding causal locus region: {args.exclude_region}")
    print()
    print(f"{'Node':<45}{'CarrierFreq':<13}{'BgFreq':<10}{'Enrichment':<10}")
    for node, cf, bf, e in results[:args.top]:
        print(f"{node:<45}{cf:<13.2f}{bf:<10.2f}{e:<10.2f}")


if __name__ == "__main__":
    main()
