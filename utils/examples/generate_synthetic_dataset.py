#!/usr/bin/env python3
"""
Generates a synthetic haplograph dataset (nodes.csv, edges.csv,
phenotypes.csv) for hackathon use, in case the real pipeline output
isn't ready in time or HPC access is unavailable during the event.

The cluster-size distribution is deliberately singleton-heavy (most
clusters carried by just 1-2 individuals, a few carried by many) to
match what real MMseqs2 clustering on population WGS data looks like
-- so code written against this dataset behaves the same way against
the real one.

Usage:
    python3 generate_synthetic_dataset.py --outdir ./synthetic_data \
        --individuals 200 --chromosomes chr21,chr22 --blocks-per-chr 50
"""
import argparse
import csv
import itertools
import os
import random
from collections import defaultdict

ANCESTRIES = ["EUR", "AFR", "EAS", "SAS", "AMR"]


def make_individuals(n):
    # HGxxxxx / NAxxxxx style IDs, matching real 1000G naming
    ids = []
    for i in range(n):
        prefix = "HG" if i % 3 else "NA"
        ids.append(f"{prefix}{10000 + i:05d}")
    return ids


def make_ancestry(individuals, rng):
    return {ind: rng.choice(ANCESTRIES) for ind in individuals}


def make_clusters_for_block(individuals, rng):
    """
    Assigns each individual (x2 haplotypes) to a cluster, with a
    singleton-heavy distribution: most clusters have 1-2 members, a
    few "common" clusters absorb a large fraction of the population --
    same shape as real MMseqs2 output on population data.
    """
    haplotypes = [f"{ind}_hap{h}" for ind in individuals for h in (0, 1)]
    rng.shuffle(haplotypes)

    clusters = defaultdict(list)
    cluster_id = itertools.count(1)

    i = 0
    n = len(haplotypes)
    while i < n:
        # 70% chance: small/singleton cluster (1-2 members)
        # 30% chance: common cluster (absorbs a larger random chunk)
        if rng.random() < 0.7:
            size = min(rng.choice([1, 1, 1, 2]), n - i)
        else:
            size = min(rng.randint(5, max(6, n // 15)), n - i)

        cid = next(cluster_id)
        for hap in haplotypes[i:i + size]:
            clusters[cid].append(hap.split("_hap")[0])  # back to individual id
        i += size

    return clusters  # cluster_id -> [individual, individual, ...] (may repeat if both haps in same cluster)


def build_dataset(individuals, chromosomes, blocks_per_chr, rng):
    node_to_inds = {}  # node_id -> set(individuals)
    individual_nodes = defaultdict(set)

    for chr_name in chromosomes:
        pos = 1_000_000
        for b in range(blocks_per_chr):
            start = pos
            end = pos + rng.randint(200_000, 800_000)
            pos = end + 1
            region = f"{chr_name}_{start}-{end}"

            clusters = make_clusters_for_block(individuals, rng)

            for cid, members in clusters.items():
                node_id = f"{region}_cluster{cid}"
                inds = set(members)
                node_to_inds[node_id] = inds
                for ind in inds:
                    individual_nodes[ind].add(node_id)

    return node_to_inds, individual_nodes


def write_nodes_csv(path, individuals, node_to_inds):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["id", "high_dim_edge"] + individuals)
        for node_id in sorted(node_to_inds):
            region = node_id.split("_cluster")[0]
            inds = node_to_inds[node_id]
            row = [node_id, region]
            row.extend(1 if ind in inds else 0 for ind in individuals)
            writer.writerow(row)


def write_edges_csv(path, individual_nodes, min_support_node_to_inds, min_edge_weight=1):
    edge_counts = defaultdict(int)
    for ind, nodes in individual_nodes.items():
        for a, b in itertools.combinations(sorted(nodes), 2):
            edge_counts[(a, b)] += 1

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["source", "target", "weight"])
        for (a, b), w in sorted(edge_counts.items()):
            if w >= min_edge_weight:
                writer.writerow([a, b, w])


def write_phenotypes_csv(path, ancestry):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["individual_id", "phenotype", "value", "source"])
        for ind, pop in ancestry.items():
            writer.writerow([ind, "ancestry", pop, "1000G_panel"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="./synthetic_data")
    parser.add_argument("--individuals", type=int, default=200)
    parser.add_argument("--chromosomes", default="chr21,chr22")
    parser.add_argument("--blocks-per-chr", type=int, default=50)
    parser.add_argument("--min-edge-weight", type=int, default=2)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    individuals = make_individuals(args.individuals)
    ancestry = make_ancestry(individuals, rng)
    chromosomes = args.chromosomes.split(",")

    node_to_inds, individual_nodes = build_dataset(
        individuals, chromosomes, args.blocks_per_chr, rng
    )

    nodes_path = os.path.join(args.outdir, "nodes.csv")
    edges_path = os.path.join(args.outdir, "edges.csv")
    phenotypes_path = os.path.join(args.outdir, "phenotypes.csv")

    write_nodes_csv(nodes_path, individuals, node_to_inds)
    write_edges_csv(edges_path, individual_nodes, node_to_inds, args.min_edge_weight)
    write_phenotypes_csv(phenotypes_path, ancestry)

    print(f"Individuals: {len(individuals)}")
    print(f"Chromosomes: {chromosomes}")
    print(f"Nodes: {len(node_to_inds)}")
    print(f"Wrote: {nodes_path}")
    print(f"Wrote: {edges_path}")
    print(f"Wrote: {phenotypes_path}")


if __name__ == "__main__":
    main()
