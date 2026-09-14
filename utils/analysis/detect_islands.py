#!/usr/bin/env python3
"""
Detects "extended haplotype islands" -- clusters of interlocking
high-lift edges that stay within a bounded total genomic span.

A single high-lift pair is just a strong dyad; a connected CLUSTER of
several interlocking high-lift edges (like chr21's ~39.0-39.9Mb
region, found manually earlier) suggests a single extended haplotype
that haploblock boundaries didn't fully resolve into one block.

Usage:
    python3 detect_islands.py <edges_lift_above_threshold.csv> [min_nodes] [chromosome_label] [max_span_bp] [min_density]

min_density (default 0.3): an island must have at least this fraction
of all possible pairwise edges (n_edges / C(n_nodes,2)) to be reported.
This is the REAL discriminator, not the span cap. Testing showed the
span cap alone doesn't work: the largest "island" found scaled almost
exactly with whatever cap was set (e.g. ~500kb at a 500kb cap, ~2Mb at
a 2Mb cap) -- a daisy-chain of individually-local hops just gets
truncated at whatever ceiling is imposed, rather than actually
stopping at a natural boundary. A genuine extended haplotype should
have most of its member blocks correlated with EACH OTHER (dense),
not just their immediate chain-neighbors (sparse): the trusted chr21
example (4 nodes, 3 edges) has 50% density; the chain artifacts found
during testing had ~3-4% density despite passing any span cap.
"""
import csv
import sys
from collections import defaultdict


class SpanAwareUnionFind:
    """Union-Find where each root tracks the [min_pos, max_pos] span of
    its component. union() refuses to merge if the resulting combined
    span would exceed max_span_bp."""

    def __init__(self, max_span_bp):
        self.parent = {}
        self.span = {}  # root -> (min_pos, max_pos)
        self.max_span_bp = max_span_bp

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def add_node(self, x, pos):
        root = self.find(x)
        if root not in self.span:
            self.span[root] = (pos, pos)

    def try_union(self, a, b, pos_a, pos_b):
        """Attempts to union a and b. Returns True if merged, False if
        refused because the combined span would exceed max_span_bp."""
        self.add_node(a, pos_a)
        self.add_node(b, pos_b)
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return True  # already in the same component

        min_a, max_a = self.span[ra]
        min_b, max_b = self.span[rb]
        combined_min = min(min_a, min_b)
        combined_max = max(max_a, max_b)

        if combined_max - combined_min > self.max_span_bp:
            return False

        self.parent[ra] = rb
        self.span[rb] = (combined_min, combined_max)
        del self.span[ra]
        return True


def parse_position(node_id):
    parts = node_id.split("_")
    chrom = parts[0]
    region = parts[1]
    start, end = region.split("-")
    midpoint = (int(start) + int(end)) / 2
    return chrom, midpoint


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    edges_path = sys.argv[1]
    min_nodes = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    chrom_label = sys.argv[3] if len(sys.argv) > 3 else None
    max_span_bp = int(sys.argv[4]) if len(sys.argv) > 4 else 500_000
    min_density = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3

    uf = SpanAwareUnionFind(max_span_bp)
    edge_list = []
    node_positions = {}

    with open(edges_path) as f:
        rows = list(csv.DictReader(f))

    for row in rows:
        a, b = row["source"], row["target"]
        for node in (a, b):
            if node not in node_positions:
                node_positions[node] = parse_position(node)

    # process edges sorted by genomic position, so islands grow outward
    # from their densest local core rather than in input-file order
    def edge_sort_key(row):
        _, pa = node_positions[row["source"]]
        _, pb = node_positions[row["target"]]
        return min(pa, pb)

    rows.sort(key=edge_sort_key)

    for row in rows:
        a, b = row["source"], row["target"]
        lift = float(row["lift"])
        chrom_a, pos_a = node_positions[a]
        chrom_b, pos_b = node_positions[b]
        if chrom_a == chrom_b:
            uf.try_union(a, b, pos_a, pos_b)
        edge_list.append((a, b, lift))

    components = defaultdict(list)
    for node in node_positions:
        components[uf.find(node)].append(node)

    component_edges = defaultdict(list)
    for a, b, lift in edge_list:
        root_a, root_b = uf.find(a), uf.find(b)
        if root_a == root_b:
            component_edges[root_a].append((a, b, lift))

    islands = []
    for root, nodes in components.items():
        if len(nodes) < min_nodes:
            continue
        edges_here = component_edges[root]
        if not edges_here:
            continue

        max_possible_edges = len(nodes) * (len(nodes) - 1) / 2
        density = len(edges_here) / max_possible_edges if max_possible_edges > 0 else 0
        if density < min_density:
            continue

        chroms = {node_positions[n][0] for n in nodes}
        chrom_str = "MIXED:" + ",".join(sorted(chroms)) if len(chroms) > 1 else next(iter(chroms))
        positions = [node_positions[n][1] for n in nodes]
        lifts = [e[2] for e in edges_here]
        islands.append({
            "chromosome": chrom_str,
            "n_nodes": len(nodes),
            "n_edges": len(edges_here),
            "density": round(density, 3),
            "span_start": int(min(positions)),
            "span_end": int(max(positions)),
            "span_kb": round((max(positions) - min(positions)) / 1000, 1),
            "max_lift": round(max(lifts), 2),
            "mean_lift": round(sum(lifts) / len(lifts), 2),
        })

    islands.sort(key=lambda x: x["n_nodes"], reverse=True)

    writer = csv.DictWriter(
        sys.stdout,
        fieldnames=["island_id", "chromosome", "n_nodes", "n_edges", "density",
                    "span_start", "span_end", "span_kb", "max_lift", "mean_lift"],
    )
    writer.writeheader()
    for i, island in enumerate(islands, 1):
        island["island_id"] = i
        writer.writerow(island)


if __name__ == "__main__":
    main()
