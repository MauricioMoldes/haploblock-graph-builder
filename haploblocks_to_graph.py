#!/usr/bin/env python3

import os
import glob
import csv
import itertools
from collections import defaultdict

ROOT = os.getenv("HAPLOBLOCK_ROOT", "/data")
OUTPUT = os.getenv("OUTPUT_DIR", "/results")

os.makedirs(OUTPUT, exist_ok=True)

nodes_path = os.path.join(OUTPUT, "nodes.csv")
edges_path = os.path.join(OUTPUT, "edges.csv")



############################################
# collapse haplotypes → individual
############################################

def extract_individual(sample):

    return sample.split("_")[0]


############################################
# find all haploblock files
############################################

def find_blocks():

    blocks = []

    for chr_dir in sorted(glob.glob(f"{ROOT}/chr*")):

        chr_name = os.path.basename(chr_dir)

        indiv_files = glob.glob(f"{chr_dir}/individual_hashes_*.tsv")

        for f in indiv_files:

            region = f.split("individual_hashes_")[1].replace(".tsv","")

            cluster_file = f"{chr_dir}/cluster_hashes_{region}.tsv"

            if os.path.exists(cluster_file):

                blocks.append((chr_name, region, f, cluster_file))

    return blocks


############################################
# load cluster map
############################################

def load_clusters(cluster_file):

    mapping = {}

    with open(cluster_file) as f:

        for line in f:

            parts = line.strip().split()

            if len(parts) != 2:
                continue

            cluster_id, hash_val = parts

            mapping[hash_val] = cluster_id

    return mapping


############################################
# parse individuals
############################################

def parse_block(indiv_file, cluster_map, chr_name, region):

    individual_clusters = defaultdict(set)

    with open(indiv_file) as f:

        for line in f:

            sample, hash_string = line.strip().split()

            individual = extract_individual(sample)

            cluster_hash = hash_string[-20:]

            cluster_id = cluster_map.get(cluster_hash)

            if cluster_id:

                node = f"{chr_name}_{region}_cluster{cluster_id}"

                individual_clusters[individual].add(node)

    return individual_clusters



############################################
# main
############################################

def main():

    print("Scanning haploblocks...")

    blocks = find_blocks()

    print("blocks found:", len(blocks))

    individual_clusters = defaultdict(set)

    ############################################
    # aggregate across genome
    ############################################

    for chr_name, region, indiv_file, cluster_file in blocks:

        print("processing", chr_name, region)

        cluster_map = load_clusters(cluster_file)

        block_data = parse_block(indiv_file, cluster_map, chr_name, region)

        for ind, nodes in block_data.items():

            individual_clusters[ind].update(nodes)

    ############################################
    # build node list
    ############################################

    print("building node table")

    individuals = sorted(individual_clusters.keys())

    all_nodes = set()

    for nodes in individual_clusters.values():

        all_nodes.update(nodes)

    all_nodes = sorted(all_nodes)

    ############################################
    # write nodes.csv
    ############################################

    with open("nodes_path.csv","w") as f:

        writer = csv.writer(f)

        header = ["id","high_dim_edge"] + individuals

        writer.writerow(header)

        for node in all_nodes:

            block = node.split("_cluster")[0]

            vec = []

            for ind in individuals:

                vec.append(1 if node in individual_clusters[ind] else 0)

            writer.writerow([node, block] + vec)

    ############################################
    # build edges
    ############################################

    print("building edges")

    edges = set()

    for ind, nodes in individual_clusters.items():

        for a, b in itertools.combinations(sorted(nodes), 2):

            edges.add((a, b))

    ############################################
    # write edges.csv
    ############################################

    with open("edges_path.csv","w") as f:

        writer = csv.writer(f)

        writer.writerow(["source","target"])

        for e in edges:

            writer.writerow(e)

    print("done")

############################################
# entry point
############################################

if __name__ == "__main__":
    main()
