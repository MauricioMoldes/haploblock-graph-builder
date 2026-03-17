#!/usr/bin/env python3

import os
import glob
import csv
import itertools
from multiprocessing import Pool, cpu_count
from collections import defaultdict

ROOT = os.getenv("HAPLOBLOCK_ROOT", "/data")
OUTPUT = os.getenv("OUTPUT_DIR", "/results")
NPROC = int(os.getenv("NPROC") or cpu_count())

TMP = os.getenv("TMPDIR", "/tmp")
TMP = os.path.join(TMP, "haploblock_graph")

os.makedirs(OUTPUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

############################################
# sample → individual
############################################

def extract_individual(sample):
    return sample.split("_")[0]


############################################
# find haploblocks
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
# load cluster mapping
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
# process block (parallel)
############################################

def process_block(args):

    chr_name, region, indiv_file, cluster_file = args

    cluster_map = load_clusters(cluster_file)

    out_file = os.path.join(TMP, f"{chr_name}_{region}.tsv")

    with open(indiv_file) as f, open(out_file, "w") as out:

        writer = csv.writer(out, delimiter="\t")

        for line in f:

            sample, hash_string = line.strip().split()

            individual = extract_individual(sample)

            cluster_hash = hash_string[-20:]

            cluster_id = cluster_map.get(cluster_hash)

            if cluster_id:

                node = f"{chr_name}_{region}_cluster{cluster_id}"

                writer.writerow([individual, node])

    return out_file


############################################
# merge nodes
############################################

def merge_nodes(block_files):

    individual_nodes = defaultdict(set)
    node_to_inds = defaultdict(set)

    for f in block_files:

        with open(f) as fh:

            for ind, node in csv.reader(fh, delimiter="\t"):

                individual_nodes[ind].add(node)
                node_to_inds[node].add(ind)

    individuals = sorted(individual_nodes.keys())
    nodes = sorted(node_to_inds.keys())

    nodes_path = os.path.join(OUTPUT, "nodes.csv")

    with open(nodes_path, "w") as out:

        writer = csv.writer(out)

        writer.writerow(["id","high_dim_edge"] + individuals)

        for node in nodes:

            block = node.split("_cluster")[0]

            inds = node_to_inds[node]

            row = [node, block]

            row.extend(1 if i in inds else 0 for i in individuals)

            writer.writerow(row)

    return individual_nodes

############################################
# build weighted edges (parallel, disk-based)
############################################

def process_edge_chunk(args):

    chunk_id, individuals, individual_nodes = args

    edge_counts = defaultdict(int)

    for ind in individuals:
        nodes = sorted(individual_nodes[ind])

        for a, b in itertools.combinations(nodes, 2):
            edge_counts[(a, b)] += 1

    out_file = os.path.join(TMP, f"edges_part_{chunk_id}.tsv")

    with open(out_file, "w") as out:
        for (a, b), count in edge_counts.items():
            out.write(f"{a}\t{b}\t{count}\n")

    return out_file


def build_edges(individual_nodes):

    print("Building weighted edges in parallel...")

    individuals = list(individual_nodes.keys())

    # Split individuals into chunks
    chunk_size = max(1, len(individuals) // NPROC)

    chunks = [
        individuals[i:i + chunk_size]
        for i in range(0, len(individuals), chunk_size)
    ]

    args = [
        (i, chunk, individual_nodes)
        for i, chunk in enumerate(chunks)
    ]

    ########################################
    # Parallel edge generation
    ########################################

    with Pool(NPROC) as pool:
        part_files = pool.map(process_edge_chunk, args)

    ########################################
    # Merge partial edge files
    ########################################

    print("Merging edge files...")

    global_counts = defaultdict(int)

    for pf in part_files:
        with open(pf) as f:
            for line in f:
                a, b, count = line.strip().split("\t")
                global_counts[(a, b)] += int(count)

    ########################################
    # Write final edges.csv
    ########################################

    edges_path = os.path.join(OUTPUT, "edges.csv")

    print("Writing final edges...")

    with open(edges_path, "w") as out:
        writer = csv.writer(out)
        writer.writerow(["source", "target", "weight"])

        for (a, b), count in global_counts.items():
            writer.writerow([a, b, count])


############################################
# main
############################################

def main():

    print("Scanning haploblocks...")

    blocks = find_blocks()

    print("Blocks found:", len(blocks))

    ########################################
    # parallel parsing
    ########################################

    print("Processing with", NPROC, "CPUs")

    with Pool(NPROC) as pool:

        block_files = pool.map(process_block, blocks)

    ########################################
    # merge nodes
    ########################################

    print("Merging node data")

    individual_nodes = merge_nodes(block_files)

    ########################################
    # build edges
    ########################################

    print("Generating edges")

    build_edges(individual_nodes)

    print("Done.")


############################################

if __name__ == "__main__":
    main()
