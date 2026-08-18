import os
import re
import glob
import csv
import itertools
from multiprocessing import Pool, cpu_count
from collections import defaultdict

ROOT = os.getenv("HAPLOBLOCK_ROOT", "/data")
OUTPUT = os.getenv("OUTPUT_DIR", "/results")
NPROC = int(os.getenv("NPROC") or cpu_count())

# Minimum number of individuals that must carry a cluster for it to be
# eligible for edge generation. Clusters below this are almost always
# MMseqs2 singleton clusters (private to one haplotype) — they can
# never produce an edge with weight > 1, contribute no population-level
# co-occurrence signal, and are the dominant term in the combinatorial
# blowup (most of an individual's nodes are private/rare clusters).
# Does NOT affect nodes.csv, which still reports the full feature matrix.
MIN_CLUSTER_SUPPORT = int(os.getenv("MIN_CLUSTER_SUPPORT", "2"))

TMP = os.getenv("TMPDIR", "/tmp")
TMP = os.path.join(TMP, "haploblock_graph")

os.makedirs(OUTPUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

# chr1_99874387-100330740_cluster.tsv -> ("chr1", "99874387-100330740")
BLOCK_FILENAME_RE = re.compile(r"^(chr[0-9XYM]+)_(\d+-\d+)_cluster\.tsv$")

############################################
# sample_hap → individual
############################################

def extract_individual(sample):
    return sample.split("_")[0]


############################################
# find haploblocks
############################################

def find_blocks():

    blocks = []

    for f in sorted(glob.glob(f"{ROOT}/chr*/clusters/*_cluster.tsv")):

        fname = os.path.basename(f)

        m = BLOCK_FILENAME_RE.match(fname)

        if not m:
            continue

        chr_name, region = m.group(1), m.group(2)

        blocks.append((chr_name, region, f))

    return blocks


############################################
# process block (parallel)
############################################
#
# Input is an MMseqs2 "createtsv" style cluster file:
#   representative_member    cluster_member
# with the representative repeated once per member of its cluster.
# We assign a stable numeric cluster id per representative (sorted,
# so re-runs are deterministic) and emit individual -> node rows.

def process_block(args):

    chr_name, region, cluster_file = args

    out_file = os.path.join(TMP, f"{chr_name}_{region}.tsv")

    representatives = set()
    rows = []

    with open(cluster_file) as f:

        for line in f:

            parts = line.split()

            if len(parts) != 2:
                continue

            representative, member = parts

            representatives.add(representative)
            rows.append((representative, member))

    cluster_id_map = {
        rep: cluster_id
        for cluster_id, rep in enumerate(sorted(representatives), start=1)
    }

    with open(out_file, "w") as out:

        writer = csv.writer(out, delimiter="\t")

        for representative, member in rows:

            individual = extract_individual(member)

            cluster_id = cluster_id_map[representative]

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

    return individual_nodes, node_to_inds


############################################
# filter low-support nodes before edge generation
############################################
#
# nodes.csv (already written) keeps the full feature matrix untouched.
# This only affects what goes into build_edges: nodes carried by fewer
# than MIN_CLUSTER_SUPPORT individuals are dropped from each
# individual's node set, since they can never form a co-occurrence
# edge with weight > 1 and dominate the pairwise combination count.

def filter_low_support_nodes(individual_nodes, node_to_inds, min_support):

    keep = {node for node, inds in node_to_inds.items() if len(inds) >= min_support}

    total_nodes = len(node_to_inds)
    kept_nodes = len(keep)

    print(
        f"Node support filter (min_support={min_support}): "
        f"keeping {kept_nodes}/{total_nodes} nodes "
        f"({total_nodes - kept_nodes} dropped)"
    )

    filtered = {
        ind: (nodes & keep)
        for ind, nodes in individual_nodes.items()
    }

    return filtered


############################################
# build weighted edges (parallel, disk-based)
############################################

def process_edge_chunk(args):

    # `chunk_nodes` is a dict containing ONLY the individuals this worker
    # needs (individual -> set of nodes), not the full population dict.
    chunk_id, chunk_nodes = args

    edge_counts = defaultdict(int)

    for nodes in chunk_nodes.values():

        for a, b in itertools.combinations(sorted(nodes), 2):
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
        (i, {ind: individual_nodes[ind] for ind in chunk})
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

    individual_nodes, node_to_inds = merge_nodes(block_files)

    ########################################
    # filter low-support nodes (see filter_low_support_nodes docstring)
    ########################################

    individual_nodes_for_edges = filter_low_support_nodes(
        individual_nodes, node_to_inds, MIN_CLUSTER_SUPPORT
    )

    ########################################
    # build edges
    ########################################

    print("Generating edges")

    build_edges(individual_nodes_for_edges)

    print("Done.")


############################################

if __name__ == "__main__":
    main()
