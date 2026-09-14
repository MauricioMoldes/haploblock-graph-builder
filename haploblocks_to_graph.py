import os
import re
import glob
import csv
import pickle
import statistics
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

# Default: only generate edges between nodes on the SAME chromosome.
# This is what makes it safe to process several chromosomes in one
# run (CHR_FILTER=chr21,chr22,chr20,...) -- cost stays the sum of each
# chromosome's own C(k_chr,2), instead of one C(k_total,2) across all
# of them combined. Set CROSS_CHROMOSOME_EDGES=1 to allow trans edges
# again for a deliberately small, specific set of chromosomes -- not
# recommended for more than 2-3 chromosomes at once given the cost
# demonstrated on the full genome-wide attempt.
CROSS_CHROMOSOME_EDGES = os.getenv("CROSS_CHROMOSOME_EDGES", "0") == "1"


def _chromosome_of(node_id):
    return node_id.split("_", 1)[0]

# If set, skip find_blocks/process_block/merge_nodes entirely and load
# the filtered individual_nodes dict from CHECKPOINT_PATH instead. Use
# this to retry/tune build_edges without re-parsing all cluster.tsv
# files and rewriting nodes.csv every time.
RESUME_FROM_CHECKPOINT = os.getenv("RESUME_FROM_CHECKPOINT", "0") == "1"
# Skip edge GENERATION and jump straight to merging existing
# edges_part_*.tsv files -- for retrying after the merge step itself
# fails, without redoing the expensive generation phase. Independent
# of RESUME_FROM_CHECKPOINT (which skips node-building/parsing instead).
RESUME_EDGE_PARTS = os.getenv("RESUME_EDGE_PARTS", "0") == "1"
CHECKPOINT_PATH = os.path.join(
    os.getenv("TMPDIR", "/tmp"), "haploblock_graph", "individual_nodes_filtered.pkl"
)

TMP = os.getenv("TMPDIR", "/tmp")
TMP = os.path.join(TMP, "haploblock_graph")

os.makedirs(OUTPUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

# chr1_99874387-100330740_cluster.tsv -> ("chr1", "99874387-100330740")
BLOCK_FILENAME_RE = re.compile(r"^(chr[0-9XYM]+)_(\d+-\d+)_cluster\.tsv$")

# Optional: restrict a run to specific chromosomes, e.g. CHR_FILTER=chr21
# or CHR_FILTER=chr21,chr22. Unset (default) processes every chromosome.
# For scoping a fast hackathon-sized subset without touching the main
# genome-wide run.
CHR_FILTER = os.getenv("CHR_FILTER")
CHR_FILTER = set(CHR_FILTER.split(",")) if CHR_FILTER else None

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

        if CHR_FILTER and chr_name not in CHR_FILTER:
            continue

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
# symmetric minor-frequency filter before edge generation
############################################
#
# nodes.csv (already written) keeps the full feature matrix untouched.
# This only affects what goes into build_edges.
#
# Rather than two independently-chosen MIN/MAX thresholds, this uses
# the same convention as minor allele frequency (MAF) filtering in
# population genetics: what matters is how many individuals are on
# the SMALLER side of the split -- present vs. absent -- not which
# side that is. A cluster present in 2400/2548 individuals is exactly
# as statistically weak as one present in 148/2548 (148 people on the
# informative side either way). A separate, independently-set ceiling
# (e.g. "MAX=2000") has no such principled basis and is typically far
# stricter on one side than the other for no good reason.
#
# min_support is therefore applied to min(support, N - support), not
# to support directly. This still correctly excludes the two extremes
# -- true singletons/near-private clusters, and true near-fixed/
# invariant clusters (zero variance, unambiguously uninformative) --
# while keeping real polymorphic signal near either tail that a naive
# asymmetric ceiling would otherwise discard.

def filter_low_support_nodes(individual_nodes, node_to_inds, min_support):

    total_individuals = len(individual_nodes)

    def is_informative(inds):
        support = len(inds)
        minor_side = min(support, total_individuals - support)
        return minor_side >= min_support

    keep = {node for node, inds in node_to_inds.items() if is_informative(inds)}

    total_nodes = len(node_to_inds)
    kept_nodes = len(keep)

    print(
        f"Symmetric minor-frequency filter (min_support={min_support}, "
        f"N={total_individuals}): keeping {kept_nodes}/{total_nodes} nodes "
        f"({total_nodes - kept_nodes} dropped)"
    )

    filtered = {
        ind: (nodes & keep)
        for ind, nodes in individual_nodes.items()
    }

    return filtered


############################################
# diagnostics: report per-individual node count (k) before the
# combinatorial step, since C(k,2) is what actually determines runtime
############################################
#
# IMPORTANT: total_pairs must mirror process_edge_chunk's actual logic
# (same-chromosome-only by default) or this estimate is meaningless --
# a flat C(total_k,2) across an individual's ENTIRE node set massively
# overstates the real workload once more than one chromosome is in
# play, since it's exactly the unrestricted, disproven computation
# this restriction exists to avoid. Measured case: an 8-chromosome
# batch reported total_k=11,602 (flat) vs. real per-chromosome-summed
# work that is roughly 9x smaller -- the flat estimate would have
# reported combinatorics equivalent to having no chromosome
# restriction at all.

def report_edge_workload(individual_nodes):

    ks = [len(nodes) for nodes in individual_nodes.values()]

    if not ks:
        print("No individuals to build edges for.", flush=True)
        return

    print(
        f"Per-individual TOTAL node count (k) after filtering: "
        f"min={min(ks)} median={statistics.median(ks)} "
        f"max={max(ks)} mean={statistics.mean(ks):.1f}",
        flush=True,
    )

    if CROSS_CHROMOSOME_EDGES:
        total_pairs = sum(k * (k - 1) // 2 for k in ks)
        print(
            f"CROSS_CHROMOSOME_EDGES=1: estimated total pairwise edge "
            f"emissions (sum of C(k,2), unrestricted): {total_pairs:,}",
            flush=True,
        )
        return

    # Same-chromosome-only (default): sum C(k_chr, 2) per chromosome,
    # per individual -- matches process_edge_chunk exactly.
    total_pairs = 0
    per_chr_k = defaultdict(list)

    for nodes in individual_nodes.values():
        by_chr = defaultdict(int)
        for node in nodes:
            by_chr[_chromosome_of(node)] += 1
        for chr_name, k_chr in by_chr.items():
            total_pairs += k_chr * (k_chr - 1) // 2
            per_chr_k[chr_name].append(k_chr)

    print(
        f"Estimated total pairwise edge emissions "
        f"(sum of C(k_chr,2), same-chromosome-only): {total_pairs:,}",
        flush=True,
    )
    print("Per-chromosome mean k (nodes/individual, within that chromosome only):", flush=True)
    for chr_name in sorted(per_chr_k):
        vals = per_chr_k[chr_name]
        print(f"  {chr_name}: mean={statistics.mean(vals):.1f} max={max(vals)}", flush=True)


############################################
# build weighted edges (parallel, disk-based)
############################################

def process_edge_chunk(args):

    # `chunk_nodes` is a dict containing ONLY the individuals this worker
    # needs (individual -> set of nodes), not the full population dict.
    chunk_id, chunk_nodes = args

    edge_counts = defaultdict(int)

    for nodes in chunk_nodes.values():

        if CROSS_CHROMOSOME_EDGES:
            groups = [sorted(nodes)]
        else:
            by_chr = defaultdict(list)
            for node in nodes:
                by_chr[_chromosome_of(node)].append(node)
            groups = [sorted(g) for g in by_chr.values()]

        for group in groups:
            for a, b in itertools.combinations(group, 2):
                edge_counts[(a, b)] += 1

    out_file = os.path.join(TMP, f"edges_part_{chunk_id}.tsv")

    with open(out_file, "w") as out:
        for (a, b), count in edge_counts.items():
            out.write(f"{a}\t{b}\t{count}\n")

    return out_file


def build_edges(individual_nodes):

    # If a prior attempt already generated edges_part_*.tsv but died
    # during the (single-threaded, unparallelized) merge step, this
    # skips regenerating them -- the expensive part -- and jumps
    # straight to merging what's already on disk. Set
    # RESUME_EDGE_PARTS=1 for a retry after an OOM/crash specifically
    # in "Merging edge files..." (not after an earlier failure).
    if RESUME_EDGE_PARTS:

        part_files = sorted(glob.glob(os.path.join(TMP, "edges_part_*.tsv")))

        print(
            f"RESUME_EDGE_PARTS=1: skipping generation, found "
            f"{len(part_files)} existing edges_part_*.tsv files in {TMP}",
            flush=True,
        )

        if not part_files:
            print(
                "WARNING: no edges_part_*.tsv files found -- nothing to "
                "resume. Falling back to normal generation.",
                flush=True,
            )
        else:
            return _merge_and_write_edges(part_files)

    print("Building weighted edges in parallel...", flush=True)

    individuals = list(individual_nodes.keys())

    # Chunk size is NOT derived from NPROC. The old formula
    # (len(individuals) // NPROC) meant chunk size grew with individual
    # count independent of per-individual k, so at large scale (big
    # chromosomes, many chromosomes combined) each of NPROC workers
    # held a complete in-memory dict for MANY individuals'
    # combinatorics simultaneously -- peak memory is
    # NPROC * (chunk_size * per-individual-dict-size), and that grows
    # unboundedly with chromosome/genome scope even though NPROC and
    # total memory don't. Measured failure: a 16-chromosome run with
    # mean_k=41,228 and chunk_size=13 individuals/worker was OOM-killed
    # during generation, before any merge step. EDGE_CHUNK_SIZE
    # decouples chunk size from NPROC entirely -- small, fixed chunks
    # bound peak per-worker memory regardless of scope; NPROC still
    # controls how many chunks run concurrently, so total task count
    # (not per-task memory) grows with scale instead.
    chunk_size = int(os.getenv("EDGE_CHUNK_SIZE", "1"))
    chunk_size = max(1, chunk_size)

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

    return _merge_and_write_edges(part_files)


def _merge_and_write_edges(part_files):

    ########################################
    # Merge partial edge files
    ########################################

    print("Merging edge files...", flush=True)

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

    print("Writing final edges...", flush=True)

    with open(edges_path, "w") as out:
        writer = csv.writer(out)
        writer.writerow(["source", "target", "weight"])

        for (a, b), count in global_counts.items():
            writer.writerow([a, b, count])


############################################
# main
############################################

def main():

    if RESUME_FROM_CHECKPOINT:

        print(f"RESUME_FROM_CHECKPOINT=1: loading {CHECKPOINT_PATH}", flush=True)

        with open(CHECKPOINT_PATH, "rb") as f:
            individual_nodes_for_edges = pickle.load(f)

        print(
            f"Loaded checkpoint: {len(individual_nodes_for_edges)} individuals",
            flush=True,
        )

    else:

        print("Scanning haploblocks...", flush=True)

        blocks = find_blocks()

        print("Blocks found:", len(blocks), flush=True)

        ########################################
        # parallel parsing
        ########################################

        print("Processing with", NPROC, "CPUs", flush=True)

        with Pool(NPROC) as pool:

            block_files = pool.map(process_block, blocks)

        ########################################
        # merge nodes
        ########################################

        print("Merging node data", flush=True)

        individual_nodes, node_to_inds = merge_nodes(block_files)

        ########################################
        # filter low-support nodes (see filter_low_support_nodes docstring)
        ########################################

        individual_nodes_for_edges = filter_low_support_nodes(
            individual_nodes, node_to_inds, MIN_CLUSTER_SUPPORT
        )

        ########################################
        # checkpoint, so a failed/retried build_edges never has to
        # redo node parsing + merging + filtering again
        ########################################

        with open(CHECKPOINT_PATH, "wb") as f:
            pickle.dump(individual_nodes_for_edges, f, protocol=pickle.HIGHEST_PROTOCOL)

        print(f"Checkpoint written: {CHECKPOINT_PATH}", flush=True)

    ########################################
    # report workload BEFORE the expensive step, so a walltime kill
    # still leaves us a diagnosis in the log
    ########################################

    report_edge_workload(individual_nodes_for_edges)

    ########################################
    # build edges
    ########################################

    print("Generating edges", flush=True)

    build_edges(individual_nodes_for_edges)

    print("Done.", flush=True)


############################################

if __name__ == "__main__":
    main()
