# haploblock-graph-builder

![Haploblock co-occurrence graph](https://haploblocks.org/figures/haploblock_co_occurence_graph.png)

Builds a **haploblock cluster co-occurrence graph** from MMseqs2 haploblock
clustering results.

Nodes represent **haploblock clusters** and edges represent **co-occurrence
across individuals**.

The pipeline is designed for **population-scale genomic datasets** and
supports efficient execution on HPC clusters.

A validated, published output of this pipeline (22 human autosomes,
1000 Genomes, N=2,548) is available at
[data.haploblocks.org/haplograph/1000G/](https://data.haploblocks.org/haplograph/1000G/).

---

# Overview

This project converts haploblock clustering results into a **graph
representation of genomic background structure**.

The resulting graph encodes how **haploblock clusters co-occur across
individuals**, enabling downstream analyses such as:

* genotype–phenotype association
* genomic background modeling
* network-based feature extraction
* machine learning on genomic hashes

The implementation follows the framework described in:

**Kubica et al. (2025)**
*Decoding Complex Genotype-Phenotype Interactions by Discretizing the Genome*

---

# Biological Context & Design Decisions

## Haploblocks and Individuals

Haploblocks are contiguous regions of the genome where variants are
inherited together.

The input data is the direct output of MMseqs2 clustering on haplotype
sequences within each haploblock. Sample/haplotype identifiers are
collapsed into per-individual representations by extracting the
individual identifier from each sample name.

This ensures the graph represents **individual genomic backgrounds**
rather than individual haplotypes.

---

## Node Definition

Each node represents a **haploblock cluster**.

Nodes are labeled:

```
<chromosome>_<region>_cluster<cluster_id>
```

Example:

```
chr21_14215892-14284114_cluster1
```

`region` is the real coordinate range for that haploblock, taken directly
from the input filename. `cluster_id` is assigned by sorting the unique
MMseqs2 cluster representatives within that region and numbering them
1, 2, 3... — deterministic across re-runs, but not comparable to any
external cluster numbering.

---

## High-Dimensional Feature Vectors

`nodes.csv` contains a binary vector describing cluster presence across
individuals, for **every** node found — this file is never filtered by
population frequency.

Structure:

```
id, high_dim_edge, individual_1, individual_2, ...
```

Example:

```
chr21_14215892-14284114_cluster1,chr21_14215892-14284114,1,0,1,0
```

Where:

```
1 = individual carries this haploblock cluster
0 = cluster absent in that individual
```

This representation forms a **high-dimensional feature matrix** suitable
for graph machine learning, clustering, synthetic node expansion, and
genotype–phenotype modeling.

---

## Edge Construction

Edges represent **co-occurrence of haploblock clusters within the same
individual**, subject to two filters applied before combinatorics run
(neither filter touches `nodes.csv`, only what feeds `edges.csv`):

### 1. Same-chromosome-only, by default

For each individual, clusters are grouped by chromosome, and pairwise
combinations are only generated **within** a chromosome:

```
clusters_on_chr21 = {cluster A, cluster B, cluster C}
clusters_on_chr22 = {cluster D, cluster E}

edges generated: A-B, A-C, B-C, D-E
NOT generated:   A-D, A-E, B-D, ... (cross-chromosome)
```

This is a deliberate, measured trade-off, not an oversight: a genuine
genome-wide attempt with unrestricted (trans-chromosome) combinations
produced a single worker-chunk output file of 738GB and was killed
after exceeding a 72-hour walltime limit twice. See *Computational
Complexity* below. Set `CROSS_CHROMOSOME_EDGES=1` to re-enable
cross-chromosome edges — only recommended for a small, deliberately
scoped set of chromosomes (2-3 at most), never genome-wide.

### 2. Symmetric minor-frequency filter (`MIN_CLUSTER_SUPPORT`)

Before combinations are generated, a cluster is dropped unless the
**smaller** of {individuals carrying it, individuals not carrying it}
is at least `MIN_CLUSTER_SUPPORT`. This mirrors minor allele frequency
(MAF) filtering convention: a cluster present in 2,400 of 2,548
individuals is exactly as statistically weak as one present in 148 —
both have only 148 people on the informative side. This single
threshold replaces what earlier versions handled as two independent,
asymmetric MIN/MAX values.

This filter exists because a cluster's edge weight is capped at the
smaller side's population count regardless — clusters near either
extreme (near-private or near-universal) can only ever produce
low-information edges, while dominating the pairwise combination count.

**Known artifact:** even after filtering, every chromosome tested
independently exhibits a single node connected to ~99.98-99.99% of
all other nodes in its graph. This is a population-frequency artifact,
not real linkage — validated against real 1000 Genomes ancestry
labels, it corresponds to the out-of-Africa serial founder effect
(the "hub" cluster is the ancestral/reference haplotype fixed during
out-of-Africa bottlenecks). **Normalize by lift before interpreting
edge weight** — see `utils/analysis/compute_lift_and_top_edges.sh`.

---

# Pipeline Architecture

```
scan haploblocks (chr*/clusters/*_cluster.tsv)
        │
parallel block parsing (MMseqs2 cluster.tsv -> individual/node pairs)
        │
intermediate membership files ($TMPDIR)
        │
merge per-individual cluster sets
        │
generate nodes.csv (full, unfiltered feature matrix)
        │
symmetric minor-frequency filter (MIN_CLUSTER_SUPPORT)
        │
checkpoint filtered node data (pickle, for cheap edge-step retries)
        │
generate edges.csv (same-chromosome-only, by default)
```

This design avoids storing the entire genome graph in memory during
parsing, and avoids re-parsing all input files just to retry the edge
step with different thresholds or after a crash.

---

# Parallelization Strategy

## Block-Level Parallelism

Each haploblock region is processed independently, via Python
multiprocessing (`process_block`). Because haploblocks are independent,
parsing scales nearly linearly with CPU count.

## Edge-Level Parallelism

Individuals are chunked across workers. **Chunk size is decoupled from
worker count** (`EDGE_CHUNK_SIZE`, default 1 individual/task) — an
earlier design tied chunk size to `NPROC`, which meant each worker held
a complete in-memory dictionary for many individuals' combinatorics
simultaneously; at large scale (16-chromosome run, mean k≈41,000) this
was OOM-killed *during generation*, not the merge step, with roughly
192 workers × a multi-individual chunk each blowing past available
memory concurrently. Setting `EDGE_CHUNK_SIZE=1` bounds peak per-worker
memory regardless of scope; verified to produce byte-identical output
at any chunk size, since it only affects task granularity, not the
underlying combinatorics.

## HPC Execution

Typical execution uses 32–200 CPUs. Temporary files, including the
checkpoint pickle, are written to `$TMPDIR` to avoid network filesystem
bottlenecks:

```
INPUT      → shared filesystem
TEMP FILES → node-local scratch ($TMPDIR)
OUTPUT     → shared filesystem
```

See `utils/pbs/` for real, tested PBS job templates covering a single
scoped chromosome, an 8-chromosome batch, checkpoint/edge-part resume
after a failure, and the full-genome run.

---

# Input Data

Expected directory structure — note the nested `clusters/` subfolder,
and that chromosome is encoded in the filename itself:

```
<HAPLOBLOCK_ROOT>/
├── chr1/
│   └── clusters/
│       ├── chr1_<start>-<end>_cluster.tsv
│       ├── chr1_<start2>-<end2>_cluster.tsv
├── chr2/
│   └── clusters/
│       ├── chr2_<start>-<end>_cluster.tsv
```

## `*_cluster.tsv` files

Standard MMseqs2 `createtsv` cluster output — two columns, tab-separated,
**no header**: representative sequence ID, then cluster member ID (the
representative is repeated once per member, including itself):

```
HG00097_chr1_region_99874387-100330740_hap0    HG00097_chr1_region_99874387-100330740_hap0
HG00119_chr1_region_99874387-100330740_hap1    HG00119_chr1_region_99874387-100330740_hap1
HG00119_chr1_region_99874387-100330740_hap1    NA19774_chr1_region_99874387-100330740_hap0
```

Cluster membership is derived directly from which rows share a
representative — no separate hash-mapping file or hash-suffix matching
is used.

---

# Output Files

## nodes.csv

Full, unfiltered feature matrix.

```
id,high_dim_edge,IND1,IND2,IND3,...
```

| column        | meaning                 |
| ------------- | ----------------------- |
| id            | node identifier         |
| high_dim_edge | haploblock region       |
| IND*          | binary cluster presence |

## edges.csv

Filtered (`MIN_CLUSTER_SUPPORT`), same-chromosome-only (by default) graph
edges, **with weight**:

```
source,target,weight
chr21_14215892-14284114_cluster1,chr21_14284114-14356253_cluster219,7
```

`weight` = number of individuals carrying both `source` and `target`.

---

# Environment Variables

| Variable | Default | Purpose |
| --- | --- | --- |
| `HAPLOBLOCK_ROOT` | `/data` | Root of the `chr*/clusters/*_cluster.tsv` input tree |
| `OUTPUT_DIR` | `/results` | Where `nodes.csv`/`edges.csv` are written |
| `TMPDIR` | `/tmp` | Scratch space for intermediate files, edge parts, and the checkpoint pickle |
| `NPROC` | CPU count | Worker pool size |
| `MIN_CLUSTER_SUPPORT` | `2` | Symmetric minor-frequency threshold (see *Edge Construction*) |
| `EDGE_CHUNK_SIZE` | `1` | Individuals per edge-generation task; decoupled from `NPROC` for memory safety at scale (see *Parallelization Strategy*) |
| `CHR_FILTER` | unset (all) | Comma-separated chromosome list to scope a run, e.g. `chr21` or `chr21,chr22` |
| `CROSS_CHROMOSOME_EDGES` | `0` | Set to `1` to re-enable trans-chromosome edges (not recommended beyond 2-3 chromosomes) |
| `RESUME_FROM_CHECKPOINT` | `0` | Set to `1` to skip node-building/parsing and reload the last filtered checkpoint |
| `RESUME_EDGE_PARTS` | `0` | Set to `1` to skip edge *generation* and jump straight to merging existing `edges_part_*.tsv` — for retrying after a merge-step failure without redoing generation |

---

# Computational Complexity

Let:

```
H = number of haploblocks (39,077 genome-wide, confirmed against real data)
N = number of individuals (2,548, confirmed against real data)
k = surviving nodes per individual, after MIN_CLUSTER_SUPPORT filtering
```

Parsing complexity: `O(H)`.

Edge construction complexity, per chromosome, with the same-chromosome
restriction: `O(N * k_chr²)` where `k_chr` is the per-individual node
count *for that chromosome only* — **not** across the whole genome.

**Genome-wide, unrestricted combinatorics is not tractable at this
scale.** This was tested directly: a single worker's output chunk
reached 738GB before the job was killed at the 72-hour walltime limit
(twice). No amount of node-support filtering alone resolves this — even
keeping only the ~2,239 most population-wide clusters genome-wide still
implied ~6.3 billion pairs, because per-individual `k` genome-wide is
inherently in the tens of thousands. The same-chromosome restriction is
what actually makes the combinatorics tractable, by bounding `k` to a
single chromosome's worth of blocks instead of the whole genome.

`mean k` also does **not** scale as a simple quadratic function of
chromosome size in practice — measured across all 22 autosomes, the
relationship is close to linear, and gene-dense/high-recombination
chromosomes (e.g. chr19) deviate from a pure size-based prediction.

---

# Performance Expectations (measured, not theoretical)

Real production run, all 22 human autosomes, 192 CPUs,
`MIN_CLUSTER_SUPPORT=25`:

| Metric | Result |
| --- | --- |
| Chromosomes | 22 (all autosomes; chrX/chrY not present in upstream clustering data) |
| Individuals | 2,548 |
| Largest single chromosome (chr1, 248Mb) | ~31 minutes for full analysis toolkit pass |
| 8-chromosome batch (chr13,16-22) | 4h53m end-to-end |
| Remaining 16-chromosome batch | ~2 days 22 hours end-to-end (checkpoint-resumed after one OOM on the first attempt) |
| Combined output (compressed) | ~26GB across all 22 chromosomes |

Memory: 2TB was sufficient for the full 16-chromosome run with
`EDGE_CHUNK_SIZE=1`; the same run previously failed at 2TB with the old
NPROC-derived chunking. 800GB was sufficient for the 8-chromosome batch.

**Genome-wide is not a target configuration for this pipeline** — see
*Computational Complexity* above.

---

# Reproducibility

To ensure deterministic output:

* cluster IDs are assigned via sorted order of MMseqs2 representatives
* nodes are sorted
* edges are sorted
* input blocks processed deterministically

Running the pipeline multiple times with the same input and thresholds
yields identical graph files.

---

# Design Principles

### 1. Portability

Uses only the Python standard library. No external dependencies.

This was a deliberate decision point, not a default: a sparse-matrix
(`scipy`) rewrite of edge generation was considered as an alternative
to the same-chromosome restriction, since it would compute co-occurrence
via vectorized linear algebra rather than a Python-level combinatorial
loop. The same-chromosome restriction, the symmetric support filter,
and decoupling `EDGE_CHUNK_SIZE` from `NPROC` proved sufficient without
adding the dependency, so the stdlib-only constraint was kept.

### 2. HPC Compatibility

Supports execution inside Apptainer, Docker, and PBS clusters.

### 3. Memory Efficiency

Intermediate files, the post-filter checkpoint pickle, and small
per-task edge chunks prevent memory overload and avoid re-doing
expensive work when only part of the pipeline needs to be retried.

### 4. Graph ML Compatibility

Output format directly supports graph neural networks, node embedding,
and network analysis.

---

# `utils/` — Auxiliary Tooling

Scripts for analyzing, publishing, and demonstrating the graph, kept
separate from the core pipeline (`haploblocks_to_graph.py`). All are
tested against synthetic data with known-correct expected output before
being run against real data.

### `utils/analysis/`

Post-processing for interpreting the graph correctly:

* `compute_lift_and_top_edges.sh` — computes lift (PMI-style
  normalization) for every edge; raw edge weight is dominated by the
  mega-hub artifact (see *Edge Construction*) and should not be used
  directly.
* `bin_position_heatmap.sh` — bins edges by genomic position for a
  Hi-C-style adjacency heatmap (local LD vs. hub artifact, visually).
* `detect_islands.py` — detects extended-haplotype "islands": dense,
  span-bounded clusters of high-lift edges. Uses a span-aware Union-Find
  **and** a minimum-density requirement — a distance cap alone is not
  sufficient, since a long chain of individually-local edges can
  daisy-chain via transitivity into an arbitrarily large, spurious
  component regardless of the cap chosen.
* `check_individual_outliers.py` — flags individuals who are
  consistent node-count outliers across *multiple* chromosomes at once
  (more likely a technical/QC signal than biology).
* `run_toolkit_all_chr.sh` / `run_toolkit_full16.sh` — wrappers running
  the above across a full chromosome set, with a per-chromosome density
  comparison summary.

### `utils/distribution/`

For splitting a combined multi-chromosome run into per-chromosome files
and packaging for publication:

* `split_by_chromosome.sh` — splits `nodes.csv`/`edges.csv` by
  chromosome (safe because edges are same-chromosome-only by default).
* `split_and_summarize_batch8.sh` — same, plus per-chromosome summary
  stats (node/edge counts, degree distribution) in one pass.
* `consolidate_all22.sh` — joins summary/degree/island stats from
  multiple split runs into one cross-chromosome comparison table.
* `compress_for_distribution.sh` — stages and gzips (using `pigz` if
  available) a full chromosome set plus analysis outputs into a
  publish-ready directory tree.
* `convert_1000g_panel.sh` — converts the IGSR/1000 Genomes sample
  panel file into this project's phenotype schema.

### `utils/plotting/`

R (ggplot2) scripts for the figures referenced above, plus the `awk`
data-prep scripts that feed them from raw CSV output.

### `utils/examples/`

* `usecase1_snp_distribution.py` — given a node (SNP's cluster), reports
  its carrier frequency by ancestry/phenotype group.
* `usecase2_background_enrichment.py` — given a list of rare-variant
  carriers, finds which background clusters are enriched relative to
  non-carriers.
* `generate_synthetic_dataset.py` — generates a synthetic
  `nodes.csv`/`edges.csv`/`phenotypes.csv` with a realistic
  singleton-heavy cluster-size distribution, for testing without real
  data.

### `utils/pbs/`

Real, tested PBS job templates: a single scoped chromosome, an
8-chromosome batch, a full-genome run, and the distribution
splitting/compression jobs. Useful as starting points — adjust paths,
resource requests, and `CHR_FILTER` for your own run.

---

# Future Extensions

Potential extensions include:

* cluster frequency statistics
* haploblock-specific graph layers
* population stratification analysis (worth noting: same-chromosome-only
  edges incidentally reduce, but do not eliminate, population-structure
  confounding — see the mega-hub / ancestry-validation discussion above)
* phenotype association modeling
* graph neural networks on genomic background
* targeted (non-exhaustive) trans-chromosome queries: one-vs-rest
  co-occurrence between a specific locus and the rest of the genome,
  which is linear rather than quadratic and does not require the full
  all-vs-all trans graph
* a reference-cluster-assignment scheme (e.g. `mmseqs search` against a
  fixed representative database, rather than de novo clustering per
  cohort) to make cluster IDs stable and comparable across independently
  processed cohorts (e.g. combining with UK Biobank/All of Us data,
  where raw genotypes typically cannot leave the source enclave and only
  aggregate summary statistics can be pooled)

---

# License

Released under the same license as the associated BioHackathon pipeline.
