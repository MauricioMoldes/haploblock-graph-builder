# haploblock-graph-builder

Builds a **genome-wide haploblock cluster co-occurrence graph** from phased haploblock hashes.

Nodes represent **haploblock clusters** and edges represent **co-occurrence across individuals**.

The pipeline is designed for **population-scale genomic datasets** and supports efficient execution on HPC clusters.

---

# Overview

This project converts haploblock clustering results into a **graph representation of genomic background structure**.

The resulting graph encodes how **haploblock clusters co-occur across individuals**, enabling downstream analyses such as:

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

Haploblocks are contiguous regions of the genome where variants are inherited together.

The input data contains **haplotype-level hashes**. These are collapsed into **per-individual representations** by extracting the individual identifier from each sample name.

This ensures the graph represents **individual genomic backgrounds** rather than individual haplotypes.

---

## Node Definition

Each node represents a **haploblock cluster**.

Nodes are labeled:

```
<chromosome>_<region>_cluster<cluster_id>
```

Example:

```
chr6_153_cluster7
```

This labeling ensures traceability to:

* chromosome
* haploblock region
* cluster identity

Nodes therefore correspond to **distinct genomic background states within haploblocks**.

---

## High-Dimensional Feature Vectors

The node table contains a binary vector describing cluster presence across individuals.

Structure:

```
node_id, haploblock_region, individual_1, individual_2, ...
```

Example:

```
chr6_153_cluster7,chr6_153,1,0,1,0
```

Where:

```
1 = individual carries this haploblock cluster
0 = cluster absent in that individual
```

This representation forms a **high-dimensional feature matrix** suitable for:

* graph machine learning
* clustering
* synthetic node expansion
* genotype–phenotype modeling

---

## Edge Construction

Edges represent **co-occurrence of haploblock clusters within the same individual**.

For each individual:

```
clusters = {all haploblock clusters observed}

edges = all pairwise combinations(clusters, 2)
```

Example:

Individual carries:

```
chr1_20_cluster3
chr4_91_cluster5
chr9_12_cluster1
```

Edges generated:

```
chr1_20_cluster3 — chr4_91_cluster5
chr1_20_cluster3 — chr9_12_cluster1
chr4_91_cluster5 — chr9_12_cluster1
```

This produces a **genome-wide co-occurrence network** capturing relationships between genomic background states.

---

# Pipeline Architecture

The pipeline follows a **map–reduce design** optimized for HPC environments.

```
scan haploblocks
        │
parallel block parsing
        │
intermediate membership files
        │
merge per-individual cluster sets
        │
generate nodes.csv
        │
generate edges.csv
```

This design avoids storing the entire genome graph in memory during parsing.

---

# Parallelization Strategy

## Block-Level Parallelism

Each haploblock region is processed independently.

Parallelization is performed using Python multiprocessing:

```
process_block(chr, region)
```

Each worker:

1. Loads cluster mappings
2. Parses haploblock hashes
3. Emits individual–cluster assignments

Because haploblocks are independent, this scales nearly linearly with CPU count.

---

## HPC Execution

Typical execution uses:

```
32–200 CPUs
```

Temporary files are written to:

```
$TMPDIR
```

to avoid network filesystem bottlenecks.

Recommended HPC pattern:

```
INPUT      → shared filesystem
TEMP FILES → node-local scratch ($TMPDIR)
OUTPUT     → shared filesystem
```

---

# Input Data

Expected directory structure:

```
haploblocks/
├── chr1/
│   ├── individual_hashes_<region>.tsv
│   ├── cluster_hashes_<region>.tsv
│
├── chr2/
│   ├── individual_hashes_<region>.tsv
│   ├── cluster_hashes_<region>.tsv
```

---

## individual_hashes files

Format:

```
sample_id    haploblock_hash
```

Example:

```
HG00123_A    010110101010...
HG00123_B    110101010111...
```

The last **20 bits** encode the cluster identifier.

---

## cluster_hashes files

Format:

```
cluster_id    hash_suffix
```

Example:

```
1    10100101101001010101
2    11100010101010100111
```

Used to map haplotype hashes to cluster identifiers.

---

# Output Files

## nodes.csv

Contains the node feature matrix.

```
id,high_dim_edge,IND1,IND2,IND3,...
```

Columns:

| column        | meaning                 |
| ------------- | ----------------------- |
| id            | node identifier         |
| high_dim_edge | haploblock region       |
| IND*          | binary cluster presence |

---

## edges.csv

Contains graph edges.

```
source,target
chr1_15_cluster3,chr4_91_cluster2
chr6_10_cluster5,chr12_77_cluster1
```

Edges represent **cluster co-occurrence within individuals**.

---

# Computational Complexity

Let:

```
H = number of haploblocks
N = number of individuals
C = clusters per individual
```

Parsing complexity:

```
O(H)
```

Edge construction complexity:

```
O(N * C²)
```

For typical datasets:

```
~2000 haploblocks
~2500 individuals
```

The pipeline remains tractable on modern HPC systems.

---

# Performance Expectations

Typical runtime:

| CPUs | runtime        |
| ---- | -------------- |
| 8    | ~1–2 hours     |
| 32   | ~20–40 minutes |
| 200  | ~5–10 minutes  |

Memory usage:

```
10–50 GB depending on dataset
```

---

# Reproducibility

To ensure deterministic output:

* nodes are sorted
* edges are sorted
* input blocks processed deterministically

Running the pipeline multiple times yields identical graph files.

---

# Design Principles

The implementation follows several principles:

### 1. Portability

Uses only the Python standard library.

No external dependencies.

---

### 2. HPC Compatibility

Supports execution inside:

* Apptainer
* Docker
* PBS / SLURM clusters

---

### 3. Memory Efficiency

Intermediate files prevent memory overload for large datasets.

---

### 4. Graph ML Compatibility

Output format directly supports:

* graph neural networks
* node embedding
* network analysis

---

# Future Extensions

Potential extensions include:

* cluster frequency statistics
* haploblock-specific graph layers
* population stratification analysis
* phenotype association modeling
* graph neural networks on genomic background

---

# License

Released under the same license as the associated BioHackathon pipeline.

---


