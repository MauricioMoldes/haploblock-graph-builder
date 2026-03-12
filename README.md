# haploblock-graph-builder
Builds a genome-wide haploblock cluster co-occurrence graph from phased haploblock hashes.
Nodes represent haploblock clusters and edges represent co-occurrence across individuals.

## Biological Context & Design Decisions

### Haploblocks and Individuals
Haploblocks are contiguous regions of the genome where variants are inherited together. In this project, we collapse per-haplotype data into **per-individual representations**, reflecting the combined genetic information from both haplotypes. This approach ensures that downstream graph analyses capture individual-level variation rather than treating haplotypes independently.

### Node Definition
Each node in the graph corresponds to a unique haploblock cluster observed in at least one individual. Nodes are labeled as `<haploblock_region>_cluster<cluster_id>` to maintain traceability back to their genomic region and cluster identity.

### High-Dimensional Feature Vectors
The `high_dim_edge` column represents the haploblock region itself, while the binary vector for each individual indicates whether they carry the given haploblock cluster. This **binary encoding of cluster presence/absence** provides a compact, high-dimensional representation suitable for synthetic node expansion and graph-based analyses.

### Edge Construction
Edges are constructed between nodes that co-occur within the same individual. By iterating over all individuals and generating all pairwise node combinations, the graph encodes **co-occurrence relationships** that reflect real genomic linkage patterns across haploblocks.

### Design Choices
1. **Genome-wide aggregation**: All chromosomes are processed, and blocks are aggregated per individual to build a complete graph.
2. **File-driven parsing**: The script reads `individual_hashes_*.tsv` and `cluster_hashes_*.tsv` files, maintaining flexibility to process blocks in any chromosome order.
3. **Standard Python**: No external dependencies are required, making the script portable across systems with Python 3.x.
4. **Extensibility**: Additional TSV files such as `haploblock_hashes.tsv`, `haploblock_boundaries_chr*.tsv`, and `variant_counts.tsv` can be integrated for future analyses without modifying core logic.
5. **Synthetic Node Expansion Ready**: The output format (`nodes.csv` with binary vectors and `edges.csv`) is fully compatible with downstream pipelines that expand nodes based on graph topology and high-dimensional features.

This design ensures that the graph faithfully represents both the **structural organization of haploblocks** and the **individual-level genetic variation**, forming a robust foundation for subsequent computational analyses.

