# clustrX: High-Performance Graph-Based Sequence Clustering

[![Version](https://img.shields.io/badge/version-0.1.5-blue.svg)](https://test.pypi.org/project/clustrX/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**clustrX** is a high-performance Python tool designed for ultra-fast clustering of protein and nucleotide sequences using the **Leiden community detection algorithm**. 

By modeling sequence homology as a weighted mathematical network, `clustrX` overcomes the limitations of traditional "greedy incremental" methods (like CD-HIT) and "single-linkage" network tools (like SiLiX), providing superior precision and biological coherence in complex evolutionary scenarios.

---

## 🔬 Scientific Core: Leiden vs. Single-Linkage

Traditional graph-based clustering often relies on *Single-Linkage*, where a single similarity hit can bridge two unrelated groups. This leads to **fragmentation** or the creation of **giant, chimeric clusters** in the presence of:
*   **Domain Bridges**: Unrelated protein families sharing a single common domain.
*   **Chimeric Sequences**: Proteins with multiple, unrelated functional modules.
*   **Highly Divergent Families**: Low-identity homologs that are missed by greedy methods.

**clustrX** addresses this by using the **Leiden Algorithm**, which maximizes network *modularity*. It identifies densely connected communities rather than just simple links, allowing it to:
1.  **Resolve Complex Architectures**: Correctly separate families even when they share a few domain-bridge hits.
2.  **Ensure Taxonomic Purity**: Minimize the creation of taxonomically inconsistent "mixed" clusters.
3.  **High Sensitivity**: Group divergent sequences (like AMPs, Toxins, and small peptides) that standard methods often break apart.

---

## 🚀 Technical Advantages

*   **Ultra-High Speed (Polars Engine)**: All hit parsing and filtering is powered by `Polars`, a Rust-based DataFrame library. This allows `clustrX` to process millions of hits in seconds with minimal RAM usage.
*   **Scalability (C-based igraph)**: Network operations are executed via `python-igraph`, enabling the analysis of massive homology networks where `NetworkX` would fail.
*   **Scientific Filters**:
    *   **Dynamic Coverage**: Automatically adjusts the alignment length requirement based on sequence size ($C_{threshold} = \text{clamp} \left( 1.0 - \frac{\min(L_{max}, 500)}{600} , 0.4, 0.8 \right)$).
    *   **Bitscore, E-value, and Identity**: Multi-parameter edge filtering to refine the network before clustering.
*   **Agnostic Input Parsing**: Supports BLAST tabular (`-outfmt 6`), HMMER `domtblout`, or any other tool via custom column mapping.

---

## 🛠️ Installation

Currently available on **TestPyPI**. You can install it using `pip`:

```bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ clustrX
```

### Dependencies
*   Python ≥ 3.8
*   `polars`, `python-igraph`, `psutil`, `numpy`

---

## 📖 Quick Start

### Basic Clustering (BLAST)
```bash
clustrx -i hits.blast -f sequences.fasta --format blast --outdir results
```

### Advanced Scientific Clustering (HMMER + Dynamic Coverage)
```bash
clustrx -i hits.domtblout -f sequences.fasta --format domhmmer --coverage dynamic --outdir results_hmmer
```

### Output Parameters
`clustrX` produces a clean output structure ready for downstream analysis:
*   **`clusters/`**: TXT files containing the IDs of sequences in each family.
*   **`fasta_files/`**: Automatically extracted FASTA files for each cluster.
*   **`alignments/`**: (Optional) Multiple Sequence Alignments (MSAs) orchestrated via **MAFFT**.

---

## 🧬 Validation Scenarios
The scientific robustness of `clustrX` has been validated against Gold Standard datasets (Pfam Seed), extreme complexity cases (AMPs/Toxins), and ground-truth numerical simulations, consistently outperforming standard methods in Precisión, Recall, and Jaccard Index.

---

## 📝 Citation
If you use **clustrX** in your research, please cite:
> Benítez-Prián, M. (2024). clustrX: High-performance graph clustering for sequence similarity networks.

## 👤 Author
**Mario Benítez-Prián**  
[mario.benitezprian@gmail.com](mailto:mario.benitezprian@gmail.com) | [GitHub](https://github.com/mario-benitez-prian)
