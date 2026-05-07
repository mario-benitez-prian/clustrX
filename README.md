# clustrX: Highly Robust and Sensitive Protein Clustering

[![Version](https://img.shields.io/badge/version-0.1.7-blue.svg)](https://pypi.org/project/clustrX/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**clustrX** is a high-performance framework designed to transform sequence similarity search results into biologically coherent protein families. By modeling homology as a weighted mathematical network and applying the **Leiden community detection algorithm**, `clustrX` provides a sensitive and robust solution for clustering sequences, especially in complex scenarios involving remote homology and short peptides.

---

## 🚀 Key Features

*   **Leiden Community Detection**: Beyond simple links, `clustrX` identifies densely connected communities, ensuring high internal cohesion and preventing artificial family merging (e.g., due to domain bridges).
*   **Agnostic Input**: Works with results from **BLAST**, **Diamond**, **MMseqs2**, and **HMMER**. Or others using the custom input option.
*   **Dynamic Coverage Filter**: Our recommended approach to handle sequences of varying lengths to obtain the most reliable and biologically sound results.
*   **Ultra-Fast Performance**: Powered by `Polars` (Rust-based) for data processing and `igraph` (C-based) for network analysis.
*   **Integrated Workflow**: From similarity hits to Multiple Sequence Alignments (MSAs) in a single command.

---

## 📦 Installation

You can install `clustrX` using two main methods. Note the difference in dependency management:

### Option A: Via Conda (Recommended)
This is the easiest way as it automatically installs all external dependencies, including **MAFFT** for alignments.
```bash
conda install -c bioconda clustrx
```

### Option B: Via Pip
If you prefer `pip`, remember that you must **install MAFFT manually** on your system if you plan to use the `--mafft` option.
```bash
pip install clustrX
```

---

## ⚙️ Input Formats & Requirements

`clustrX` is designed to be a post-processing layer. It requires two main inputs:
1.  **Similarity Hits**: A tabular file (BLAST-like or HMMER).
2.  **Sequences**: A FASTA file containing the sequences referenced in the hits.

### Using BLAST
`clustrX` works natively with the default tabular output of BLAST (`-outfmt 6`).
```bash
blastp -query sequences.fasta -db database -out hits.tsv -outfmt 6
```

### Using Diamond or MMseqs2
If you use these tools, you **must** ensure the output is in **BLAST tabular format (outfmt 6)**:

*   **Diamond**:
    ```bash
    diamond blastp -q query.fasta -d db.dmnd -o hits.tsv --outfmt 6
    ```
*   **MMseqs2**:
    ```bash
    mmseqs easy-search query.fasta target.fasta hits.tsv tmp --format-mode 0
    ```

### Using HMMER
HMMER outputs require specific flags depending on the filtering level you need:

*   **`domtblout` (Recommended)**: Use the `--domtblout` flag in `hmmsearch` or `phmmer`. This format provides alignment coordinates, which are **required** for using the **Dynamic Coverage** filter.
    ```bash
    hmmsearch --domtblout hits.domtblout profile.hmm database.fasta
    ```
*   **`tblout`**: Use the `--tblout` flag. Note that this format lacks coordinate information; therefore, **Dynamic Coverage cannot be applied** (only E-value and Bitscore filters will be used).
    ```bash
    hmmsearch --tblout hits.tblout profile.hmm database.fasta
    ```

---

## 🧬 The Power of Dynamic Coverage

We strongly recommend using the **Dynamic Coverage** mode (`--coverage dynamic`) for most scientific applications. For more information about this, please, read the paper.

Standard clustering methods often use fixed thresholds that fail to resolve relationships between sequences of very different sizes. Our dynamic filter uses a **hyperbolic decay function** (calibrated with a 50-residue scale factor) that:
1.  Increases stringency for **short peptides** (up to 0.8 coverage) to filter out statistical noise.
2.  Gradually relaxes for **larger proteins** (down to 0.4 coverage) to maximize sensitivity in detecting remote homology.

---

## 🛠️ Workflow & Usage

The `clustrX` pipeline follows a clear 3-step logic:

1.  **Filter**: Hits are filtered based on E-value, Bitscore, and (recommended) Dynamic Coverage.
2.  **Cluster**: A similarity network is built where edges are weighted by Bitscore, then partitioned using Leiden algorithm.
3.  **Output**: Results are exported. **Note: Fasta generation and alignments are optional.**

### Example: Recommended Scientific Run
```bash
clustrx -i hits.tsv -f sequences.fasta --coverage dynamic --write-fasta --mafft --outdir results_full
```
*   `--write-fasta`: (Optional) Creates a FASTA file for each generated cluster.
*   `--mafft`: (Optional) Automatically performs Multiple Sequence Alignment for each cluster.

---

## 💡 Use Cases

*   **Protein Family Discovery**: Organizing large proteomes into evolutionarily related groups.
*   **Short Peptide Classification**: Specifically tuned for the discovery of **Antimicrobial Peptides (AMPs)**, toxins, signaling peptides or others.
*   **Remote Homology Exploration**: Identifying relationships in the "twilight zone" (identity < 30%) where traditional greedy methods fragment families.
*   **Domain-Aware Clustering**: Using HMMER `domtblout` inputs to cluster sequences based on specific functional domains.

---

## 📝 Citation
If you use **clustrX** in your research, please cite:
> Benítez-Prián, M. & San Mauro, D. (2026). clustrX: Highly Robust and Sensitive Protein Clustering Using Similarity Networks and Leiden Community Detection.

## 👤 Authors
**Mario Benítez-Prián** & **Diego San Mauro**

Contact: [mario.benitezprian@gmail.com](mailto:mario.benitezprian@gmail.com) | [GitHub](https://github.com/mario-benitez-prian)
