# clustRX

A lightweight Python tool to cluster sequences from BLAST/HMMER outputs and generate FASTA files for each cluster. Designed to quickly organize sequences into connected groups based on similarity.

## Description

`clustRX` takes two text files containing query and hit IDs from a homology search (BLAST, HMMER, etc.) and a FASTA file with all corresponding sequences. It uses graph-based clustering (NetworkX) to identify connected components, then outputs both text files with cluster IDs and FASTA files for each cluster. This allows for easy downstream analysis such as alignments, gene family studies, or metagenomic binning. The tool is flexible and works with outputs from any homology search software.

## Getting Started

### Dependencies

* Python ≥3.8
* `networkx` Python package
* Works on Linux, MacOS, Windows with Python installed

### Installing

* Install via PyPI (once published):
```bash
pip install clustRX
```

### Executing Program

* Run the program using the CLI, providing a BLAST tabular file or HMMER tblout as input:

```bash
clustrx -i hits.txt --format blast -f all_sequences.fasta --min-cluster-size 2 --outdir my_output
```

* Output example:

``` bash
my_output/
├── clusters/
│ ├── cluster_1.txt
│ └── cluster_2.txt
└── fasta_files/
├── cluster_1.fasta
└── cluster_2.fasta
```

* Step-by-step:

1. Prepare `hits.txt` (standard BLAST tabular `-outfmt 6` or HMMER `-tblout`).
2. Prepare `all_sequences.fasta` containing all sequences referenced in the hits.
3. Run the `clustrx` command with required arguments.
4. Check `my_output/clusters` for ID lists and `my_output/fasta_files` for FASTA sequences.

## Help

* For usage instructions:

```bash
clustrx -h
```

## Author information

Mario Benítez-Prián  
(https://github.com/mario-benitez-prian)  
mario.benitezprian@gmail.com  

## Version History

* 0.1.0
    * Initial release: basic clustering functionality, CLI, FASTA output

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Potential Applications

1. **Gene family clustering across multiple species** – group homologous genes for downstream analysis.  
2. **Protein isoform or variant grouping** – identify clusters of highly similar sequences.  
3. **Metagenomic sequence binning** – group contigs or ORFs from environmental datasets.  
4. **Preprocessing for multiple sequence alignments** – create clusters ready for alignment tools.  
5. **Homology-based candidate gene discovery** – focus on connected sequences within gene families.  
6. **Educational / teaching utility** – demonstrate graph-based clustering in a reproducible way.







