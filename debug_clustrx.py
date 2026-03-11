import sys
import os
sys.path.append(os.getcwd())
from clustrx.clustrx import read_hits, read_fasta, build_clusters, write_clusters

print("Starting debug script...")

hits_file = "examples/hits.txt"
fasta_file = "examples/all_seqs.fasta"
outdir = "examples/output_debug"

print(f"Reading hits from {hits_file}")
edges = read_hits(hits_file, format='blast')
print(f"Edges found: {len(edges)}")
print(f"Sample edges: {edges[:5]}")

print(f"Reading fasta from {fasta_file}")
sequences = read_fasta(fasta_file)
print(f"Sequences found: {len(sequences)}")
print(f"Sample sequence IDs: {list(sequences.keys())[:5]}")

print("Building clusters...")
components = build_clusters(edges, min_size=2)
print(f"Components found: {len(components)}")
print(f"Sample components: {components[:2]}")

print(f"Writing to {outdir}...")
write_clusters(components, sequences, os.path.join(outdir, "clusters"), os.path.join(outdir, "fasta_files"))
print("Done.")
