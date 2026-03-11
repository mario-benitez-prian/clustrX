import os
import subprocess
from pathlib import Path
import pytest

def test_cli_integration(tmp_path):
    # Setup files
    hits = tmp_path / "hits.txt"
    hits.write_text("# BLASTP\n"
                    "A\tB\t100.0\t50\t0\t0\t1\t50\t1\t50\t1e-10\t100.0\n"
                    "B\tC\t100.0\t50\t0\t0\t1\t50\t1\t50\t1e-10\t100.0\n"
                    "D\tE\t100.0\t50\t0\t0\t1\t50\t1\t50\t1e-10\t100.0\n")
    
    fasta = tmp_path / "seqs.fasta"
    fasta.write_text(">A\nAAA\n>B\nBBB\n>C\nCCC\n>D\nDDD\n>E\nEEE\n")
    
    outdir = tmp_path / "output"
    
    import sys
    # Run CLI
    result = subprocess.run([
        sys.executable, "-m", "clustrx.cli",
        "-i", str(hits),
        "-f", str(fasta),
        "--outdir", str(outdir),
        "--format", "blast",
        "-min", "2"
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    assert "2 clusters written" in result.stdout
    
    # Verify files
    assert (outdir / "clusters" / "cluster_1.txt").exists()
    assert (outdir / "fasta_files" / "cluster_1.fasta").exists()
    
    # Verify cluster contents exist
    cluster_files = list((outdir / "clusters").glob("cluster_*.txt"))
    assert len(cluster_files) == 2
    
    clusters = []
    for cf in cluster_files:
        with open(cf) as f:
            clusters.append(set(f.read().splitlines()))
    
    assert {"A", "B", "C"} in clusters
    assert {"D", "E"} in clusters
