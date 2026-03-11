import pytest
from clustrx.clustrx import read_hits, read_fasta, build_clusters
import os

def test_read_hits_blast(tmp_path):
    hits_file = tmp_path / "hits.txt"
    hits_file.write_text("# BLASTP\n"
                         "seq1\tseq2\t99.0\t50\t0\t0\t1\t50\t1\t50\t1e-20\t200.0\n"
                         "seq2\tseq3\t90.0\t50\t0\t0\t1\t50\t1\t50\t1e-2\t90.0\n")
    
    edges, weights, names = read_hits(str(hits_file), format='blast')
    assert len(edges) == 2
    assert len(names) == 3

def test_read_hits_hmmer(tmp_path):
    hits_file = tmp_path / "hits.tblout"
    # Simplified HMMER tblout format
    content = """# tblout format
# target name          accession  query name           accession    E-value  score  bias
#------------------- ---------- -------------------- ---------- --------- ------ -----
target1              -          query1               -            1e-10   100.0   0.1
target2              -          query1               -            1e-5     50.0   0.1
"""
    hits_file.write_text(content.strip())
    
    edges, weights, names = read_hits(str(hits_file), format='hmmer')
    assert len(edges) == 2
    assert len(names) == 3

def test_read_fasta(tmp_path):
    fasta_file = tmp_path / "seqs.fasta"
    fasta_file.write_text(">seq1 first seq\nATGC\n>seq2\nGATT\n")
    
    sequences = read_fasta(str(fasta_file))
    assert sequences == {"seq1": "ATGC", "seq2": "GATT"}

import numpy as np

def test_build_clusters():
    names = ['A', 'B', 'C', 'D', 'E']
    edges = np.array([[0, 1], [1, 2], [3, 4]], dtype=np.int32)
    weights = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    data = (edges, weights, names)
    
    # Min size 2
    clusters = build_clusters(data, min_size=2)
    assert len(clusters) == 2
    assert sorted(clusters[0]) == ['A', 'B', 'C'] or sorted(clusters[1]) == ['A', 'B', 'C']
    assert sorted(clusters[0]) == ['D', 'E'] or sorted(clusters[1]) == ['D', 'E']
    
    # Min size 4
    clusters = build_clusters(data, min_size=4)
    assert len(clusters) == 0

def test_complex_graph():
    # Multi-edge component (1-2-3 loop, 4-5)
    names = ['1', '2', '3', '4', '5']
    edges = np.array([[0, 1], [1, 2], [2, 0], [3, 4]], dtype=np.int32)
    weights = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    data = (edges, weights, names)
    
    clusters = build_clusters(data, min_size=1)
    assert len(clusters) == 2
    assert set(clusters[0]) | set(clusters[1]) == {'1', '2', '3', '4', '5'}
