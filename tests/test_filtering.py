import pytest
from clustrx.clustrx import read_hits
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"

def test_blast_filtering():
    blast_file = DATA_DIR / "sample_blast.out"
    
    # No filtering
    edges, _, _ = read_hits(blast_file, format='blast')
    assert len(edges) == 3
    
    # E-value filter
    edges, _, names = read_hits(blast_file, format='blast', evalue=1e-20)
    assert len(edges) == 1
    assert tuple(names[i] for i in edges[0]) == ("SEQ1", "SEQ2")
    
    # Bitscore filter
    edges, _, names = read_hits(blast_file, format='blast', bitscore=150.0)
    assert len(edges) == 1
    assert tuple(names[i] for i in edges[0]) == ("SEQ1", "SEQ2")

    # Both filters
    edges, _, _ = read_hits(blast_file, format='blast', evalue=1.0, bitscore=80.0)
    assert len(edges) == 3

def test_hmmer_tblout_filtering():
    tblout_file = DATA_DIR / "sample_hmmer.tblout"
    
    # No filtering
    edges, _, _ = read_hits(tblout_file, format='hmmer')
    assert len(edges) == 3
    
    # E-value filter
    edges, _, names = read_hits(tblout_file, format='hmmer', evalue=1e-20)
    assert len(edges) == 1
    assert tuple(names[i] for i in edges[0]) == ("SEQ1", "SEQ2") or tuple(names[i] for i in edges[0]) == ("SEQ2", "SEQ1")
    
    # Bitscore filter
    edges, _, names = read_hits(tblout_file, format='hmmer', bitscore=150.0)
    assert len(edges) == 1
    assert tuple(names[i] for i in edges[0]) == ("SEQ1", "SEQ2") or tuple(names[i] for i in edges[0]) == ("SEQ2", "SEQ1")


