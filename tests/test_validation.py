import pytest
from clustrx.clustrx import read_hits
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"

def test_blast_file_as_hmmer():
    """Ensure passing a BLAST format file requesting HMMER format throws ValueError."""
    blast_file = DATA_DIR / "sample_blast.out"
    with pytest.raises(ValueError, match="does not appear to be a valid HMMER output"):
        read_hits(blast_file, format='hmmer')

def test_hmmer_file_as_blast():
    """Ensure passing a HMMER format file requesting BLAST format throws ValueError."""
    hmmer_file = DATA_DIR / "sample_hmmer.tblout"
    with pytest.raises(ValueError, match="does not appear to be a valid BLAST tabular file"):
        read_hits(hmmer_file, format='blast')

def test_invalid_custom_format(tmp_path):
    """Ensure passing a file with less columns than the requested custom format throws ValueError."""
    custom_file = tmp_path / "custom.tsv"
    # Only 2 columns
    custom_file.write_text("SEQ1\tSEQ2\n")
    
    # But requesting columns 0,1,2,3 (which require at least 4 cols)
    custom_cols = {'q': 0, 't': 1, 'bitscore': 2, 'evalue': 3}
    with pytest.raises(ValueError, match="fewer columns than requested"):
        read_hits(str(custom_file), format='custom', custom_cols=custom_cols)
