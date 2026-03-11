import os
import shutil
import subprocess
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).parent.absolute()
project_root = str(BASE_DIR.parent)
CONDA_PYTHON = "/home/mario/miniconda3/envs/amp_miner_env/bin/python3"
TEST_DATA_DIR = BASE_DIR / "data"

# Test input files
FASTA = TEST_DATA_DIR / "sample.fasta"
BLAST_HITS = TEST_DATA_DIR / "sample_blast.out"
HMMER_HITS = TEST_DATA_DIR / "sample_hmmer.tblout"
RESULTS_DIR = BASE_DIR / "results"

def setup_dummy_data():
    TEST_DATA_DIR.mkdir(parents=True, exist_ok=True)
    
    with open(FASTA, "w") as f:
        f.write(">SEQ1\nACDEFGH\n>SEQ2\nACDEFGH\n>SEQ3\nACDEFGH\n>SEQ4\nACDEFGH\n")
            
    with open(BLAST_HITS, "w") as f:
        # Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        f.write("SEQ1\tSEQ2\t100.0\t50\t0\t0\t1\t50\t1\t50\t1e-10\t100.0\n")
        f.write("SEQ1\tSEQ3\t90.0\t50\t0\t0\t1\t50\t1\t50\t1e-8\t90.0\n")
        f.write("SEQ3\tSEQ4\t95.0\t50\t0\t0\t1\t50\t1\t50\t1e-9\t95.0\n")
        
    with open(HMMER_HITS, "w") as f:
        # Standard HMMER tblout format with irregular spaces
        f.write("# tblout format\n")
        f.write("# comment line 2\n")
        f.write("SEQ2        -          SEQ1                 -    1e-10   100.0   0.0   1.2e-10  100.0   0.0   1.0   1   0   0   1   1   1   1 - \n")
        f.write("SEQ3      -         SEQ1                  -    1e-8    90.0    0.0   1.2e-8   90.0    0.0   1.0   1   0   0   1   1   1   1 - \n")
        f.write("SEQ4             -  SEQ3                  -    1e-9    95.0    0.0   1.2e-9   95.0    0.0   1.0   1   0   0   1   1   1   1 - \n")

def run_cmd(cmd):
    env = os.environ.copy()
    env["PYTHONPATH"] = project_root + os.pathsep + env.get("PYTHONPATH", "")
    print(f">> Executing: {cmd}")
    ret = subprocess.run(cmd, shell=True, env=env)
    return ret.returncode

def parse_clustrx(outdir):
    clusters = []
    p = Path(outdir) / "clusters"
    if not p.exists(): return []
    for cf in p.glob("cluster_*.txt"):
        with open(cf) as f:
            members = sorted([line.strip() for line in f if line.strip()])
            if members: clusters.append(members)
    return sorted(clusters)

def run_test():
    setup_dummy_data()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Run using Standard BLAST format
    out_blast = RESULTS_DIR / "out_blast_hmmer"
    if out_blast.exists(): shutil.rmtree(out_blast)
    cmd_blast = f"{CONDA_PYTHON} -m clustrx.cli -i {BLAST_HITS} -f {FASTA} -fmt blast -e 1e-5 --seed 42 --min-cluster-size 1 --outdir {out_blast}"
    run_cmd(cmd_blast)
    
    # Run using Native HMMER format
    out_hmmer = RESULTS_DIR / "out_hmmer"
    if out_hmmer.exists(): shutil.rmtree(out_hmmer)
    cmd_hmmer = f"{CONDA_PYTHON} -m clustrx.cli -i {HMMER_HITS} -f {FASTA} -fmt hmmer -e 1e-5 --seed 42 --min-cluster-size 1 --outdir {out_hmmer}"
    run_cmd(cmd_hmmer)
    
    # Compare
    clusters_blast = parse_clustrx(out_blast)
    clusters_hmmer = parse_clustrx(out_hmmer)
    
    print("\n--- TEST RESULTS ---")
    print(f"Standard BLAST clusters: {clusters_blast}")
    print(f"Native HMMER clusters:   {clusters_hmmer}")
    
    if clusters_blast == clusters_hmmer and len(clusters_blast) > 0:
        print("\n✅ SUCCESS: Native HMMER parser produces identical mathematical outputs as standard BLAST parser!")
    else:
        print("\n❌ FAILED: Clusters do not match or are empty.")

if __name__ == "__main__":
    run_test()
