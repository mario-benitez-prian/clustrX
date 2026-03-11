import os
import shutil
import subprocess
import time
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).parent.absolute()
project_root = str(BASE_DIR.parent)
CONDA_PYTHON = "/home/mario/miniconda3/envs/amp_miner_env/bin/python3"
TEST_DATA_DIR = BASE_DIR / "data"

# Test input files (Ensure these exist or we generate small dummies)
FASTA = TEST_DATA_DIR / "sample.fasta"
BLAST_HITS = TEST_DATA_DIR / "sample_blast.out"
RESULTS_DIR = BASE_DIR / "results"

def setup_dummy_data():
    TEST_DATA_DIR.mkdir(parents=True, exist_ok=True)
    
    with open(FASTA, "w") as f:
        f.write(">SEQ1\nACDEFGH\n>SEQ2\nACDEFGH\n>SEQ3\nACDEFGH\n>SEQ4\nACDEFGH\n")
            
    with open(BLAST_HITS, "w") as f:
        # Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        f.write("# BLASTP\n")
        f.write("SEQ1\tSEQ2\t100.0\t50\t0\t0\t1\t50\t1\t50\t1e-10\t100.0\n")
        f.write("SEQ1\tSEQ3\t90.0\t50\t0\t0\t1\t50\t1\t50\t1e-8\t90.0\n")
        f.write("SEQ3\tSEQ4\t95.0\t50\t0\t0\t1\t50\t1\t50\t1e-9\t95.0\n")

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
    out_blast = RESULTS_DIR / "out_blast"
    if out_blast.exists(): shutil.rmtree(out_blast)
    cmd_blast = f"{CONDA_PYTHON} -m clustrx.cli -i {BLAST_HITS} -f {FASTA} -fmt blast -e 1e-5 --seed 42 --min-cluster-size 1 --outdir {out_blast}"
    run_cmd(cmd_blast)
    
    # Run using Custom format masquerading as BLAST
    # column_indices: q=0, t=1, pident=2, length=3, evalue=10, bitscore=11
    out_custom = RESULTS_DIR / "out_custom"
    if out_custom.exists(): shutil.rmtree(out_custom)
    cmd_custom = f"{CONDA_PYTHON} -m clustrx.cli -i {BLAST_HITS} -f {FASTA} --col-query 0 --col-target 1 --col-pident 2 --col-length 3 --col-evalue 10 --col-bitscore 11 -e 1e-5 --seed 42 --min-cluster-size 1 --outdir {out_custom}"
    run_cmd(cmd_custom)
    
    # Compare
    clusters_blast = parse_clustrx(out_blast)
    clusters_custom = parse_clustrx(out_custom)
    
    print("\n--- TEST RESULTS ---")
    print(f"Standard BLAST clusters: {clusters_blast}")
    print(f"Custom format clusters:  {clusters_custom}")
    
    if clusters_blast == clusters_custom and len(clusters_blast) > 0:
        print("\n✅ SUCCESS: Custom mapping parser produces identical mathematical outputs as standard C parser!")
    else:
        print("\n❌ FAILED: Clusters do not match or are empty.")

if __name__ == "__main__":
    run_test()
