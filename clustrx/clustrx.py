import os
import subprocess
import time
import gc
import psutil
import igraph as ig
import polars as pl
import numpy as np
import math
import random
from pathlib import Path

def get_ram_usage():
    """
    Returns the current RSS (Resident Set Size) memory usage of the process in GB.
    Used for monitoring performance during large-scale sequence processing.
    """
    process = psutil.Process()
    return process.memory_info().rss / (1024**3)

def get_fasta_info(fasta_file):
    """
    Selective I/O: Extracts sequence names and lengths from FASTA headers.
    
    CRITICAL FOR PERFORMANCE: This function NEVER loads the full sequences into memory.
    It only reads the file line-by-line to identify headers ('>') and increment 
    the byte/character count for lengths. This allows processing millions of 
    sequences with minimal RAM footprint.
    
    Returns:
        names (list): Ordered list of sequence IDs found in the FASTA.
        lengths (dict): Mapping of sequence ID -> length (in characters).
    """
    names = []
    lengths = {}
    with open(fasta_file) as f:
        curr_id = None
        curr_len = 0
        for line in f:
            if line.startswith(">"):
                if curr_id:
                    lengths[curr_id] = curr_len
                # Extract the first word after '>' as the ID
                curr_id = line[1:].split()[0]
                names.append(curr_id)
                curr_len = 0
            else:
                # Accumulate length, stripping potential whitespace/newlines
                curr_len += len(line.strip())
        if curr_id:
            lengths[curr_id] = curr_len
    return names, lengths

def validate_format(file_path, fmt, custom_cols=None):
    """
    Stric validation of the input hits file against the requested format.
    Ensures that technical errors (e.g., using a BLAST file with --fmt hmmer)
    are caught before processing begins.
    """
    try:
        with open(file_path, 'r') as f:
            # Check first 20 lines for identifying headers or structure
            lines = [f.readline() for _ in range(20)]
            lines = [l for l in lines if l]
            
            if not lines:
                raise ValueError(f"Input file '{file_path}' is empty.")

            if fmt == 'hmmer' or fmt == 'domhmmer' or fmt == 'tblhmmer':
                # HMMER files always start with identifying comments
                if not any(l.startswith("# tblout") or l.startswith("# target name") or l.startswith("# domtblout") for l in lines):
                    raise ValueError(f"File '{file_path}' does not appear to be a valid HMMER output. Did you specify '--fmt hmmer' by mistake?")
                    
            elif fmt == 'blast':
                is_blast = False
                for l in lines:
                    # Look for BLAST headers (tab-separated) or typical 10+ columns structure
                    if l.startswith("# BLAST") or l.startswith("# Iteration"): 
                        is_blast = True
                        break
                    if not l.startswith("#") and len(l.split('\t')) >= 10:
                        is_blast = True
                        break
                if not is_blast:
                    raise ValueError(f"File '{file_path}' does not appear to be a valid BLAST tabular file (-outfmt 6).")
                    
            elif fmt == 'custom':
                is_tabular = False
                max_idx = max(custom_cols.values()) if custom_cols else 0
                for l in lines:
                    if not l.startswith("#") and '\t' in l:
                        cols = l.split('\t')
                        if len(cols) > max_idx:
                            is_tabular = True
                            break
                if not is_tabular:
                    raise ValueError(f"File '{file_path}' does not appear to be a valid tab-separated file, or has fewer columns than requested (max index: {max_idx}).")
    except Exception as e:
        if isinstance(e, ValueError):
            raise e
        pass

def detect_format(file_path):
    """
    Heuristic-based format detection.
    Scans the beginning of the file to determine if it's BLAST, HMMER (tblout), 
    or HMMER (domtblout).
    """
    try:
        with open(file_path, 'r') as f:
            for _ in range(30): # Scan first 30 lines
                line = f.readline()
                if not line: break
                
                # Keywords for domtblout (alignment coordinates)
                if any(x in line for x in ["# domtblout", "hmm coord", "ali coord", "env coord"]):
                    return 'domhmmer'
                
                # Keywords for tblout (no coordinates)
                if any(x in line for x in ["# tblout", "--- best 1 domain ---"]):
                    return 'tblhmmer'
                
                # Keywords for BLAST
                if line.startswith("# BLASTP") or line.startswith("# Iteration"): return 'blast'
                if not line.startswith("#") and len(line.split('\t')) >= 10: return 'blast'
    except: pass
    return 'blast' # Default to BLAST if uncertain

def read_hits(file_path, format=None, evalue=None, bitscore=None, pident=None, coverage=None, lengths=None, pident_override=90.0, mapping_df=None, lengths_df=None, custom_cols=None):
    """
    Optimized hit reading and filtering using Polars (DataFrame library).
    
    This function implements:
    1. Parsing different bioinformatics formats (BLAST, HMMER, MMseqs...) into a unified internal schema.
    2. Applying hard filters for E-value, Bitscore, Identity, and Coverage.
    3. Calculating Dynamic Coverage: Ensuring hits represent meaningful alignments relative 
       to the full sequence length (fetched from FASTA info).
    4. Deduplicating hits: Keeping only the strongest hit (max Bitscore) between two sequences.
    
    Returns:
        u_v (ndarray): 2-column array of vertex ID pairs (edges).
        weights (ndarray): 1-column array of bitscores (edge weights).
        v_list (list): Map of vertex IDs back to sequence names.
    """
def _detect_needed_columns(format, evalue, pident, coverage, custom_cols):
    """
    Determines required columns and their names based on the input format and filters.
    
    Returns:
        needed_cols (list): List of column indices to read.
        col_names (dict): Map of indices to internal semantic names.
        schema (dict): Polars schema overrides for these columns.
    """
    if format == 'domhmmer':
        needed_cols = [0, 3, 13] # target, query, score
        if evalue is not None: needed_cols.append(6)
        if coverage is not None or pident is not None:
            needed_cols.extend([15, 16, 17, 18]) # coord columns
        col_names = {0: 't', 3: 'q', 13: 'bitscore', 6: 'evalue', 15: 'hmm_f', 16: 'hmm_t', 17: 'ali_f', 18: 'ali_t'}
        schema = None
    elif format in ['hmmer', 'tblhmmer']:
        needed_cols = [0, 2, 5] # target, query, score
        if evalue is not None: needed_cols.append(4)
        col_names = {0: 't', 2: 'q', 5: 'bitscore', 4: 'evalue'}
        schema = None 
    elif format == 'custom':
        if not custom_cols or 'q' not in custom_cols or 't' not in custom_cols or 'bitscore' not in custom_cols:
            raise ValueError("Custom format requires: q, t, bitscore columns.")
            
        needed_cols = [custom_cols['q'], custom_cols['t'], custom_cols['bitscore']]
        col_names = {custom_cols['q']: 'q', custom_cols['t']: 't', custom_cols['bitscore']: 'bitscore'}
        
        if evalue is not None and 'evalue' in custom_cols:
            needed_cols.append(custom_cols['evalue'])
            col_names[custom_cols['evalue']] = 'evalue'
        if (pident is not None or coverage is not None) and 'pident' in custom_cols:
            needed_cols.append(custom_cols['pident'])
            col_names[custom_cols['pident']] = 'pident'
        if coverage is not None and 'length' in custom_cols:
            needed_cols.append(custom_cols['length'])
            col_names[custom_cols['length']] = 'length'
            
        needed_cols = list(set(needed_cols))
        schema = {f'column_{i+1}': pl.String if i in [custom_cols['q'], custom_cols['t']] 
                  else (pl.Int32 if 'length' in custom_cols and i == custom_cols['length'] else pl.Float64)
                  for i in needed_cols}
    else:
        # BLAST (Default)
        needed_cols = [0, 1, 11]
        if evalue is not None: needed_cols.append(10)
        if pident is not None or coverage is not None: needed_cols.append(2)
        if coverage is not None: needed_cols.append(3)
        
        col_names = {0: 'q', 1: 't', 11: 'bitscore'}
        if evalue is not None: col_names[10] = 'evalue'
        if 2 in needed_cols: col_names[2] = 'pident'
        if 3 in needed_cols: col_names[3] = 'length'
        schema = {f'column_{i+1}': pl.String if i in [0, 1] else (pl.Int32 if i == 3 else pl.Float64) for i in needed_cols}

    return needed_cols, col_names, schema

def _parse_hmmer_hits(file_path, format, needed_cols, col_names):
    """Parses HMMER-specific variable-whitespace formats into a DataFrame."""
    raw_df = pl.read_csv(file_path, separator='\n', has_header=False, comment_prefix='#', ignore_errors=True)
    if len(raw_df) == 0: return None
    
    col = raw_df.columns[0]
    parsed_df = raw_df.select(
        pl.col(col).str.strip_chars().str.replace_all(r"\s+", "\t").str.split("\t").alias("cols")
    )
    
    exprs = []
    for idx in needed_cols:
        name = col_names[idx]
        if name in ['q', 't']:
            exprs.append(pl.col("cols").list.get(idx).alias(name))
        else:
            exprs.append(pl.col("cols").list.get(idx).cast(pl.Float64, strict=False).alias(name))
            
    df = parsed_df.select(exprs).drop_nulls(subset=['q', 't', 'bitscore'])

    if format == 'domhmmer' and 'hmm_f' in df.columns:
        df = df.with_columns([
            (pl.col('hmm_t').cast(pl.Float32) - pl.col('hmm_f').cast(pl.Float32)).abs().alias('length'),
            pl.lit(0.0).alias('pident')
        ])
    return df

def _apply_scientific_filters(df, lengths_df, evalue, bitscore, pident, coverage, pident_override):
    """
    Applies the core scientific filtering logic using Polars Lazy-ready operations.
    Includes E-value, Bitscore, Identity, and the Dynamic Coverage formula.
    """
    # 1. Base Filters
    mask = pl.col('q') != pl.col('t')
    if evalue is not None and 'evalue' in df.columns: 
        mask = mask & (pl.col('evalue') <= evalue)
    if bitscore is not None: 
        mask = mask & (pl.col('bitscore') >= bitscore)
    if pident is not None and 'pident' in df.columns:
        mask = mask & (pl.col('pident') >= float(pident))
    
    df = df.filter(mask)

    # 2. Coverage & Dynamic Logic
    if lengths_df is not None and 'length' in df.columns:
        df = df.join(lengths_df, left_on='q', right_on='name', how='left').rename({'len': 'l_q'})
        df = df.join(lengths_df, left_on='t', right_on='name', how='left').rename({'len': 'l_t'})
        df = df.with_columns([pl.col('l_q').fill_null(1), pl.col('l_t').fill_null(1)])
        
        m_lens = pl.max_horizontal('l_q', 'l_t')
        c_lens = m_lens.clip(upper_bound=500).cast(pl.Float32)
        
        if coverage is not None:
            # UMBERAL DINÁMICO: t_cov = 1.0 - (L_max / 600) clamped [0.4, 0.8]
            t_cov = (pl.lit(1.0) - (c_lens / 600.0)).clip(0.4, 0.8) if coverage == 'dynamic' else pl.lit(float(coverage), dtype=pl.Float32)
            cov_mask = (pl.col('length') / m_lens) >= t_cov
            
            if pident_override is not None and 'pident' in df.columns:
                rescue_mask = (pl.col('pident') >= pident_override) & (m_lens >= 20)
                cov_mask = cov_mask | rescue_mask
            
            df = df.filter(cov_mask)
            
    return df

def read_hits(file_path, format=None, lengths_df=None, mapping_df=None, 
              evalue=None, bitscore=None, pident=None, coverage=None, 
              pident_override=90.0, custom_cols=None):
    """
    Parses a hits file and returns edge list with weights.
    Efficiently handles BLAST, HMMER, and custom tabular formats.
    """
    if format is None:
        format = detect_format(file_path)
    else:
        validate_format(file_path, format, custom_cols)
    
    print(f"  [RAM] Starting Data Processing ({format}): {get_ram_usage():.2f} GB", flush=True)
    t0 = time.perf_counter()

    # 1. Detection
    needed_cols, col_names, schema = _detect_needed_columns(format, evalue, pident, coverage, custom_cols)

    # 2. Parsing
    if format in ['blast', 'custom']:
        df = pl.read_csv(file_path, has_header=False, comment_prefix='#', ignore_errors=True,
                        columns=[f'column_{i+1}' for i in needed_cols],
                        schema_overrides=schema, separator='\t')
        df = df.rename({f'column_{i+1}': name for i, name in col_names.items()})
    else:
        df = _parse_hmmer_hits(file_path, format, needed_cols, col_names)
        if df is None: return [], [], []

    # 3. Scientific Filtering
    # Early type-cast sequence names based on reference (Enums optimization)
    q_type = pl.String
    if mapping_df is not None: q_type = mapping_df.schema['name']
    elif lengths_df is not None: q_type = lengths_df.schema['name']

    df = df.with_columns([
        pl.col('q').cast(q_type, strict=False),
        pl.col('t').cast(q_type, strict=False)
    ]).drop_nulls(subset=['q', 't'])

    df = _apply_scientific_filters(df, lengths_df, evalue, bitscore, pident, coverage, pident_override)

    # 4. Community Preparation
    df = df.group_by(['q', 't']).agg(pl.col('bitscore').max())

    if mapping_df is not None:
        df = df.join(mapping_df.select(['name', 'id']), left_on='q', right_on='name')
        df = df.join(mapping_df.select(['name', 'id']), left_on='t', right_on='name', suffix='_t')
        u_v, weights = df.select(['id', 'id_t']).to_numpy(), df.get_column('bitscore').cast(pl.Float64).to_numpy()
        v_list = mapping_df.sort("id").get_column("name").to_list()
    else:
        all_verts = pl.concat([df['q'], df['t']]).unique().sort()
        v_list = all_verts.to_list()
        v_map = {name: i for i, name in enumerate(v_list)}
        u_v = df.with_columns([pl.col('q').replace(v_map).cast(pl.Int32), pl.col('t').replace(v_map).cast(pl.Int32)]).select(['q', 't']).to_numpy()
        weights = df.get_column('bitscore').cast(pl.Float64).to_numpy()

    print(f"  [PHASE] Data Processing took {time.perf_counter()-t0:.4f}s. Valid Edges: {len(u_v)}", flush=True)
    del df; gc.collect()
    return u_v, weights, v_list

def build_clusters(data, min_size=1, seed=None, resolution=1.0):
    """
    Applies the Leiden community detection algorithm to identify sequence clusters.
    
    Logic:
    1. Creates an undirected graph (via igraph).
    2. Weights edges using Bitscores or similar score from sequence similarity search tools.
    3. Maximizes Modularity with a specific 'resolution' parameter to control cluster size.
    
    Resolution:
    - 1.0: Standard behavior.
    - > 1.0: Produces more, smaller clusters (higher stringency).
    - < 1.0: Produces fewer, larger clusters.
    """
    print(f"  [RAM] Start build_clusters: {get_ram_usage():.2f} GB", flush=True)
    edges, weights, names = data
    if len(edges) == 0: return []

    g = ig.Graph(len(names), edges, directed=False)
    del edges; gc.collect()

    if seed is not None: random.seed(seed)
    partition = g.community_leiden(
        weights=weights, 
        objective_function='modularity', 
        resolution_parameter=resolution
    )
    
    membership = partition.membership
    del weights; del partition; gc.collect()

    # Extract clusters and filter by minimum size
    df = pl.DataFrame({"name": names, "cluster": membership})
    del names; del membership; del g; gc.collect()

    clusters_df = df.group_by("cluster").agg(
        pl.col("name").sort()
    ).filter(pl.col("name").list.len() >= min_size)
    
    return clusters_df.get_column("name").to_list()

def read_fasta(fasta_file):
    sequences = {}
    with open(fasta_file) as f:
        seq_id, seq_lines = None, []
        for line in f:
            if line.startswith(">"):
                if seq_id: sequences[seq_id] = "".join(seq_lines)
                seq_id = line[1:].split()[0]
                seq_lines = []
            else: seq_lines.append(line.strip())
        if seq_id: sequences[seq_id] = "".join(seq_lines)
    return sequences

def write_clusters(components, sequences, outdir_clusters, outdir_fastas, outdir_alignments=None):
    os.makedirs(outdir_clusters, exist_ok=True)
    os.makedirs(outdir_fastas, exist_ok=True)
    if outdir_alignments: os.makedirs(outdir_alignments, exist_ok=True)
    for i, comp in enumerate(components, start=1):
        with open(Path(outdir_clusters) / f"cluster_{i}.txt", "w") as cf:
            cf.write("\n".join(comp) + "\n")
        fasta_file = Path(outdir_fastas) / f"cluster_{i}.fasta"
        with open(fasta_file, "w") as ff:
            found = 0
            for seq_id in comp:
                seq = sequences.get(seq_id)
                if seq:
                    ff.write(f">{seq_id}\n{seq}\n"); found += 1
        if outdir_alignments and found > 1:
            try:
                cmd = ["mafft", "--auto", str(fasta_file)]
                with open(Path(outdir_alignments) / f"cluster_{i}.fasta", "w") as out:
                    subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=True)
            except: pass
