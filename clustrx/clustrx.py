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
    """Returns current RSS memory usage in GB."""
    process = psutil.Process()
    return process.memory_info().rss / (1024**3)

def get_fasta_info(fasta_file):
    """
    Extracts sequence names and lengths from FASTA headers without loading sequences into memory.
    Returns: names (list), lengths (dict).
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
                curr_id = line[1:].split()[0]
                names.append(curr_id)
                curr_len = 0
            else:
                curr_len += len(line.strip())
        if curr_id:
            lengths[curr_id] = curr_len
    return names, lengths

def validate_format(file_path, fmt, custom_cols=None):
    """Validates that the file contents match the explicitly requested format."""
    try:
        with open(file_path, 'r') as f:
            lines = [f.readline() for _ in range(20)]
            lines = [l for l in lines if l]
            
            if not lines:
                raise ValueError(f"Input file '{file_path}' is empty.")

            if fmt == 'hmmer':
                if not any(l.startswith("# tblout") or l.startswith("# target name") for l in lines):
                    raise ValueError(f"File '{file_path}' does not appear to be a valid HMMER tblout file. Did you specify '--fmt hmmer' by mistake?")
                    
            elif fmt == 'blast':
                is_blast = False
                for l in lines:
                    if l.startswith("# BLAST") or l.startswith("# Iteration"): 
                        is_blast = True
                        break
                    if not l.startswith("#") and len(l.split('\t')) >= 10:
                        is_blast = True
                        break
                if not is_blast:
                    raise ValueError(f"File '{file_path}' does not appear to be a valid BLAST tabular file (-outfmt 6). Did you specify '--fmt blast' by mistake?")
                    
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
                    raise ValueError(f"File '{file_path}' does not appear to be a valid tab-separated file, or it has fewer columns than your maximum requested index ({max_idx}).")
    except Exception as e:
        if isinstance(e, ValueError):
            raise e
        pass

def detect_format(file_path):
    """Detects if a file is BLAST (tabular) or HMMER (tblout)."""
    try:
        with open(file_path, 'r') as f:
            for _ in range(20): # Check first 20 lines
                line = f.readline()
                if not line: break
                if line.startswith("# tblout"): return 'hmmer'
                if line.startswith("# BLASTP") or line.startswith("# Iteration"): return 'blast'
                if not line.startswith("#") and len(line.split('\t')) >= 10: return 'blast'
    except: pass
    return 'blast' # Default

def read_hits(file_path, format=None, evalue=None, bitscore=None, pident=None, coverage=None, lengths=None, pident_override=90.0, mapping_df=None, lengths_df=None, custom_cols=None):
    """
    Optimized hit reading using Selective I/O.
    Only loads the columns strictly necessary for requested filters.
    """
    if format is None:
        format = detect_format(file_path)
    else:
        validate_format(file_path, format, custom_cols)
    
    print(f"  [RAM] Starting Data Processing ({format}): {get_ram_usage():.2f} GB", flush=True)
    t0 = time.perf_counter()

    # Define column indices and schema based on format
    # BLAST: 0:q, 1:t, 2:pident, 3:length, 10:evalue, 11:bitscore
    # HMMER: 0:target, 2:query, 4:evalue, 5:score
    if format == 'hmmer':
        needed_cols = [0, 2, 5] # target, query, score
        if evalue is not None: needed_cols.append(4)
        
        col_names = {0: 't', 2: 'q', 5: 'bitscore'}
        if evalue is not None: col_names[4] = 'evalue'
        
        # We will parse as a single text column first
        schema = None 
    elif format == 'custom':
        if custom_cols is None or 'q' not in custom_cols or 't' not in custom_cols or 'bitscore' not in custom_cols:
            raise ValueError("For custom format, you must specify at least query, target, and bitscore columns (--col-query, --col-target, --col-bitscore).")
            
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
            
        needed_cols = list(set(needed_cols)) # Deduplicate
        
        schema = {}
        for i in needed_cols:
            if i in [custom_cols['q'], custom_cols['t']]:
                schema[f'column_{i+1}'] = pl.String
            elif 'length' in custom_cols and i == custom_cols['length']:
                schema[f'column_{i+1}'] = pl.Int32
            else:
                schema[f'column_{i+1}'] = pl.Float64
    else:
        # BLAST
        needed_cols = [0, 1, 11] # q, t, bitscore
        if evalue is not None: needed_cols.append(10)
        if pident is not None or coverage is not None: needed_cols.append(2) # pident
        if coverage is not None: needed_cols.append(3) # length
        
        col_names = {0: 'q', 1: 't', 11: 'bitscore'}
        if evalue is not None: col_names[10] = 'evalue'
        if 2 in needed_cols: col_names[2] = 'pident'
        if 3 in needed_cols: col_names[3] = 'length'

        schema = {f'column_{i+1}': pl.String if i in [0, 1] else pl.Float64 for i in needed_cols}
        if 3 in needed_cols: schema['column_4'] = pl.Int32

    if format in ['blast', 'custom']:
        # Use native read_csv with 'columns' optimization for tabs
        read_csv_kwargs = {
            'has_header': False, 
            'comment_prefix': '#', 
            'ignore_errors': True,
            'columns': [f'column_{i+1}' for i in needed_cols],
            'schema_overrides': schema,
            'separator': '\t'
        }
        df = pl.read_csv(file_path, **read_csv_kwargs)
        # Rename columns to standard internal names
        df = df.rename({f'column_{i+1}': name for i, name in col_names.items()})

    elif format == 'hmmer':
        # Native Polars String Manipulation for variable whitespace
        # 1. Read entire file as single column string ignoring comments
        raw_df = pl.read_csv(file_path, separator='\n', has_header=False, comment_prefix='#', ignore_errors=True)
        if len(raw_df) == 0:
            return [], [], []
            
        col = raw_df.columns[0]
        # 2. Strip edge whitespace, replace multi-space with tab, and split into list
        parsed_df = raw_df.select(
            pl.col(col).str.strip_chars().str.replace_all(r"\s+", "\t").str.split("\t").alias("cols")
        )
        
        # 3. Extract the requested indexed fields
        exprs = []
        for idx in needed_cols:
            name = col_names[idx]
            if name in ['q', 't']:
                exprs.append(pl.col("cols").list.get(idx).alias(name))
            else:
                exprs.append(pl.col("cols").list.get(idx).cast(pl.Float64, strict=False).alias(name))
                
        df = parsed_df.select(exprs).drop_nulls(subset=['q', 't', 'bitscore'])

    # Determine Name type (Enum if zero-string, else String)
    q_type = pl.String
    if mapping_df is not None:
        q_type = mapping_df.schema['name']
    elif lengths_df is not None:
        q_type = lengths_df.schema['name']

    # Cast to Enum but with strict=False to convert unknown IDs to null.
    # Then drop nulls to automatically filter out hits not present in mapping/FASTA.
    df = df.with_columns([
        pl.col('q').cast(q_type, strict=False),
        pl.col('t').cast(q_type, strict=False)
    ]).drop_nulls(subset=['q', 't'])

    # 1. Base Filters (No Self-loops)
    mask = pl.col('q') != pl.col('t')
    if evalue is not None and 'evalue' in df.columns: 
        mask = mask & (pl.col('evalue') <= evalue)
    if bitscore is not None: 
        mask = mask & (pl.col('bitscore') >= bitscore)
    if pident is not None and 'pident' in df.columns:
        mask = mask & (pl.col('pident') >= float(pident))
    
    df = df.filter(mask)

    # 2. Coverage & Dynamic Filters (Only if necessary columns exist)
    if lengths_df is not None and 'pident' in df.columns and 'length' in df.columns:
        df = df.join(lengths_df, left_on='q', right_on='name', how='left').rename({'len': 'l_q'})
        df = df.join(lengths_df, left_on='t', right_on='name', how='left').rename({'len': 'l_t'})
        df = df.with_columns([pl.col('l_q').fill_null(1), pl.col('l_t').fill_null(1)])
        
        m_lens = pl.max_horizontal('l_q', 'l_t')
        c_lens = m_lens.clip(upper_bound=500).cast(pl.Float32)
        
        cb_mask = pl.lit(True)

        if coverage is not None:
            t_cov = (pl.lit(1.0) - (c_lens / 600.0)).clip(0.4, 0.8) if coverage == 'dynamic' else pl.lit(float(coverage), dtype=pl.Float32)
            cov_mask = (pl.col('length') / m_lens) >= t_cov
            if pident_override is not None:
                rescue_mask = (pl.col('pident') >= pident_override) & (m_lens >= 20)
                cov_mask = cov_mask | rescue_mask
            cb_mask = cb_mask & cov_mask
        
        df = df.filter(cb_mask)

    # 3. Deduplication
    df = df.group_by(['q', 't']).agg(pl.col('bitscore').max())

    # 4. ID Mapping
    if mapping_df is not None:
        # Zero-String: Extract IDs directly
        df = df.join(mapping_df.select(['name', 'id']), left_on='q', right_on='name')
        df = df.join(mapping_df.select(['name', 'id']), left_on='t', right_on='name', suffix='_t')
        
        u_v = df.select(['id', 'id_t']).to_numpy()
        weights = df.get_column('bitscore').cast(pl.Float64).to_numpy()
        v_list = mapping_df.sort("id").get_column("name").to_list()
    else:
        # Legacy Discovery Mode
        print("  [LEGACY] Performing vertex discovery...", flush=True)
        all_verts = pl.concat([df['q'], df['t']]).unique().sort()
        v_list = all_verts.to_list()
        v_map = {name: i for i, name in enumerate(v_list)}
        
        u_v = df.with_columns([
            pl.col('q').replace(v_map).cast(pl.Int32),
            pl.col('t').replace(v_map).cast(pl.Int32)
        ]).select(['q', 't']).to_numpy()
        weights = df.get_column('bitscore').cast(pl.Float64).to_numpy()

    t1 = time.perf_counter()
    print(f"  [PHASE] Data Processing took {t1-t0:.4f}s", flush=True)
    print(f"  [STATS] Total Hits: {len(u_v)}", flush=True)
    
    del df
    gc.collect()
    return u_v, weights, v_list

def build_clusters(data, min_size=1, seed=None, resolution=1.0):
    """Build graph and extract clusters efficiently using Leiden algorithm."""
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
