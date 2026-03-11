import argparse
import polars as pl
from pathlib import Path
from .clustrx import read_hits, get_fasta_info, build_clusters, write_clusters

def main():
    parser = argparse.ArgumentParser(
        description=(
        "clustRX: Cluster sequences from BLAST/HMMER outputs and generate FASTA groups.\n\n"
        "Author: Mario Benítez-Prián | Please cite clustRX if used in your research.\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True, help="Input hits file (BLAST tabular or HMMER tblout)")
    parser.add_argument("-fmt", "--format", choices=["blast", "hmmer", "custom"], help="Input format (optional, auto-detected if not provided)")
    parser.add_argument("-f", "--fasta", required=True, help="FASTA file with all sequences")
    parser.add_argument("-min", "--min-cluster-size", type=int, default=2, help="Minimum cluster size to output (default=2).")
    parser.add_argument("-e", "--evalue", type=float, help="E-value threshold for filtering hits (default: no filter)")
    parser.add_argument("-b", "--bitscore", type=float, help="Bitscore threshold for filtering hits (default: no filter)")
    parser.add_argument("-pi", "--pidentity", type=float, help="Minimum percentage identity (0-100) (default: no filter)")
    parser.add_argument("-c", "--coverage", help="Alignment coverage filter (0.0-1.0 or 'dynamic'). Percentage of the longest sequence covered by the alignment.")
    parser.add_argument("--id-override", type=float, default=90.0, help="Identity threshold to override coverage filter (default: 90.0).")
    parser.add_argument("--seed", type=int, help="Seed for reproducibility (default: None)")
    parser.add_argument("--resolution", type=float, default=1.0, help="Resolution parameter for Leiden algorithm (default: 1.0)")
    parser.add_argument("--outdir", default="clustrx_output", help="Output directory")
    parser.add_argument("--mafft", action="store_true", help="Automatically run MAFFT on generated clusters (requires MAFFT installed)")
    
    # Custom format groups
    custom_group = parser.add_argument_group('Custom Format Options', 'Specify 0-indexed column numbers for custom tabular inputs.')
    custom_group.add_argument("--col-query", type=int, help="Column index for query ID")
    custom_group.add_argument("--col-target", type=int, help="Column index for target ID")
    custom_group.add_argument("--col-bitscore", type=int, help="Column index for alignment score (bitscore)")
    custom_group.add_argument("--col-evalue", type=int, help="Column index for E-value")
    custom_group.add_argument("--col-pident", type=int, help="Column index for percentage identity")
    custom_group.add_argument("--col-length", type=int, help="Column index for query length (used in dynamic coverage)")
    args = parser.parse_args()

    if args.mafft:
        import shutil
        if not shutil.which("mafft"):
            print("Error: MAFFT is not installed or not found in system PATH. Please install it (e.g., 'sudo apt install mafft' or via Bioconda) to use the --mafft flag.")
            return

    # ZERO-STRING: Read only names and lengths (very light)
    names, lengths = get_fasta_info(args.fasta)
    
    # Create an Enum for all sequence names
    name_enum = pl.Enum(names)
    
    # Create the ID mapping and lengths DataFrames for Polars
    mapping_df = pl.DataFrame({
        "name": names,
        "id": list(range(len(names)))
    }).with_columns([
        pl.col("name").cast(name_enum),
        pl.col("id").cast(pl.Int32)
    ])
    
    lengths_df = pl.DataFrame({
        'name': names,
        'len': [lengths.get(n, 0) for n in names]
    }).with_columns([
        pl.col('name').cast(name_enum),
        pl.col('len').cast(pl.Int32)
    ])
    
    # Process coverage arg
    cov_val = None
    if args.coverage:
        if args.coverage.lower() == 'dynamic':
            cov_val = 'dynamic'
        else:
            try:
                cov_val = float(args.coverage)
            except ValueError:
                print(f"Error: coverage must be a float or 'dynamic'. Got '{args.coverage}'")
                return

    # Process custom cols
    custom_cols = {}
    if args.col_query is not None: custom_cols['q'] = args.col_query
    if args.col_target is not None: custom_cols['t'] = args.col_target
    if args.col_bitscore is not None: custom_cols['bitscore'] = args.col_bitscore
    if args.col_evalue is not None: custom_cols['evalue'] = args.col_evalue
    if args.col_pident is not None: custom_cols['pident'] = args.col_pident
    if args.col_length is not None: custom_cols['length'] = args.col_length
    
    input_format = args.format
    if custom_cols and input_format is None:
        input_format = 'custom'

    # Call read_hits with mapping and lengths
    u_v, weights, v_list = read_hits(
        args.input, format=input_format, evalue=args.evalue, bitscore=args.bitscore, 
        pident=args.pidentity, coverage=cov_val,
        mapping_df=mapping_df, lengths_df=lengths_df,
        pident_override=args.id_override,
        custom_cols=custom_cols if input_format == 'custom' else None
    )

    # Build clusters
    components = build_clusters(
        (u_v, weights, v_list), 
        min_size=args.min_cluster_size, 
        seed=args.seed,
        resolution=args.resolution
    )

    # Write outputs (Need full sequences for FASTA generation)
    # Lazy sequence loading only for final writing
    from .clustrx import read_fasta
    sequences = read_fasta(args.fasta)
    
    out_clusters = Path(args.outdir) / "clusters"
    out_fastas = Path(args.outdir) / "fasta_files"
    out_alignments = Path(args.outdir) / "alignments" if args.mafft else None
    
    write_clusters(components, sequences, out_clusters, out_fastas, out_alignments)

    print(f"Done. {len(components)} clusters written to {args.outdir}/")

if __name__ == "__main__":
    main()
