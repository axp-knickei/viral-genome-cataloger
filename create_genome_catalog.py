#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
create_genome_catalog.py

A bioinformatic pipeline to create a dereplicated viral genome catalog from a
collection of genome assemblies.

This script automates the process of:
1. Aggregating individual FASTA files.
2. Performing all-vs-all Average Nucleotide Identity (ANI) comparison using skani.
3. Reformatting skani's output for clustering.
4. Performing greedy clustering to identify viral Operational Taxonomic Units (vOTUs).
5. Generating a final, non-redundant FASTA file of representative genomes.
"""

import argparse
import os
import subprocess
import sys
import time
import platform
import resource
import glob

def run_command(command, error_message):
    """Executes a shell command and checks for errors."""
    print(f"--> Running: {' '.join(command)}")
    try:
        # Using capture_output=True to hide stdout/stderr unless there's an error
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"🔥🔥🔥 ERROR: {error_message}", file=sys.stderr)
        print(f"Command failed: {' '.join(command)}", file=sys.stderr)
        print(f"Stderr: {e.stderr}", file=sys.stderr)
        sys.exit(1)

def perform_clustering(params):
    """
    Reads ANI results and performs greedy clustering to group sequences into species.
    This logic is adapted from the original notebook's clustering cell.
    """
    def parse_seqs(path):
        """Helper function to parse a FASTA file."""
        with open(path, 'r') as handle:
            try:
                id = next(handle).split()[0][1:]
            except StopIteration:
                return
            seq = ''
            for line in handle:
                if line.startswith('>'):
                    yield id, seq
                    id = line.split()[0][1:]
                    seq = ''
                else:
                    seq += line.rstrip()
            yield id, seq

    def max_mem_usage():
        """Calculates peak memory usage. Note: Works on Linux/macOS."""
        if platform.system() == 'Linux':
            return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 2)
        else: # macOS uses bytes, Linux uses kilobytes
            return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2)

    start_time = time.time()

    # 1. Read sequences and sort by length (descending)
    print("--> Reading and sorting sequences by length...")
    seqs_dict = {id: len(seq) for id, seq in parse_seqs(params['fna'])}
    seqs_sorted = [x[0] for x in sorted(seqs_dict.items(), key=lambda item: item[1], reverse=True)]
    print(f"    {len(seqs_sorted)} sequences loaded for clustering.")

    # 2. Store graph edges from the ANI file
    print("\n--> Building similarity graph from ANI results...")
    edges = {x: [] for x in seqs_sorted}
    num_edges = 0
    with open(params['ani'], 'r') as handle:
        for line in handle:
            try:
                qname, tname, _, ani, qcov, tcov = line.split()
                if qname == tname or qname not in edges or tname not in edges:
                    continue
                if float(ani) >= params['min_ani'] and float(qcov) >= params['min_qcov'] and float(tcov) >= params['min_tcov']:
                     edges[qname].append(tname)
                     edges[tname].append(qname) # Make the graph undirected
                     num_edges += 1
            except ValueError:
                continue
    print(f"    {num_edges} edges retained that meet thresholds.")

    # 3. Perform greedy clustering
    print("\n--> Performing greedy clustering...")
    clust_to_seqs = {}
    seq_to_clust = {}
    for seq_id in seqs_sorted:
        if seq_id not in seq_to_clust:
            clust_to_seqs[seq_id] = [seq_id]
            seq_to_clust[seq_id] = seq_id
            # Check neighbors and neighbors of neighbors (for better clustering)
            to_visit = list(edges[seq_id])
            visited = {seq_id}
            while to_visit:
                mem_id = to_visit.pop(0)
                if mem_id not in visited:
                    visited.add(mem_id)
                    if mem_id not in seq_to_clust:
                        clust_to_seqs[seq_id].append(mem_id)
                        seq_to_clust[mem_id] = seq_id
                        # Add its neighbors to the list to visit
                        to_visit.extend(edges[mem_id])

    print(f"    Found {len(clust_to_seqs)} total clusters.")
    print(f"    (Time elapsed: {time.time() - start_time:.2f} seconds, Peak RAM: {max_mem_usage()} GB)")

    # 4. Write output file
    print("\n--> Writing cluster results...")
    with open(params['out'], 'w') as out_handle:
        out_handle.write("representative\tmembers\n")
        for seq_id, mem_ids in clust_to_seqs.items():
            out_handle.write(f"{seq_id}\t{','.join(mem_ids)}\n")
    print(f"✅ Clustering complete. Results saved to '{params['out']}'.")


def main():
    """Main function to orchestrate the pipeline."""
    parser = argparse.ArgumentParser(
        description="A pipeline to create a dereplicated viral genome catalog.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Directory containing input FASTA files (e.g., *.fa, *.fasta)."
    )
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="Directory to store all output files."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=4,
        help="Number of threads to use for skani (default: 4)."
    )
    parser.add_argument(
        "--min_ani",
        type=float,
        default=95.0,
        help="Minimum ANI percentage for clustering (default: 95.0)."
    )
    parser.add_argument(
        "--min_tcov",
        type=float,
        default=85.0,
        help="Minimum target coverage percentage for clustering (default: 85.0)."
    )
    parser.add_argument(
        "--prefix",
        default="catalog",
        help="Prefix for output file names (default: 'catalog')."
    )

    args = parser.parse_args()

    # --- 1. Setup ---
    print("🚀 Starting Genome Catalog Creation Pipeline 🚀")
    os.makedirs(args.output_dir, exist_ok=True)

    # --- 2. Aggregate FASTA files ---
    print("\n--- Step 1: Aggregating FASTA files ---")
    all_fasta_path = os.path.join(args.output_dir, "all_genomes.fa")
    fasta_files = glob.glob(os.path.join(args.input_dir, '*.fa')) + \
                  glob.glob(os.path.join(args.input_dir, '*.fasta')) + \
                  glob.glob(os.path.join(args.input_dir, '*.fna'))

    if not fasta_files:
        print(f"🔥🔥🔥 ERROR: No FASTA files (*.fa, *.fasta, *.fna) found in '{args.input_dir}'.", file=sys.stderr)
        sys.exit(1)

    with open(all_fasta_path, 'w') as outfile:
        for fname in fasta_files:
            with open(fname, 'r') as infile:
                outfile.write(infile.read())
    print(f"✅ Combined {len(fasta_files)} files into '{all_fasta_path}'.")

    # --- 3. Run skani ---
    print("\n--- Step 2: Calculating all-vs-all ANI with skani ---")
    skani_output_path = os.path.join(args.output_dir, "skani_results.txt")
    skani_command = [
        "skani", "triangle", all_fasta_path,
        "-o", skani_output_path,
        "-t", str(args.threads),
        "-i",        # Treat each contig as an individual genome
        "-m", "200", # Recommended for small genomes
        "--slow",    # More accurate alignment fraction (AF)
        "-E",        # Output as an edge list
        "--faster-small", # Speed-up for large datasets
        "-s", "90"   # Screen for approximate ANI >90%
    ]
    run_command(skani_command, "skani execution failed.")
    print(f"✅ skani results saved to '{skani_output_path}'.")

    # --- 4. Reformat skani Output ---
    print("\n--- Step 3: Reformatting skani output for clustering ---")
    reformatted_ani_path = os.path.join(args.output_dir, "ani_formatted.txt")
    awk_command = [
        "awk",
        'BEGIN{OFS="\t"} NR>1{print $1, $2, "1", $3, $4*100, $5*100}',
        skani_output_path
    ]
    reformatted_output = run_command(awk_command, "awk reformatting failed.")
    with open(reformatted_ani_path, 'w') as f:
        f.write(reformatted_output)
    print(f"✅ Formatted ANI data saved to '{reformatted_ani_path}'.")

    # --- 5. Cluster Sequences ---
    print("\n--- Step 4: Clustering genomes into vOTUs ---")
    cluster_path = os.path.join(args.output_dir, f"{args.prefix}_clusters.tsv")
    clustering_params = {
        'fna': all_fasta_path,
        'ani': reformatted_ani_path,
        'out': cluster_path,
        'min_ani': args.min_ani,
        'min_qcov': 0.0, # As per original script
        'min_tcov': args.min_tcov,
    }
    perform_clustering(clustering_params)

    # --- 6. Create Final Catalog ---
    print("\n--- Step 5: Generating final representative catalog ---")
    representative_ids_path = os.path.join(args.output_dir, "representative_ids.txt")
    # Use NR>1 to skip header, then print first column
    awk_command_extract = ["awk", 'NR>1{print $1}', cluster_path]
    ids_output = run_command(awk_command_extract, "Extracting representative IDs failed.")
    with open(representative_ids_path, 'w') as f:
        f.write(ids_output)

    final_fasta_path = os.path.join(args.output_dir, f"{args.prefix}_vOTU_catalog.fasta")
    seqkit_command = [
        "seqkit", "grep",
        "-f", representative_ids_path,
        all_fasta_path,
        "-o", final_fasta_path
    ]
    run_command(seqkit_command, "seqkit grep failed.")
    print(f"✅ Final dereplicated catalog saved to '{final_fasta_path}'.")

    # --- 7. Final Verification ---
    final_count = run_command(["grep", "-c", ">", final_fasta_path], "Counting sequences failed.")
    print(f"\nTotal representative sequences in the final catalog: {final_count.strip()}")
    print("\n🎉 Workflow Complete! 🎉")


if __name__ == "__main__":
    main()