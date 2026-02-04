"""Core pipeline logic for viral genome catalog creation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

resource = None
if sys.platform != "win32":  # pragma: no cover - Windows
    import resource  # type: ignore  # noqa: F401


@dataclass(frozen=True)
class ClusteringParams:
    """Parameters for greedy star clustering."""

    fna: Path
    ani: Path
    out: Path
    min_ani: float
    min_qcov: float
    min_tcov: float


def check_external_dependencies() -> None:
    """Ensure required external tools are available."""

    missing = [tool for tool in ("skani", "seqkit") if shutil.which(tool) is None]
    if missing:
        hint = (
            "Install skani from https://github.com/bluenote-1577/skani and "
            "seqkit from https://bioinf.shenwei.me/seqkit/. "
            "Ensure they are on your PATH."
        )
        raise RuntimeError(
            f"Missing external dependencies: {', '.join(missing)}. {hint}"
        )


def run_command(
    command: Sequence[str],
    error_message: str,
    stdout_path: Optional[Path] = None,
) -> None:
    """Execute a command while streaming output to file when requested."""

    print(f"--> Running: {' '.join(command)}")
    try:
        if stdout_path:
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            with stdout_path.open("w", encoding="utf-8") as handle:
                subprocess.run(command, check=True, stdout=handle)
        else:
            subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"{error_message} Command failed: {' '.join(command)}"
        ) from exc


def parse_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (sequence_id, sequence) tuples from a FASTA file."""

    seq_id: Optional[str] = None
    seq_chunks: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(seq_chunks)
                seq_id = line.split()[0][1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if seq_id is not None:
            yield seq_id, "".join(seq_chunks)


def read_fasta_lengths(path: Path) -> Dict[str, int]:
    """Return a mapping of sequence ID to sequence length."""

    return {seq_id: len(seq) for seq_id, seq in parse_fasta(path)}


def max_mem_usage_gb() -> Optional[float]:
    """Return peak memory usage in GB when supported."""

    if resource is None:
        return None
    if platform.system() == "Linux":
        return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 2)
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2)


def aggregate_fastas(input_dir: Path, output_path: Path) -> List[Path]:
    """Combine FASTA files from input_dir into output_path."""

    patterns = ("*.fa", "*.fasta", "*.fna")
    fasta_files = sorted(
        {path for pattern in patterns for path in input_dir.glob(pattern)}
    )
    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA files (*.fa, *.fasta, *.fna) found in '{input_dir}'."
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as outfile:
        for fasta_path in fasta_files:
            outfile.write(fasta_path.read_text(encoding="utf-8"))
    return fasta_files


def format_skani_output(input_path: Path, output_path: Path) -> None:
    """Reformat skani edge list output into clustering-ready columns."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with input_path.open("r", encoding="utf-8") as infile, output_path.open(
        "w", encoding="utf-8"
    ) as outfile:
        header = infile.readline()
        if not header:
            return
        for line in infile:
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            qname, tname, ani, qcov, tcov = parts[:5]
            try:
                ani_val = float(ani)
                qcov_val = float(qcov) * 100
                tcov_val = float(tcov) * 100
            except ValueError:
                continue
            outfile.write(
                f"{qname}\t{tname}\t1\t{ani_val:.2f}\t{qcov_val:.2f}\t{tcov_val:.2f}\n"
            )


def load_ani_edges(
    path: Path,
    min_ani: float,
    min_qcov: float,
    min_tcov: float,
    valid_ids: Iterable[str],
) -> Dict[str, List[str]]:
    """Load ANI edges that meet thresholds."""

    valid_set = set(valid_ids)
    edges: Dict[str, List[str]] = {seq_id: [] for seq_id in valid_set}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 6:
                continue
            qname, tname, _, ani, qcov, tcov = parts[:6]
            if qname == tname or qname not in valid_set or tname not in valid_set:
                continue
            try:
                if (
                    float(ani) >= min_ani
                    and float(qcov) >= min_qcov
                    and float(tcov) >= min_tcov
                ):
                    edges[qname].append(tname)
                    edges[tname].append(qname)
            except ValueError:
                continue
    return edges


def greedy_star_clustering(
    seq_lengths: Dict[str, int],
    edges: Dict[str, List[str]],
) -> Dict[str, List[str]]:
    """Cluster sequences using greedy representative (star topology).

    Assumes edges only represent direct ANI + coverage matches. Clusters are
    built by choosing the longest unassigned genome as representative and
    assigning only genomes directly connected to that representative.
    """

    seqs_sorted = sorted(
        seq_lengths.keys(), key=lambda seq: (-seq_lengths[seq], seq)
    )
    unassigned = set(seqs_sorted)
    clusters: Dict[str, List[str]] = {}

    for seq_id in seqs_sorted:
        if seq_id not in unassigned:
            continue
        unassigned.remove(seq_id)
        neighbor_list = sorted(
            (neighbor for neighbor in edges.get(seq_id, []) if neighbor in unassigned),
            key=lambda neighbor: (-seq_lengths[neighbor], neighbor),
        )
        cluster_members = [seq_id]
        for neighbor in neighbor_list:
            if neighbor in unassigned:
                unassigned.remove(neighbor)
                cluster_members.append(neighbor)
        clusters[seq_id] = cluster_members
    return clusters


def write_clusters(path: Path, clusters: Dict[str, List[str]]) -> None:
    """Write cluster assignments to TSV file."""

    with path.open("w", encoding="utf-8") as out_handle:
        out_handle.write("representative\tmembers\n")
        for rep_id, members in clusters.items():
            out_handle.write(f"{rep_id}\t{','.join(members)}\n")


def write_representative_ids(path: Path, clusters: Dict[str, List[str]]) -> None:
    """Write representative IDs to file."""

    with path.open("w", encoding="utf-8") as handle:
        for rep_id in clusters:
            handle.write(f"{rep_id}\n")


def count_fasta_records(path: Path) -> int:
    """Count FASTA records in a file."""

    return sum(1 for _ in parse_fasta(path))


def perform_clustering(params: ClusteringParams) -> Dict[str, List[str]]:
    """Perform greedy star clustering and write output."""

    start_time = time.time()
    print("--> Reading and sorting sequences by length...")
    seq_lengths = read_fasta_lengths(params.fna)
    print(f"    {len(seq_lengths)} sequences loaded for clustering.")

    print("\n--> Building similarity graph from ANI results...")
    edges = load_ani_edges(
        params.ani,
        min_ani=params.min_ani,
        min_qcov=params.min_qcov,
        min_tcov=params.min_tcov,
        valid_ids=seq_lengths.keys(),
    )
    print(
        f"    {sum(len(neigh) for neigh in edges.values()) // 2} edges retained "
        "that meet thresholds."
    )

    print("\n--> Performing greedy clustering...")
    clusters = greedy_star_clustering(seq_lengths, edges)
    elapsed = time.time() - start_time
    peak_mem = max_mem_usage_gb()
    if peak_mem is None:
        print(f"    Found {len(clusters)} total clusters. (Time elapsed: {elapsed:.2f} seconds)")
    else:
        print(
            f"    Found {len(clusters)} total clusters. (Time elapsed: {elapsed:.2f} seconds, "
            f"Peak RAM: {peak_mem} GB)"
        )

    print("\n--> Writing cluster results...")
    write_clusters(params.out, clusters)
    print(f"✅ Clustering complete. Results saved to '{params.out}'.")
    return clusters


def run_pipeline(
    input_dir: Path,
    output_dir: Path,
    threads: int,
    min_ani: float,
    min_tcov: float,
    prefix: str,
) -> None:
    """Run the full genome catalog pipeline."""

    check_external_dependencies()
    print("🚀 Starting Genome Catalog Creation Pipeline 🚀")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n--- Step 1: Aggregating FASTA files ---")
    all_fasta_path = output_dir / "all_genomes.fa"
    fasta_files = aggregate_fastas(input_dir, all_fasta_path)
    print(f"✅ Combined {len(fasta_files)} files into '{all_fasta_path}'.")

    print("\n--- Step 2: Calculating all-vs-all ANI with skani ---")
    skani_output_path = output_dir / "skani_results.txt"
    skani_command = [
        "skani",
        "triangle",
        str(all_fasta_path),
        "-o",
        str(skani_output_path),
        "-t",
        str(threads),
        "-i",
        "-m",
        "200",
        "--slow",
        "-E",
        "--faster-small",
        "-s",
        "90",
    ]
    run_command(skani_command, "skani execution failed.")
    print(f"✅ skani results saved to '{skani_output_path}'.")

    print("\n--- Step 3: Reformatting skani output for clustering ---")
    reformatted_ani_path = output_dir / "ani_formatted.txt"
    format_skani_output(skani_output_path, reformatted_ani_path)
    print(f"✅ Formatted ANI data saved to '{reformatted_ani_path}'.")

    print("\n--- Step 4: Clustering genomes into vOTUs ---")
    cluster_path = output_dir / f"{prefix}_clusters.tsv"
    clustering_params = ClusteringParams(
        fna=all_fasta_path,
        ani=reformatted_ani_path,
        out=cluster_path,
        min_ani=min_ani,
        min_qcov=0.0,
        min_tcov=min_tcov,
    )
    clusters = perform_clustering(clustering_params)

    print("\n--- Step 5: Generating final representative catalog ---")
    representative_ids_path = output_dir / "representative_ids.txt"
    write_representative_ids(representative_ids_path, clusters)

    final_fasta_path = output_dir / f"{prefix}_vOTU_catalog.fasta"
    seqkit_command = [
        "seqkit",
        "grep",
        "-f",
        str(representative_ids_path),
        str(all_fasta_path),
        "-o",
        str(final_fasta_path),
    ]
    run_command(seqkit_command, "seqkit grep failed.")
    print(f"✅ Final dereplicated catalog saved to '{final_fasta_path}'.")

    final_count = count_fasta_records(final_fasta_path)
    print(f"\nTotal representative sequences in the final catalog: {final_count}")
    print("\n🎉 Workflow Complete! 🎉")
