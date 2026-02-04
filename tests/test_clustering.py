from pathlib import Path

from viral_cataloger.pipeline import greedy_star_clustering, load_ani_edges, read_fasta_lengths


def test_greedy_star_clustering_avoids_chaining() -> None:
    fasta_path = Path("tests/data/sample.fasta")
    ani_path = Path("tests/data/ani_formatted.txt")
    lengths = read_fasta_lengths(fasta_path)
    edges = load_ani_edges(
        ani_path,
        min_ani=95.0,
        min_qcov=85.0,
        min_tcov=85.0,
        valid_ids=lengths.keys(),
    )
    clusters = greedy_star_clustering(lengths, edges)

    assert "genome_a" in clusters
    assert clusters["genome_a"] == ["genome_a", "genome_b"]
    assert clusters["genome_c"] == ["genome_c"]
