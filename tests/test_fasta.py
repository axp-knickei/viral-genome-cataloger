from pathlib import Path

from viral_cataloger.pipeline import format_skani_output, parse_fasta, read_fasta_lengths


def test_parse_fasta_sequences() -> None:
    fasta_path = Path("tests/data/sample.fasta")
    records = list(parse_fasta(fasta_path))
    assert [record[0] for record in records] == [
        "genome_a",
        "genome_b",
        "genome_c",
    ]
    assert records[0][1].startswith("ATGCGTACGTA")


def test_read_fasta_lengths() -> None:
    fasta_path = Path("tests/data/sample.fasta")
    lengths = read_fasta_lengths(fasta_path)
    assert lengths["genome_a"] > lengths["genome_b"] > lengths["genome_c"]


def test_format_skani_output(tmp_path: Path) -> None:
    skani_path = Path("tests/data/skani_results.txt")
    output_path = tmp_path / "formatted.tsv"
    format_skani_output(skani_path, output_path)
    lines = output_path.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0].startswith("genome_a\tgenome_b\t1\t96.50\t90.00\t90.00")
