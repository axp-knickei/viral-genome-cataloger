"""Command-line entry point for viral-cataloger."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from viral_cataloger.pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""

    parser = argparse.ArgumentParser(
        description="A pipeline to create a dereplicated viral genome catalog.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_dir",
        required=True,
        help="Directory containing input FASTA files (e.g., *.fa, *.fasta).",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        required=True,
        help="Directory to store all output files.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use for skani (default: 4).",
    )
    parser.add_argument(
        "--min_ani",
        type=float,
        default=95.0,
        help="Minimum ANI percentage for clustering (default: 95.0).",
    )
    parser.add_argument(
        "--min_tcov",
        type=float,
        default=85.0,
        help="Minimum target coverage percentage for clustering (default: 85.0).",
    )
    parser.add_argument(
        "--prefix",
        default="catalog",
        help="Prefix for output file names (default: 'catalog').",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the viral catalog pipeline."""

    parser = build_parser()
    args = parser.parse_args(argv)
    run_pipeline(
        input_dir=Path(args.input_dir),
        output_dir=Path(args.output_dir),
        threads=args.threads,
        min_ani=args.min_ani,
        min_tcov=args.min_tcov,
        prefix=args.prefix,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
