"""Slurm script generator for viral-cataloger."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Sequence


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for Slurm generation."""
    parser = argparse.ArgumentParser(
        description="Generate a Slurm batch script for viral-cataloger.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Slurm Options ---
    slurm = parser.add_argument_group("Slurm Configuration")
    slurm.add_argument(
        "--job-name", "-J", default="viral_cat", help="Name of the job."
    )
    slurm.add_argument(
        "--time", default="01:00:00", help="Wall time limit (HH:MM:SS)."
    )
    slurm.add_argument(
        "--mem", default="8G", help="Memory allocation (e.g., 8G, 16000M)."
    )
    slurm.add_argument(
        "--partition", "-p", help="Partition to submit the job to."
    )
    slurm.add_argument(
        "--account", "-A", help="Account/Project name to charge."
    )
    slurm.add_argument(
        "--cpus-per-task",
        type=int,
        default=4,
        help="Number of CPU cores. Also sets --threads for the pipeline.",
    )
    slurm.add_argument(
        "--mail-user", help="Email address for notifications."
    )
    slurm.add_argument(
        "--mail-type",
        default="FAIL",
        help="Notify on state change (BEGIN, END, FAIL, ALL).",
    )
    slurm.add_argument(
        "--sbatch-file",
        default="submit_cataloger.sbatch",
        help="Name of the generated sbatch file.",
    )
    slurm.add_argument(
        "--log-output",
        default="%x_%j.out",
        help="Filename pattern for stdout/stderr (%%x=job_name, %%j=job_id).",
    )

    # --- Pipeline Options ---
    pipeline = parser.add_argument_group("Pipeline Configuration")
    pipeline.add_argument(
        "-i", "--input_dir", required=True, help="Input directory (passed to viral-cataloger)."
    )
    pipeline.add_argument(
        "-o", "--output_dir", required=True, help="Output directory (passed to viral-cataloger)."
    )
    pipeline.add_argument(
        "--min_ani", type=float, default=95.0, help="Minimum ANI (passed to viral-cataloger)."
    )
    pipeline.add_argument(
        "--min_tcov", type=float, default=85.0, help="Minimum coverage (passed to viral-cataloger)."
    )
    pipeline.add_argument(
        "--prefix", default="catalog", help="Output prefix (passed to viral-cataloger)."
    )

    return parser


def generate_sbatch(args: argparse.Namespace) -> str:
    """Construct the content of the sbatch file."""
    lines = ["#!/bin/bash", ""]

    # Slurm Directives
    lines.append(f"#SBATCH --job-name={args.job_name}")
    lines.append(f"#SBATCH --time={args.time}")
    lines.append(f"#SBATCH --mem={args.mem}")
    lines.append(f"#SBATCH --cpus-per-task={args.cpus_per_task}")
    lines.append(f"#SBATCH --output={args.log_output}")

    if args.partition:
        lines.append(f"#SBATCH --partition={args.partition}")
    if args.account:
        lines.append(f"#SBATCH --account={args.account}")
    if args.mail_user:
        lines.append(f"#SBATCH --mail-user={args.mail_user}")
        lines.append(f"#SBATCH --mail-type={args.mail_type}")

    lines.append("")
    lines.append("# Exit on error")
    lines.append("set -e")
    lines.append("")
    lines.append("echo \"Starting viral-cataloger job on $(hostname)\"")
    lines.append("date")
    lines.append("")

    # Construct command
    cmd = [
        "viral-cataloger",
        "--input_dir", f'"{args.input_dir}"',
        "--output_dir", f'"{args.output_dir}"',
        "--threads", str(args.cpus_per_task),
        "--min_ani", str(args.min_ani),
        "--min_tcov", str(args.min_tcov),
        "--prefix", f'"{args.prefix}"',
    ]

    lines.append(" ".join(cmd))
    lines.append("")
    lines.append("echo \"Job complete.\"")
    lines.append("date")

    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    """Run the Slurm generator."""
    parser = build_parser()
    args = parser.parse_args(argv)

    content = generate_sbatch(args)

    output_path = Path(args.sbatch_file)
    try:
        output_path.write_text(content, encoding="utf-8")
        print(f"✅ Generated Slurm script: {output_path}")
        print(f"   Submit with: sbatch {output_path}")
    except IOError as e:
        print(f"❌ Error writing file: {e}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
