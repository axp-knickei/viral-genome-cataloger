import sys
from unittest.mock import patch
from viral_cataloger.slurm import build_parser, generate_sbatch, main

def test_generate_sbatch_content():
    parser = build_parser()
    args = parser.parse_args([
        "--input_dir", "data/in",
        "--output_dir", "data/out",
        "--time", "02:00:00",
        "--cpus-per-task", "16",
        "--job-name", "test_job"
    ])

    content = generate_sbatch(args)

    assert "#SBATCH --job-name=test_job" in content
    assert "#SBATCH --time=02:00:00" in content
    assert "#SBATCH --cpus-per-task=16" in content
    assert 'viral-cataloger' in content
    assert '--threads 16' in content
    assert '--input_dir "data/in"' in content

def test_slurm_main(tmp_path):
    sbatch_file = tmp_path / "test.sbatch"
    argv = [
        "--input_dir", "in",
        "--output_dir", "out",
        "--sbatch-file", str(sbatch_file)
    ]

    # We pass argv directly to main so sys.argv patch isn't strictly necessary
    # but good for safety if main used sys.argv defaults.
    # The main() signature is main(argv: Sequence[str] | None = None).
    ret = main(argv)
    assert ret == 0

    assert sbatch_file.exists()
    content = sbatch_file.read_text(encoding="utf-8")
    assert "#SBATCH" in content
    assert '--input_dir "in"' in content
