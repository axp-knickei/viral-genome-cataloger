# Viral Genome Cataloging Pipeline (voter-cataloger)

[![Python Version](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

The **Viral Genome Cataloging Pipeline** is a specialized bioinformatic tool designed to dereplicate viral genome sequences and establish non-redundant catalogs of **viral Operational Taxonomic Units (vOTUs)**. 

In viral metagenomics, researchers often recover thousands of sequences that may represent the same biological entity. This pipeline automates the process of grouping these sequences into species-level clusters using high-resolution **Average Nucleotide Identity (ANI)** and a deterministic **Greedy Representative (Star Topology)** clustering algorithm.

### Key Features
*   **High Performance:** Leverages `skani` for ultra-fast ANI calculations.
*   **Scientific Rigor:** Implements star-topology clustering to prevent "chaining" artifacts, ensuring each vOTU is defined by direct similarity to a high-quality representative.
*   **Accessibility:** Designed for researchers who prefer a straightforward command-line interface with minimal setup.
*   **Portability:** Refactored for modern Python environments with comprehensive dependency checking.

---

## 🔬 Scientific Methodology

### 1. Species Demarcation (vOTUs)
The pipeline defaults to a **95% ANI** over **85% Alignment Fraction (Coverage)** threshold. This is the widely accepted international standard for demarcating viral species (vOTUs) in metagenomic datasets (see [Roux et al., 2019](https://doi.org/10.1038/s41564-019-0357-z)).

### 2. Clustering Strategy
Unlike transitive clustering (which can merge distinct species via intermediate "hybrid" sequences), this pipeline uses a **Greedy Representative** approach:
1.  Genomes are sorted by length (longest to shortest).
2.  The longest genome becomes the "Seed" or **Representative** for the first vOTU.
3.  Any other genome matching the ANI and coverage thresholds *directly* against this seed is assigned to its cluster.
4.  The process repeats for the remaining unassigned sequences.

This ensures that every member of a vOTU is directly similar to the longest (and presumably most complete) genome in that group.

---

## ⚙️ Prerequisites & Dependencies

The pipeline requires two external bioinformatic tools to be installed on your system. 

### 1. Required External Tools
| Tool | Purpose | Source |
| :--- | :--- | :--- |
| **skani** | Fast ANI calculation | [GitHub](https://github.com/bluenote-1577/skani) |
| **SeqKit** | FASTA manipulation | [Bioinf.shenwei.me](https://bioinf.shenwei.me/seqkit/) |

#### Quick Installation (Linux/macOS)
If you have `conda` or `mamba` installed, we recommend:
```bash
conda install -c bioconda skani seqkit
```
Alternatively, you can download the binaries manually and place them in your `PATH` (e.g., `/usr/local/bin/`).

---

## 📦 Installation Guide

1. **Clone the Repository**
   ```bash
   git clone https://github.com/axp-knickei/viral-genome-cataloger.git
   cd viral-genome-cataloger
   ```

2. **Install the Pipeline**
   We recommend using a virtual environment. You can install the package using `pip` or `uv`:
   ```bash
   # Using standard pip
   pip install .

   # Or using uv (recommended for speed)
   uv pip install .
   ```

---

## 🚀 Execution Guide

### Preparing Data
Place all your viral genome sequences (FASTA format) into a single directory. The files should end in `.fa`, `.fna`, or `.fasta`.

### Running the Pipeline
Once installed, you can run the pipeline using the `viral-cataloger` command:

```bash
viral-cataloger \
  --input_dir ./my_genomes/ \
  --output_dir ./vOTU_results/ \
  --threads 8 \
  --min_ani 95.0 \
  --min_tcov 85.0
```

### HPC Job Submission (Slurm)
For users on High-Performance Computing (HPC) clusters using Slurm, a helper utility `viral-cataloger-slurm` generates a submission script automatically:

```bash
viral-cataloger-slurm \
  --job-name my_vOTU_job \
  --time 02:00:00 \
  --cpus-per-task 16 \
  --mem 32G \
  --input_dir ./my_genomes/ \
  --output_dir ./vOTU_results/
```

This command creates a file `submit_cataloger.sbatch` which can be submitted via:
```bash
sbatch submit_cataloger.sbatch
```

### Command-Line Arguments
| Argument | Shorthand | Default | Description |
| :--- | :--- | :--- | :--- |
| `--input_dir` | `-i` | (Required) | Directory containing input FASTA files. |
| `--output_dir` | `-o` | (Required) | Directory to store results. |
| `--threads` | `-t` | `4` | CPU cores to use for `skani`. |
| `--min_ani` | - | `95.0` | Minimum % ANI for clustering. |
| `--min_tcov` | - | `85.0` | Minimum % alignment coverage for the target genome. |
| `--prefix` | - | `catalog` | Prefix for the output filenames. |

---

## 💻 Hardware Requirements & Performance

The pipeline is highly optimized and can be executed on most modern laptops or high-performance computing (HPC) nodes. Below are estimated requirements based on the number of input sequences (assuming an average viral genome size of ~50kb).

| Number of Sequences | Recommended CPU Cores | Recommended RAM | Estimated Run Time* |
| :--- | :--- | :--- | :--- |
| **< 1,000** | 2 - 4 cores | 4 GB | < 2 minutes |
| **1,000 - 10,000** | 8 - 16 cores | 8 - 16 GB | 5 - 15 minutes |
| **10,000 - 50,000** | 32+ cores | 32 - 64 GB | 30 - 60 minutes |

*\*Run times are estimates based on standard SSD storage and reflect the combined duration of ANI calculation and clustering.*

### Performance Notes
*   **CPU Scaling:** `skani` (the engine behind the ANI calculation) scales linearly with thread count. Increasing the `--threads` parameter is the most effective way to reduce run time.
*   **Memory Efficiency:** The pipeline is designed to be memory-efficient. Even for large datasets, the peak RAM usage is primarily determined by the number of unique sequences and the density of the similarity graph.

---

## 📁 Interpreting Results

The pipeline generates the following files in your specified `--output_dir`:

| File | Description |
| :--- | :--- |
| **`catalog_vOTU_catalog.fasta`** | **Primary Output.** The dereplicated catalog containing one representative sequence per vOTU. |
| `catalog_clusters.tsv` | A table mapping each representative to its cluster members. |
| `skani_results.txt` | The raw all-vs-all ANI comparison data. |
| `all_genomes.fa` | A concatenated file of all input sequences used for the analysis. |

---

## 📜 Citations

If you use this pipeline in your research, please cite the underlying tools:

*   **skani:** Ondov, B.D., et al. (2023). Fast and accurate average nucleotide identity estimation with skani. *Nat Biotechnol*. [doi:10.1038/s41587-023-01774-x](https://doi.org/10.1038/s41587-023-01774-x)
*   **SeqKit:** Shen, W., et al. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. *PLoS ONE*. [doi:10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962)
*   **vOTU Standards:** Roux, S., et al. (2019). Minimum Information about an Uncultivated Virus Genome (MIUViG). *Nat Biotechnol*. [doi:10.1038/s41564-019-0357-z](https://doi.org/10.1038/s41564-019-0357-z)

---

## License
Distributed under the MIT License. See `LICENSE` for more information.