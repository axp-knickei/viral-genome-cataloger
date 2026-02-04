# Viral Genome Cataloging Pipeline

A reproducible pipeline to dereplicate a collection of viral genome sequences and create a non-redundant genome catalog of viral Operational Taxonomic Units (vOTUs). The workflow uses `skani` for fast Average Nucleotide Identity (ANI) calculations and a greedy representative (star topology) clustering strategy to define vOTUs.

---

## ⚙️ Dependencies

The pipeline requires the following command-line tools to be installed and available on your `PATH`:

1. **Python 3.9+**
2. **skani**: For calculating ANI.
   * Installation (Linux):
     ```bash
     wget https://github.com/bluenote-1577/skani/releases/download/latest/skani -O skani
     chmod +x skani
     sudo mv skani /usr/local/bin/
     ```
3. **SeqKit**: A toolkit for FASTA/Q file manipulation.
   * Installation (Linux):
     ```bash
     wget https://github.com/shenwei356/seqkit/releases/download/v2.8.0/seqkit_linux_amd64.tar.gz
     tar -xzf seqkit_linux_amd64.tar.gz
     chmod +x seqkit
     sudo mv seqkit /usr/local/bin/
     ```

---

## 📦 Installation

Clone the repository and install the package locally:

```bash
git clone <your-repository-url>
cd <your-repository-name>
pip install .
```

---

## 🚀 Usage

Prepare a directory of FASTA files (`.fa`, `.fna`, `.fasta`) and run the CLI:

```bash
viral-cataloger \
  -i input_genomes/ \
  -o results/ \
  -t 8
```

You can also invoke the legacy entry point:

```bash
python create_genome_catalog.py -i input_genomes/ -o results/
```

### Command-Line Arguments

| Argument | Shorthand | Description | Default |
|---|---|---|---|
| `--input_dir` | `-i` | **(Required)** Directory containing input FASTA files. | |
| `--output_dir`| `-o` | **(Required)** Directory to store all output files. | |
| `--threads` | `-t` | Number of threads to use for `skani`. | `4` |
| `--min_ani` | | Minimum ANI percentage for clustering. | `95.0` |
| `--min_tcov` | | Minimum target coverage percentage for clustering. | `85.0` |
| `--prefix` | | Prefix for the final output files. | `catalog` |

---

## 📁 Output Files

After a successful run, the output directory (e.g., `results/`) will contain:

* `all_genomes.fa`: A single FASTA file with all input genomes combined.
* `skani_results.txt`: The raw all-vs-all comparison output from `skani`.
* `ani_formatted.txt`: The reformatted `skani` output, ready for clustering.
* `catalog_clusters.tsv`: A tab-separated file with two columns: `representative` and `members` (comma-separated list of IDs).
* `catalog_vOTU_catalog.fasta`: The final dereplicated genome catalog containing only the representative sequences for each vOTU.

---

## 🔬 Clustering Methodology (Greedy Representative / Star Topology)

The clustering algorithm is deterministic and avoids transitive chaining by using a greedy representative (star topology) strategy:

1. Sort genomes by descending length (ties broken by ID).
2. Select the longest unassigned genome as the cluster representative.
3. Assign only genomes that directly meet ANI and coverage thresholds relative to that representative.
4. Repeat until all genomes are assigned.

This approach ensures cluster assignments are stable across runs and reflect direct similarity to the representative genome rather than indirect links through intermediate genomes.
