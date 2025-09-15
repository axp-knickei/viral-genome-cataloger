# Viral Genome Cataloging Pipeline

A simple and reproducible pipeline to dereplicate a collection of viral genome sequences and create a non-redundant genome catalog of viral Operational Taxonomic Units (vOTUs).

This workflow uses `skani` for fast and accurate Average Nucleotide Identity (ANI) calculations and a greedy clustering algorithm to group genomes based on a species-level threshold (typically ≥95% ANI).

---

## ⚙️ Dependencies

Before you begin, you must have the following command-line tools installed and available in your system's `PATH`.

1.  **Python 3**:
    * The script is written in Python 3. It uses standard libraries, so no special packages are needed.

2.  **skani**: For calculating ANI.
    * **Installation (Linux):**
        ```bash
        wget [https://github.com/bluenote-1577/skani/releases/download/latest/skani](https://github.com/bluenote-1577/skani/releases/download/latest/skani) -O skani
        chmod +x skani
        sudo mv skani /usr/local/bin/
        ```

3.  **SeqKit**: A toolkit for FASTA/Q file manipulation.
    * **Installation (Linux):**
        ```bash
        wget [https://github.com/shenwei356/seqkit/releases/download/v2.8.0/seqkit_linux_amd64.tar.gz](https://github.com/shenwei356/seqkit/releases/download/v2.8.0/seqkit_linux_amd64.tar.gz)
        tar -xzf seqkit_linux_amd64.tar.gz
        chmod +x seqkit
        sudo mv seqkit /usr/local/bin/
        ```

---

## 🚀 Usage

1.  **Clone the Repository**
    ```bash
    git clone <your-repository-url>
    cd <your-repository-name>
    ```

2.  **Prepare Your Input Data**
    * Create a directory (e.g., `input_genomes`) and place all your individual viral genome FASTA files inside it. The files should have a `.fa`, `.fna`, or `.fasta` extension.
    ```
    .
    ├── create_genome_catalog.py
    ├── README.md
    └── input_genomes/
        ├── virus1.fa
        ├── virus2.fa
        └── ...
    ```

3.  **Run the Pipeline**
    * Execute the script, providing the path to your input directory and a desired output directory.

    * **Basic Example:**
        ```bash
        python create_genome_catalog.py \
          -i input_genomes/ \
          -o results/ \
          -t 8
        ```
    * This command will process all genomes in `input_genomes/`, use 8 threads for `skani`, and save all results in a new `results/` directory.

---

### Command-Line Arguments

| Argument | Shorthand | Description | Default |
|---|---|---|---|
| `--input_dir` | `-i` | **(Required)** Directory containing input FASTA files. | |
| `--output_dir`| `-o` | **(Required)** Directory to store all output files. | |
| `--threads` | `-t` | Number of threads to use for `skani`. | `4` |
| `--min_ani` | | Minimum ANI percentage for clustering. | `95.0` |
| `--min_tcov` | | Minimum target coverage percentage for clustering. | `85.0` |
| `--prefix` | | Prefix for the final output files. | `catalog`|

---

## 📁 Output Files

After a successful run, the specified output directory (e.g., `results/`) will contain:

* `all_genomes.fa`: A single FASTA file with all input genomes combined.
* `skani_results.txt`: The raw all-vs-all comparison output from `skani`.
* `ani_formatted.txt`: The reformatted `skani` output, ready for clustering.
* **`catalog_clusters.tsv`**: The main clustering result. A tab-separated file with two columns: `representative` and `members` (a comma-separated list of all genome IDs in that cluster).
* **`catalog_vOTU_catalog.fasta`**: The final, dereplicated genome catalog containing only the representative sequences for each vOTU.

This final FASTA file is your non-redundant set of genomes, ready for downstream analyses.