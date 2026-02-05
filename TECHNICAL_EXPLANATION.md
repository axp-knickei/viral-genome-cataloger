# Technical Explanation of `viral-genome-cataloger`

## 1. Purpose
The `viral-genome-cataloger` is a specialized bioinformatic pipeline designed to **dereplicate** viral genome assemblies. Its primary goal is to transform a redundant collection of viral sequences (often from metagenomics) into a **non-redundant catalog of Viral Operational Taxonomic Units (vOTUs)**.

It implements widely accepted community standards for viral species demarcation (typically 95% Average Nucleotide Identity (ANI) over 85% coverage).

## 2. Architecture
The system is designed as a lightweight Python wrapper that orchestrates high-performance external binaries.

*   **Language:** Python 3.9+ (Standard Library only, ensuring high portability).
*   **Core Engine:**
    *   **Orchestration:** `src/viral_cataloger/pipeline.py` manages the workflow.
    *   **ANI Calculation:** **skani** (written in Rust) is used for ultra-fast, alignment-free ANI estimation.
    *   **Sequence Manipulation:** **SeqKit** (written in Go) is used for efficient FASTA parsing and extraction.
*   **Design Philosophy:** "Minimal Dependencies, Maximum Speed." By avoiding heavy Python libraries like `pandas` or `numpy` and relying on optimized binaries, the pipeline aims to handle tens of thousands of genomes with minimal memory overhead.

## 3. Workflow & Data Flow

The pipeline executes a linear 5-step process:

### Step 1: Input Aggregation
*   **Input:** Directory of individual FASTA files (`*.fa`, `*.fasta`, `*.fna`).
*   **Operation:** Concatenates all valid input files into a single master file (`all_genomes.fa`).
*   **Rationale:** Batching sequences allows `skani` to build a single index, significantly optimizing the all-vs-all comparison step.

### Step 2: High-Resolution ANI Calculation
*   **Tool:** `skani triangle`
*   **Key Parameters:**
    *   `-m 200`: Sets marker size to 200bp (optimized for viral genomes, which are smaller than bacteria).
    *   `--slow`: Forces rigorous alignment-based Average Fraction (AF) calculation, essential for accurate viral taxonomy.
    *   `-s 90`: Optimizes speed by only reporting pairs with >90% identity (as the clustering threshold is usually 95%).
*   **Output:** A raw edge list (`skani_results.txt`) containing ANI, Query Coverage, and Reference Coverage for pairs.

### Step 3: Data Parsing & Normalization
*   **Operation:** Python parses the raw `skani` output.
*   **Logic:**
    *   Converts fractional values (0.95) to percentages (95.0).
    *   Filters out self-matches.
    *   Formats data into a standardized tab-delimited edge list (`ani_formatted.txt`).

### Step 4: Greedy Star-Topology Clustering (The Core Algorithm)
This is the critical scientific step. The pipeline uses a deterministic **Greedy Representative** algorithm to define vOTUs.

**Algorithm Logic:**
1.  **Sort:** All genomes are ranked by length (longest to shortest). The assumption is that longer genomes are more "complete" and make better representatives.
2.  **Iterate:** The algorithm walks down the sorted list.
3.  **Seed Selection:** If the current genome is **unassigned**, it becomes the **Seed (Representative)** of a new vOTU.
4.  **Recruitment:** The algorithm looks at all **unassigned** neighbors of this Seed (based on the pre-computed ANI edges).
    *   If a neighbor meets the thresholds (default: ≥95% ANI, ≥85% Coverage), it is recruited into the cluster.
    *   **Crucially:** Recruitment is based strictly on similarity to the **Seed**, not to other members.
5.  **Finalize:** Once recruited, a genome is marked as "assigned" and cannot be picked again.

**Why Star Topology?**
This method prevents **"Chaining"** (or transitive clustering), a common error where Sequence A is similar to B, and B is similar to C, so A, B, and C are grouped together—even if A and C are totally different. In Star Topology, A and C are only grouped if they *both* match the Seed directly.

### Step 5: Catalog Extraction
*   **Input:** The list of selected "Seed" sequences.
*   **Tool:** `seqkit grep`
*   **Operation:** Extracts the full nucleotide sequences of the representatives from `all_genomes.fa`.
*   **Output:** A final FASTA file (`*_vOTU_catalog.fasta`) containing exactly one sequence per vOTU.

## 4. Key Files
*   `src/viral_cataloger/pipeline.py`: Contains the `greedy_star_clustering` function and `run_command` wrappers.
*   `src/viral_cataloger/main.py`: Handles CLI arguments (`--min_ani`, `--min_tcov`, etc.).
