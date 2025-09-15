# Pipeline Methodology and Tools

This document provides a detailed, step-by-step breakdown of the bioinformatics tools, algorithms, and parameters used in the `create_genome_catalog.py` script.

---

## 🧬 Step 1: Genome Aggregation

* **Tool**: Standard Python I/O (`glob` module)
* **Purpose**: To combine all individual genome FASTA files from the input directory into a single, large FASTA file.
* **Why this is performed**: Performing an "all-vs-all" comparison requires a single file containing all sequences. This step prepares the data for efficient processing by `skani`.
* **Parameters**: None. This is a hardcoded data preparation step.

---

## 🔬 Step 2: All-vs-All ANI Calculation

* **Tool**: **skani**
* **Purpose**: To calculate the **Average Nucleotide Identity (ANI)** for every possible pair of genomes in the aggregated file. ANI is a measure of genetic similarity and is the primary metric used for clustering.
* **Why this is performed**: This step generates the raw similarity data needed to determine which genomes are closely related enough to be considered part of the same species-level group (vOTU).
* **Parameters & Rationale**:
    * `skani triangle`: The command used for an efficient all-vs-all comparison.
    * `-i`: **(Critical)** This tells `skani` to treat each individual FASTA sequence (contig) as a separate entity for comparison. This is essential for analyzing metagenomic assemblies where each sequence is a potential genome or genome fragment.
    * `-m 200`: Sets the number of minimizers (short, representative k-mers) used for the initial sketch. A lower value like 200 is recommended for smaller viral genomes to ensure enough markers are used for a reliable comparison.
    * `--slow`: This flag enables a more accurate alignment-based calculation of the alignment fraction (AF) or coverage. This accuracy is important for the clustering step to avoid false groupings based on incomplete alignments.
    * `-E`: Sets the output to an "edge list" format (`Ref-ID`, `Query-ID`, `ANI`, etc.), which is a simple table that is easy to parse in downstream steps.
    * `--faster-small`: An optimization flag recommended for large datasets of small genomes (like viruses) to speed up the computation.
    * `-s 90`: This is a pre-screening threshold. `skani` first performs a very fast, approximate comparison and only proceeds with the full, slower alignment for pairs that are likely above 90% ANI. This dramatically reduces runtime without affecting the accuracy of high-identity pairs.
* **Credit & Citation**:
    * **GitHub**: [bluenote-1577/skani](https://github.com/bluenote-1577/skani)
    * **Publication**: Ondov, B.D., Starrett, G.J., Sappington, A. et al. Fast and accurate average nucleotide identity estimation with `skani`. *Nat Biotechnol* **41**, 1733–1736 (2023). [https://doi.org/10.1038/s41587-023-01774-x](https://doi.org/10.1038/s41587-023-01774-x)

---

## ✍️ Step 3: Reformatting `skani` Output

* **Tool**: **awk**
* **Purpose**: To convert the output from `skani` into the specific 6-column format required by the clustering script.
* **Why this is performed**: This is a data-wrangling step to ensure compatibility between different components of the pipeline. The clustering script expects coverage values as percentages (0-100), while `skani` outputs them as fractions (0-1). `awk` handles this conversion and column reordering.
* **Parameters & Rationale**:
    * `'BEGIN{OFS="\t"} NR>1{print $1, $2, "1", $3, $4*100, $5*100}'`: This command tells `awk` to:
        * `BEGIN{OFS="\t"}`: Set the output field separator to a tab.
        * `NR>1`: Skip the first line (the header).
        * `{print ...}`: For all other lines, print the columns in the desired order, multiplying the coverage values in columns 4 and 5 by 100.

---

## 🧩 Step 4: Greedy Clustering

* **Tool**: Custom Python Function (within `create_genome_catalog.py`)
* **Purpose**: To group genomes into clusters that represent vOTUs based on the calculated ANI and coverage values.
* **Why this is performed**: This is the core step where the raw similarity data is used to define species-level boundaries and identify a single representative member for each group, thus dereplicating the dataset.
* **Process**: The script uses a **greedy clustering algorithm**:
    1.  All genomes are sorted by length in descending order.
    2.  The script iterates through the sorted list. The longest genome that has not yet been assigned to a cluster becomes the **representative** of a new cluster.
    3.  All other unclustered genomes that meet the ANI and coverage thresholds with this representative are assigned to its cluster.
    4.  This process repeats until every genome is in a cluster.
* **Parameters & Rationale**:
    * `--min_ani 95.0`: This sets the minimum ANI for two genomes to be in the same cluster to **95%**. This is the widely accepted standard threshold for demarcating viral species.
    * `--min_tcov 85.0`: This sets the minimum target coverage to **85%**. This ensures that the 95% ANI similarity extends across at least 85% of the target genome's length, preventing false positives from short, highly conserved regions (e.g., a single gene).

---

## ✅ Step 5: Final Catalog Generation

* **Tool**: **SeqKit**
* **Purpose**: To extract the full FASTA sequences of only the representative genomes identified during the clustering step.
* **Why this is performed**: This final step creates the main output of the pipeline: a clean, non-redundant FASTA file containing one high-quality sequence for each vOTU, ready for downstream analysis.
* **Parameters & Rationale**:
    * `seqkit grep`: The command used to retrieve sequences based on their IDs.
    * `-f representative_ids.txt`: This flag tells `seqkit` to read the list of sequence IDs to extract from the specified file.
    * `-o <output.fasta>`: Specifies the path for the final output FASTA file.
* **Credit & Citation**:
    * **GitHub**: [shenwei356/seqkit](https://github.com/shenwei356/seqkit)
    * **Publication**: Shen, W., Le, S., Li, Y. et al. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. *PLoS ONE* **11(10)**: e0163962 (2016). [https://doi.org/10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962)