# Viro-Flow & Cataloger: The Complete Viral Discovery Workflow

This document explains the intentional relationship between the **Viral Discovery Pipeline (Viro-Flow)** and the **Viral Genome Cataloger**. Together, they form a complete ecosystem for transforming raw metagenomic data into species-level viral catalogs.

## Executive Summary

*   **Viral Discovery Pipeline (Viro-Flow)** is a comprehensive, end-to-end **metagenomics pipeline**. It takes raw sequencing data (FASTQ) and performs every step necessary to find viruses, including cleaning, assembly, identification, and abundance profiling.
*   **Viral Genome Cataloger** is a specialized **clustering tool**. It takes already-assembled viral genomes (FASTA) and groups them into non-redundant species catalogs (vOTUs) using a specific algorithm to prevent clustering artifacts.

---

## 1. Purpose and Scope

| Feature | Viral Discovery Pipeline (Viro-Flow) | Viral Genome Cataloger |
| :--- | :--- | :--- |
| **Primary Goal** | **Discovery**: Turn raw sequencing reads into a viral community abundance matrix. | **Dereplication**: Turn a collection of viral genomes into a non-redundant catalog. |
| **Input** | Raw sequencing data (**FASTQ** files). | Assembled genome sequences (**FASTA** files). |
| **Output** | Quality reports, viral contigs, and a sample-by-virus abundance matrix (TPM). | A single FASTA file of representative viral species (vOTUs) and a cluster map. |
| **Use Case** | You have raw data from a sequencer and want to know "What viruses are in my samples?" | You have viral genomes from multiple studies/samples and want to merge them into a single database. |

---

## 2. Methodological Differences

The most significant difference lies in how they handle **Clustering** (grouping similar viruses into species).

### Viral Discovery Pipeline (Viro-Flow)

*   **Method:** Uses a "connected components" approach.
*   **Tools:** Uses **CD-HIT-EST** for initial dereplication (99% identity) followed by **FastANI** for species-level clustering.
*   **Algorithm:** It calculates Average Nucleotide Identity (ANI). If Genome A matches Genome B, they are linked. It uses `networkx` to find connected components in this graph.
*   **Pros/Cons:** Standard approach, but susceptible to "chaining" (where A matches B, and B matches C, so A and C are grouped even if they are not similar).

### Viral Genome Cataloger

*   **Method:** Uses a "Greedy Star-Topology" approach.
*   **Tools:** Uses **skani** (a newer, faster alternative to FastANI).
*   **Algorithm:**
    1.  Sorts all genomes by length.
    2.  Picks the longest as a "centroid."
    3.  Recruits only genomes that match the centroid directly.
    4.  Removes them from the pool and repeats.
*   **Pros/Cons:** Explicitly designed to prevent the "chaining" artifacts mentioned above, ensuring all members of a cluster are similar to the representative.

---

## 3. Workflow and Architecture

### Viral Discovery Pipeline (Viro-Flow)

*   **Complexity:** High. It orchestrates a massive suite of bioinformatics tools including **fastp** (QC), **BWA** (host removal), **MEGAHIT** (assembly), **VirSorter2** & **GeNoMad** (detection), and **CheckV** (quality).
*   **Implementation:** Hybrid. It provides traditional Bash scripts (`viral_pipeline.sh`) for academic HPC users and a modern Python CLI (`viral-pipeline`) for developers.
*   **Environment:** Heavy reliance on **Conda** environments or **Docker/Singularity** containers due to the complex dependency tree.

### Viral Genome Cataloger

*   **Complexity:** Low. It focuses on doing one task efficiently.
*   **Implementation:** Pure Python package (`viral_cataloger`) that wraps external binaries.
*   **Dependencies:** Lightweight. Requires only **skani** and **SeqKit** installed externally, plus standard Python libraries.

---

## 4. User Experience

*   **Viro-Flow** is designed as a "set and forget" pipeline. You point it at a directory of raw data, and it runs for hours or days, handling logic for QC, assembly, and mapping automatically. It includes extensive config files to tune parameters for memory, threads, and databases.
*   **Viral Genome Cataloger** is a "utility" tool. It runs quickly (minutes to an hour). The user provides a folder of FASTA files, and it outputs the cleaned catalog. It has fewer configuration options, focused specifically on ANI thresholds and coverage.

---

## Recommendation

*   **Use `viral-discovery-pipeline`** if you are starting a new project with raw sequencing data and need to perform the heavy lifting of assembly and viral identification.
*   **Use `viral-genome-cataloger`** if you have already run the discovery pipeline on multiple different datasets and now need to combine all your results into a single, high-quality reference database without duplicates.

---

## 🔄 Workflow Integration

The combined workflow follows a "Discovery then Dereplication" logic:

1.  **Mining (Viro-Flow)**: Processes raw FASTQ files to identify potential viral sequences per sample.
2.  **Librarianship (Cataloger)**: Aggregates discoveries from all samples and organizes them into a non-redundant catalog of Viral Operational Taxonomic Units (vOTUs).

```mermaid
graph LR
    A[Raw Reads] --> B(Viro-Flow)
    B --> C[Sample-Specific Viruses]
    C --> D(Genome Cataloger)
    D --> E[Final vOTU Catalog]
```

### 🛠 Step 1: Mining with Viro-Flow

Run the `viral_pipeline.sh` for each of your sequencing samples. This stage handles quality control, host removal, assembly, and viral identification.

*   **Input:** Raw sequencing reads (`.fastq` files).
*   **Action:** It cleans the data, removes host DNA, assembles reads into contigs, and identifies which contigs are viral.
*   **Output:** For every sample, it produces a FASTA file containing identified viral genomes, named something like `Sample1_final_viruses.fna`.
*   **Location:** Found in the individual analysis directory for each sample.

### 🛠 Step 2: The Handover

Collect the `.fna` files from all your Viro-Flow analysis directories and place them into a single input folder for the cataloger.

```bash
# Example of gathering outputs
mkdir all_discovered_viruses
cp *_analysis/*_final_viruses.fna all_discovered_viruses/
```

### 🛠 Step 3: Cataloging with Genome Cataloger

Run the `viral-cataloger` on your gathered sequences. It uses `skani` for high-resolution ANI estimation and a Greedy Star-Topology algorithm to prevent "chaining" artifacts.

*   **Scientific Standard:** It defaults to a **95% ANI** and **85% Coverage** threshold, which is the international standard for demarcating viral species.
*   **Primary Output:** `catalog_vOTU_catalog.fasta` (one representative per species).

---

## 🧠 Scientific Rationale

This two-repo approach was designed to solve specific challenges in viral metagenomics:

*   **Data Integrity**: Viro-Flow uses `CheckV` to assess the quality of every discovered contig before it enters the catalog.
*   **Clustering Precision**: While Viro-Flow offers a multi-sample script, the **Genome Cataloger** is the preferred method for building the final reference database. Its "star-topology" ensures that every member of a vOTU is directly similar to a high-quality, long representative genome, rather than being grouped via a weak chain of intermediates.
*   **Performance**: By offloading the final clustering to `skani` within the Cataloger, you achieve much higher speeds than traditional BLAST-based methods.
