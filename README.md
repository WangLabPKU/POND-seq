# Analysis Pipeline for RNA-Binding Protein Interacting Transcripts Based on POND-seq data

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

👤 **Author:** Gang Xie, PhD Candidate  
🏫 **Affiliation:** PKU-THU-NIBS Joint Graduate Program, Academy for Advanced Interdisciplinary Studies, Peking University, Beijing, China  
✉️ **Email:** gangx1e@stu.pku.edu.cn  
📅 **Date:** May 8th, 2026  
✅ **Version:** 1.0  

## Overview
This repository contains the custom analysis pipeline and scripts used in the manuscript: *"Longitudinal monitoring of cytoplasmic RBP-RNA interactions and transcriptome in living cells by engineered protein nanocages"*. 

It provides a comprehensive workflow for analyzing transcripts interacting with RNA-binding proteins (RBPs) utilizing our novel omics technology, **POND-seq**. The pipeline covers the main process from raw high-throughput sequencing data preprocessing and alignment to downstream bioinformatics and statistical analyses.

## System Requirements
- **OS:** Linux (Tested on Ubuntu 20.04 / CentOS 7)
- **Memory:** Minimum 64GB RAM is recommended for full dataset processing.
- **Languages:** R (>= 4.2.0), Python (>= 3.9)

### Key Dependencies
- **Command-line Tools:** `trim_galore`, `STAR`, `subread` (featureCounts), `samtools`, `bowtie`, `bedtools`, `umi_tools`
- **Python Packages:** `pandas`, `numpy`, `matplotlib`, `seaborn`
- **R Packages:** `DESeq2`, `clusterProfiler`, `Mfuzz`

> **Note:** We highly recommend using [Miniconda](https://docs.conda.io/en/latest/miniconda.html) to manage your software environments and dependencies.

## Data Availability
The raw sequencing data (FASTQ format) and processed count matrices generated in this study have been deposited in the NCBI Gene Expression Omnibus (GEO) and are publicly accessible under accession number: **[GSE293919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293919)**.

## Usage & Repository Structure
The scripts provided in this repository are intended as a reference for reproducing the analysis described in our study. The repository is organized as follows:

*   **`scripts/`**: Contains example shell scripts for upstream sequence processing, quality control, and read alignment.
*   **`plot/`**: Contains scripts for downstream analyses, including target enrichment evaluation, comparative transcriptomics, and functional profiling (e.g., GO/KEGG analysis).

**Additional Information:**
The code provided here serves as a representative framework. Additional custom scripts or specific intermediate data processing codes are available upon reasonable request via email. If you have any questions, encounter bugs, or need assistance running the pipeline, please feel free to open an issue or reach out directly.

## Citation
If you find this code, the POND-seq methodology, or our datasets useful for your research, please cite our paper:

> #Hu LF, #Xie G, Wu YX, Li YX, Wan ZL, Mi L, Wang JZ, *Wang Y. *"Longitudinal monitoring of cytoplasmic RBP-RNA interactions and transcriptome in living cells by engineered protein nanocages."* **Molecular Cell**, 2026 (Accepted).

## License
This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)** - see the [LICENSE](https://www.gnu.org/licenses/gpl-3.0) page for details.