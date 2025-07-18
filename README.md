# Oral Microbiome Metagenomic Analysis

A pipeline for characterizing the oral microbial community via shotgun metagenomics:  
- **Taxonomic profiling** (differential abundance, boxplots)  
- **Beta‑diversity analyses** (PCoA)  
- **Functional annotation** (AMR genes, mobile genetic elements)  
- **Statistical modeling** (MaAsLin2, ANCOM-BC)

---

## Prerequisites

- R (≥4.0) with packages: phyloseq, MaAsLin2, ANCOMBC, ggplot2, vegan  
- Python 3 with: pandas, seaborn (for any custom scripts)  
- Data:  
  - `matrix_allranks_conteo.tsv`  
  - `metadata.xlsx`

---

## Repository structure

├── barplot_boxplot/ # Rmd and PNGs for differential abundance
├── lefse_pcoa/ # PCoA scripts and results
├── rgi_ensamblajes/ # Functional annotation & MaAsLin2 outputs
├── rgi_reads_heatmap/ # Raw counts and heatmap scripts
├── *.docx # Drafts, methodology, final reports
└── README.md # ← you are here

---

## Usage

1. Open the appropriate R Markdown (`.Rmd`) in RStudio.  
2. Knit to HTML or PDF to reproduce figures and tables.  
3. Review output in each subfolder (e.g. `barplot_boxplot/`, `lefse_pcoa/`, etc.).  

---
