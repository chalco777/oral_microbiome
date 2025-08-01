# Oral Microbiome Metagenomic Analysis

A pipeline for characterizing the oral microbiome in health-versus-caries cohorts via shotgun metagenomics:

* **Host read removal & QC:** filter human reads against T2T-CHM13v2.0 (BWA-MEM2 + SAMtools).
* **Taxonomic profiling:** (Kraken2 → Bracken)/Sylph → Pavian for absolute counts and relative abundances.
* **Beta-diversity:** PCoA on species matrices to assess community-level shifts.
* **Assembly & binning:** [nf-core/mag](https://nf-co.re/mag/3.3.0/) for contigs, bins, viral identification, and QC summaries.
* **Resistome analysis:** RGI/HTSeq-count → RPKM normalization → MaAsLin2 & ANCOM-BC2 for differential ARGs.

## Prerequisites

- R (≥4.0) with packages: phyloseq, MaAsLin2, ANCOMBC, ggplot2, vegan  
- Python 3 with: pandas, seaborn (for any custom scripts)  
- Data:  
  - `matrix_allranks_conteo.tsv`  
  - `metadata.xlsx`


## Repository structure

```
├── barplot_boxplot/         # Rmd and PNGs for differential abundance
├── lefse_pcoa/              # PCoA scripts and results
├── rgi_ensamblajes/         # Functional annotation & MaAsLin2 outputs
├── rgi_reads_heatmap/       # Raw counts and heatmap scripts
├── software_scripts/        # Bash and Nextflow workflows
└── README.md                # ← you are here
```

---

## Recommended Analysis Sequence

1. **Input:** Get paired-end Illumina reads (R1/R2) per sample.

2. **Host read removal (keep only non-human pairs)**
   - Align to T2T-CHM13v2.0 with BWA-MEM2 → convert SAM→BAM, sort by name.
   - Remove reads mapped to human; keep pairs where both mates are unmapped and convert back to FASTQ for downstream use.
   - *Hint*: `samtools view -f 12` selects read pairs where both mates are unmapped.

3. **Taxonomic profiling**
   - Kraken2 (Standard DB) → Bracken (e.g., `-r 150 -l S`), as well as Sylph, to estimate species-level abundances.
   - Use Pavian to export absolute counts and relative abundances across taxonomic ranks (produce matrices such as `matrix_allranks_conteo.tsv`).

4. **Differential abundance & beta diversity**
   - Rarefaction of full-rank count matrices for LEfSe; perform a separate species-level rarefaction for PCoA.
   - LEfSe example:
     ```bash
     lefse_format_input.py lefse_fullranks.tsv species_reads_fullrank.in -c 2 -s -1 -u 1 -o 1000000
     lefse_run.py species_reads_fullrank.in species_bracken_fullrank.res
     lefse_plot_res.py species_bracken_fullrank.res species_bracken_ranksall.png --dpi 1000
     ```
   - PCoA: compute on the (rarefied) species matrix; log/center as appropriate.

5. **Assembly & binning (nf-core/mag)**
   - Assembly-only ([nf-core/mag v3.0.3](https://nf-co.re/mag/3.0.3/)): run to generate contigs and basic metrics (e.g., MEGAHIT/SPAdes, QUAST); skip downstream classification/annotation steps.
     ```bash
     nextflow run nf-core/mag --input samplesheet.csv --outdir . -profile docker \
       --cleanup_workdir true --skip_clipping true \
       --skip_gtdbtk true --skip_prodigal true --skip_prokka true \
       --skip_metaeuk true --skip_ancient_damagecorrection true --keep_phix .
     ```
   - Binning & QC ([nf-core/mag v3.3.0](https://nf-co.re/mag/3.3.0/)): enable viral identification (geNomad), binning (MetaBAT2/MaxBin2/CONCOCT), DAS Tool refinement, CheckM2/BUSCO for QC, GUNC for contamination, and MultiQC summaries.
     ```bash
     nextflow run nf-core/mag --input samplesheet.csv --assembly_input samplesheetassembly.csv \
       --outdir . -profile docker --skip_prodigal true --skip_prokka true --skip_metaeuk true \
       --run_virus_identification true --refine_bins_dastool true \
       --binqc_tool checkm2 --binqc_tool busco --run_gunc true
     ```

6. **ARGs (resistance genes): quantification & differential analysis**
   - Detection: run RGI on scaffolds (local DIAMOND mode) and on raw reads (rgi bwt), using CARD.
   - Counting per gene: sort/index BAM (SAMtools); add read groups (Picard AddOrReplaceReadGroups); mark duplicates (GATK MarkDuplicatesSpark); convert RGI calls to GTF; count with HTSeq-count (e.g., `--format=bam --order=pos --stranded=no --minaqual=20 --mode=union --nonunique=none --secondary-alignments=ignore --supplementary-alignments=score`).
   - Normalization & DA: compute RPKM; test differential abundance with MaAsLin2 (fixed effect: status Healthy vs. Caries; prevalence ≥0.2; TSS; log; linear model) and ANCOM-BC2 for compositional bias and LFCs.
   - Visualization: heatmaps/boxplots in R (ggplot2).


## Scripts description

### 1. `barplot_boxplot/`
- **Differential_taxonomic_analysis.Rmd**
  - **Role:** Performs differential taxonomic analysis to identify species with significant abundance differences between groups (e.g., caries vs. control). Generates boxplots and barplots for visualization.
  - **Sequence:** Run after generating abundance tables with Kraken2-Bracken and Sylph. Use the output plots and tables for downstream interpretation.

### 2. `lefse_pcoa/`
- **PCoA.Rmd**
  - **Role:** Conducts Principal Coordinates Analysis (PCoA) on species abundance data obtained through Kraken2-Bracken and Sylph to visualize sample clustering and diversity patterns.
  - **Sequence:** Run after generating the abundance matrix with the taxonomic profiling tools and have at hand the metadata. Use the resulting plots to interpret community structure.
- **pavian_to_lefse_preprocessing.R**
  - **Role:** Preprocesses data for LEfSe analysis, converting outputs from Pavian into the required format.
  - **Sequence:** Run before LEfSe script in `software_scripts/lefse.nf`  if using Pavian outputs.

### 3. `rgi_ensamblajes/`
- **assemblys_maaslin_ancom.Rmd**
  - **Role:** Integrates assembly data with MaAsLin2 and ANCOM-BC for differential abundance analysis at the gene or class level of antibiotic resistance genes detected with RGI.
  - **Sequence:** Run after assembling metagenomes and generating count tables from RGI. Use the results for downstream statistical analysis and visualization.

### 4. `rgi_reads_heatmap/`
- **rgi.Rmd**
  - **Role:** Analyzes raw read counts for resistance genes according to the RGI output and generates heatmaps to visualize gene abundance across samples.
  - **Sequence:** Run after obtaining read count tables. Use the heatmaps for interpretation of resistance gene distribution.


## Notes

- Ensure all required input files (abundance matrices, metadata, count tables) are prepared before running each script.
- Output plots and tables are saved in their respective folders for reporting and interpretation.
- For further details, please refer to comments within each RMD script or contact the project maintainer.

## Script Usage

1. Open the appropriate R Markdown (`.Rmd`) in RStudio.  
2. Knit to HTML or PDF to reproduce figures and tables.  
3. Review output in each subfolder (e.g. `barplot_boxplot/`, `lefse_pcoa/`, etc.).  

---
