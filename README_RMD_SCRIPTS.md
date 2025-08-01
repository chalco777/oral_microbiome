# Microbiome Oral Project

This README provides an overview of the R Markdown (`.Rmd`) scripts found in the `Microbioma oral` folder and its subfolders. Each script is described with its role and the recommended sequence for analysis.

## Folder Structure and RMD Scripts

### 1. `barplot_boxplot/`
- **Differential_taxonomic_analysis.Rmd**
  - **Role:** Performs differential taxonomic analysis to identify species with significant abundance differences between groups (e.g., caries vs. control). Generates boxplots and barplots for visualization.
  - **Sequence:** Run after preparing abundance and metadata tables. Use the output plots and tables for downstream interpretation.

### 2. `lefse_pcoa/`
- **PCoA.Rmd**
  - **Role:** Conducts Principal Coordinates Analysis (PCoA) on species abundance data to visualize sample clustering and diversity patterns.
  - **Sequence:** Run after generating the abundance matrix and metadata. Use the resulting plots to interpret community structure.
- **pavian_to_lefse_preprocessing.R** (not RMD, but relevant)
  - **Role:** Preprocesses data for LEfSe analysis, converting outputs from Pavian or other tools into the required format.
  - **Sequence:** Run before LEfSe or differential analysis scripts if using Pavian outputs.

### 3. `rgi_ensamblajes/`
- **assemblys_maaslin_ancom.Rmd**
  - **Role:** Integrates assembly data with MaAsLin2 and ANCOM-BC for differential abundance analysis at the gene or class level.
  - **Sequence:** Run after assembling metagenomes and generating count tables. Use the results for downstream statistical analysis and visualization.

### 4. `rgi_reads_heatmap/`
- **rgi.Rmd**
  - **Role:** Analyzes raw read counts for resistance genes and generates heatmaps to visualize gene abundance across samples.
  - **Sequence:** Run after obtaining read count tables. Use the heatmaps for interpretation of resistance gene distribution.

## Recommended Analysis Sequence
1. **Data Preprocessing:**
   - Use `pavian_to_lefse_preprocessing.R` (if needed) to format data for downstream analysis.
2. **Taxonomic and Diversity Analysis:**
   - Run `Differential_taxonomic_analysis.Rmd` for taxonomic comparisons.
   - Run `PCoA.Rmd` for diversity and clustering analysis.
3. **Functional and Resistance Analysis:**
   - Run `assemblys_maaslin_ancom.Rmd` for gene/class-level differential analysis.
   - Run `rgi.Rmd` for resistance gene abundance and visualization.

## Notes
- Ensure all required input files (abundance matrices, metadata, count tables) are prepared before running each script.
- Output plots and tables are saved in their respective folders for reporting and interpretation.

For further details, refer to comments within each RMD script or contact the project maintainer.
