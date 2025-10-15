---
title: "rgi_assembly_DA"
output: html_document
date: "2024-10-29"
---

This notebook does the following: loads metagenomic ARG count data per assembler, normalizes it (RPKM/TPM), assigns sample statuses, performs differential abundance analysis by gene and drug class using Maaslin2 and ANCOM-BC, and generates heatmaps and boxplots to visualize the differences.

## Setting working directory, libraries and importing read count files
 
```{r}
setwd("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/metmon/rgi_heatmap_boxplot")
library(tidyverse)
mega<-read_tsv("megahit_read_counts.tsv")
spa<-read_tsv("spades_read_counts.tsv")

lib<-read_table("reads_per_sample.tsv",   # Space delimiter
          col_names = c("Sample", "Size"))
lib_processed <- lib %>%
  mutate(sample = str_extract(Sample, "^SL\\d+")) %>%  # Extract 'SL' followed by numbers
  group_by(sample) %>% 
  summarise(final=sum(Size)) %>% mutate(final=final/4)
mega
elim<-c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
```
## NORMALIZATION BY rpkm

```{r}

spa2<-spa %>% filter(!gene_name%in% elim) %>% mutate(length=end-start+1) %>% inner_join(lib_processed, by="sample") %>% 
  mutate(fraction=conteo/(final/1000000)) %>% 
  mutate(rpkm=fraction/(length/1000))

mega2<-mega %>% filter(!gene_name%in% elim) %>% mutate(length=end-start+1) %>% inner_join(lib_processed, by="sample") %>% 
  mutate(fraction=conteo/(final/1000000)) %>% 
  mutate(rpkm=fraction/(length/1000))

```
## NORMALIZATION BY TPM

```{r}

spa3<-spa %>% filter(!gene_name%in% elim) %>% mutate(length=end-start+1) %>%
  group_by(sample) %>%
  mutate(RPK = conteo / (length/1000),     # Step 2: Calculate reads per length
     scaling_factor = sum(RPK) / 1e6,  # Sum RPK and scale to millions
     TPM = RPK / scaling_factor) %>%   # Step 3: Calculate TPM
  ungroup()

mega3<-mega %>% filter(!gene_name%in% elim) %>% mutate(length=end-start+1) %>%
  group_by(sample) %>%
  mutate(RPK = conteo / (length/1000),     # Step 2: Calculate reads per length
     scaling_factor = sum(RPK) / 1e6,  # Sum RPK and scale to millions
     TPM = RPK / scaling_factor) %>%   # Step 3: Calculate TPM
  ungroup()
  
status <- c("caries_free","caries_active", "caries_active", "caries_active", "caries_active", 
      "caries_active", "caries_active", "caries_free", "caries_free", "caries_free", 
      "caries_free", "caries_free", "caries_active", "caries_free", "caries_free", 
      "caries_active", "caries_free", "caries_active")

# Replace "free" with "health"
status <- gsub("free", "health", status)
unique_samples <- unique(mega2$sample)
sample_status <- data.frame(sample = unique_samples, status = status)
sample_status <- sample_status %>%
  mutate(status = case_when(
  status == "caries_health" ~ "Health",
  status == "caries_active" ~ "Caries",
  TRUE ~ status
  ))
# Assign row names of `sample_status` from the `sample` column
rownames(sample_status) <- sample_status$sample
# Remove the 'sample' column as it is now in the row names
sample_status$sample <- NULL
```

## RPKM
### MEGAHIT BY GENE AND APPLY MAASLIN2

```{r}
wide<-mega2 %>% group_by(sample, gene_name) %>% 
  summarise(conteo_final=sum(rpkm)) %>% 
  pivot_wider(names_from = sample, values_from = conteo_final, values_fill = list(conteo_final=0))
#wide<-wide %>% filter(startsWith(gene, "tet"))
megahit_df <- as.data.frame(t(wide[,-1])) # Transpose excluding gene name column
colnames(megahit_df) <- wide$gene_name     # Assign gene names as column names
       # Remove row names
colnames(megahit_df) <- gsub(" ", "_", colnames(megahit_df))

#BiocManager::install("Maaslin2")
library(Maaslin2)
fit_data=Maaslin2(input_data = megahit_df,
          input_metadata = sample_status,
            output = "Maaslin2_output_megahit_genes_rpkm",       # Output folder for results
          fixed_effects = c("status"), 
          min_prevalence = 0.2,
          normalization="TSS",
          transform = "LOG",
          analysis_method = "LM")

```

### BY CATEGORY
```{r}
mega_cat <- mega2 %>%
  separate_rows(drug_class, sep = ";") %>%
  mutate(drug_class = str_trim(drug_class)) %>% group_by(sample,drug_class ) %>% summarise(count=sum(rpkm)) %>%
  mutate(class=drug_class) %>%
  select(-drug_class)

wide<-mega_cat %>% 
  pivot_wider(names_from = sample, values_from = count, values_fill = list(count=0))
#wide<-wide %>% filter(startsWith(gene, "tet"))
megahit_class <- as.data.frame(t(wide[,-1])) # Transpose excluding gene name column
colnames(megahit_class) <- wide$class     # Assign gene names as column names
       # Remove row names
colnames(megahit_class) <- gsub(" ", "_", colnames(megahit_class))

#BiocManager::install("Maaslin2")
library(Maaslin2)
fit_data=Maaslin2(input_data = megahit_class,
          input_metadata = sample_status,
            output = "Maaslin2_output_megahit_class_rpkm",       # Output folder for results
          fixed_effects = c("status"), 
          min_prevalence = 0.2,
          normalization="TSS",
          transform = "LOG",
          analysis_method = "LM",
          heatmap_first_n = 5)
```

### SPADES BY GENE AND APPLY MAASLIN2

```{r}
wide<-spa2 %>% group_by(sample, gene_name) %>% 
  summarise(conteo_final=sum(rpkm)) %>% 
  pivot_wider(names_from = sample, values_from = conteo_final, values_fill = list(conteo_final=0))
#wide<-wide %>% filter(startsWith(gene, "tet"))
spa_df <- as.data.frame(t(wide[,-1])) # Transpose excluding gene name column
colnames(spa_df) <- wide$gene_name     # Assign gene names as column names
       # Remove row names
colnames(spa_df) <- gsub(" ", "_", colnames(spa_df))

#BiocManager::install("Maaslin2")
library(Maaslin2)
fit_data=Maaslin2(input_data = spa_df,
          input_metadata = sample_status,
            output = "Maaslin2_output_spa_genes_rpkm",       # Output folder for results
          fixed_effects = c("status"), 
          min_prevalence = 0.2,
          normalization="TSS",
          transform = "LOG",
          analysis_method = "LM")
```

### SPADES BY CATEGORY
```{r}
spa_cat <- spa2 %>%
  separate_rows(drug_class, sep = ";") %>%
  mutate(drug_class = str_trim(drug_class)) %>% group_by(sample,drug_class ) %>% summarise(count=sum(rpkm)) %>%
  mutate(class=drug_class) %>%
  select(-drug_class)

wide<-spa_cat %>% 
  pivot_wider(names_from = sample, values_from = count, values_fill = list(count=0))
#wide<-wide %>% filter(startsWith(gene, "tet"))
spa_class <- as.data.frame(t(wide[,-1])) # Transpose excluding gene name column
colnames(spa_class) <- wide$class     # Assign gene names as column names
       # Remove row names
colnames(spa_class) <- gsub(" ", "_", colnames(spa_class))

#BiocManager::install("Maaslin2")
library(Maaslin2)
fit_data=Maaslin2(input_data = spa_class,
          input_metadata = sample_status,
            output = "Maaslin2_output_spa_class_rpkm",       # Output folder for results
          fixed_effects = c("status"), 
          min_prevalence = 0.2,
          normalization="TSS",
          transform = "LOG",
          analysis_method = "LM",
          heatmap_first_n = 5)
```

## TPM: DOES NOT APPLY, FINAL ANALYSES WERE DONE USING RPKM
### MEGAHIT BY GENE AND APPLY MAASLIN2
```{r}
wide<-mega3 %>% group_by(sample, gene_name) %>% 
  summarise(conteo_final=sum(TPM)) %>% 
  pivot_wider(names_from = sample, values_from = conteo_final, values_fill = list(conteo_final=0))
#wide<-wide %>% filter(startsWith(gene, "tet"))
megahit_df <- as.data.frame(t(wide[,-1])) # Transpose excluding gene name column
colnames(megahit_df) <- wide$gene_name     # Assign gene names as column names
       # Remove row names
colnames(megahit_df) <- gsub(" ", "_", colnames(megahit_df))

#BiocManager::install("Maaslin2")
library(Maaslin2)
fit_data=Maaslin2(input_data = megahit_df,
          input_metadata = sample_status,
            output = "Maaslin2_output_megahit_genes_tpm",       # Output folder for results
          fixed_effects = c("status"), 
          min_prevalence = 0.2,
          normalization="TSS",
          transform = "LOG",
          analysis_method = "LM")

```
### BY CATEGORY
```{r}
mega_cat <- mega3 %>%
  separate_rows(drug_class, sep = ";") %>%
  mutate(drug_class = str_trim(drug_class)) %>% group_by(sample,drug_class ) %>% summarise(count=sum(TPM)) %>%
  mutate(class=drug_class) %>%
  select(-drug_class)

wide<-mega_cat %>% 
  pivot_wider(names_from = sample, values_from = count, values_fill = list(count=0))
#wide<-wide %>% filter(startsWith(gene, "tet"))
megahit_class <- as.data.frame(t(wide[,-1])) # Transpose excluding gene name column
colnames(megahit_class) <- wide$class     # Assign gene names as column names
       # Remove row names
colnames(megahit_class) <- gsub(" ", "_", colnames(megahit_class))

#BiocManager::install("Maaslin2")
library(Maaslin2)
fit_data=Maaslin2(input_data = megahit_class,
          input_metadata = sample_status,
            output = "Maaslin2_output_megahit_class_tpm",       # Output folder for results
          fixed_effects = c("status"), 
          min_prevalence = 0.2,
          normalization="TSS",
          transform = "LOG",
          analysis_method = "LM",
          heatmap_first_n = 5)
```

### SPADES BY GENE AND APPLY MAASLIN2

```{r}
wide<-spa3 %>% group_by(sample, gene_name) %>% 
  summarise(conteo_final=sum(TPM)) %>% 
  pivot_wider(names_from = sample, values_from = conteo_final, values_fill = list(conteo_final=0))
#wide<-wide %>% filter(startsWith(gene, "tet"))
spa_df <- as.data.frame(t(wide[,-1])) # Transpose excluding gene name column
colnames(spa_df) <- wide$gene_name     # Assign gene names as column names
       # Remove row names
colnames(spa_df) <- gsub(" ", "_", colnames(spa_df))

#BiocManager::install("Maaslin2")
library(Maaslin2)
fit_data=Maaslin2(input_data = spa_df,
          input_metadata = sample_status,
            output = "Maaslin2_output_spa_genes_tpm",       # Output folder for results
          fixed_effects = c("status"), 
          min_prevalence = 0.2,
          normalization="TSS",
          transform = "LOG",
          analysis_method = "LM")
```

### BY CATEGORY
```{r}
spa_cat <- spa3 %>%
  separate_rows(drug_class, sep = ";") %>%
  mutate(drug_class = str_trim(drug_class)) %>% group_by(sample,drug_class ) %>% summarise(count=sum(TPM)) %>%
  mutate(class=drug_class) %>%
  select(-drug_class)

wide<-spa_cat %>% 
  pivot_wider(names_from = sample, values_from = count, values_fill = list(count=0))
#wide<-wide %>% filter(startsWith(gene, "tet"))
spa_class <- as.data.frame(t(wide[,-1])) # Transpose excluding gene name column
colnames(spa_class) <- wide$class     # Assign gene names as column names
       # Remove row names
colnames(spa_class) <- gsub(" ", "_", colnames(spa_class))

#BiocManager::install("Maaslin2")
library(Maaslin2)
fit_data=Maaslin2(input_data = spa_class,
          input_metadata = sample_status,
            output = "Maaslin2_output_spa_class_tpm",       # Output folder for results
          fixed_effects = c("status"), 
          min_prevalence = 0.2,
          normalization="TSS",
          transform = "LOG",
          analysis_method = "LM",
          heatmap_first_n = 5)
```

## ANCOM-BC

TEST (Do not run)
```{r}
# Load devtools
library(devtools)

# Install the bugfix version of ANCOMBC
#install_github("FrederickHuangLin/ANCOMBC", ref = "bugfix")
library(ANCOMBC)

data(atlas1006, package = "microbiome")
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(atlas1006)

# subset to baseline
tse = tse[, tse$time == 0]

# Re-code the bmi group
tse$bmi = recode(tse$bmi_group,
         obese = "obese",
         severeobese = "obese",
         morbidobese = "obese")
# Subset to lean, overweight, and obese subjects
tse = tse[, tse$bmi %in% c("lean", "overweight", "obese")]
tse$bmi = factor(tse$bmi, levels = c("obese", "overweight", "lean"))
# You can verify the change by checking:
# levels(sample_data(tse)$bmi)

# Create the region variable
tse$region = recode(as.character(tse$nationality),
          Scandinavia = "NE", UKIE = "NE", SouthEurope = "SE", 
          CentralEurope = "CE", EasternEurope = "EE",
          .missing = "unknown")

# Discard "EE" as it contains only 1 subject
# Discard subjects with missing values of region
tse = tse[, ! tse$region %in% c("EE", "unknown")]

print(tse)
```

BELOW, NOTE WHERE THE FIRST DATA FRAME USED TO DETERMINE IF IT IS RPKM COMES FROM

### FOR MEGAHIT GENES 
FROM MEGAHIT RPKM/TPM. Use megahit_df

```{r}
# Combine `megahit_df` with `sample_status`
megahit_df_2 <- megahit_df %>%
  mutate(sample = rownames(.))
rownames(megahit_df_2)<-NULL
# Transpose the data frame
otu_mat <- megahit_df_2 %>%
  column_to_rownames(var = "sample") %>%
  t() %>%
  as.data.frame()
sample_status$sample <- rownames(sample_status)

sample_status <- sample_status %>%
  arrange(match(sample, colnames(otu_mat)))
library(TreeSummarizedExperiment)
library(S4Vectors)
# Create the TSE object
megahit_tse <- TreeSummarizedExperiment(
  assays = list(counts = as.matrix(otu_mat)),
  colData = DataFrame(sample_status)
)
megahit_tse$status <- factor(megahit_tse$status, levels = c("Health", "Caries"))

print(megahit_tse)
out_megahit_genes <- ancombc2(
  data = megahit_tse,
  assay_name = "counts",
  fix_formula = "status",
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.20,
  lib_cut = 0,
  group = "status",
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.25,
  verbose = TRUE
)
res=out_megahit_genes$res
ancom_matrix=out_megahit_genes$bias_correct_log_table

library(dplyr)
library(tidyr)
library(ggplot2)

# Filter relevant data for statusHealth
df_status = res %>%
  dplyr::select(taxon, contains("status")) 

# Prepare data for log fold-change and color values
df_fig_status1 = df_status %>%
  dplyr::filter(diff_statusCaries == TRUE) %>%
  dplyr::mutate(lfc_value = ifelse(diff_statusCaries == TRUE, 
                   round(lfc_statusCaries, 2), 0)) %>%
  dplyr::select(taxon, lfc_value)

df_fig_status2 = df_status %>%
  dplyr::filter(diff_statusCaries == TRUE) %>%
  dplyr::mutate(color_value = ifelse(passed_ss_statusCaries == TRUE & diff_statusCaries == TRUE, 
                   "aquamarine3", "black")) %>%
  dplyr::select(taxon, color_value)

# Join log fold-change and color data
df_fig_status = df_fig_status1 %>%
  dplyr::left_join(df_fig_status2, by = "taxon")

# Adjust limits and mid-point for heatmap color
lo = floor(min(df_fig_status$lfc_value, na.rm = TRUE))
up = ceiling(max(df_fig_status$lfc_value, na.rm = TRUE))
mid = (lo + up) / 2

# Plot the heatmap
G<-ggplot(df_fig_status, aes(x = "Caries vs Health", y = taxon, fill = lfc_value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
             na.value = "white", midpoint = mid, limit = c(lo, up),
             name = "Log Fold Change") +
  geom_text(aes(label = lfc_value, color = color_value), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log Fold Changes as Compared to Healthy Subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
```

### FOR MEGAHIT CLASSES
FROM MEGAHIT RPKM/TPM. Use megahit_class
```{r}
# Combine `megahit_df` with `sample_status`
megahit_df_class <- megahit_class %>%
  mutate(sample = rownames(.))
rownames(megahit_df_class)<-NULL
# Transpose the data frame
otu_mat <- megahit_df_class %>%
  column_to_rownames(var = "sample") %>%
  t() %>%
  as.data.frame()
sample_status$sample <- rownames(sample_status)

sample_status <- sample_status %>%
  arrange(match(sample, colnames(otu_mat)))

library(TreeSummarizedExperiment)
library(S4Vectors)

# Create the TSE object
megahit_tse <- TreeSummarizedExperiment(
  assays = list(counts = as.matrix(otu_mat)),
  colData = DataFrame(sample_status)
)
megahit_tse$status <- factor(megahit_tse$status, levels = c("Health", "Caries"))

print(megahit_tse)
out_megahit_classes <- ancombc2(
  data = megahit_tse,
  assay_name = "counts",
  fix_formula = "status",
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.20,
  lib_cut = 0,
  group = "status",
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.25,
  verbose = TRUE
)
res=out_megahit_classes$res
ancom_matrix=out_megahit_classes$bias_correct_log_table

library(dplyr)
library(tidyr)
library(ggplot2)

# Filter relevant data for statusHealth
df_status = res %>%
  dplyr::select(taxon, contains("status")) 

# Prepare data for log fold-change and color values
df_fig_status1 = df_status %>%
  dplyr::filter(diff_statusCaries == TRUE) %>%
  dplyr::mutate(lfc_value = ifelse(diff_statusCaries == TRUE, 
                   round(lfc_statusCaries, 2), 0)) %>%
  dplyr::select(taxon, lfc_value)

df_fig_status2 = df_status %>%
  dplyr::filter(diff_statusCaries == TRUE) %>%
  dplyr::mutate(color_value = ifelse(passed_ss_statusCaries == TRUE & diff_statusCaries == TRUE, 
                   "aquamarine3", "black")) %>%
  dplyr::select(taxon, color_value)

# Join log fold-change and color data
df_fig_status = df_fig_status1 %>%
  dplyr::left_join(df_fig_status2, by = "taxon")

# Adjust limits and mid-point for heatmap color
lo = floor(min(df_fig_status$lfc_value, na.rm = TRUE))
up = ceiling(max(df_fig_status$lfc_value, na.rm = TRUE))
mid = (lo + up) / 2

# Plot the heatmap
G<-ggplot(df_fig_status, aes(x = "Caries vs Health", y = taxon, fill = lfc_value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
             na.value = "white", midpoint = mid, limit = c(lo, up),
             name = "Log Fold Change") +
  geom_text(aes(label = lfc_value, color = color_value), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log Fold Changes as Compared to Healthy Subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
```

### FOR GENES (SPADES)
RPKM/TPM
```{r}
# Combine `megahit_df` with `sample_status`
spa_df_class <- spa_df %>%
  mutate(sample = rownames(.))
rownames(spa_df_class)<-NULL
# Transpose the data frame
otu_mat <- spa_df_class %>%
  column_to_rownames(var = "sample") %>%
  t() %>%
  as.data.frame()
sample_status$sample <- rownames(sample_status)

sample_status <- sample_status %>%
  arrange(match(sample, colnames(otu_mat)))

library(TreeSummarizedExperiment)
library(S4Vectors)

# Create the TSE object
megahit_tse <- TreeSummarizedExperiment(
  assays = list(counts = as.matrix(otu_mat)),
  colData = DataFrame(sample_status)
)
megahit_tse$status <- factor(megahit_tse$status, levels = c("Caries", "Health"))

print(megahit_tse)
out_megahit_classes <- ancombc2(
  data = megahit_tse,
  assay_name = "counts",
  fix_formula = "status",
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.20,
  lib_cut = 0,
  group = "status",
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.25,
  verbose = TRUE
)
res=out_megahit_classes$res

library(dplyr)
library(tidyr)
library(ggplot2)

# Filter relevant data for statusHealth
df_status = res %>%
  dplyr::select(taxon, contains("status")) 

# Prepare data for log fold-change and color values
df_fig_status1 = df_status %>%
  dplyr::filter(diff_statusHealth == TRUE) %>%
  dplyr::mutate(lfc_value = ifelse(diff_statusHealth == TRUE, 
                   round(lfc_statusHealth, 2), 0)) %>%
  dplyr::select(taxon, lfc_value)

df_fig_status2 = df_status %>%
  dplyr::filter(diff_statusHealth == TRUE) %>%
  dplyr::mutate(color_value = ifelse(passed_ss_statusHealth == TRUE & diff_statusHealth == TRUE, 
                   "aquamarine3", "black")) %>%
  dplyr::select(taxon, color_value)

# Join log fold-change and color data
df_fig_status = df_fig_status1 %>%
  dplyr::left_join(df_fig_status2, by = "taxon")

# Adjust limits and mid-point for heatmap color
lo = floor(min(df_fig_status$lfc_value, na.rm = TRUE))
up = ceiling(max(df_fig_status$lfc_value, na.rm = TRUE))
mid = (lo + up) / 2

# Plot the heatmap
G<-ggplot(df_fig_status, aes(x = "Health vs Caries", y = taxon, fill = lfc_value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
             na.value = "white", midpoint = mid, limit = c(lo, up),
             name = "Log Fold Change") +
  geom_text(aes(label = lfc_value, color = color_value), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log Fold Changes for Health vs Caries") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
```
### FOR CLASSES (SPADES)
RPKM/TPM
```{r}
# Combine `megahit_df` with `sample_status`
spas_df_class <- spa_class %>%
  mutate(sample = rownames(.))
rownames(spas_df_class)<-NULL
# Transpose the data frame
otu_mat <- spas_df_class %>%
  column_to_rownames(var = "sample") %>%
  t() %>%
  as.data.frame()
sample_status$sample <- rownames(sample_status)

sample_status <- sample_status %>%
  arrange(match(sample, colnames(otu_mat)))

library(TreeSummarizedExperiment)
library(S4Vectors)

# Create the TSE object
megahit_tse <- TreeSummarizedExperiment(
  assays = list(counts = as.matrix(otu_mat)),
  colData = DataFrame(sample_status)
)
megahit_tse$status <- factor(megahit_tse$status, levels = c("Caries", "Health"))

print(megahit_tse)
out_megahit_classes <- ancombc2(
  data = megahit_tse,
  assay_name = "counts",
  fix_formula = "status",
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.20,
  lib_cut = 0,
  group = "status",
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.25,
  verbose = TRUE
)
res=out_megahit_classes$res

library(dplyr)
library(tidyr)
library(ggplot2)

# Filter relevant data for statusHealth
df_status = res %>%
  dplyr::select(taxon, contains("status")) 

# Prepare data for log fold-change and color values
df_fig_status1 = df_status %>%
  dplyr::filter(diff_statusHealth == TRUE) %>%
  dplyr::mutate(lfc_value = ifelse(diff_statusHealth == TRUE, 
                   round(lfc_statusHealth, 2), 0)) %>%
  dplyr::select(taxon, lfc_value)

df_fig_status2 = df_status %>%
  dplyr::filter(diff_statusHealth == TRUE) %>%
  dplyr::mutate(color_value = ifelse(passed_ss_statusHealth == TRUE & diff_statusHealth == TRUE, 
                   "aquamarine3", "black")) %>%
  dplyr::select(taxon, color_value)

# Join log fold-change and color data
df_fig_status = df_fig_status1 %>%
  dplyr::left_join(df_fig_status2, by = "taxon")

# Adjust limits and mid-point for heatmap color
lo = floor(min(df_fig_status$lfc_value, na.rm = TRUE))
up = ceiling(max(df_fig_status$lfc_value, na.rm = TRUE))
mid = (lo + up) / 2

# Plot the heatmap
G<-ggplot(df_fig_status, aes(x = "Health vs Caries", y = taxon, fill = lfc_value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
             na.value = "white", midpoint = mid, limit = c(lo, up),
             name = "Log Fold Change") +
  geom_text(aes(label = lfc_value, color = color_value), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log Fold Changes for Health vs Caries") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
```

## Figures ANCOMBC
Process: First run megahit genes and then run the two chunks below, then run megahit classes and then run the two chunks below.

MEGAHIT: TEST BOXPLOTS RPKM/TPM FOR ANCOMBC (DATA NORMALIZED BY BIAS AND LOG).

```{r}
# Define the taxa to plot
# Required libraries
library(ggplot2)
library(stringr)

library(ggrepel)
library(dplyr)
# Transpose the matrix to have taxa in columns and samples in rows
ancom_matrix_t <- as.data.frame(t(ancom_matrix))
ancom_matrix_t$sample <- rownames(ancom_matrix_t)
# Join the status of each sample
data <- ancom_matrix_t %>%
  left_join(sample_status, by = "sample")  # Make sure sample_status has a `sample` column with sample names

# Create a loop to iterate over each taxon and generate the plots
for (taxon in df_fig_status$taxon) {
  
  # Create the data frame for the current taxon
  data_taxon <- data %>%
  dplyr::select(all_of(df_fig_status$taxon), sample, status)
  
  # Create the plot
  gg <- ggplot(data_taxon, aes(x = status, y = .data[[taxon]], fill = status, label = sample)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +  # Boxplot without outliers and narrower
  geom_jitter(width = 0.1, size = 2, alpha = 0.7, color = "black") +  # Scatter points
  geom_text_repel(aes(label = sample), size = 3, color = "black", max.overlaps = 10) +  # Labels on points
  scale_fill_manual(values = c("Health" = "#1f78b4", "Caries" = "#33a02c")) +  # Custom colors
 labs(title = paste(paste0(toupper(substring(gsub("_", " ", taxon), 1, 1)), 
           substring(gsub("_", " ", taxon), 2)), "Normalised Count Between Health and Caries"),
     x = "Status",
     y = paste(taxon, "RPKM"))+
  theme_bw(base_size = 15) +  # Minimalist theme style for publications
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )
  
  # Save the plot to a file with a defined name pattern
  ggsave(filename = paste0("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/metmon/rgi_heatmap_boxplot/test/megahit_norm_classes_", taxon, "_rpkm.png"), plot = gg, width = 8, height = 6, bg = "white")
}
```
TEST BOXPLOTS FOR RAW ANCOM (DOES NOT USE ANCOM NORMALIZED TABLE)

```{r}
# Required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)

# Create a loop to iterate over each taxon and generate the plots
for (taxon in df_fig_status$taxon) {
  ###HERE CHANGE TO MEGAHIT_DF IF COMPARING GENES
  # Create the data frame for the current taxon
  data <- megahit_class %>%
  dplyr::select(!!sym(taxon)) %>%
  dplyr::mutate(sample = sample_status$sample, status = sample_status$status)
  
  # Create the plot
  gg <- ggplot(data, aes(x = status, y = !!sym(taxon), fill = status, label = sample)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +  # Boxplot without outliers and narrower
  geom_jitter(width = 0.1, size = 2, alpha = 0.7, color = "black") +  # Scatter points
  geom_text_repel(aes(label = sample), size = 3, color = "black", max.overlaps = 10) +  # Labels on points
  scale_fill_manual(values = c("Health" = "#1f78b4", "Caries" = "#33a02c")) +  # Custom colors
 labs(title = paste(paste0(toupper(substring(gsub("_", " ", taxon), 1, 1)), 
           substring(gsub("_", " ", taxon), 2)), "Count Between Health and Caries"),
     x = "Status",
     y = paste(taxon, "RPKM")) +
  theme_minimal(base_size = 15) +  # Minimalist theme style for publications
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )
  
  # Save the plot to a file with a defined name pattern
  ggsave(filename = paste0("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/metmon/rgi_heatmap_boxplot/test/megahit_raw_classes_", taxon, "_rpkm.png"), plot = gg, width = 8, height = 6, bg = "white")
}
```

