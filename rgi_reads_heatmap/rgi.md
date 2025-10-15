---
title: "rgi"
output: html_document
date: "2024-10-20"
---

This script performs a complete read-based ARG analysis by: importing raw read counts, annotating samples by caries status, normalizing gene counts using log transformation, TPM, TPKM, and RPKM methods, visualizing gene and drug class profiles with hierarchical clustering and k-means heatmaps, applying subsampling for read depth normalization, and conducting differential abundance testing with Maaslin2 for both gene and drug class levels.

## Setup and read count results from RGI preprocessing

```{r}
library(tidyverse)
setwd("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/metmon/rgi_heatmap_boxplot")
reads_raw<-read_tsv("raw_reads_count.tsv.txt")

reads<-reads_raw %>% separate(col=`SAMPLE ARO Term`, into = c("sample","gene"),sep =" ", extra="merge") %>% select(c(1,2,"All Mapped Reads")) %>% #26 drug class, 2 is the gene
  mutate(across(c(3),as.numeric))

status <- c("caries_free", "caries_active", "caries_active", "caries_active", "caries_active", 
            "caries_active", "caries_active", "caries_free", "caries_free", "caries_free", 
            "caries_free", "caries_free", "caries_active", "caries_free", "caries_free", 
            "caries_active", "caries_free", "caries_active")

# Replace "free" with "health"
status <- gsub("free", "health", status)
unique_samples <- unique(reads$sample)
sample_status <- data.frame(sample = unique_samples, status = status)
sample_status <- sample_status %>%
  mutate(status = case_when(
    status == "caries_health" ~ "Health",
    status == "caries_active" ~ "Caries",
    TRUE ~ status
  ))
reads <- reads %>%
  left_join(sample_status, by = "sample") %>% mutate(sample = paste(status, sample, sep = "_")) %>%
  select(-status)  # Optional: remove the `status` column if it's no longer needed
```


```{r}
wide<-reads %>% 
  pivot_wider(names_from = sample, values_from = `All Mapped Reads`, values_fill = list(`All Mapped Reads`=0))

#wide<-wide %>% filter(startsWith(gene, "OXA"))
reads_mat<-as.matrix(wide[,-1])
rownames(reads_mat)<-wide$gene

reads_mat <- log1p(reads_mat)

```

## By Individual genes

In theory, normalization follows from here. I understand that I am normalizing not based on the total, but based on the mean of each gene in each row, so that the resulting heatmap is interpretable.

### Hierarchical clustering and Heatmap

```{r}
reads_scaled <- t(scale(t(reads_mat)))
#reads_scaled<-reads_mat
dis_genes<-dist(reads_scaled, method ="euclidean")
dis_samples<-dist(t(reads_scaled), method = "euclidean")

hc_genes <- hclust(dis_genes, method = "average")
hc_samples <- hclust(dis_samples, method = "average")

library(pheatmap)
library(RColorBrewer)

# Create the heatmap
g<-pheatmap(reads_scaled,
         cluster_rows = hc_genes,
         cluster_cols = hc_samples,
         show_rownames = TRUE,  # Hide gene names if there are too many
         fontsize_col = 8,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of Number of reads per gene vs Samples (Hierarchical clustering)")

```

### Now apply k-means

```{r}

reads_mat_sc <- t(scale(t(reads_mat)))
mat_transposed <- t(reads_mat_sc)
kmeans_result <- kmeans(mat_transposed, centers = 2)
#obtain clusters
cluster_order <- order(kmeans_result$cluster)
# Reorder the transposed matrix according to clusters
mat_clustered <- mat_transposed[cluster_order, ]
# Create an annotation for the clusters (metadata of which cluster each sample belongs to)
annotation_col <- data.frame(Cluster = factor(kmeans_result$cluster[cluster_order]))
rownames(annotation_col) <- rownames(mat_clustered)

# Create the heatmap
h<-pheatmap(mat_clustered,
         annotation_row = annotation_col,
         show_colnames = TRUE,  # Hide gene names if there are too many
         fontsize_row = 8,
         color = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100),
         main = "Heatmap of Genes vs Samples (Grouped by K-means)")
```

## By Gene Categories

```{r}
reads<-reads_raw %>% separate(col=`SAMPLE ARO Term`, into = c("sample","gene"),sep =" ", extra="merge") %>% select(c(1,26,"All Mapped Reads")) %>% #26 drug class, 2 is the gene

  mutate(across(c(3),as.numeric))

status <- c("caries_free", "caries_active", "caries_active", "caries_active", "caries_active", 
            "caries_active", "caries_active", "caries_free", "caries_free", "caries_free", 
            "caries_free", "caries_free", "caries_active", "caries_free", "caries_free", 
            "caries_active", "caries_free", "caries_active")

# Replace "free" with "health"
status <- gsub("free", "health", status)
unique_samples <- unique(reads$sample)
sample_status <- data.frame(sample = unique_samples, status = status)
sample_status <- sample_status %>%
  mutate(status = case_when(
    status == "caries_health" ~ "Health",
    status == "caries_active" ~ "Caries",
    TRUE ~ status
  ))
reads <- reads %>%
  left_join(sample_status, by = "sample") %>% mutate(sample = paste(status, sample, sep = "_")) %>%
  select(-status)  # Optional: remove the `status` column if it's no longer needed
reads_expanded <- reads %>%
  separate_rows(`Drug Class`, sep = ";") %>% group_by(sample,`Drug Class`) %>% summarise(count=sum(`All Mapped Reads`)) %>% 
  mutate(class=`Drug Class`) %>% 
  select(-`Drug Class`)

# New version of tidyr
wide <- reads_expanded %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = list(count = 0))

#wide<-wide %>% filter(startsWith(gene, "OXA"))
reads_mat<-as.matrix(wide[,-1])
rownames(reads_mat)<-wide$class

reads_mat <- log1p(reads_mat)
```

### Hierarchical clustering
```{r}
reads_scaled <- t(scale(t(reads_mat)))
#reads_scaled<-reads_mat
dis_genes<-dist(reads_scaled, method ="euclidean")
dis_samples<-dist(t(reads_scaled), method = "euclidean")

hc_genes <- hclust(dis_genes, method = "average")
hc_samples <- hclust(dis_samples, method = "average")

library(pheatmap)
library(RColorBrewer)

# Create the heatmap
g<-pheatmap(reads_scaled,
         cluster_rows = hc_genes,
         cluster_cols = hc_samples,
         show_rownames = TRUE,  # Hide gene names if there are too many
  fontsize_col = 8, fontsize_row = 6, # Adjust row (genes) and column (samples or antibiotics) font size
  angle_col = 45,
         
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of Category vs Samples (Hierarchical clustering)")

```

### k-means
```{r}

# Scale the remaining rows
reads_mat_scaled <- t(scale(t(reads_mat)))


# Continue with the analysis
mat_transposed <- t(reads_mat_scaled)
kmeans_result <- kmeans(mat_transposed, centers = 2)
#obtain clusters
cluster_order <- order(kmeans_result$cluster)
# Reorder the transposed matrix according to clusters
mat_clustered <- mat_transposed[cluster_order, ]
# Create an annotation for the clusters (metadata of which cluster each sample belongs to)
annotation_col <- data.frame(Cluster = factor(kmeans_result$cluster[cluster_order]))
rownames(annotation_col) <- rownames(mat_clustered)

# Create the heatmap
h <- pheatmap(
  mat_clustered,
  annotation_row = annotation_col,
  show_colnames = TRUE,  # Show column names
  fontsize_row = 8,  # Adjust row (genes) font size
  fontsize_col = 8,  # Adjust column (samples or antibiotics) font size
  angle_col = 45,    # Tilt column labels to 45 degrees
  color = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100),
  main = "Heatmap of Categories vs Samples (Grouped by K-means)"
)
```

Categories normalizing by TPKM and then log1p, result for clustering and k-means directly

```{r}
reads<-reads_raw %>% separate(col=`SAMPLE ARO Term`, into = c("sample","gene"),sep =" ", extra="merge") %>% select(c(1,2,24,26,"All Mapped Reads")) %>% #26 drug class, 2 is the gene
  mutate(across(c(3,5),as.numeric))
reads_norm<-reads %>% group_by(sample) %>%
  mutate(RPK = `All Mapped Reads` / (`Reference Length`*1000),     # Step 2: Calculate reads per length
         scaling_factor = sum(RPK) / 1e6,  # Sum RPK and scale to millions
         TPM = RPK / scaling_factor) %>%   # Step 3: Calculate TPM
  ungroup() %>% select(c(1,2,4,8))


status <- c("caries_free", "caries_active", "caries_active", "caries_active", "caries_active", 
            "caries_active", "caries_active", "caries_free", "caries_free", "caries_free", 
            "caries_free", "caries_free", "caries_active", "caries_free", "caries_free", 
            "caries_active", "caries_free", "caries_active")

# Replace "free" with "health"
status <- gsub("free", "health", status)
unique_samples <- unique(reads$sample)
sample_status <- data.frame(sample = unique_samples, status = status)
sample_status <- sample_status %>%
  mutate(status = case_when(
    status == "caries_health" ~ "Health",
    status == "caries_active" ~ "Caries",
    TRUE ~ status
  ))
##ADD SAMPLE_STATUS TO READS_NORM
reads_expanded <- reads_norm %>%
  separate_rows(`Drug Class`, sep = ";") %>% group_by(sample,`Drug Class`) %>% summarise(count=sum(tpk)) %>% 
  mutate(class=`Drug Class`) %>% 
  select(-`Drug Class`)

# New version of tidyr
wide <- reads_expanded %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = list(count = 0))

#wide<-wide %>% filter(startsWith(gene, "OXA"))
reads_mat<-as.matrix(wide[,-1])
rownames(reads_mat)<-wide$class

reads_mat <- log1p(reads_mat)

```

## By Individual genes 

TPKM

```{r}
reads<-reads_raw %>% separate(col=`SAMPLE ARO Term`, into = c("sample","gene"),sep =" ", extra="merge") %>% select(c(1,2,24,26,"All Mapped Reads")) %>% #26 drug class, 2 is the gene
  mutate(across(c(3,5),as.numeric))
reads_norm<-reads %>% mutate(rpk = `All Mapped Reads` / `Reference Length`*1000) %>% mutate(tpk=rpk/sum(rpk)*1000000) %>% select(c(1,2,7))

status <- c("caries_free", "caries_active", "caries_active", "caries_active", "caries_active", 
            "caries_active", "caries_active", "caries_free", "caries_free", "caries_free", 
            "caries_free", "caries_free", "caries_active", "caries_free", "caries_free", 
            "caries_active", "caries_free", "caries_active")

# Replace "free" with "health"
status <- gsub("free", "health", status)
unique_samples <- unique(reads$sample)
sample_status <- data.frame(sample = unique_samples, status = status)
sample_status <- sample_status %>%
  mutate(status = case_when(
    status == "caries_health" ~ "Health",
    status == "caries_active" ~ "Caries",
    TRUE ~ status
  ))
reads <- reads_norm %>%
  left_join(sample_status, by = "sample") %>% mutate(sample = paste(status, sample, sep = "_")) %>%
  select(-status)  # Optional: remove the `status` column if it's no longer needed

wide<-reads %>% 
  pivot_wider(names_from = sample, values_from = tpk, values_fill = list(tpk=0))

wide<-wide %>% filter(startsWith(gene, "OXA"))
reads_mat<-as.matrix(wide[,-1])
rownames(reads_mat)<-wide$gene

reads_mat <- log1p(reads_mat)
```
## Subsampling of table
```{r}
library(readr)
library(vegan)

lib<-read_table("reads_per_sample.tsv",   # Space delimiter
                  col_names = c("Size", "Sample"))
reads<-reads_raw %>% separate(col=`SAMPLE ARO Term`, into = c("sample","gene"),sep =" ", extra="merge") %>% select(c(1,2,24,26,"All Mapped Reads")) %>% #26 drug class, 2 is the gene
  mutate(across(c(3,5),as.numeric))

lib_processed <- lib %>%
  mutate(sample = str_extract(Sample, "^SL\\d+")) %>%  # Extract 'SL' followed by numbers
  group_by(sample) %>% 
  summarise(final=sum(Size))%>% filter(final >= 520000)
samples_to_keep <- lib_processed$sample
# Filter the 'reads' dataframe to include only the selected samples
reads_filtered <- reads %>%
  filter(sample %in% samples_to_keep)
#reads_expanded <- reads_filtered %>%
  #separate_rows(`Drug Class`, sep = ";\\s*")

reads_summarized <- reads_filtered %>%
  group_by(sample, gene) %>%
  summarise(total_reads = sum(`All Mapped Reads`), .groups = 'drop')

# Convert the data into a matrix with samples as rows and antibiotic classes as columns
reads_matrix <- reads_summarized %>%
  pivot_wider(names_from = gene, values_from = total_reads, values_fill = 0)

# Convert 'sample' into row names
reads_counts <- reads_matrix %>%
  column_to_rownames(var = "sample")
# Calculate the total reads per sample
total_counts_per_sample <- rowSums(reads_counts)

# Find the minimum total reads among the samples
min_counts <- min(total_counts_per_sample)

# Find the minimum total reads among the samples
rarefied_counts <- rrarefy(reads_counts, sample = min_counts)

rarefied_df <- as.data.frame(rarefied_counts) %>%
  rownames_to_column(var = "sample")

# Convert to long format
rarefied_long <- rarefied_df %>%
  pivot_longer(
    cols = -sample,
    names_to = "gene",
    values_to = "rarefied_counts"
  )

reads_with_rarefied <- reads %>%
  left_join(rarefied_long, by = c("sample", "gene")) %>%
  filter(sample %in% samples_to_keep)
```

Next step after subsampling

```{r}
reads_norm<-reads_with_rarefied %>% mutate(rpk = `rarefied_counts` / `Reference Length`*1000) %>% mutate(tpk=rpk/sum(rpk)*1000000) %>% select(c(1,2,4,8))
# reads_norm<-reads_with_rarefied %>% mutate(tpk = `rarefied_counts` / `Reference Length`*1000) %>% select(c(1,2,4,7))

status <- c("caries_active", "caries_active", "caries_active", "caries_active", 
            "caries_active", "caries_active", "caries_free", "caries_free", "caries_free", 
            "caries_free", "caries_free", "caries_active", "caries_free", "caries_free", 
            "caries_active", "caries_free", "caries_active")

# Replace "free" with "health"
status <- gsub("free", "health", status)
unique_samples <- unique(reads_with_rarefied$sample)
sample_status <- data.frame(sample = unique_samples, status = status)
sample_status <- sample_status %>%
  mutate(status = case_when(
    status == "caries_health" ~ "Health",
    status == "caries_active" ~ "Caries",
    TRUE ~ status
  ))
reads <- reads_norm %>%
  left_join(sample_status, by = "sample") %>% mutate(sample = paste(status, sample, sep = "_")) %>%
  select(-status)  # Optional: remove the `status` column if it's no longer needed
reads_expanded <- reads %>%
  separate_rows(`Drug Class`, sep = ";") %>% group_by(sample,`Drug Class`) %>% summarise(count=sum(tpk)) %>% 
  mutate(class=`Drug Class`) %>% 
  select(-`Drug Class`)

# New version of tidyr
wide <- reads_expanded %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = list(count = 0))

#wide<-wide %>% filter(startsWith(gene, "OXA"))
reads_mat<-as.matrix(wide[,-1])
rownames(reads_mat)<-wide$class

reads_mat <- log1p(reads_mat)
```


## Transform to RPKM (AS IN LUIS' PAPER)
```{r}
library(readr)
library(vegan)


reads<-reads_raw %>% separate(col=`SAMPLE ARO Term`, into = c("sample","gene"),sep =" ", extra="merge") %>% select(c(1,2,24,26,"All Mapped Reads")) %>% #26 drug class, 2 is the gene
  mutate(across(c(3,5),as.numeric))
lib<-read_table("reads_per_sample.tsv",   # Space delimiter
                  col_names = c("Sample", "Size"))
lib_processed <- lib %>%
  mutate(sample = str_extract(Sample, "^SL\\d+")) %>%  # Extract 'SL' followed by numbers
  group_by(sample) %>% 
  summarise(final=sum(Size)) %>% mutate(final=final/4) %>% inner_join(reads, by="sample")%>%  rename(read_count = `All Mapped Reads`) %>% rename(gene_length=`Reference Length` ) %>% 
  rename(drug_class=`Drug Class` )
# divided by four because I had used wc -l to get the number of lines

###RPKM
reads_rpkm<-lib_processed %>% mutate(fraction=read_count/(final/1000000)) %>% 
  mutate(rpkm=fraction/(gene_length/1000)) %>% select(c(1,3,5,8))

status <- c("caries_free","caries_active", "caries_active", "caries_active", "caries_active", 
            "caries_active", "caries_active", "caries_free", "caries_free", "caries_free", 
            "caries_free", "caries_free", "caries_active", "caries_free", "caries_free", 
            "caries_active", "caries_free", "caries_active")

# Replace "free" with "health"
status <- gsub("free", "health", status)
unique_samples <- unique(reads$sample)
sample_status <- data.frame(sample = unique_samples, status = status)
sample_status <- sample_status %>%
  mutate(status = case_when(
    status == "caries_health" ~ "Health",
    status == "caries_active" ~ "Caries",
    TRUE ~ status
  ))
reads <- reads_rpkm %>%
  left_join(sample_status, by = "sample") %>% mutate(sample = paste(status, sample, sep = "_")) %>%
  select(-status)

# ###BY DRUG CLASS
reads_expanded <- reads %>%
  separate_rows(drug_class, sep = ";") %>% group_by(sample,drug_class ) %>% summarise(count=sum(rpkm)) %>%
  mutate(class=drug_class) %>%
  select(-drug_class)

wide<-reads_expanded %>% 
   pivot_wider(names_from = sample, values_from = count, values_fill = list(count=0))
reads_mat<-as.matrix(wide[,-c(1)])
rownames(reads_mat)<-wide$class
###BY GENE
# wide<-reads %>% 
#   pivot_wider(names_from = sample, values_from = rpkm, values_fill = list(rpkm=0))
# wide<-wide %>% filter(startsWith(gene, "tet"))
#reads_mat<-as.matrix(wide[,-c(1,2)])
#rownames(reads_mat)<-wide$gene

reads_mat <- log1p(reads_mat)

```

## Differential ARG analysis
### TPM FOR MAASLIN2
```{r}
reads<-reads_raw %>% separate(col=`SAMPLE ARO Term`, into = c("sample","gene"),sep =" ", extra="merge") %>% select(c(1,2,24,26,"All Mapped Reads")) %>% #26 drug class, 2 is the gene
  mutate(across(c(3,5),as.numeric))
reads_norm<-reads %>% group_by(sample) %>%
  mutate(RPK = `All Mapped Reads` / (`Reference Length`*1000),     # Step 2: Calculate reads per length
         scaling_factor = sum(RPK) / 1e6,  # Sum RPK and scale to millions
         TPM = RPK / scaling_factor) %>%   # Step 3: Calculate TPM
  ungroup() %>% select(c(1,2,4,8))


status <- c("caries_free", "caries_active", "caries_active", "caries_active", "caries_active", 
            "caries_active", "caries_active", "caries_free", "caries_free", "caries_free", 
            "caries_free", "caries_free", "caries_active", "caries_free", "caries_free", 
            "caries_active", "caries_free", "caries_active")

# Replace "free" with "health"
status <- gsub("free", "health", status)
unique_samples <- unique(reads$sample)
sample_status <- data.frame(sample = unique_samples, status = status)
sample_status <- sample_status %>%
  mutate(status = case_when(
    status == "caries_health" ~ "Health",
    status == "caries_active" ~ "Caries",
    TRUE ~ status
  ))
# Assign the row names of `sample_status` from the `sample` column
rownames(sample_status) <- sample_status$sample
# Remove the 'sample' column since it's now in the row names
sample_status$sample <- NULL


####FROM READS
wide_genes<-reads_norm %>% select(-`Drug Class`) %>% 
  pivot_wider(names_from = sample, values_from = TPM, values_fill = list(TPM=0))
df_genes <- as.data.frame(t(wide_genes[,-1])) # Transpose excluding the gene name column
colnames(df_genes) <- wide_genes$gene    # Assign gene names as column names
           # Remove row names
colnames(df_genes) <- gsub(" ", "_", colnames(df_genes))

####FROM CATEGORIES
reads_class <- reads_norm %>%
  separate_rows(`Drug Class`, sep = ";") %>% 
mutate(`Drug Class`= str_trim(`Drug Class`))%>% group_by(sample,`Drug Class`) %>% summarise(count=sum(TPM)) %>% 
  mutate(class=`Drug Class`) %>% 
  select(-`Drug Class`)

# New version of tidyr
wide_class <- reads_class %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = list(count = 0))

df_class <- as.data.frame(t(wide_class[,-1])) # Transpose excluding the gene name column
colnames(df_class) <- wide_class$class    # Assign gene names as column names
           # Remove row names
colnames(df_class) <- gsub(" ", "_", colnames(df_class))


```
### RUN MAASLIN2
```{r}
library(Maaslin2)
fit_data=Maaslin2(input_data = df_genes,
                  input_metadata = sample_status,
                      output = "Maaslin2_rawreads_genes_tpm",       # Output folder for results
                  fixed_effects = c("status"), 
                  min_prevalence = 0.2,
                  normalization="TSS",
                  transform = "LOG",
                  analysis_method = "LM")
library(Maaslin2)
fit_data=Maaslin2(input_data = df_class,
                  input_metadata = sample_status,
                      output = "Maaslin2_rawreads_class_tpm",       # Output folder for results
                  fixed_effects = c("status"), 
                  min_prevalence = 0.2,
                  normalization="TSS",
                  transform = "LOG",
                  analysis_method = "LM")
```


 
