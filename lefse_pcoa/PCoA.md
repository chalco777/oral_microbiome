---
title: "PCoA"
author: "Adrián Chalco"
date: "2025-01-12"
output: html_document
---

In this script Principal Coordinates Analysis (PCoA) on rarefied abundance data is performed, integrating metadata to visualize sample clustering by caries status, severity and others, and computing correlations between individual species and PCoA axes. The output are annotated plots and summary tables for interpretation.

```{r}

library(tidyverse)
library(vegan)
library(glue)
library(ggrepel)
library(openxlsx)
library(knitr)
library(kableExtra)
library(RColorBrewer)

```
## INPUT: SYLPH
```{r}
# df_sylph<-read_tsv("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/metmon/sylph_results/results_unknown.tsv") %>% 
#   rename(name=Contig_name) %>% 
#   mutate(sample = gsub("^(.*?)_.*", "\\1", Sample_file)) %>% # the ? makes it non-greedy and only takes up to what is before the FIRST "_", otherwise it would take everything up to the last "_"
#   select(c(16,4,15)) %>% 
#   mutate(conteo=(Sequence_abundance*909628)/100) %>% 
#   select(-c("Sequence_abundance")) %>% 
#   arrange(sample) %>% 
#   pivot_wider(names_from = sample,values_from = conteo,
#               values_fill = 0)
# df<-df_sylph
# df$name<-gsub(" ","_",df$name)
# rarefied_matrix<-df[, -c(1)]
# rarefied_matrix<-t(rarefied_matrix)
```
## INPUT: BRAKEN->PAVIAN COUNTS: MATRIX ALL_RANKS
Rarefaction
```{r}
df<-read_tsv("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/metmon/lefse/lefse_krakenout_bracken/matrix_allranks_conteo.tsv")
df<-df %>% filter(taxRank=="S")
df<-df[,-c(2,3,4,ncol(df))]
#change column name to only the ID
names(df)<-sub("_.*","",names(df))
df[,-1]<-sapply(df[,-1],as.numeric)
df[is.na(df)]<-0
#connect df name by underscore
df$name<-gsub(" ","_",df$name)
df<-as.data.frame(df)
##RAREFACTION OF BRACKEN
set.seed(100)
count_matrix <- df[, -c(1,2)] #remove the smallest one SL021
rownames(count_matrix) <- make.unique(df$name)
count_matrix_t <- t(count_matrix)
min_count <- min(rowSums(count_matrix_t))
rarefied_matrix <- rrarefy(count_matrix_t, sample = min_count)
rarefied_matrix <- as.data.frame(rarefied_matrix)

```
## PCoA
```{r}
#Alternatively: Additional filters
#rarefied_matrix <- rarefied_matrix[, which(colSums(rarefied_matrix) > 30)]
#rarefied_matrix <- rarefied_matrix[, apply(rarefied_matrix, 2, function(col) sum(col > 0) >= 10)]
#Convert to numeric matrix
rarefied_matrix<-as.matrix(rarefied_matrix)
log_rarefied_matrix<-log2(rarefied_matrix+1)
#Calculate Bray-Curtis distance
dist<-vegdist(rarefied_matrix, method = "bray")
dist<-vegdist(log_rarefied_matrix, method = "bray")

# Perform PCoA with cmdscale acting on a (dist_matrix), add correction and get eigen values
pcoa <- cmdscale(dist, k = 2, eig = TRUE, add=TRUE)
positions <- pcoa$points #extract the vectors
colnames(positions) <- c("pcoa1", "pcoa2") #assign names to coordinates
positions<-positions[order(rownames(positions)), ] #sort by row name
#calculate explained percentage
expl<-(pcoa$eig / sum(pcoa$eig))*100
exp<-format(round(expl[1:2],digits=1),nsmall=1, trim=TRUE)
labs<-c(glue("PCo 1 ({exp[1]}%)"),
  glue("PCo 2 ({exp[2]}%)"))
```
## Getting metadata
```{r}
###METADATA
metadata<-read.xlsx("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/metmon/pcoa/metadata.xlsx")
##ADD CARIES FREE IF SL021 IS CONSIDERED
status <- c("caries_active", "caries_active", "caries_active", "caries_active", 
      "caries_active", "caries_active", "caries_free", "caries_free", "caries_free", 
      "caries_free", "caries_free", "caries_active", "caries_free", "caries_free", 
      "caries_active", "caries_free", "caries_active")
# Create a data frame where rows are samples
status_df <- data.frame(
  samples = rownames(positions),  # Samples come from rownames of positions
  status = status,                # The status vector
  stringsAsFactors = FALSE
)
colnames(metadata)[1]<-'samples'
metadata_dos<-inner_join(status_df, metadata, by='samples')
metadata_dos <- metadata_dos %>% 
  mutate(grado = if_else(samples %in% c("SL038", "SL068", "SL200"), "very_severe",
       if_else(status == "caries_active", "severe", "no_caries")))
metadata_dos$edad_ni <- factor(gsub("[^0-9]", "", metadata_dos$edad_ni))


```

## PCoA of Microbial Species: Caries, Severe Caries, and Health
```{r}

# Convert the ordered PCoA results to tibble and join with df_variable
# PCoA plot with improvements and sample labels
gg1<-positions %>%
  as_tibble(rownames = "samples") %>%
  left_join(metadata_dos, by = "samples") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = grado)) +
  geom_point(size = 2.5, show.legend = TRUE)+
  #geom_jitter(size = 2.5, stroke = 0.8,width = 0.005, height = 0.005) +  # Place geom_point without size in aes()
  geom_text_repel(aes(label = samples), size = 3, max.overlaps = 2, force=2)+#tras
  # Add dashed lines at the origin (x=0, y=0)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Add ellipses for each "grado" group
  stat_ellipse(aes(group = grado),
         type = "norm",       # You can change to "t" or "euclid" as you prefer
         linetype = "solid",  # Or "dashed", etc.
         alpha = 0.5,         # Transparency for the ellipse
         size = 0.8,
         na.rm = TRUE  # Ignore groups with fewer points
  ) +          # Line thickness of the ellipse
  
    labs(x = labs[1], y = labs[2]) +  
  theme_minimal(base_size = 16) +  
  theme(
    axis.title = element_text(size = 14, face = "bold"),  
    axis.text = element_text(size = 14),  
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  
  ) +
  scale_color_manual(
    values = c("severe" = "#E41A1C", 
         "very_severe" = "#984EA3", 
         "no_caries" = "#377EB8"),
    name = "Condition",
    breaks = c("very_severe", "severe", "no_caries"),
    labels = c("Severe caries", "Caries", "Health")
  )+  
  ggtitle("PCoA of Microbial Species: Caries, Severe Caries, and Health")
#+ xlim(c(-0.105, -0.095))

# set1_colors
# [1] "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3"
# [5] "#FF7F00" "#FFFF33" "#A65628" "#F781BF"
# [9] "#999999"

ggsave("PCoA_species_3.png", plot = gg1, width = 8, height = 6, dpi = 300,bg="white")



```
## PCoA of Microbial Species: Caries vs. Health
```{r}

gg2<-positions %>%
  as_tibble(rownames = "samples") %>%
  left_join(metadata_dos, by = "samples") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = status)) +
  geom_point(size = 2.5, show.legend = TRUE)+
  #geom_jitter(size = 2.5, stroke = 0.8,width = 0.005, height = 0.005) +  # Place geom_point without size in aes()
  geom_text_repel(aes(label = samples), size = 3, max.overlaps = 2, force=2)+#tras
  # Add dashed lines at the origin (x=0, y=0)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Add ellipses for each "status" group
  stat_ellipse(aes(group = status), 
         type = "norm",       # You can change to "t" or "euclid" as you prefer
         linetype = "solid",  # Or "dashed", etc.
         alpha = 0.5,         # Transparency for the ellipse
         size = 0.8) +          # Line thickness of the ellipse
  
  labs(x = labs[1], y = labs[2]) +  
  theme_minimal(base_size = 16) +  
  theme(
    axis.title = element_text(size = 14, face = "bold"),  
    axis.text = element_text(size = 14),  
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  
  ) +
  scale_color_brewer(
    palette = "Set1",
    name = "Condition",                   # Legend title
    breaks = c("caries_active",        # Original values in 'status'
         "caries_free"),
    labels = c("Caries", "Health")     # New labels in the legend
  )+  
  ggtitle("PCoA of Microbial Species: Caries vs. Health")
#+ xlim(c(-0.105, -0.095))

ggsave("PCoA_species_2.png", plot = gg2, width = 8, height = 6, dpi = 300, bg="white")


```
## Correlation of species with PCoA Axes
```{r}
# 1. Prepare data frames of positions and check the order of the samples
positions_df <- as.data.frame(positions)
positions_df$samples <- rownames(positions_df)

# Sort the rows of 'positions_df' and 'log_rarefied_matrix' according to the sample name
positions_df <- positions_df[order(positions_df$samples), ]
log_rarefied_matrix <- log_rarefied_matrix[order(rownames(log_rarefied_matrix)), ]

# Check for name match
if(!all(rownames(log_rarefied_matrix) == positions_df$samples)){
  stop("Samples in log_rarefied_matrix and positions_df do not match in order/names.")
}

# 2. Calculate the correlation species vs. PCoA1 and PCoA2

cor_list <- lapply(colnames(log_rarefied_matrix), function(sp){
  # Abundances (log) of the species
  abund_sp <- log_rarefied_matrix[, sp]
  
  # Correlation with pcoa1
  cor_pcoa1 <- cor(abund_sp, positions_df$pcoa1, method = "kendall")
  cor_pcoa1_test <- cor.test(abund_sp, positions_df$pcoa1, method = "kendall")
  pval_pcoa1 <- cor_pcoa1_test$p.value
  
  # Correlation with pcoa2
  cor_pcoa2 <- cor(abund_sp, positions_df$pcoa2, method = "spearman")
  cor_pcoa2_test <- cor.test(abund_sp, positions_df$pcoa2, method = "spearman")
  pval_pcoa2 <- cor_pcoa2_test$p.value
  
  data.frame(
    species = sp,
    cor_pcoa1 = cor_pcoa1,
    pval_pcoa1 = pval_pcoa1,
    cor_pcoa2 = cor_pcoa2,
    pval_pcoa2 = pval_pcoa2
  )
})

# 3. Combine the results into a single data frame
cor_results <- do.call(rbind, cor_list)
# For PCoA1
cor_results_pcoa1 <- cor_results %>%
  select(species, cor_pcoa1, pval_pcoa1) %>%
  mutate(abs_cor_pcoa1 = abs(cor_pcoa1)) %>%
  arrange(desc(abs_cor_pcoa1)) %>%
  rename(
    Spearman_Correlation = cor_pcoa1,
    P_value              = pval_pcoa1,
    Absolute_Correlation = abs_cor_pcoa1
  )

# For PCoA2
cor_results_pcoa2 <- cor_results %>%
  select(species, cor_pcoa2, pval_pcoa2) %>%
  mutate(abs_cor_pcoa2 = abs(cor_pcoa2)) %>%
  arrange(desc(abs_cor_pcoa2)) %>%
  rename(
    Spearman_Correlation = cor_pcoa2,
    P_value              = pval_pcoa2,
    Absolute_Correlation = abs_cor_pcoa2
  )

cor_results_pcoa1_final <- cor_results_pcoa1 %>%
  select(-Absolute_Correlation)  # used for display/export

cor_results_pcoa2_final <- cor_results_pcoa2 %>%
  select(-Absolute_Correlation)

pcoa1_html <- cor_results_pcoa1_final %>%
  head(10) %>%
    rename_with(~ gsub("_", " ", .x)) %>%  # Replace "_" with space

  kbl(
    format  = "html",
    caption = "Spearman correlation with PCoA1 (TOP 10)",
    # Optional: change column width, alignment, etc. here
  ) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

pcoa2_html <- cor_results_pcoa2_final %>%
  head(10) %>%
    rename_with(~ gsub("_", " ", .x)) %>%  # Replace "_" with space

  kbl(
    format  = "html",
    caption = "Spearman correlation with PCoA2 (TOP 10)"
  ) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

# Save the HTML tables if desired
save_kable(pcoa1_html, "pcoa1_top10.html")
save_kable(pcoa2_html, "pcoa2_top10.html")

###############################################################################
# 4. Export the COMPLETE tables in CSV (without the abs. corr. column)
###############################################################################
write.csv(cor_results_pcoa1_final, "pcoa1_full.csv", row.names = FALSE)
write.csv(cor_results_pcoa2_final, "pcoa2_full.csv", row.names = FALSE)

message("Process completed: HTML tables (TOP 10) and full CSVs (without abs. corr.) have been generated.")
```


