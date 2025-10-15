---
title: "Differential_taxonomic_analysis"
author: "Adrián Chalco"
date: "2024-12-22"
output: html_document
---

This notebook performs differential taxonomic analysis of the oral microbiome by generating barplots and boxplots of relative abundances and raw counts (before and after rarefaction), identifying the top species associated with caries and health, and preparing input files for LEfSe analysis.

## Directory setup

```{r}
library(tidyverse)
library(ggtext)
library(vegan)
setwd("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/metmon")
# Performing LefSe at all taxonomic ranks and rarefaction
res<-read.delim("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/metmon/lefse/lefse_krakenout_bracken/matrix_allranks_conteo.tsv",sep = "\t")
# RELATIVE ABUNDANCE TABLE
rela<-read_tsv("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/metmon/Stackbar_abundance/species_rel_abundance.tsv")
```

## Processing relative abundance table
```{r}
names(rela)<-sub("\\..*","",names(rela))
names(rela)[-c(1:4,ncol(rela))]<-sub("$","e",names(rela)[-c(1:4,ncol(rela))])
# from NAs to 0
rela[is.na(rela)]<-0
```

## Barplot of relative abundance
```{r}
rel_abundance<-rela %>% pivot_longer(cols=starts_with("SL"),names_to = "sample_id2", values_to = "rel_abundance") %>% 
  select(sample_id2,rel_abundance, taxRank,name) %>%
  mutate(sample_id=sub("_.*","",sample_id2),
            status=sub(".*_","",sample_id2)) %>% select(sample_id, status, everything(),-sample_id2) %>% 
  rename(taxon=name) %>% 
  arrange(sample_id) %>% 
  mutate(sample_id=factor(sample_id, 
                          levels=c("SL021", "SL085", "SL086", "SL088", "SL090", "SL103", "SL112", "SL114", "SL125","SL036", "SL037", "SL038", "SL048", "SL050", "SL068", "SL111", "SL118", "SL200")))

# FILTER BY TAXONOMIC RANK
rel_abundance_dos<-rel_abundance %>% 
  filter(taxRank=="S") %>% 
  group_by(status,sample_id,taxon)%>% summarize(rel2=sum(rel_abundance), .groups="drop") 
# SELECT TOP_20
top_20<-rel_abundance_dos %>%
  group_by(taxon) %>%
  summarize(total_abundance = median(rel2)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:20) %>%  # Select the top 20 species
  pull(taxon)

# GROUP Others AND EDIT TAXON NAMES TO DISPLAY NICELY IN THE PLOT
rel_abundance_dos_edit<-rel_abundance_dos %>%   
  mutate(taxon = ifelse(taxon %in% top_20, taxon, "Others")) %>% 
  group_by(status,sample_id,taxon)%>% summarize(rel=sum(rel2), .groups="drop") %>% # merge the Others
  mutate(taxon=str_replace(taxon,
                           "^(?!Others$).*","*\\0*"))  %>% 
  mutate(taxon = reorder(taxon, -rel)) %>%
  mutate(taxon = factor(taxon, levels = c(levels(taxon)[levels(taxon) != "Others"], "Others")))

# check that the sum is 100 for each taxon
rel_abundance_dos_edit %>% group_by(sample_id) %>% summarise(total=sum(rel))

custom_palette <- c("#3C8A9B", "#E17597", "#D4AE2D", "#E884B9", "#54A453", "#BD6066", "#FF990A", 
                    "#469F6C", "#FFCF20", "#526E9F", "#E97422", "#94539E", "#9B445D", 
                    "#C08EA9", "#999999","#AF6729", "#BF6357", "#747B78", "#FAF632", "#E41A1C", "gray")
# Plot
rel_abundance_dos_edit %>%
  group_by(status, sample_id) %>%
  mutate(abundancia_Others = sum(rel[taxon == "Others"])) %>%  # Calculate the abundance of "Others"
  ungroup() %>%
  ggplot(aes(x = reorder(sample_id, abundancia_Others), y = rel, fill = taxon))+ ## order according to abundance of Others
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = custom_palette, name = NULL) + # remove the legend title for fill colors
  #+  scale_x_discrete(breaks=c(), labels=c(<br>*italic*<br>)) change bar names on x axis and use ggtext (<br> is line break)
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1),
        legend.text=element_markdown(),
        legend.key.size=unit(10,"pt"))+
  labs(y = "Relative abundance (%)", x = NULL)+
  facet_wrap(~ status,scales = "free_x", nrow = 1,     labeller = as_labeller(c("active" = "Caries", "free" = "Health"))
) # each facet with its own adjusted X axis and both in the same row

ggsave(
  filename = "relative_abundance_plot.png", # File name
  width = 12, # Width in inches
  height = 6, # Height in inches
  dpi = 300,   # Resolution in DPI
  bg="white"

  )
```


## Processing table with raw counts
```{r}
str(res)
table(res$taxRank)
res<-res[,-c(2,3,4,ncol(res))]
dim(res)
# convert to numeric
res[,-1]<-sapply(res[,-1],as.numeric)
# change NA to 0.0
res[is.na(res)]<-0
sapply(lapply(res,is.na),sum)
# change column name to only the ID
names(res)<-sub("_.*","",names(res))

```

## Rarefaction of counts to minimum sample size to then obtain the lefse input file:
Obtain lefse_fullranks_counts.tsv file from dataframe res
Lefse directory: 
```{r}
# connect names by underscore
res$name<-gsub(" ","_",res$name)

# another option: rename_with(~ sub("_.*", "", .))
## RAREFACTION
set.seed(100)
count_matrix <- res[, -1]
rownames(count_matrix) <- make.unique(res$name)
count_matrix_t <- t(count_matrix)
min_count <- min(rowSums(count_matrix_t))
rarefied_matrix <- rrarefy(count_matrix_t, sample = min_count)
rarefied_matrix <- t(rarefied_matrix)
rarefied_df <- as.data.frame(rarefied_matrix)
status<-c("status", "caries_free", "caries_active", "caries_active", "caries_active", "caries_active", 
          "caries_active", "caries_active", "caries_free", "caries_free", "caries_free", "caries_free", 
          "caries_free", "caries_active", "caries_free", "caries_free", "caries_active", "caries_free", 
          "caries_active")
status<-as.data.frame(t(status))
names(status)<-colnames(res)

# Add row names as a new column
rarefied_df$name <- rownames(rarefied_df)

# Reorder columns so that 'name' is first
rarefied_df <- rarefied_df[, c("name", colnames(rarefied_df)[-ncol(rarefied_df)])]

df<-rbind(status,rarefied_df)
dim(df)
any(is.na(df))
sum(sapply(df,is.infinite))

write.table(df,file="lefse_fullranks_counts.tsv",sep="\t",row.names=FALSE,quote = FALSE)
```


## Boxplot of relative abundance before rarefaction

```{r}
especies_caries <- c(
  "Leptotrichia sp. HMT-225",
  "Leptotrichia trevisanii",
  "Leptotrichia hofstadii",
  "Leptotrichia wadei",
  "Leptotrichia hongkongensis",
  "Schaalia sp. HMT-172",
  "Streptococcus mutans",
  "Campylobacter rectus",
  "Selenomonas sputigena",
  "Leptotrichia sp. oral taxon 847"
)
rela<-rela %>% filter(name %in% especies_caries) %>% as.data.frame()
rela<-rela[,-c(2:4,23)]
names(rela)<-c("name",sub("(^[^_]*).*","\\1",names(rela)[-1]))
data_long <- rela %>%
  pivot_longer(cols = -name, names_to = "Sample", values_to = "Value") %>%
  left_join(status %>% pivot_longer(cols = -name, names_to = "Sample", values_to = "Status"), by = "Sample")
# Check the result to make sure the transformation is correct
str(data_long)
medianas <- data_long %>%
  group_by(name.x) %>%
  summarise(median_value = median(Value, na.rm = TRUE))
data_long <- data_long %>%
  mutate(name.x = factor(name.x, levels = medianas$name.x[order(medianas$median_value)]))

# Suppose `data_long_abundance` is the name of the data frame containing the relative abundance
ggplot(data_long, aes(x = name.x, y = Value, fill = Status)) +
  geom_boxplot(alpha = 0.7, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA) +  # Boxplot without visible outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), 
              size = 1.5, alpha = 0.8, aes(color = Status))  +  # Points instead of jitter, with border
   scale_fill_manual(
    values = c("caries_free" = "#66C2A5", "caries_active" = "#FC8D62"),
    labels = c("caries_free" = "Health", "caries_active" = "Caries")  # Modify legend labels
  ) +
  scale_color_manual(
    values = c("caries_free" = "#66C2A5", "caries_active" = "#FC8D62"),
    labels = c("caries_free" = "Health", "caries_active" = "Caries")  # Modify legend labels
  )+
  labs(x = "Species", y = "Relative Abundance", fill = "Status", color = "Status",
       title = "Relative Abundance of Bacterial Species Associated with Caries") +  # Labels and title
  theme_minimal(base_size = 14) +  # Minimalist theme with base text size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),  # Rotate x labels
    axis.text.y = element_text(size = 10),  # y text size
    axis.title = element_text(size = 14, face = "bold"),  # Axis title size and bold
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered title
    legend.title = element_text(size = 14),  # Legend title size
    legend.text = element_text(size = 12)  # Legend text size
  ) +
  scale_y_log10(breaks = c(0,0.00005,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05)*100,
                labels = paste0(c(0,0.00005,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05)*100,"%"))

ggsave("boxplot_abundance_diffspecies_caries_beforerare.png", width = 10, height = 6, dpi = 300, bg="white")



```

## Boxplot of absolute count before rarefaction

```{r}
### ABSOLUTE COUNT
res2<-res %>% filter(name %in% gsub("\\s","_",especies_caries))
data_long <- res2%>% 
  pivot_longer(cols = -name, names_to = "Sample", values_to = "Value") %>%
  left_join(status %>% pivot_longer(cols = -name, names_to = "Sample", values_to = "Status"), by = "Sample")
# Check the result to make sure the transformation is correct
str(data_long)
data_long$Value <- as.numeric(data_long$Value)
# group by median
# Calculate the median of each group
medianas <- data_long %>%
  group_by(name.x) %>%
  summarise(median_value = median(Value, na.rm = TRUE))

data_long <- data_long %>%
  mutate(name.x = factor(name.x, levels = medianas$name.x[order(medianas$median_value)]))

ggplot(data_long, aes(x = name.x, y = Value, fill = Status)) +
  geom_boxplot(alpha = 0.7, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA) +  # Adjustment to avoid overlap
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), 
              size = 1.5, alpha = 0.8, aes(color = Status)) +  # Jitter adjustment
    scale_fill_manual(
    values = c("caries_free" = "#56B4E9", "caries_active" = "#D55E00"),
    labels = c("caries_free" = "Health", "caries_active" = "Caries")  # Modify legend labels
  ) +
  scale_color_manual(
    values = c("caries_free" = "#56B4E9", "caries_active" ="#D55E00"),
    labels = c("caries_free" = "Health", "caries_active" = "Caries")  # Modify legend labels
  )+
  labs(x = "Species", y = "Read Count", fill = "Status", color = "Status",
       title = "Read Count of Bacterial Species Associated with Caries") +  # Labels and title
  theme_minimal(base_size = 14) +  # Minimalist theme with base text size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),  # Rotate x labels
    axis.text.y = element_text(size = 12),  # y text size
    axis.title = element_text(size = 14, face = "bold"),  # Axis title size and bold
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered title
    legend.title = element_text(size = 14),  # Legend title size
    legend.text = element_text(size = 12)  # Legend text size
  ) +# Ensure the legend is displayed correctly
  scale_y_log10(
    breaks = c(1, 10, 100, 200,500,2000, 5000,20000,50000),
    labels = c(0, 10, 100,200,500,2000, 5000,20000,50000)
  )+  scale_x_discrete(labels = function(x) gsub("_", " ", x))
ggsave("boxplot_count_diffspecies_caries_final.png", width = 10, height = 6, dpi = 300, bg="white")

```

## Boxplot of absolute count after rarefaction (LefSe)

```{r}
rownames(rarefied_df)<-NULL
### ABSOLUTE COUNT
rarefied_df2<-rarefied_df %>% filter(name %in% gsub("\\s","_",especies_caries))
data_long <- rarefied_df2%>%
  pivot_longer(cols = -name, names_to = "Sample", values_to = "Value") %>%
  left_join(status %>% pivot_longer(cols = -name, names_to = "Sample", values_to = "Status"), by = "Sample")
# Check the result to make sure the transformation is correct
str(data_long)
data_long$Value <- as.numeric(data_long$Value)
# group by median
# Calculate the median of each group
medianas <- data_long %>%
  group_by(name.x) %>%
  summarise(median_value = median(Value, na.rm = TRUE))

data_long <- data_long %>%
  mutate(name.x = factor(name.x, levels = medianas$name.x[order(medianas$median_value)]))

ggplot(data_long, aes(x = name.x, y = Value, fill = Status)) +
  geom_boxplot(alpha = 0.7, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA) +  # Adjustment to avoid overlap
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), 
              size = 1.5, alpha = 0.8, aes(color = Status)) +  # Jitter adjustment
    scale_fill_manual(
    values = c("caries_free" = "#56B4E9", "caries_active" = "#D55E00"),
    labels = c("caries_free" = "Health", "caries_active" = "Caries")  # Modify legend labels
  ) +
  scale_color_manual(
    values = c("caries_free" = "#56B4E9", "caries_active" ="#D55E00"),
    labels = c("caries_free" = "Health", "caries_active" = "Caries")  # Modify legend labels
  )+
  labs(x = "Species", y = "Read Count", fill = "Status", color = "Status",
       title = "Read Count of Bacterial Species Associated with Caries (after rarefaction)") +  # Labels and title
  theme_minimal(base_size = 14) +  # Minimalist theme with base text size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),  # Rotate x labels
    axis.text.y = element_text(size = 12),  # y text size
    axis.title = element_text(size = 14, face = "bold"),  # Axis title size and bold
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered title
    legend.title = element_text(size = 14),  # Legend title size
    legend.text = element_text(size = 12)  # Legend text size
  ) +# Ensure the legend is displayed correctly
  scale_y_log10(
    breaks = c(1, 10, 100, 200,500,2000, 5000,20000,50000),
    labels = c(0, 10, 100,200,500,2000, 5000,20000,50000)
  )
ggsave("boxplot_count_differential_species_after_raref.png", width = 10, height = 6, dpi = 300, bg="white")


```
## Boxplot of relative abundance after rarefaction (LefSe)

```{r}
totals <- as.numeric(df[df$name == "root", -1]) # Exclude the 'name' column to get only the totals
sample_names <- names(rarefied_df2)[-1] # Exclude the 'name' column
# Create a data frame to store relative abundances
relative_abundance <- rarefied_df2
# Calculate relative abundance (as percentage)
for (sample in sample_names) {
  relative_abundance[[sample]] <- (rarefied_df2[[sample]] / totals[which(sample_names == sample)]) * 100
}
# View the result
print(relative_abundance)
data_long <- relative_abundance %>%
  pivot_longer(cols = -name, names_to = "Sample", values_to = "Value") %>%
  left_join(status %>% pivot_longer(cols = -name, names_to = "Sample", values_to = "Status"), by = "Sample")
# Check the result to make sure the transformation is correct
str(data_long)
medianas <- data_long %>%
  group_by(name.x) %>%
  summarise(median_value = median(Value, na.rm = TRUE))
data_long <- data_long %>%
  mutate(name.x = factor(name.x, levels = medianas$name.x[order(medianas$median_value)]))

# Suppose `data_long_abundance` is the name of the data frame containing the relative abundance
ggplot(data_long, aes(x = name.x, y = Value, fill = Status)) +
  geom_boxplot(alpha = 0.7, color = "black", position = position_dodge(width = 0.75), outlier.shape = NA) +  # Boxplot without visible outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), 
              size = 1.5, alpha = 0.8, aes(color = Status))  +  # Points instead of jitter, with border
   scale_fill_manual(
    values = c("caries_free" = "#66C2A5", "caries_active" = "#FC8D62"),
    labels = c("caries_free" = "Health", "caries_active" = "Caries")  # Modify legend labels
  ) +
  scale_color_manual(
    values = c("caries_free" = "#66C2A5", "caries_active" = "#FC8D62"),
    labels = c("caries_free" = "Health", "caries_active" = "Caries")  # Modify legend labels
  )+
  labs(x = "Species", y = "Relative Abundance", fill = "Status", color = "Status",
       title = "Relative Abundance of Bacterial Species Associated with Caries (after rarefaction)") +  # Labels and title
  theme_minimal(base_size = 14) +  # Minimalist theme with base text size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),  # Rotate x labels
    axis.text.y = element_text(size = 10),  # y text size
    axis.title = element_text(size = 14, face = "bold"),  # Axis title size and bold
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered title
    legend.title = element_text(size = 14),  # Legend title size
    legend.text = element_text(size = 12)  # Legend text size
  ) +
  scale_y_log10(breaks = c(0,0.00005,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05)*100,
                labels = paste0(c(0,0.00005,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05)*100,"%"))+  scale_x_discrete(labels = function(x) gsub("_", " ", x))

ggsave("boxplot_abundance_diffspecies_caries_afterrerare.png", width = 10, height = 6, dpi = 300, bg="white")




```

