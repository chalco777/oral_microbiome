# This script processes a taxonomic count matrix coming from Pavian for LEfSe analysis.
# It cleans the data, rarefies counts to the minimum sample depth,
# and formats the output for LEfSe, including sample status.

library(tidyverse)
library(vegan)

# Set working directory
setwd("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/metmon/lefse/lefse_krakenout_bracken")

# Load count matrix
counts <- read.delim("matrix_allranks_conteo.tsv", sep = "\t")
str(counts)
table(counts$taxRank)

# Remove unnecessary columns (2, 3, 4, and last column)
counts <- counts[, -c(2, 3, 4, ncol(counts))]
dim(counts)

# Convert all columns except the first to numeric
counts[,-1] <- sapply(counts[,-1], as.numeric)

# Replace NA values with 0
counts[is.na(counts)] <- 0
sapply(lapply(counts, is.na), sum)

# Optional: If working with fractions, divide all numeric columns by 100
# counts[sapply(counts, is.numeric)] <- sapply(counts[sapply(counts, is.numeric)], function(x) x / 100)

# Replace spaces with underscores in taxon names
counts$name <- gsub(" ", "_", counts$name)

# Rename columns to keep only the ID before the first underscore
names(counts) <- sub("_.*", "", names(counts))

# Save pre-rarefaction counts
pre_rarefaction_counts <- counts

# Rarefaction
count_matrix <- counts[, -1]
rownames(count_matrix) <- make.unique(counts$name)
count_matrix_t <- t(count_matrix)
min_count <- min(rowSums(count_matrix_t))
rarefied_matrix <- rrarefy(count_matrix_t, sample = min_count)
rarefied_matrix <- t(rarefied_matrix)
rarefied_df <- as.data.frame(rarefied_matrix)

# Define sample status (first row for LEfSe format)
status <- c("status", "caries_free", "caries_active", "caries_active", "caries_active", "caries_active", 
            "caries_active", "caries_active", "caries_free", "caries_free", "caries_free", "caries_free", 
            "caries_free", "caries_active", "caries_free", "caries_free", "caries_active", "caries_free", 
            "caries_active")
status <- as.data.frame(t(status))
names(status) <- colnames(counts)

# Add row names as a new column
rarefied_df$name <- rownames(rarefied_df)

# Reorder columns so 'name' is first
rarefied_df <- rarefied_df[, c("name", colnames(rarefied_df)[-ncol(rarefied_df)])]

# Combine status and rarefied data
final_df <- rbind(status, rarefied_df)
dim(final_df)
any(is.na(final_df))
sum(sapply(final_df, is.infinite))

# Write output for LEfSe
write.table(final_df, file = "lefse_fullranks_counts.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
