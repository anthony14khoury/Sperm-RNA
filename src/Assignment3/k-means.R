# K-Means Clustering

library(tidyverse)
library(M3Drop)
library(ggplot2)
library(cluster)


# File Directories
data_dir = file.path("Raw Data", "SRP053246")
data_file = file.path("CleanData", "CleanedData.tsv")
metadata_file = file.path(data_dir, "metadata_SRP053246.tsv")

# Load Data
metadata <- readr::read_tsv(metadata_file)
expression_df <- readr::read_tsv(data_file)%>%
  tibble::column_to_rownames("Gene")

# Ensuring Data and Metadata are in the Same Order
expression_df <- expression_df %>% 
  dplyr::select(metadata$refinebio_accession_code)
all.equal(colnames(expression_df), metadata$refinebio_accession_code)


# 5000 most variable genes
var_genes <- read.csv("CleanData/5000VariableGenes.csv")

# Removing Gene Names
genes_filtered <- subset(var_genes, select = -c(Gene))


# Find optimum number of clusters using Average Silhouette Width
sil <- rep(0, 20)
for (i in 2:20) {
  k_1to20 <- kmeans(genes_filtered, centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k_1to20$cluster, dist(genes_filtered))
  sil[i] <- mean(ss[, 3])
}

# Clustering Results
cat("Optimal Number ff Clusters:", which.max(sil), "\n")
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)


