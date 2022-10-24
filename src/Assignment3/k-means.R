# K-Means Clustering

library(tidyverse)
library(M3Drop)
library(ggplot2)
library(cluster)
library(pheatmap)


# ---------------------------------
print("File Directories")
data_dir = file.path("Raw Data", "SRP053246")
data_file = file.path("CleanData", "CleanedData.tsv")
metadata_file = file.path(data_dir, "metadata_SRP053246.tsv")


# ---------------------------------
print("Loading the Data")
metadata <- readr::read_tsv(metadata_file)
expression_df <- readr::read_tsv(data_file)%>%
  tibble::column_to_rownames("Gene")


# ---------------------------------
print("Ensuring Data and Metadata are in the Same Order")
expression_df <- expression_df %>% 
  dplyr::select(metadata$refinebio_accession_code)
all.equal(colnames(expression_df), metadata$refinebio_accession_code)


# ---------------------------------
print("Reading in the 5000 most variable genes")
var_genes <- read.csv("CleanData/5000VariableGenes.csv")
genes_filtered <- subset(var_genes, select = -c(Gene)) # Removing Gene Names


# ---------------------------------
# Find optimum number of clusters using Average Silhouette Width
# Uncomment below to run it, but previously calculated optimal is 2 clusters
# sil <- rep(0, 20)
# for (i in 2:20) {
#   k_1to20 <- kmeans(genes_filtered, centers = i, nstart = 25, iter.max = 20)
#   ss <- silhouette(k_1to20$cluster, dist(genes_filtered))
#   sil[i] <- mean(ss[, 3])
# }

# Clustering Results
# cat("Optimal Number ff Clusters:", which.max(sil), "\n")
# plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
# abline(v = which.max(sil), lty = 2)


# ---------------------------------
print("Clustering the Data")
clusters = 2
set.seed(20)
kClust  <- kmeans(genes_filtered, centers = clusters, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster


# ---------------------------------
print("Heatmap")
heatmap <- pheatmap(
  kClust,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  #annotation_col = annotation_df,
  main = "Sperm RNA Complex Heatmap",
  colorRampPalette(c(
    "deepskyblue",
    "black",
    "yellow"
  ))(25
  ),
  scale = "row"
)


