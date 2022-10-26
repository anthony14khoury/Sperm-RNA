# K-Means Clustering

library(tidyverse)
library(M3Drop)
library(ggplot2)
library(cluster)
library(pheatmap)
library(gplots)
library(NbClust)


# ---------------------------------
print("File Directories")
data_dir = file.path("Raw Data", "SRP053246")
data_file = file.path("CleanData", "CleanedData.tsv")
metadata_file = file.path(data_dir, "metadata_SRP053246.tsv")


# ---------------------------------
print("Loading the Data")
metadata = readr::read_tsv(metadata_file)
expression_df = readr::read_tsv(data_file)
View(expression_df)

# get gene expression data into matrix
expMatrix <- multiClust::input_file(data_file)

# extract only certain num of most variable genes
numGenes = 10

d <- multiClust::probe_ranking(data_file, 
                               numGenes, 
                               probe_num_selection = "Fixed_Probe_Num",
                               expMatrix, 
                               method = "SD_Rank")

df <- scale(d)
# ---------------------------------
print("Reading in the X most variable genes")
set.seed(123)
k = kmeans(df, centers = 2, nstart = 25)
kclust = k$cluster
kclust = as.data.frame(kclust)


# ---------------------------------
print("Find optimum number of clusters using Average Silhouette Width")
# sil = rep(0, 20)
# for (i in 2:20) {
#  k_1to20 = kmeans(df, centers = i, nstart = 25, iter.max = 20)
#  ss = silhouette(k_1to20$cluster, dist(df))
#  sil[i] = mean(ss[, 3])
# }

# Clustering Results
# cat("Optimal Number ff Clusters:", which.max(sil), "\n")
# plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
# abline(v = which.max(sil), lty = 2)


# Create Annotations
# annotation_df <- metadata %>%
#   dplyr::mutate(
#     mutation = dplyr::case_when(
#       startsWith(refinebio_title, "S") ~ "S",
#       startsWith(refinebio_title, "H") ~ "H",
#       TRUE ~ "unknown"
#     )
#   ) %>%
#   dplyr::select(
#     refinebio_accession_code,
#     mutation,
#     refinebio_title
#   ) %>%
#   tibble::column_to_rownames("refinebio_accession_code")
# 
# set.seed(123)
# k = kmeans(df, 2)
# pheatmap((as.matrix(df)[order(k$cluster),]), 
#         scale="row",
#         color=colorRampPalette(c("blue", "white"))(2),
#         show_rownames=F,
#         show_colnames=F,
#         cluster_cols=T, 
#         cluster_rows=T, 
#         clustering_method="complete",
#         )



# num_genes = 10
# df = data.matrix(expression_df, rownames.force = NA)
# mads = apply(df, 1, mad)
# df = df[rev(order(mads))[1:num_genes],]
# df = sweep(df, 1, apply(df, 1, median, na.rm=T))
