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



# ---------------------------------
print("Preprocessing the data") 
numGenes = 10
expMatrix <- multiClust::input_file(data_file)
d <- multiClust::probe_ranking(data_file, 
                               numGenes, 
                               probe_num_selection = "Fixed_Probe_Num",
                               expMatrix, 
                               method = "SD_Rank")
df = scale(d)
df = t(df)


# ---------------------------------
print("Clustering the Data with K-Means")
set.seed(123)
k = kmeans(df, centers = 2, nstart = 25)
kclust = k$cluster
kclust = data.frame(kclust)


# ---------------------------------
print("Find optimum number of clusters using Average Silhouette Width")
sil = rep(0, 20)
for (i in 2:20) {
  k_1to20 = kmeans(df, centers = i, nstart = 25, iter.max = 20)
  ss = silhouette(k_1to20$cluster, dist(df))
  sil[i] = mean(ss[, 3])
}


# ---------------------------------
print("Clustering Results")
cat("Optimal Number ff Clusters:", which.max(sil), "\n")
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)


# ---------------------------------
print("Heatmap")
set.seed(123)
k = kmeans(df, 2)
pheatmap((as.matrix(df)[order(k$cluster),]),
         scale="row",
         color=colorRampPalette(c("blue", "white"))(2),
         show_rownames=F,
         show_colnames=F,
         cluster_cols=T,
         cluster_rows=T,
         clustering_method="complete",
)
