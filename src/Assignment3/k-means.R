# K-Means Clustering

library(tidyverse)
library(M3Drop)
library(ggplot2)
library(cluster)
library(pheatmap)
library(gplots)


# ---------------------------------
print("File Directories")
data_dir = file.path("Raw Data", "SRP053246")
data_file = file.path("CleanData", "CleanedData.tsv")
metadata_file = file.path(data_dir, "metadata_SRP053246.tsv")


# ---------------------------------
print("Loading the Data")
metadata = readr::read_tsv(metadata_file)
expression_df = readr::read_tsv(data_file)%>%
  tibble::column_to_rownames("Gene")


# ---------------------------------
print("Ensuring Data and Metadata are in the Same Order")
expression_df = expression_df %>% 
  dplyr::select(metadata$refinebio_accession_code)
all.equal(colnames(expression_df), metadata$refinebio_accession_code)


# ---------------------------------
print("Reading in the X most variable genes")
num_genes = 1000
df = data.matrix(expression_df, rownames.force = NA)
mads = apply(df, 1, mad)
df = df[rev(order(mads))[1:num_genes],]
df = sweep(df, 1, apply(df, 1, median, na.rm=T))


# Find optimum number of clusters using Average Silhouette Width
sil = rep(0, 20)
for (i in 2:20) {
 k_1to20 = kmeans(df, centers = i, nstart = 25, iter.max = 20)
 ss = silhouette(k_1to20$cluster, dist(df))
 sil[i] = mean(ss[, 3])
}

# Clustering Results
cat("Optimal Number ff Clusters:", which.max(sil), "\n")
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)

k = kmeans(df, 2)
pheatmap(as.matrix(df)[order(k$cluster),],)



# Create Annotations
annotation_df <- metadata %>%
  dplyr::mutate(
    mutation = dplyr::case_when(
      startsWith(refinebio_title, "S") ~ "S",
      startsWith(refinebio_title, "H") ~ "H",
      TRUE ~ "unknown"
    )
  ) %>%
  dplyr::select(
    refinebio_accession_code,
    mutation,
    refinebio_title
  ) %>%
  tibble::column_to_rownames("refinebio_accession_code")
#col <- metadata[,1]

# Clustering
set.seed(123)
kClust = kmeans(df, 2)
kMeans = cbind(df, kClust$cluster)
dimension = dim(kMeans)[2]
o = order(kMeans[, 73])
kMeans<- kMeans[o, ]

# Create a HEATMAP
pheatmap(
  kMeans,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = annotation_df,
  main = "Sperm RNA Complex Heatmap",
  colorRampPalette(c(
    "deepskyblue",
    "black",
    "yellow"
  ))(25
  ),
  scale = "row"
)

pheatmap(kMeans[,1:72], cluster_rows = F, cluster_cols = F, col= annotation_df, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)







heatmap(kClust, scale="row")

mtscaled <- as.matrix(scale(kCluster))
heatmap(mtscaled, Colv=F, scale='none')
#kMeans <- cbind(genes_filtered, kCluster)
# o <- order(kMeans[,73])
# kMeans <- kMeans[o,]
# heatmap <- pheatmap( 
#               kMeans[,1:73], 
#               cluster_rows = F, 
#               cluster_cols = F, 
#               row= col, 
#               legend=FALSE, 
#               show_rownames=FALSE, 
#               show_colnames=FALSE)




# heatmap <- pheatmap(kMeans[,1:73], 
#                     cluster_rows = TRUE,
#                     cluster_cols = TRUE,
#                     show_rownames = FALSE, 
#                     annotation_col = annotation_df, 
#                     show_rownames=FALSE,
#                     scale = "row",
#                     main = "K-Means Heatmap")













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
# print("Clustering the Data")
# clusters = 3
# set.seed(20)
# kClust  <- kmeans(genes_filtered, centers = clusters, nstart = 1000, iter.max = 20)
# kClusters <- kClust$cluster
# kClusters <- as.data.frame(kClusters)
# KClusters <- as.matrix(sapply(kClusters, as.numeric))  
# 
# # ---------------------------------
# print("Heatmap")
# kClusters <- data.matrix(kClusters)
# image(kClusters)
# #heatmap.2(kClusters)


