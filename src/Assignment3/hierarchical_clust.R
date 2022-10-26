# Hierarchical Clustering

# library
library(tidyverse)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install("multiClust")

# get gene expression data into matrix
expMatrix <- multiClust::input_file("CleanData/CleanedData.tsv")
View(expMatrix)

# extract only certain num of most variable genes
numGenes = 5000
mostVarExpMatrix <- multiClust::probe_ranking("CleanData/CleanedData.tsv", 
                                              numGenes, 
                                              probe_num_selection = "Fixed_Probe_Num",
                                              expMatrix, 
                                              method = "SD_Rank")
View(mostVarExpMatrix)

# get number of clusters object (necessary for cluster analysis function)
numClusters <- multiClust::number_clusters(mostVarExpMatrix, Fixed = NULL, gap_statistic = TRUE)

# run the actual clustering analysis function
multiClust::cluster_analysis(mostVarExpMatrix, 
                 cluster_type = "HClust", 
                 seed = NULL,
                 distance = "euclidean", 
                 linkage_type = "ward.D2",
                 gene_distance = "correlation", 
                 numClusters, 
                 data_name= "Male Fertility using RNA-Seq Data",
                 probe_rank = "SD_Rank", 
                 probe_num_selection = "Fixed_Probe_Num",
                 cluster_num_selection = "Gap_Statistic")

