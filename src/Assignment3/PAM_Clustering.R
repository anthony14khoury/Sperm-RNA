install.packages(c("cluster", "factoextra", "NbClust"))
BiocManager::install("multiClust")


library(cluster)
library(factoextra)
library(NbClust)


# get gene expression data into matrix
expMatrix <- multiClust::input_file("CleanData/CleanedData.tsv")

# extract only certain num of most variable genes
numGenes = 10000

d <- multiClust::probe_ranking("CleanData/CleanedData.tsv", 
                                              numGenes, 
                                              probe_num_selection = "Fixed_Probe_Num",
                                              expMatrix, 
                                              method = "SD_Rank")
# Scale data
df <- scale(d)
#transpose data
tdf<-t(df)

# Silhouette method
clusters <- fviz_nbclust(tdf, pam, method = "silhouette")+
  labs(subtitle = "Silhouette method")

k <- which.max(clusters[1][["data"]][["y"]])
k

#PAM clustering
pamResult <- pam(tdf, k = k)
pamResult
result <- pamResult$clustering

Patient_Clusters = data.frame(result)

#Write to CSV
write.csv(Patient_Clusters,"src/Assignment3/PAM_CSVs/10000.csv")

#plotting
pm <- eclust(tdf,FUNcluster="pam", k=k ,hc_metric = "euclidean")



