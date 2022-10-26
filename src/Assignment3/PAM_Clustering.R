install.packages(c("cluster", "factoextra", "NbClust"))
BiocManager::install("multiClust")


library(cluster)
library(factoextra)
library(NbClust)


# get gene expression data into matrix
expMatrix <- multiClust::input_file("CleanData/CleanedData.tsv")

# extract only certain num of most variable genes
numGenes = 100

d <- multiClust::probe_ranking("CleanData/CleanedData.tsv", 
                                              numGenes, 
                                              probe_num_selection = "Fixed_Probe_Num",
                                              expMatrix, 
                                              method = "SD_Rank")
# Scale data
df <- scale(d)

# Silhouette method
clusters <- fviz_nbclust(df, pam, method = "silhouette")+
  labs(subtitle = "Silhouette method")

k <- which.max(clusters[1][["data"]][["y"]])

pamResult <- pam(df, k = k)
pamResult

pm <- eclust(df,FUNcluster="pam", k=k,hc_metric = "euclidean")

fviz_cluster(pamResult, 
             palette =c("#007892","#D9455F"),
             ellipse.type ="euclid",
             repel =TRUE,
             ggtheme =theme_minimal())

#change