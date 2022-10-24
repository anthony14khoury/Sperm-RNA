data_dir = file.path("UF/data", "SRP053246")
data_file = file.path(data_dir, "5000VariableGenes.csv")

df = read.csv(data_file, header=TRUE, sep=",")
View(df)
#Data is already the 5000 most variable genes, sorted and prepared using Excel.

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ConsensusClusterPlus")

library(ConsensusClusterPlus)

d = data.matrix(df, rownames.force=NA)
d = sweep(d,1, apply(d,1,median,na.rm=T))

results = ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf",writeTable=TRUE)

results[[2]][["consensusMatrix"]][1:5,1:5]
results[[2]][["consensusTree"]]
results[[2]][["consensusClass"]][1:5]


icl = calcICL(results,title=title,plot="png")
icl[["clusterConsensus"]]
icl[["itemConsensus"]][1:5,]

