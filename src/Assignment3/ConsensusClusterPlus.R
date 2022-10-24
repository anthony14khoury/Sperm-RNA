data_dir = file.path("UF/data", "SRP053246")
data_file = file.path(data_dir, "CleanedData.tsv")

df = readr::read_tsv(data_file)
View(df)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ConsensusClusterPlus")

library(ConsensusClusterPlus)

num_genes = 10
d = data.matrix(df, rownames.force=NA)
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:num_genes],]
d = sweep(d,1, apply(d,1,median,na.rm=T))
results = ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",title=as.character(num_genes),distance="pearson",seed=1262118388.71279,plot="pdf",writeTable=TRUE)

num_genes = 100
d = data.matrix(df, rownames.force=NA)
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:num_genes],]
d = sweep(d,1, apply(d,1,median,na.rm=T))
results = ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",title=as.character(num_genes),distance="pearson",seed=1262118388.71279,plot="pdf",writeTable=TRUE)

num_genes = 1000
d = data.matrix(df, rownames.force=NA)
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:num_genes],]
d = sweep(d,1, apply(d,1,median,na.rm=T))
results = ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",title=as.character(num_genes),distance="pearson",seed=1262118388.71279,plot="pdf",writeTable=TRUE)


num_genes = 5000
d = data.matrix(df, rownames.force=NA)
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:num_genes],]
d = sweep(d,1, apply(d,1,median,na.rm=T))
results = ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",title=as.character(num_genes),distance="pearson",seed=1262118388.71279,plot="pdf",writeTable=TRUE)


num_genes = 10000
d = data.matrix(df, rownames.force=NA)
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:num_genes],]
d = sweep(d,1, apply(d,1,median,na.rm=T))
results = ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",title=as.character(num_genes),distance="pearson",seed=1262118388.71279,plot="pdf",writeTable=TRUE)



results[[2]][["consensusMatrix"]][1:5,1:5]
results[[2]][["consensusTree"]]
results[[2]][["consensusClass"]][1:5]


icl = calcICL(results,title=title,plot="png")
icl[["clusterConsensus"]]
icl[["itemConsensus"]][1:5,]

