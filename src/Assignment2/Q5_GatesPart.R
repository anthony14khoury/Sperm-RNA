# Q5 piece 1/4

# Cluster Profile and Gene Ontology

organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# necessary imports
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# load differential expressed genes into diffExpDF from tsv from #3 
diffExpDF <- readr::read_tsv("CleanData/diff_expr_results.tsv")

# convert back to ENSEMBLE notation for gseGO function
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
toConvert <- diffExpDF$Gene
ens_g_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values= toConvert,mart= mart)

# extract only log2FoldChange column (gene labels are lost but order is preserved)
l2fc <- diffExpDF$log2FoldChange

# re-attach gene labels but in ENSEMBLE 
names(l2fc) <- diffExpDF$Gene

# remove bad values
l2fc <- na.omit(l2fc)

# sort decreasing (necessary for cluster profiler)
l2fc = sort(l2fc, decreasing = TRUE)

# gene set enrichment analysis
gse <- gseGO(geneList=l2fc, 
             ont ="BP", 
             keyType = "HUGO",
             OrgDb = organism, 
             pAdjustMethod = "none")