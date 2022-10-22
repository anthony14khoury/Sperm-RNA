# Install Libraries & Packages
# ---------------------------

# Install Tidyverse Package
install.packages("tidyverse")

# Install M3Drop - Highly Expressed Genes
BiocManager::install("M3Drop")

# Install BiocManager
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.15")

# Install BiocManager Package
if (!("org.Hs.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}
if (!("umap" %in% installed.packages())) {
  BiocManager::install("umap", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  BiocManager::install("EnhancedVolcano")
}
if (!("biomaRt" %in% installed.packages())) {
  BiocManager::install("biomaRt")
}
if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2")
}
if (!("apeglm" %in% installed.packages())) {
  BiocManager::install("apeglm")
}
if (!require("BiocManager", quietly = TRUE)){
  BiocManager::install("ComplexHeatmap")
}
if (!require("pheatmap", quietly = TRUE)){
  BiocManager::install("pheatmap")
}

# Install gProfiler2 and jsonlite
install.packages("gprofiler2")
install.packages('jsonlite', version = '1.8.1') #DO NOT COMPILE JSONLITE

