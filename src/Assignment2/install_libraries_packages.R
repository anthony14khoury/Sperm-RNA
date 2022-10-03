# Install Libraries & Packages
# ---------------------------

# Install Tidyverse Package
install.packages("tidyverse")

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
