# Install Libraries & Packages
# ---------------------------

# Install BiocManager Package
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.15")

# Install more BiocManager Packages
if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

# Install Tidyverse Package
install.packages("tidyverse")

# Install BiocManager Packages
BiocManager::install("biomaRt")
