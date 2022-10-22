# K-Means Clustering

library(tidyverse)
library(M3Drop)
library(ggplot2)


# File Directories
data_dir = file.path("Raw Data", "SRP053246")
data_file = file.path("CleanData", "CleanedData.tsv")
metadata_file = file.path(data_dir, "metadata_SRP053246.tsv")

# Load Data
metadata <- readr::read_tsv(metadata_file)
expression_df <- readr::read_tsv(data_file)%>%
  tibble::column_to_rownames("Gene")

# Ensuring Data and Metadata are in the Same Order
expression_df <- expression_df %>% 
  dplyr::select(metadata$refinebio_accession_code)
all.equal(colnames(expression_df), metadata$refinebio_accession_code)


# Subsetting data to the 5000 most variable genes (not 100% sure this is right)
x <- apply(expression_df, 1, IQR) # Calculate IQR
y <- expression_df[x > quantile(x, 0.835), ] # Selecting top ~5000 highly variable genes

