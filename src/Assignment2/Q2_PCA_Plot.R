# Q2: Generate a PCA plot

# Libraries
library(DESeq2)
library(magrittr)
library(M3C)

# File Directories
data_dir = file.path("Raw Data", "SRP053246")
data_file = file.path("CleanData", "CleanedData.tsv")
metadata_file = file.path(data_dir, "metadata_SRP053246.tsv")

# Load Data
metadata <- readr::read_tsv(metadata_file)
expression_df <- readr::read_tsv(data_file)


# Ensuring Data and Metadata are in the Same Order
expression_df <- expression_df %>% 
  dplyr::select(metadata$refinebio_accession_code)
all.equal(colnames(expression_df), metadata$refinebio_accession_code)


# Prepare Metadata for DESEq2
metadata <- metadata %>% 
  dplyr::mutate(refinebio_title = factor(
    refinebio_title, levels = c("S", "H")
    ),
    refinebio_disease = as.factor(refinebio_disease)
)

# Define a Minimum Counts Cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)
filtered_expression_df <- round(filtered_expression_df)


# Create a DESeqDataset
dds <- DESeqDataSetFromMatrix(
  countData = round(filtered_expression_df),colData = metadata, design = ~1,
)

#Perform DESeq2 normalization and transformation
dds_norm <- vst(dds)

# PCA Plot
plotPCA(dds_norm, intgroup = "refinebio_title")
