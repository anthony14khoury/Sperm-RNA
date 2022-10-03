# Q2: Generate a UMAP plot

# Libraries
library(DESeq2)
library(magrittr)
library(umap)
library(ggplot2)

# File Directories
data_dir = file.path("Raw Data", "SRP053246")
data_file = file.path("CleanData", "CleanedData.tsv")
metadata_file = file.path(data_dir, "metadata_SRP053246.tsv")

# Load Data
metadata <- readr::read_tsv(metadata_file)
expression_df <- readr::read_tsv(data_file)

# Make the Data in the Order of the Metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the Same Order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

# Prepare Metadata for DESEq2
metadata <- metadata %>%
  dplyr::select( # select only the columns that we will need for plotting
    refinebio_accession_code,
    refinebio_title
  ) %>%
  dplyr::mutate( # Convert the annotation variables into factors
    refinebio_title = factor(
      refinebio_title,
      levels = c("S", "H")
    )
  )

# Define a minimum counts cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)
filtered_expression_df <- round(filtered_expression_df)

# Create a DESeqDataset
dds <- DESeqDataSetFromMatrix(
  countData = filtered_expression_df,
  colData = metadata,
  design = ~1
)

# Perform DESeq2 Normalization and Transformation
dds_norm <- vst(dds)

# Perform UMAP
normalized_counts <- assay(dds_norm) %>% t()
umap_results <- umap::umap(normalized_counts)

# Plot UMAP
umap_plot_df <- data.frame(umap_results$layout) %>%
  tibble::rownames_to_column("refinebio_accession_code") %>%
  dplyr::inner_join(metadata, by = "refinebio_accession_code")

ggplot(umap_plot_df, 
       aes(x = X1, y = X2, color=refinebio_title)) + geom_point()
