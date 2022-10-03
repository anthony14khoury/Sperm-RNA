# Q4: Generate a Heatmap

# Libraries
library(ComplexHeatmap)
library(DESeq2)
library(magrittr)
library(pheatmap)

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

# Define a Minimum Counts Cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)
filtered_expression_df <- round(filtered_expression_df)

# Create a DESeqDataset
dds <- DESeqDataSetFromMatrix(
  countData = filtered_expression_df,colData = metadata, design = ~1,
)

# Perform DESeq2 normalization and transformation
dds_norm <- vst(dds)

# Choose Genes of Interest
variances <- apply(assay(dds_norm), 1, var)
upper_var <- quantile(variances, 0.75)
df_by_var <- data.frame(assay(dds_norm)) %>%
  dplyr::filter(variances > upper_var)

# Create Annotations
annotation_df <- metadata %>%
  dplyr::mutate(
    mutation = dplyr::case_when(
      startsWith(refinebio_title, "S") ~ "S",
      startsWith(refinebio_title, "H") ~ "H",
      TRUE ~ "unknown"
    )
  ) %>%
  dplyr::select(
    refinebio_accession_code,
    mutation,
    refinebio_title
  ) %>%
  tibble::column_to_rownames("refinebio_accession_code")


# Create a HEATMAP
heatmap <- pheatmap(
  df_by_var,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = annotation_df,
  main = "Sperm RNA Complex Heatmap",
  colorRampPalette(c(
    "deepskyblue",
    "black",
    "yellow"
  ))(25
  ),
  scale = "row"
)
