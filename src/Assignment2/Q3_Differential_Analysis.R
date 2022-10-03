# Q3: Differential Analysis

# Import Libraries
library(DESeq2)   # Attach the DESeq2 library
library(ggplot2)  # Attach the ggplot2 library for plotting
library(magrittr) # We will need this so we can use the pipe: %>%


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

# Set up the MetaData
metadata <- metadata %>%
  dplyr::mutate(mutation_status = dplyr::case_when(
    stringr::str_detect(refinebio_title, "S") ~ "affected",
    stringr::str_detect(refinebio_title, "H") ~ "reference"
  ))

# Make mutation starus a factor and set the levels appropriately
metadata <- metadata %>%
  dplyr::mutate(
    mutation_status = factor(mutation_status, levels = c("reference", "affected"))
  )

# Define a minimum counts cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)

# Create a DESeq2 Dataset
gene_matrix <- round(filtered_expression_df)# Round all expression count
ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metadata,
  # Supply our experimental variable to `design`
  design = ~mutation_status
)

# Run differential expression analysis
deseq_object <- DESeq(ddset)
deseq_results <- results(deseq_object)
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)

# DESeqResults -> Data Frame
deseq_df <- deseq_results %>%
  # Make into data.frame
  as.data.frame() %>%
  # Row -> Column
  tibble::rownames_to_column("Gene") %>%
  # Add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  dplyr::arrange(dplyr::desc(log2FoldChange))


# Volcano Plot
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01
)

# Print out plot here
volcano_plot

# Create a table of differentially expressed genes
readr::write_tsv(
  deseq_df,
  file.path(
    data_dir,
    "diff_expr_results.tsv"
  )
)
