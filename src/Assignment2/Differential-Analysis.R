# Import Libraries
library(DESeq2)   # Attach the DESeq2 library
library(ggplot2)  # Attach the ggplot2 library for plotting
library(magrittr) # We will need this so we can use the pipe: %>%

data_dir = file.path("Raw Data", "SRP053246")
data_file = file.path("CleanData", "CleanedData.tsv")
metadata_file = file.path(data_dir, "metadata_SRP053246.tsv")


# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)

# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

head(metadata$refinebio_title)

metadata <- metadata %>%
  # Let's get the RPL10 mutation status from this variable
  dplyr::mutate(mutation_status = dplyr::case_when(
    stringr::str_detect(refinebio_title, "S") ~ "affected",
    stringr::str_detect(refinebio_title, "H") ~ "reference"
  ))

dplyr::select(metadata, refinebio_title, mutation_status)

str(metadata$mutation_status)

metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    mutation_status = factor(mutation_status, levels = c("reference", "affected"))
  )

levels(metadata$mutation_status)

filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)

# round all expression counts
gene_matrix <- round(filtered_expression_df)

ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metadata,
  # Supply our experimental variable to `design`
  design = ~mutation_status
)

deseq_object <- DESeq(ddset)

deseq_results <- results(deseq_object)

deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)

head(deseq_results)

# this is of class DESeqResults -- we want a data frame
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))

head(deseq_df)

readr::write_tsv(
  deseq_df,
  file.path(
    data_dir,
    "diff_expr_results.tsv"
  )
)

# We'll assign this as `volcano_plot`
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

ggsave(
  plot = volcano_plot + ylim(0,10),
  file.path(data_dir, "volcano_plot2.png")
)