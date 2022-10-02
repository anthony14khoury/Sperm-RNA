# File Paths
data_dir = file.path("../../Raw Data", "SRP053246")
data_file = file.path(data_dir, "SRP053246.tsv")
metadata_file = file.path(data_dir, "metadata_SRP053246.tsv")


# Libraries
library(org.Hs.eg.db)
library(magrittr)# We will need this so we can use the pipe: %>%
library('biomaRt')

 
# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)

# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

# Bring back the "Gene" column in preparation for mapping
expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- expression_df$Gene
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
final <- merge(expression_df,G_list,by.x="Gene",by.y="ensembl_gene_id")


# Write mapped and annotated data frame to output file
readr::write_tsv(final, file.path(
  data_dir,
  "CleanedData2.tsv"
))