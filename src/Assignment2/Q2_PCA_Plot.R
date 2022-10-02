# Q2: Generate a PCA plot
  # We do not have counts data

# File Directories
data_dir = file.path("Raw Data", "SRP053246")
data_file = file.path("CleanData", "CleanedData.tsv")
metadata_file = file.path(data_dir, "metadata_SRP053246.tsv")

# Load Expression Data
expression_df <- readr::read_tsv(data_file)
