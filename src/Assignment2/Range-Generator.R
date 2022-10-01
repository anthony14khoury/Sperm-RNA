data_dir = file.path("UF/data", "SRP053246")
data_file = file.path(data_dir, "CleanedData.tsv")

range_df = readr::read_tsv(data_file)



column <- c()

for (row in 1:nrow(range_df)) {
  range = max(range_df[row , 2:73], na.rm=TRUE) - min(range_df[row, 2:73], na.rm=TRUE)
  column = c(column, range)
}

range_df$Range <- column

#Ranges.tsv was not added to the GitHub repository for the 
#sake of space, you can generate it if needed - it looks exactly like 
#CleanedData.tsv with an extra column for the range
readr::write_tsv(range_df, file.path(
  data_dir,
  "Ranges.tsv"
))

library(tidyverse)

range_df %>% ggplot(aes(x = Range)) + geom_density()