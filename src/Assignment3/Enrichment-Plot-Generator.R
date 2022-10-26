require(lattice)
require(ggplot2)

data_dir = file.path("UF/data", "SRP053246")
data_file = file.path(data_dir, "5000ClustResults.csv")

df = read.csv(data_file)

splom(~df[2:6])