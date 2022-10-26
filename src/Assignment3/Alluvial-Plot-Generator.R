install.packages("ggalluvial")
library(ggalluvial)

data_dir = file.path("UF/data", "SRP053246")
data_file = file.path(data_dir, "Consensus.csv")

df = read.csv(data_file)

ggplot(data = df,
       aes(axis1 = Sample, axis2 = Consensus10, axis3 = Consensus100,
           axis4 = Consensus1000, axis5 = Consensus5000, axis6 = Consensus10000)) +
  scale_x_discrete(limits = c("Samples", "10 Genes", "100 Genes", "1000 Genes", "5000 Genes", "10000 Genes"), expand = c(.2, .05)) +
  xlab("Cluster") +
  geom_alluvium(aes(fill = substr(Sample,1,1))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("Clustering of Samples through ConsensusClusterPlus",
          "For varying numbers of highly variable genes") + 
  guides(fill=guide_legend(title="Fertility"))

data_dir = file.path("UF/data", "SRP053246")
data_file = file.path(data_dir, "PAM_Clustering.csv")

df = read.csv(data_file)

ggplot(data = df,
       aes(axis1 = Sample, axis2 = PAM10, axis3 = PAM100,
           axis4 = PAM1000, axis5 = PAM5000, axis6 = PAM10000)) +
  scale_x_discrete(limits = c("Samples", "10 Genes", "100 Genes", "1000 Genes", "5000 Genes", "10000 Genes"), expand = c(.2, .05)) +
  xlab("Cluster") +
  geom_alluvium(aes(fill = substr(Sample,1,1))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("Clustering of Samples through PAM Clustering",
          "For varying numbers of highly variable genes") + 
  guides(fill=guide_legend(title="Fertility"))

data_dir = file.path("UF/data", "SRP053246")
data_file = file.path(data_dir, "KMeans.csv")

df = read.csv(data_file)

ggplot(data = df,
       aes(axis1 = Sample, axis2 = KMeans10, axis3 = KMeans100,
           axis4 = KMeans1000, axis5 = KMeans5000, axis6 = KMeans10000)) +
  scale_x_discrete(limits = c("Samples", "10 Genes", "100 Genes", "1000 Genes", "5000 Genes", "10000 Genes"), expand = c(.2, .05)) +
  xlab("Cluster") +
  geom_alluvium(aes(fill = substr(Sample,1,1))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("Clustering of Samples through K Means Clustering",
          "For varying numbers of highly variable genes") + 
  guides(fill=guide_legend(title="Fertility"))

data_dir = file.path("UF/data", "SRP053246")
data_file = file.path(data_dir, "HClust.csv")

df = read.csv(data_file)

ggplot(data = df,
       aes(axis1 = Sample, axis2 = Hclust10, axis3 = Hclust100,
           axis4 = Hclust1000, axis5 = Hclust5000, axis6 = Hclust10000)) +
  scale_x_discrete(limits = c("Samples", "10 Genes", "100 Genes", "1000 Genes", "5000 Genes", "10000 Genes"), expand = c(.2, .05)) +
  xlab("Cluster") +
  geom_alluvium(aes(fill = substr(Sample,1,1))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("Clustering of Samples through HClust",
          "For varying numbers of highly variable genes") + 
  guides(fill=guide_legend(title="Fertility"))


data_dir = file.path("UF/data", "SRP053246")
data_file = file.path(data_dir, "5000ClustResults.csv")

df = read.csv(data_file)

ggplot(data = df,
       aes(axis1 = Sample, axis2 = Consensus5000, axis3 = Hclust5000,
           axis4 = KMeans5000, axis5 = PAM5000)) +
  scale_x_discrete(limits = c("Samples", "ConsensusClustPlus Clustering", "HClust Clustering", "K-Means Clustering", "PAM Clustering"), expand = c(.2, .05)) +
  xlab("Cluster") +
  geom_alluvium(aes(fill = substr(Sample,1,1))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("Clustering of Samples through each clustering method",
          "For 5000 highly variable genes") + 
  guides(fill=guide_legend(title="Fertility"))

data_dir = file.path("UF/data", "SRP053246")
data_file = file.path(data_dir, "HClust.csv")




