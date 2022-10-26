# Getting Chi-Squared Values
options(warn=-1)

# ---------------------------------
print("Reading in Cluster Values")
cluster_results = read.csv("CleanData/ClusterResults/5000ClustResults.csv")
results = cluster_results$Results
consensus = cluster_results$Consensus5000
kmeans = cluster_results$KMeans5000
hclust = cluster_results$Hclust5000
pam = cluster_results$PAM5000


# ---------------------------------
print("Comparing Each Cluster to groups in Assignment 1")
chisq.test(table(hclust, pam))

