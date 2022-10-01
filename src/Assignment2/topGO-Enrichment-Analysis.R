data_dir = file.path("UF/data", "SRP053246")
data_file = file.path(data_dir, "CleanedData.tsv")
interest_file = file.path(data_dir, "diff_expr_results.tsv")
metadata_file = file.path(data_dir, "metadata_SRP053246.tsv")

BiocManager::install("topGO")
BiocManager::install("GO.db")
BiocManager::install("biomaRt")
BiocManager::install("Rgraphviz")

library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)

exp_data= read.table(data_file, header=TRUE, sep='\t')
bg_genes=as.character(exp_data[,1])

candidate_list =read.table(interest_file, header=TRUE, sep='\t')
candidate_list= as.character(candidate_list[,1])

length(bg_genes)
head(bg_genes)
length(candidate_list)
head(candidate_list)

db= useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=bg_genes, mart=db)

# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GO=unstack(go_ids[,c(1,2)])
 
# remove any candidate genes without GO annotation
keep = candidate_list %in% go_ids[,2]
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]
 
# make named factor showing which genes are of interest
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes

GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)

# define test using the classic algorithm with fisher (refer to [1] if you want to understand how the different algorithms work)
classic_fisher_result=runTest(GOdata, algorithm='classic', statistic='fisher')

# define test using the weight01 algorithm (default) with fisher
weight_fisher_result=runTest(GOdata, algorithm='weight01', statistic='fisher') 
 
# generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
allGO=usedGO(GOdata)
all_res=GenTable(GOdata, weightFisher=weight_fisher_result, orderBy='weightFisher', topNodes=length(allGO))

#performing BH correction on our p values
p.adj=round(p.adjust(all_res$weightFisher,method="BH"),digits = 4)
 
# create the file with all the statistics from GO analysis
all_res_final=cbind(all_res,p.adj)
all_res_final=all_res_final[order(all_res_final$p.adj),]
 
#get list of significant GO before multiple testing correction
results.table.p= all_res_final[which(all_res_final$weightFisher<=0.001),]
 
#get list of significant GO after multiple testing correction
results.table.bh=all_res_final[which(all_res_final$p.adj<=0.05),]
 
#save first top 50 ontolgies sorted by adjusted pvalues
write.table(all_res_final[1:50,],"summary_topGO_analysis.csv",sep=",",quote=FALSE,row.names=FALSE)
 
# PLOT the GO hierarchy plot: the enriched GO terms are colored in yellow/red according to significance level
 
pdf(file='topGOPlot_fullnames.pdf', height=12, width=12, paper='special', pointsize=18)
showSigOfNodes(GOdata, score(weight_fisher_result), useInfo = "none", sigForAll=FALSE, firstSigNodes=2,.NO.CHAR=50)
dev.off()


myterms =results.table.p$GO.ID # change it to results.table.bh$GO.ID if working with BH corrected values
mygenes = genesInTerm(GOdata, myterms)
 
var=c()
for (i in 1:length(myterms))
{
myterm=myterms[i]
mygenesforterm= mygenes[myterm][[1]]
mygenesforterm=paste(mygenesforterm, collapse=',')
var[i]=paste("GOTerm",myterm,"genes-",mygenesforterm)
}
 
write.table(var,"genetoGOmapping.txt",sep="\t",quote=F)
