#build coexpression network, we only extract top 1% sub-network
#Then identify the module and generate module-realted list
# rm(list=ls())

library(matrixStats)
library(reshape2)
library(igraph)
library(compiler)
library(reshape2)
library(ggplot2)
library(survival)

require(matrixStats)
require(data.table)
require(igraph)
require(Hmisc)
require(compiler)
require(reshape2)

setwd("{PATH}")
source('networkson_change.R') # source of functions

top_quantile=0.95 #extract the coexpressed gene pair with r value at top 5%
num_nodes=30 # min number of gene
transitivity=0.5# a cutoff to identify module in which the genes have high Connectivity
exp_th=1 # cut-off value for expressed genes

symbol_exp<-TCGA_TP@assays@data@listData[["tpm_unstrand"]]
rownames(symbol_exp)<-TCGA_TP@rowRanges@elementMetadata@listData[["gene_id"]]
colnames(symbol_exp)<-TCGA_TP@colData@listData[["barcode"]]

#symbol_exp<-t(symbol_exp)
symbol_exp<-data.matrix(symbol_exp)
mode(symbol_exp)="numeric"

#symbol_exp <- data.frame(apply(symbol_exp, 2, function(x) as.numeric(as.character(x))))
symbol_exp <- symbol_exp[which(rowMeans(symbol_exp)>=1),] #only keep the genes average TPM >=1 
#save(file="{PATH}/symbol_exp.Rdata",symbol_exp)
path_out<-"{PATH}"
setwd(path_out)

corMatrix = makeCorTable(symbol_exp, cutoff = top_quantile, mode = "spearman", self=F, debug =F)##Creat co-expression matrix (corMatrix); take longest time
#save(file="coexp_network_top5.Rdata",corMatrix)
save(file="corMatrix.Rdata",corMatrix)

corNet = makeCorNet(corMatrix)#Creat the network graph based on the co-expression matrix
save(file="corNet.Rdata",corNet)
moduleList= makeModuleList(corNet, debug = F) #Use random walk method to extract interactive module; take some time
#save(file="module_list_high_low.Rdata",moduleList)
# save(file="moduleList.Rdata",moduleList)

cytoscapematerial = annotateModulesByCC(corNet, moduleList, cutCluster = num_nodes, cutCC = transitivity, debug = F)#Extract the modue list with number of nodes (cutCluster) >30; If connectivity (cutCC) > 0.5, label as "HighCC", otherwise "LowCC".
save(file="cytoscapematerial.Rdata",cytoscapematerial)
write.node.cytoscape(cytoscapematerial$nodeTable, 'cytoNodefile_H_L.txt')#output the module number and its transitivity/Connectivity

############################################################33
# Load saved results
path_raw<-"{PATH}"
setwd(path_raw)
load("module_list_high_low.Rdata")
module_info<-as.matrix(read.csv("cytoNodefile_H_L.txt",header=T,sep="\t"))

load("{PATH}/corMatrix.Rdata")
load("{PATH}/corNet.Rdata")
load("{PATH}/moduleList.Rdata")

bg_num<-dim(symbol_exp)[1] # background gene list (universe)
View(TCGA_HPA_INFO)
View(ENS_Sym_HGNC)
# gene_lists is a list congaing up regulated and down regulated genes from DEG results
# up_marker is a list of gene symbol of interested genes list
up_marker<-as.matrix(gene_lists[["Low_Down"]])# 
up_marker<-as.matrix(gene_lists[["Low_Up"]])# 
up_marker<-as.matrix(gene_lists[["High_Down"]])# 
up_marker<-as.matrix(gene_lists[["High_Up"]])# 

TPS_Gene_Comp_S<-list_of_dataframes[["TPS_Gene_Comp_S"]] # TPS_Gene_Comp_S includes TPS correlation result with strict cut-off
# Filter the genes based on positive and negative Spearman correlation
positive_correlation_genes <- TPS_Gene_Comp_S$Row.names[TPS_Gene_Comp_S$Spearman_Correlation >= 0.5]
negative_correlation_genes <- TPS_Gene_Comp_S$Row.names[TPS_Gene_Comp_S$Spearman_Correlation <= -0.5]


up_marker<-as.matrix(TPS_Gene_Comp_S$Row.names)# 
up_marker<-as.matrix(positive_correlation_genes)# 
up_marker<-as.matrix(negative_correlation_genes)# 

# two_gene_list_overlap is a function operating hypergeometric test to find significant overlap of two gene list
result_end<-NULL
for (i in 1:dim(module_info)[1]){
  num_module<-as.numeric(module_info[i,1])
  gene_list_1<-as.matrix(moduleList[num_module][[1]])
  result_each<-two_gene_list_overlap(gene_list_1,up_marker,bg_num)
  result_each<-cbind(result_each,num_module)
  result_end<-rbind(result_end,result_each)
}
#Modules_Positive_Cor<-result_end
#Modules_Negative_Cor<-result_end

write.table(Modules_Negative_Cor,file="Modules_Negative_Cor_TCGA.txt",sep="\t",row.names=F,col.names=T,quote=F)
############################################

# Topology of network
#calculate the degree, betweeness and closeness of each node for a subnetwork
rm(list=ls())
rm(list=ls())
rm(list=ls())

setwd("code_path\\")
source('igraph_network_features.R') # is a source of function

path_raw<-"{PATH}"
setwd(path_raw)
load("coexp_network_top5.Rdata")# Load graph network
load("module_list_TCGA.Rdata") # module lists

int_gene<-as.matrix(moduleList[7][[1]])#------------change to interested module number
int_corMatrix<-corMatrix[int_gene,int_gene] # gene-gene cor matrix, it take too long

corNet<-coExpressNetwork(int_corMatrix) # construct the network
features<-networkNodeAnno(corNet) # annotaion of the modules

result_end_features <- NULL
# Hypergeometric test shows significantly overlapped module with the TPS genes.
# 13 and 67 overlapped with negatively correlated TPS genes
# 7, 8, 36 overlapped with positively correlated TPS genes
module_numbers <- c(7,8,36,13,67) 
# for each module it make calculation
for (i in module_numbers) {
  int_gene <- as.matrix(moduleList[i][[1]])
  int_corMatrix <- corMatrix[int_gene, int_gene]
  corNet <- coExpressNetwork(int_corMatrix)
  features <- networkNodeAnno(corNet)
  # add a column for the module number
  features <- cbind(features, ModuleNumber = i)
  
  result_end_features <- rbind(result_end_features, features)
}

save(file="result_end_features.Rdata",result_end_features)

write.table(result_end_features,file="result_end_features.txt",sep="\t",row.names=F,col.names=T,quote=F)
result_end_features_HPA <- merge(result_end_features,TCGA_HPA_INFO,  by = 0, all.x = TRUE)
result_end_features_HPA <- merge(result_end_features,TCGA_HPA_INFO, by.x = 'Row.names', by = 0)
################################################################################
# get information for selected modules' genes
# Load necessary libraries
library(dplyr)
library(openxlsx)

# ENS_Sym_HGNC is already a dataframe
ENS_Sym_HGNC_2 <- data.frame(Ensembl_ID = c(ENS_Sym_HGNC$V1), GeneSymbol = c(ENS_Sym_HGNC$V2))
rownames(ENS_Sym_HGNC_2) <- ENS_Sym_HGNC_2$Ensembl_ID

module_info2 <- as.data.frame(module_info)

# an empty dataframe to store the final results
final_result_df <- data.frame(node = integer(), Ensembl_ID = character(), GeneSymbol = character(), size = integer(), summary = character(), cc = numeric(), stringsAsFactors = FALSE)

# Iterate through each module
for (i in 1:nrow(module_info2)) {
  module_genes <- moduleList[[as.numeric(module_info2$node[i])]]
  module_info_row <- module_info2[i, ]
  
  # an empty dataframe to store results for the current module
  module_df <- data.frame(node = integer(), Ensembl_ID = character(), GeneSymbol = character(), size = integer(), summary = character(), cc = numeric(), stringsAsFactors = FALSE)
  
  for (gene in module_genes) {
    gene_symbol <- ENS_Sym_HGNC_2 %>% filter(Ensembl_ID == gene) %>% pull(GeneSymbol)
    
    new_row <- data.frame(
      node = module_info_row$node,
      Ensembl_ID = gene,
      GeneSymbol = ifelse(length(gene_symbol) > 0, gene_symbol, NA),
      size = module_info_row$size,
      summary = module_info_row$summary,
      cc = module_info_row$cc,
      stringsAsFactors = FALSE
    )
    
    module_df <- rbind(module_df, new_row)
  }
  
  # Concatenate module_df to final_result_df
  final_result_df <- rbind(final_result_df, module_df)
}

# Save the final dataframe to an Excel file
write.xlsx(final_result_df, "All_Modules_Info.xlsx", rowNames = FALSE)


