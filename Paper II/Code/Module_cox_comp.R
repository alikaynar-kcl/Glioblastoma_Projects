#################################################################
## COMPARISON OF MODULE GENES AND COX RESULTS ##

setwd("code_path\\")
source('igraph_network_features.R')

path_raw<-"{path}"
setwd(path_raw)

load("coexp_network_top1_TCGA.Rdata")#
load("module_list_TCGA.Rdata")

# List_DGE_UP_Com_TCGA_CGGA_693<-list(rownames(gene_DEG_TCGA_Count_D_UP),rownames(gene_DEG_CGGA_693_UP), TCGA_vs_CGGA_693_Cox[["over_up"]])
# Common_UP_genes_TCGA_CGGA_693<-Reduce(intersect, List_DGE_UP_Com_TCGA_CGGA_693)
#mode(result_end_TCGA_vs_CGGA_693$p_hyper)="numeric"
Significant_modules_TCGA<-c("157","227")
for (i in 1:length(Significant_modules_TCGA)){
  Module_file_names<-paste("features_subnetwork_TCGA_vs_CGGA_693_Cox_TCGA_module_",Significant_modules_TCGA[i],".txt",sep="")
  int_gene<-as.matrix(moduleList_TCGA[[Significant_modules_TCGA[i]]])#
  int_corMatrix<-corMatrix[int_gene,int_gene]
  match_Genes<-NULL
  corNet<-coExpressNetwork(int_corMatrix)
  features<-networkNodeAnno(corNet)
  gene_up_C<-as.matrix(TCGA_vs_CGGA_693_Cox[["over_up"]])
  for (i in 1:length(features$node)){
    if (features$node[i] %in% gene_up_C == "FALSE")
      
    {match_Genes[i]=0}
    else {
      {match_Genes[i]=1}
    }
  }
  match_Genes<-as.matrix(match_Genes)
  
  features<-cbind(features,match_Genes)
  
  write.table(features,file=Module_file_names,sep="\t",row.names=F,col.names=T,quote=F)
}



load("coexp_network__CGGA_693_T1.Rdata")#
load("module_list_CGGA_693.Rdata")

Significant_modules_CGGA_693<-c("15","23")
for (i in 1:length(Significant_modules_CGGA_693)){
  Module_file_names<-paste("features_subnetwork_TCGA_vs_CGGA_693_Cox_CGGA_693_module_",Significant_modules_CGGA_693[i],".txt",sep="")
  int_gene<-as.matrix(moduleList_CGGA_693[[Significant_modules_CGGA_693[i]]])#
  int_corMatrix<-corMatrix_CGGA_693[int_gene,int_gene]
  match_Genes<-NULL
  corNet<-coExpressNetwork(int_corMatrix)
  features<-networkNodeAnno(corNet)
  gene_up_C<-as.matrix(TCGA_vs_CGGA_693_Cox[["over_up"]])
  for (i in 1:length(features$node)){
    if (features$node[i] %in% gene_up_C == "FALSE")
      
    {match_Genes[i]=0}
    else {
      {match_Genes[i]=1}
    }
  }
  match_Genes<-as.matrix(match_Genes)
  
  features<-cbind(features,match_Genes)
  
  write.table(features,file=Module_file_names,sep="\t",row.names=F,col.names=T,quote=F)
}





load("coexp_network_top1_TCGA.Rdata")#
load("module_list_TCGA.Rdata")

# List_DGE_UP_Com_TCGA_CGGA_693<-list(rownames(gene_DEG_TCGA_Count_D_UP),rownames(gene_DEG_CGGA_693_UP), TCGA_vs_CGGA_693_Cox[["over_up"]])
# Common_UP_genes_TCGA_CGGA_693<-Reduce(intersect, List_DGE_UP_Com_TCGA_CGGA_693)
#mode(result_end_TCGA_vs_CGGA_693$p_hyper)="numeric"
Significant_modules_TCGA<-c("157","227", "105", "60", "310", "74", "121")
for (i in 1:length(Significant_modules_TCGA)){
  Module_file_names<-paste("features_subnetwork_TCGA_vs_CGGA_325_Cox_TCGA_module_",Significant_modules_TCGA[i],".txt",sep="")
  int_gene<-as.matrix(moduleList_TCGA[[Significant_modules_TCGA[i]]])#
  int_corMatrix<-corMatrix[int_gene,int_gene]
  match_Genes<-NULL
  corNet<-coExpressNetwork(int_corMatrix)
  features<-networkNodeAnno(corNet)
  gene_up_C<-as.matrix(TCGA_vs_CGGA_325_Cox[["over_up"]])
  for (i in 1:length(features$node)){
    if (features$node[i] %in% gene_up_C == "FALSE")
      
    {match_Genes[i]=0}
    else {
      {match_Genes[i]=1}
    }
  }
  match_Genes<-as.matrix(match_Genes)
  
  features<-cbind(features,match_Genes)
  
  write.table(features,file=Module_file_names,sep="\t",row.names=F,col.names=T,quote=F)
}



load("coexp_network__CGGA_325_T1.Rdata")#
load("module_list_CGGA_325.Rdata")

Significant_modules_CGGA_325<-c("15","11", "7", "4")
for (i in 1:length(Significant_modules_CGGA_325)){
  Module_file_names<-paste("features_subnetwork_TCGA_vs_CGGA_325_Cox_CGGA_325_module_",Significant_modules_CGGA_325[i],".txt",sep="")
  int_gene<-as.matrix(moduleList_CGGA_325[[Significant_modules_CGGA_325[i]]])#
  int_corMatrix<-corMatrix_CGGA_325[int_gene,int_gene]
  match_Genes<-NULL
  corNet<-coExpressNetwork(int_corMatrix)
  features<-networkNodeAnno(corNet)
  gene_up_C<-as.matrix(TCGA_vs_CGGA_325_Cox[["over_up"]])
  for (i in 1:length(features$node)){
    if (features$node[i] %in% gene_up_C == "FALSE")
      
    {match_Genes[i]=0}
    else {
      {match_Genes[i]=1}
    }
  }
  match_Genes<-as.matrix(match_Genes)
  
  features<-cbind(features,match_Genes)
  
  write.table(features,file=Module_file_names,sep="\t",row.names=F,col.names=T,quote=F)
}



load("coexp_network__CGGA_325_T1.Rdata")#
load("module_list_CGGA_325.Rdata")

Significant_modules_CGGA_325<-c("7", "4")
for (i in 1:length(Significant_modules_CGGA_325)){
  Module_file_names<-paste("features_subnetwork_CGGA_693_vs_CGGA_325_Cox_CGGA_325_module_",Significant_modules_CGGA_325[i],".txt",sep="")
  int_gene<-as.matrix(moduleList_CGGA_325[[Significant_modules_CGGA_325[i]]])#
  int_corMatrix<-corMatrix_CGGA_325[int_gene,int_gene]
  match_Genes<-NULL
  corNet<-coExpressNetwork(int_corMatrix)
  features<-networkNodeAnno(corNet)
  gene_up_C<-as.matrix(CGGA_vs_CGGA_325_Cox[["over_up"]])
  for (i in 1:length(features$node)){
    if (features$node[i] %in% gene_up_C == "FALSE")
      
    {match_Genes[i]=0}
    else {
      {match_Genes[i]=1}
    }
  }
  match_Genes<-as.matrix(match_Genes)
  
  features<-cbind(features,match_Genes)
  
  write.table(features,file=Module_file_names,sep="\t",row.names=F,col.names=T,quote=F)
}



load("coexp_network__CGGA_693_T1.Rdata")#
load("module_list_CGGA_693.Rdata")

Significant_modules_CGGA_693<-c("23", "15")
for (i in 1:length(Significant_modules_CGGA_693)){
  Module_file_names<-paste("features_subnetwork_CGGA_693_vs_CGGA_325_Cox_CGGA_693_module_",Significant_modules_CGGA_693[i],".txt",sep="")
  int_gene<-as.matrix(moduleList_CGGA_693[[Significant_modules_CGGA_693[i]]])#
  int_corMatrix<-corMatrix_CGGA_693[int_gene,int_gene]
  match_Genes<-NULL
  corNet<-coExpressNetwork(int_corMatrix)
  features<-networkNodeAnno(corNet)
  gene_up_C<-as.matrix(CGGA_vs_CGGA_325_Cox[["over_up"]])
  for (i in 1:length(features$node)){
    if (features$node[i] %in% gene_up_C == "FALSE")
      
    {match_Genes[i]=0}
    else {
      {match_Genes[i]=1}
    }
  }
  match_Genes<-as.matrix(match_Genes)
  
  features<-cbind(features,match_Genes)
  
  write.table(features,file=Module_file_names,sep="\t",row.names=F,col.names=T,quote=F)
}



#################################################################
## COMPARISON OF MODULE GENES AND DEG RESULTS ##


List_DEG_U<-list(rownames(gene_DEG_CGGA_693_UP), rownames(gene_DEG_CGGA_325_UP), rownames(gene_DEG_TCGA_Count_D_UP))

load("coexp_network_top1_TCGA.Rdata")#
load("module_list_TCGA.Rdata")

# List_DGE_UP_Com_TCGA_CGGA_693<-list(rownames(gene_DEG_TCGA_Count_D_UP),rownames(gene_DEG_CGGA_693_UP), TCGA_vs_CGGA_693_Cox[["over_up"]])
# Common_UP_genes_TCGA_CGGA_693<-Reduce(intersect, List_DGE_UP_Com_TCGA_CGGA_693)
#mode(result_end_TCGA_vs_CGGA_693$p_hyper)="numeric"
Significant_modules_TCGA<-c("157","227")
for (i in 1:length(Significant_modules_TCGA)){
  Module_file_names<-paste("features_subnetwork_TCGA_vs_CGGA_693_DEG_UP_TCGA_module_",Significant_modules_TCGA[i],".txt",sep="")
  int_gene<-as.matrix(moduleList_TCGA[[Significant_modules_TCGA[i]]])#
  int_corMatrix<-corMatrix[int_gene,int_gene]
  match_Genes<-NULL
  corNet<-coExpressNetwork(int_corMatrix)
  features<-networkNodeAnno(corNet)
  gene_up_C<-as.matrix(intersect(rownames(gene_DEG_CGGA_693_UP), rownames(gene_DEG_TCGA_Count_D_UP)))
  for (i in 1:length(features$node)){
    if (features$node[i] %in% gene_up_C == "FALSE")
      
    {match_Genes[i]=0}
    else {
      {match_Genes[i]=1}
    }
  }
  match_Genes<-as.matrix(match_Genes)
  
  features<-cbind(features,match_Genes)
  
  write.table(features,file=Module_file_names,sep="\t",row.names=F,col.names=T,quote=F)
}



load("coexp_network__CGGA_693_T1.Rdata")#
load("module_list_CGGA_693.Rdata")

Significant_modules_CGGA_693<-c("15","23")
for (i in 1:length(Significant_modules_CGGA_693)){
  Module_file_names<-paste("features_subnetwork_TCGA_vs_CGGA_693_DEG_UP_CGGA_693_module_",Significant_modules_CGGA_693[i],".txt",sep="")
  int_gene<-as.matrix(moduleList_CGGA_693[[Significant_modules_CGGA_693[i]]])#
  int_corMatrix<-corMatrix_CGGA_693[int_gene,int_gene]
  match_Genes<-NULL
  corNet<-coExpressNetwork(int_corMatrix)
  features<-networkNodeAnno(corNet)
  gene_up_C<-as.matrix(intersect(rownames(gene_DEG_CGGA_693_UP), rownames(gene_DEG_TCGA_Count_D_UP)))
  for (i in 1:length(features$node)){
    if (features$node[i] %in% gene_up_C == "FALSE")
      
    {match_Genes[i]=0}
    else {
      {match_Genes[i]=1}
    }
  }
  match_Genes<-as.matrix(match_Genes)
  
  features<-cbind(features,match_Genes)
  
  write.table(features,file=Module_file_names,sep="\t",row.names=F,col.names=T,quote=F)
}





load("coexp_network_top1_TCGA.Rdata")#
load("module_list_TCGA.Rdata")

# List_DGE_UP_Com_TCGA_CGGA_693<-list(rownames(gene_DEG_TCGA_Count_D_UP),rownames(gene_DEG_CGGA_693_UP), TCGA_vs_CGGA_693_Cox[["over_up"]])
# Common_UP_genes_TCGA_CGGA_693<-Reduce(intersect, List_DGE_UP_Com_TCGA_CGGA_693)
#mode(result_end_TCGA_vs_CGGA_693$p_hyper)="numeric"
Significant_modules_TCGA<-c("157","227", "105", "60", "310", "74", "121")
for (i in 1:length(Significant_modules_TCGA)){
  Module_file_names<-paste("features_subnetwork_TCGA_vs_CGGA_325_DEG_UP_TCGA_module_",Significant_modules_TCGA[i],".txt",sep="")
  int_gene<-as.matrix(moduleList_TCGA[[Significant_modules_TCGA[i]]])#
  int_corMatrix<-corMatrix[int_gene,int_gene]
  match_Genes<-NULL
  corNet<-coExpressNetwork(int_corMatrix)
  features<-networkNodeAnno(corNet)
  gene_up_C<-as.matrix(intersect(rownames(gene_DEG_CGGA_325_UP), rownames(gene_DEG_TCGA_Count_D_UP)))
  for (i in 1:length(features$node)){
    if (features$node[i] %in% gene_up_C == "FALSE")
      
    {match_Genes[i]=0}
    else {
      {match_Genes[i]=1}
    }
  }
  match_Genes<-as.matrix(match_Genes)
  
  features<-cbind(features,match_Genes)
  
  write.table(features,file=Module_file_names,sep="\t",row.names=F,col.names=T,quote=F)
}



load("coexp_network__CGGA_325_T1.Rdata")#
load("module_list_CGGA_325.Rdata")

Significant_modules_CGGA_325<-c("15","11", "7", "4")
for (i in 1:length(Significant_modules_CGGA_325)){
  Module_file_names<-paste("features_subnetwork_TCGA_vs_CGGA_325_DEG_UP_CGGA_325_module_",Significant_modules_CGGA_325[i],".txt",sep="")
  int_gene<-as.matrix(moduleList_CGGA_325[[Significant_modules_CGGA_325[i]]])#
  int_corMatrix<-corMatrix_CGGA_325[int_gene,int_gene]
  match_Genes<-NULL
  corNet<-coExpressNetwork(int_corMatrix)
  features<-networkNodeAnno(corNet)
  gene_up_C<-as.matrix(intersect(rownames(gene_DEG_CGGA_325_UP), rownames(gene_DEG_TCGA_Count_D_UP)))
  for (i in 1:length(features$node)){
    if (features$node[i] %in% gene_up_C == "FALSE")
      
    {match_Genes[i]=0}
    else {
      {match_Genes[i]=1}
    }
  }
  match_Genes<-as.matrix(match_Genes)
  
  features<-cbind(features,match_Genes)
  
  write.table(features,file=Module_file_names,sep="\t",row.names=F,col.names=T,quote=F)
}



load("coexp_network__CGGA_325_T1.Rdata")#
load("module_list_CGGA_325.Rdata")

Significant_modules_CGGA_325<-c("7", "4")
for (i in 1:length(Significant_modules_CGGA_325)){
  Module_file_names<-paste("features_subnetwork_CGGA_693_vs_CGGA_325_DEG_UP_CGGA_325_module_",Significant_modules_CGGA_325[i],".txt",sep="")
  int_gene<-as.matrix(moduleList_CGGA_325[[Significant_modules_CGGA_325[i]]])#
  int_corMatrix<-corMatrix_CGGA_325[int_gene,int_gene]
  match_Genes<-NULL
  corNet<-coExpressNetwork(int_corMatrix)
  features<-networkNodeAnno(corNet)
  gene_up_C<-as.matrix(intersect(rownames(gene_DEG_CGGA_693_UP), rownames(gene_DEG_CGGA_325_UP)))
  for (i in 1:length(features$node)){
    if (features$node[i] %in% gene_up_C == "FALSE")
      
    {match_Genes[i]=0}
    else {
      {match_Genes[i]=1}
    }
  }
  match_Genes<-as.matrix(match_Genes)
  
  features<-cbind(features,match_Genes)
  
  write.table(features,file=Module_file_names,sep="\t",row.names=F,col.names=T,quote=F)
}



load("coexp_network__CGGA_693_T1.Rdata")#
load("module_list_CGGA_693.Rdata")

Significant_modules_CGGA_693<-c("23", "15")
for (i in 1:length(Significant_modules_CGGA_693)){
  Module_file_names<-paste("features_subnetwork_CGGA_693_vs_CGGA_325_DEG_UP_CGGA_693_module_",Significant_modules_CGGA_693[i],".txt",sep="")
  int_gene<-as.matrix(moduleList_CGGA_693[[Significant_modules_CGGA_693[i]]])#
  int_corMatrix<-corMatrix_CGGA_693[int_gene,int_gene]
  match_Genes<-NULL
  corNet<-coExpressNetwork(int_corMatrix)
  features<-networkNodeAnno(corNet)
  gene_up_C<-as.matrix(intersect(rownames(gene_DEG_CGGA_693_UP), rownames(gene_DEG_CGGA_325_UP)))
  for (i in 1:length(features$node)){
    if (features$node[i] %in% gene_up_C == "FALSE")
      
    {match_Genes[i]=0}
    else {
      {match_Genes[i]=1}
    }
  }
  match_Genes<-as.matrix(match_Genes)
  
  features<-cbind(features,match_Genes)
  
  write.table(features,file=Module_file_names,sep="\t",row.names=F,col.names=T,quote=F)
}
