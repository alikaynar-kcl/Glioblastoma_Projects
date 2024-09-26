#Topology analysis for the modules identified by WGCNA
#Degree (intramodule connectivity) of genes in a specific module
#Clustering coefficient (cc) for genes in a specific module
#Module memership of genes in a specific module
library(WGCNA)
###############################################################################
###############################################################################
#Degrees of genes in a specific module
degree_WGCNA<-function(adjacency,moduleColors,moduleGenes){
  #input: 
  #adjacency: row and col names are gene symbols; the adjacency of whole network
  #moduleColors: its order should be mapped with the rowname or colnames (genes) of adjacency
  
  #output: 
  #degree: 4 colunms, total connectivity (kTotal), intramodular connectivity (kWithin), extra-modular connectivity (kOut), and the difference of the intra- and extra-modular connectivities (kDiff) for all genes; kTotal represents the degree of genes in the specific module.
  
  degree<-intramodularConnectivity(adjacency, as.matrix(moduleColors), scaleByMax = FALSE)
  degree<-degree[moduleGenes,]#order is same as moduleGenes
  return(degree)#mode is list
}
###############################################################################
###############################################################################
#Clustering coefficients of genes in a specific module
clustering_coefficient_WGCNA<-function(adjacency,moduleGenes){
  #adjacency can be all measured genes' or specific modules' adjacency
  int_adjacency<-adjacency[moduleGenes,moduleGenes]
  cc<-as.matrix(clusterCoef(int_adjacency))
  cc<-as.matrix(cc[moduleGenes,1])
  colnames(cc)<-"cc"
  cc<-as.data.frame(cc)
  return(cc)#mode is list
}
###############################################################################
###############################################################################
#correlation between module eigengenes and trait for meta_data
corr_module_eigengene_vs_trait<-function(MEs,trait,uni_cluster){
  #MEs: eigengenes matrix; row is sample, col is module first component
  #MEs = moduleEigengenes(exp_vst, moduleColors)$eigengenes
  #MEs = orderMEs(MEs)
  
  #trait: from meta_data; e.g. age; its order should be matched with row order of MEs
  #uni_cluster: cluster information
  nSamples<-nrow(MEs)
  moduleTraitCor = cor(MEs, trait, method="spearman",use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  
  loc<-match(uni_cluster$moduleColors,substring(rownames(moduleTraitCor),3))
  uni_cluster$moduleTraitCor=moduleTraitCor[loc]
  uni_cluster$moduleTraitPvalue=moduleTraitPvalue[loc]
  return(uni_cluster)#mode is list
}
###############################################################################
###############################################################################
#Module memberships of genes in a specific module
module_membership_WGCNA<-function(exp_vst,MEs,int_module,moduleGenes){
  #exp_vst: the vst normalization of count exp by DEseq2; row is sample, col is gene
  #MEs: eigengenes matrix; row is sample, col is module
  #MEs = moduleEigengenes(exp_vst, moduleColors)$eigengenes
  #MEs = orderMEs(MEs)
  #int_module: interested module's color
  #int_module and moduleGenes should be mapped (same module)
  
  nSamples=nrow(exp_vst)
  geneModuleMembership = as.data.frame(cor(exp_vst, MEs, method="spearman", use = "p"))#mode is list
  #geneModuleMembership: row is gene (13700), col is module (19)
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))#mode is list
  
  int_module=paste0("ME",int_module)
  MM<-as.matrix(geneModuleMembership[moduleGenes,int_module])
  MMP<-as.matrix(MMPvalue[moduleGenes,int_module])
  
  MM_matrix<-cbind(MM,MMP)
  colnames(MM_matrix)<-c("MM","MMPValue")
  rownames(MM_matrix)<-moduleGenes
  MM_matrix<-as.data.frame(MM_matrix)
  return(MM_matrix)#mode is list
}