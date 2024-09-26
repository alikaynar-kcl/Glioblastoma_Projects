two_gene_list_overlap<-function(gene_list_1,gene_list_2,bg_num){
   #note: gene_list_1 and gene_list_2 are come from bg_gene 
  over_gene<-intersect(gene_list_1,gene_list_2)
  p_hyper=1-phyper(length(over_gene)-1,length(gene_list_1),bg_num-length(gene_list_1),length(gene_list_2))
  result<-cbind(length(gene_list_1),length(gene_list_2),length(over_gene),p_hyper)
  colnames(result)<-c("deg_1","deg_2","overlaps","p_hyper")
  return(result)
}
