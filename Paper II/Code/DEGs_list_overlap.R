DEGs_list_overlap<-function(deg_up_1,deg_down_1,deg_up_2,deg_down_2,bg_gene){
 bg=dim(bg_gene)[1]
deg_all_1<-as.matrix(rbind(deg_up_1,deg_down_1))
deg_all_2<-as.matrix(rbind(deg_up_2,deg_down_2))

deg_up_1<-as.matrix(intersect(deg_up_1,bg_gene))
deg_down_1<-as.matrix(intersect(deg_down_1,bg_gene))
deg_up_2<-as.matrix(intersect(deg_up_2,bg_gene))
deg_down_2<-as.matrix(intersect(deg_down_2,bg_gene))

num_deg_11<-dim(deg_up_1)[1]+dim(deg_down_1)[1]
num_deg_22<-dim(deg_up_2)[1]+dim(deg_down_2)[1]

deg_inform<-rbind(c(num_deg_11,dim(deg_up_1)[1],dim(deg_down_1)[1]),c(num_deg_22,dim(deg_up_2)[1],dim(deg_down_2)[1]))
colnames(deg_inform)<-c("all","up","down")
rownames(deg_inform)<-c("deg_1","deg_2")

over_all<-as.matrix(intersect(deg_all_1,deg_all_2))
over_up<-as.matrix(intersect(deg_up_1,deg_up_2))
over_down<-as.matrix(intersect(deg_down_1,deg_down_2))
p_binom=1-pbinom(dim(over_up)[1]+dim(over_down)[1]-1,dim(over_all)[1],0.5)
p_hyper=1-phyper(dim(over_up)[1]+dim(over_down)[1]-1,num_deg_11,bg-num_deg_11,num_deg_22)

over_result<-t(as.matrix(c(num_deg_11,num_deg_22,dim(over_all)[1],dim(over_up)[1]+dim(over_down)[1],
                           (dim(over_up)[1]+dim(over_down)[1])/dim(over_all)[1],
                           dim(over_up)[1],dim(over_down)[1],p_binom,p_hyper)))
colnames(over_result)<-c("deg_1","deg_2","overlap","consistent","ratio","co_up","co_down","p_binom","p_hyper")
output<-list(deg_inform=deg_inform,over_result=over_result,over_up=over_up,over_down=over_down)
return(output)
}