library(data.table)
library(ggplot2)
library(ggsignif)
library(scater)


sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}



mm_only_data<- read.table("mm_only_genes_new", header=F)
EG_mm_shared_data<- read.table("EG_mm_shared_genes_new", header=F)
GRG_mm_shared_data<- read.table("GRG_mm_shared_genes_new", header=F)

mm_only_genes<- toupper(mm_only_data[,1])
EG_mm_shared_genes<- toupper(EG_mm_shared_data[,1])
GRG_mm_shared_genes<- toupper(GRG_mm_shared_data[,1])

gene_vec0<- c(mm_only_genes, EG_mm_shared_genes, GRG_mm_shared_genes)

data0<- data.frame(fread("GSE152248_AllStages_TrophoblastNuclei_datamatrix.txt", nThread=24,header=T))

rownames(data0)<- data0[,1]
data<- data0[,-1]
data<- apply(data, 2, function(x) as.numeric(as.character(x)))
rownames(data)<- rownames(data0)

cpm=apply(data,2,function(x) x/sum(x)*10^6)
rownames(cpm)<- toupper(rownames(data))

cpm_subset<- cpm[rownames(cpm) %in% toupper(gene_vec0),]

cluster_info<- read.table("GSE152248_AllStages_TrophoblastNuclei_clusters.txt", header=T)

cluster_vec<- unique(as.character(cluster_info$annot))





#####################
get_plot<- function(my_data, my_type){

log2cpm_mmOnly<- as.numeric(log2(1+rowMeans(my_data[rownames(my_data) %in% mm_only_genes,])))
log2cpm_shared<- as.numeric(log2(1+rowMeans(my_data[rownames(my_data) %in%  c(EG_mm_shared_genes, GRG_mm_shared_genes),])))

data_for_plot<- cbind.data.frame(c(log2cpm_mmOnly, log2cpm_shared), c(rep("mmOnly", length(log2cpm_mmOnly)), rep("shared", length(log2cpm_shared))) );

colnames(data_for_plot)<- c("log2CPM", "Group");
data_for_plot$Group<- factor(data_for_plot$Group, levels=c("mmOnly", "shared"));

ymax<- max(data_for_plot[,"log2CPM"])

#p<- ggplot(data_for_plot, aes(x=Group, y=log2CPM, fill=Group)) + geom_violin() + geom_signif(comparisons=list(c("mmOnly", "shared") ), map_signif_level=sigFunc, y_position=seq(ymax+0.9, ymax+2, length=3))
p<- ggplot(data_for_plot, aes(x=Group, y=log2CPM, fill=Group)) + geom_violin() + geom_signif(comparisons=list(c("mmOnly", "shared") ),  y_position=seq(ymax+0.9, ymax+2, length=3))

p<- p+ labs(y=expression(paste(Log[2], "(CPM+1)")), x="",title=paste0(my_type, " n=", ncol(my_data)));
p<- p+ theme(plot.title=element_text(size=10))
p<- p+ ylim(0, ymax+2)
p<- p+ scale_x_discrete(labels=c("mmOnly"="Mouse\nonly", "shared"="Mouse-human\noverlap"));

p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=20), axis.title=element_text(size=20) )

return(p);
}
#################


t_list<- vector('list', length(cluster_vec));

for(i in 1:length(cluster_vec)){

    cpm_sub<- cpm_subset[, colnames(cpm_subset) %in% cluster_info[cluster_info$annot== cluster_vec[i],"cluster"]]
    t_list[[i]]<- get_plot(cpm_sub, cluster_vec[i]);
   }

mat_vec<- rep(0, 3*ceiling( length(cluster_vec) /3) )
mat_vec[1:length(cluster_vec)]<- 1:length(cluster_vec);

mat<- matrix(mat_vec, ncol=3, byrow=T)


pdf("figS9b_uniq_shared_genes_GSE152248_trophoblast.pdf", width=15, height=20);
print(multiplot(plotlist=t_list, layout=mat) );
dev.off();






