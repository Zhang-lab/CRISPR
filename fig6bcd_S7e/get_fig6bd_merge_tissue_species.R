#2021-02-26 rank merged celltypes based on median
library(gplots);
library(data.table)
library(ggplot2)
library(ggsignif)
library(scater)

library(RColorBrewer)
library(RPMG)

library(tidyr)

sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}

metadata2<- read.table("meta_ss2.txt", sep="\t", comment.char="", header=T)
count_data2<- read.table("raw_data_ss2.txt", sep="\t", header=T, row.names=1, comment.char="")
count_data2<- data.frame(count_data2);

count_data1<- fread("raw_data_10x.txt", header=T, nThread=24);
count_data1<- data.frame(count_data1);

rownames(count_data1)<- count_data1[,1]
count_data1<- count_data1[,-1]

metadata1<- read.table("meta_10x.txt", sep="\t", comment.char="", header=T, stringsAsFactor=F)
#Change tag name compatile to colnames of count table
rownames(metadata2)<- paste0("X",gsub("#", ".",rownames(metadata2) ) )

metadata<- rbind(metadata1, metadata2)
count_data<- cbind(count_data1, count_data2);

mm_only_data<- read.table("mm_only_genes_new", header=F)
EG_mm_shared_data<- read.table("EG_mm_shared_genes_new", header=F)
GRG_mm_shared_data<- read.table("GRG_mm_shared_genes_new", header=F)

mm_only_genes<- toupper(mm_only_data[,1])
EG_mm_shared_genes<- toupper(EG_mm_shared_data[,1])
GRG_mm_shared_genes<- toupper(GRG_mm_shared_data[,1])

#Get alias
EG_mm_shared_genes[grep("VIRMA", EG_mm_shared_genes)]<- "KIAA1429"
mm_only_genes[grep("ARMH3", mm_only_genes)]<- "C10orf76"
mm_only_genes[grep("SUPT3", mm_only_genes)]<- "SUPT3H"

gene_vec0<- c(mm_only_genes, EG_mm_shared_genes, GRG_mm_shared_genes)

genes<- sapply(strsplit(rownames(count_data), "_"), head, 1)
count_data<- apply(count_data, 2, function(x) as.numeric(as.character(x)))

cpm=apply(count_data,2,function(x) x/sum(x)*10^6)

rownames(cpm)<- genes;

gene_vec<- gene_vec0[gene_vec0 %in% genes]
cpm_subset<- cpm[match(gene_vec, genes),]

locations<- unique(metadata[,"location"])

#Merge celltypes
metadata[metadata[,"annotation"] %in% c("dNK1", "dNK2", "dNK3", "dNK p"),"annotation"]<- "dNK"
metadata[metadata[,"annotation"] %in% c("fFB1", "fFB2"),"annotation"]<- "fFB"
metadata[metadata[,"annotation"] %in% c("DC1", "DC2"),"annotation"]<- "DC"
metadata[metadata[,"annotation"] %in% c("dM1", "dM2","dM3"),"annotation"]<- "dM"
metadata[metadata[,"annotation"] %in% c("dP1", "dP2"),"annotation"]<- "dP"
metadata[metadata[,"annotation"] %in% c("dS1", "dS2", "dS3"),"annotation"]<- "dS"
metadata[metadata[,"annotation"] %in% c("Endo L", "Endo (m)", "Endo (f)"),"annotation"]<- "Endo"
metadata[metadata[,"annotation"] %in% c("Epi1", "Epi2"),"annotation"]<- "Epi"
metadata[metadata[,"annotation"] %in% c("NK CD16-", "NK CD16+"),"annotation"]<- "NK"

#Combine all tissues
celltypes<- unique(metadata[,"annotation"])
mean_table<- data.frame(matrix(nrow=length(gene_vec),  ncol=length(celltypes)));
for(j in 1:length(celltypes)){
    cells<- colnames(cpm_subset)[colnames(cpm_subset) %in% rownames(metadata[metadata[,"annotation"]== celltypes[j], ] )]
print(paste0("cell num:",length(cells) ));
     cpm_celltype<- cpm_subset[, cells]
     if(length(cells)>1){
        mean_table[,j]<-  rowMeans(cpm_celltype)
      }else{
        mean_table[,j]<- cpm_celltype
       }
    }
rownames(mean_table)<- gene_vec
colnames(mean_table)<- celltypes
ec<- 1;
log2cpm_table<- log2(ec+mean_table);
log2cpm_data<- apply(log2cpm_table, 2, function(x) as.numeric(as.character(x)))
rownames(log2cpm_data)<- rownames(log2cpm_table);

#Split into mm_only & mm_hs_shared
log2cpm_data_mmOnly<- log2cpm_data[rownames(log2cpm_data) %in% mm_only_genes,]
log2cpm_data_mmHs<- log2cpm_data[!(rownames(log2cpm_data) %in% mm_only_genes),]

log2cpm_data_mmHs_longer<- pivot_longer(data.frame(log2cpm_data_mmHs), cols=1:ncol(log2cpm_data_mmHs), names_to="Celltype", values_to="log2cpm")

median_celltype<- apply(log2cpm_data_mmHs, 2, function(x) median(x))
celltype_order<- rownames(data.frame(median_celltype[order(-median_celltype)]))

col_vec<- c(alpha(rainbow(3)[1], 0.5), rep("lightgrey", 17) )

p1<- ggplot(log2cpm_data_mmHs_longer, aes(x=factor(Celltype, levels=celltype_order), y=log2cpm))
p1<- p1+ geom_boxplot(fill=col_vec)

p1<- p1+ theme(axis.text.x = element_text(angle = 45, hjust=1))
p1<- p1+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#boxplot_title<- expression(paste("Mouse-human overlapping genes (n=20)\nexpression in Vento-Tormo",italic("et al."), " " )  );
p1<- p1+ labs(x="Cell types", y=expression(paste(Log[2], "CPM",sep="")) )

p1<- p1+ theme(legend.position="none", plot.title=element_text(hjust=0.5, vjust=3, size=13, face="bold"), axis.text=element_text(size=14), axis.title=element_text(size=20) )

pdf("fig6b_boxplot_across_celltype_combineTissue_mergeCelltype_mmHs.pdf", width=8, height=7);
p1;
dev.off()














hs_mm_shared_genes<- c(EG_mm_shared_genes, GRG_mm_shared_genes)

target_celltypes<- c("VCT", "EVT", "SCT")
k_vec<- match(target_celltypes, colnames(log2cpm_data))
target_celltypes<- c("VCT (CTB)", "EVT", "SCT (STB)");

p_list<- vector('list', length(k_vec))
p_count<- 0;
#for(k in 1:ncol(log2cpm_data)){
for(s in 1:length(k_vec) ){

     k<- k_vec[s]
     mm_only_vec<- as.vector(log2cpm_data[rownames(log2cpm_data) %in% mm_only_genes,k])
     hs_mm_shared_vec<- as.vector(log2cpm_data[rownames(log2cpm_data) %in%hs_mm_shared_genes,k])

     group_data<- data.frame(cbind(mm_only_vec, rep("mm_only", length(mm_only_vec))) )
     colnames(group_data)<- c("log2CPM", "Group");

     hs_mm_shared_data<- data.frame(cbind(hs_mm_shared_vec, rep("hs_mm_shared", length(hs_mm_shared_vec))))
     colnames(hs_mm_shared_data)<- c("log2CPM", "Group");

     group_data<- rbind(group_data, hs_mm_shared_data)
     group_data[,"log2CPM"]<- as.numeric(as.character(group_data[,"log2CPM"]) )

     ymax<- max(group_data[,"log2CPM"])
     n_cell<- nrow(metadata[metadata[,"annotation"]== celltypes[k],])
print(paste0("ymax=", ymax))
      p<- ggplot(group_data, aes(x=Group, y=log2CPM, fill=Group)) + geom_violin() + geom_signif(comparisons=list(c("mm_only", "hs_mm_shared") ), map_signif_level=sigFunc, y_position=seq(ymax+0.9, ymax+2, length=3))

      p<- p+ labs(y=expression(paste(Log[2], "(CPM+1)")), x="",title=target_celltypes[s]);
      p<- p+ theme(plot.title=element_text(size=10))
      p<- p+ ylim(0, ymax+3)
      p<- p+ scale_x_discrete(labels=c("mm_only"="Mouse\nonly", "hs_mm_shared"="Mouse-human\noverlap"));

      p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

      p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=33, face="bold"), axis.text=element_text(size=25), axis.title=element_text(size=25) )

      p_count<- p_count +1;
      p_list[[p_count]]<- p;
     }

pdf(paste0("fig6d_violin_EMTAB78_01_geneGroup_mergeCelltype.pdf"), width=16, height=5);
print(multiplot(plotlist=p_list, cols=3) )
dev.off();





