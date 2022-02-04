library(ggplot2)
library(data.table)

data<- read.table("umap.txt",sep="\t", header=T, row.names=1, comment.char="")

count_data<- fread("raw_data_10x.txt", header=T, nThread=24)
count_data<- data.frame(count_data);
metadata<- read.table("meta_10x.txt", sep="\t", comment.char="", header=T, stringsAsFactor=F)

count_data2<- fread("raw_data_ss2.txt", header=T, nThread=24)
count_data2<- data.frame(count_data2);
metadata2<- read.table("meta_ss2.txt", sep="\t", comment.char="", header=T, stringsAsFactor=F)

print("finish reading files");

gene_list<- count_data[,1]
gene_list2<- count_data2[,1]

count_data<- count_data[,-1]
count_data2<- count_data2[,-1]

count_data<- apply(count_data, 2, function(x) as.numeric(as.character(x)))
count_data2<- apply(count_data2, 2, function(x) as.numeric(as.character(x)))

cpm=apply(count_data,2,function(x) x/sum(x)*10^6)
cpm2=apply(count_data2,2,function(x) x/sum(x)*10^6)

rownames(cpm)<- gene_list
rownames(cpm2)<- gene_list2

print("finish processing data");

EG_mm_shared_data<- read.table("EG_mm_shared_genes_new", header=F)
GRG_mm_shared_data<-read.table("GRG_mm_shared_genes_new", header=F)
EG_mm_shared_genes<- toupper(EG_mm_shared_data[,1])
GRG_mm_shared_genes<- toupper(GRG_mm_shared_data[,1])

#Get alias
EG_mm_shared_genes[grep("VIRMA", EG_mm_shared_genes)]<- "KIAA1429"

gene_vec0<- c(EG_mm_shared_genes, GRG_mm_shared_genes)
genes<- sapply(strsplit(gene_list,  "_"), head, 1)
gene_vec<- gene_vec0[gene_vec0 %in% genes]
cpm_combined<- cbind(cpm, cpm2)
rownames(cpm_combined)<- gene_list
cpm_genes<- cpm_combined[which(genes %in% gene_vec),]
t_cpm_genes<- t(cpm_genes)
ec<- 1;
t_log10cpm_genes<- log10(ec+ t_cpm_genes);
rownames(t_log10cpm_genes)<- gsub("^X", "",rownames(t_log10cpm_genes) )
rownames(data)<- gsub("#", ".",rownames(data) )
data1<- cbind(data,t_log10cpm_genes[rownames(data),])
colnames(data1)<- sapply(strsplit(colnames(data1), "_"), head, 1)

Metadata<- rbind(metadata2, metadata);
rownames(Metadata)<- gsub("#", ".", rownames(Metadata))
#Merge celltypes
Metadata[Metadata[,"annotation"] %in% c("dNK1", "dNK2", "dNK3", "dNK p"),"annotation"]<- "dNK"
Metadata[Metadata[,"annotation"] %in% c("fFB1", "fFB2"),"annotation"]<- "fFB"
Metadata[Metadata[,"annotation"] %in% c("DC1", "DC2"),"annotation"]<- "DC"
Metadata[Metadata[,"annotation"] %in% c("dM1", "dM2","dM3"),"annotation"]<- "dM"
Metadata[Metadata[,"annotation"] %in% c("dP1", "dP2"),"annotation"]<- "dP"
Metadata[Metadata[,"annotation"] %in% c("dS1", "dS2", "dS3"),"annotation"]<- "dS"
Metadata[Metadata[,"annotation"] %in% c("Endo L", "Endo (m)", "Endo (f)"),"annotation"]<- "Endo"
Metadata[Metadata[,"annotation"] %in% c("Epi1", "Epi2"),"annotation"]<- "Epi"
Metadata[Metadata[,"annotation"] %in% c("NK CD16-", "NK CD16+"),"annotation"]<- "NK"

#UMAP for merged celltypes
Dat1<- cbind(data1[,1:2], Metadata[match(rownames(Metadata), rownames(data1) ),"annotation"])
colnames(Dat1)[3]<- "annotation";
celltypes<- unique(Metadata[,"annotation"])
celltypes<- celltypes[order(celltypes)]

celltypes<- c("DC", "dM", "dNK", "dP", "dS", "EVT", "Endo", "Epi", "fFB", "Granulocytes", "HB", "ILC3", "MO", "NK", "Plasma", "SCT", "Tcells", "VCT");

Dat1[,"annotation"]<- factor(Dat1[,"annotation"], levels=celltypes);

q<- ggplot(Dat1, aes(x=UMAP1, y=UMAP2, color=factor(annotation, levels=celltypes)), size=0.0000001) + geom_point()

q<- q+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title=element_text(hjust=0.5))
q<- q + labs(colour="Cell type")

pdf("fig6c_UMAP_visualize_celltypes.pdf", width=8);
q
dev.off();

########### For UMAP of each gene

genes_for_interest<- c("TRAF2", "SMG9", "KIAA1429", "WRAP53", "TRUB2", "NHLRC2", "SQLE", "TIMMDC1", "CRLS1", "DHODH");
for(j in 1:length(genes_for_interest) ){

    i<- which(colnames(data1)== genes_for_interest[j])

    gene_plot<- ggplot(data1, aes(UMAP1, UMAP2)) + geom_point(aes(colour = data1[,i]), size=0.0001) + scale_colour_gradient(low="lightgrey", high="brown", name= expression(Log[2](CPM+1)  )   ) + ggtitle(colnames(data1)[i] ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),  plot.title=element_text(hjust=0.5, size=40, face="bold"), axis.text=element_text(size=15), axis.title=element_text(size=21) )

    pdf(paste0("UMAP_visuzlize_", colnames(data1)[i],".pdf"), width=8);
    print(gene_plot)
    dev.off();
   }



