library(tidyr)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(scater)

EG_list<- c("GATA2", "SKP2", "ARID3A", "TEAD1", "TCAF1", "ARID5B", "ZFHX3", "TBL1X", "CTNNB1");

GRG_list<- c("GCM1", "TET2", "AJUBA", "NR6A1",  "NR2F2") 

celltypes<- c("CTB", "EPI", "PrE");


fpkm_data<- read.table("GSE136447_555-samples-fpkm.txt",sep="\t", header=T, row.names=1)
fpkm_data1<- fpkm_data[grep("chr",fpkm_data[,"Reference"]),]

gene_data<- read.table("EG_vs_invivo_genes", header=F, sep="\t")
gene_list<- gene_data[,1]

metadata<- read.table("41586_2019_1875_MOESM10_ESM.txt",sep="\t", row.names=1, header=T)

metadata1<- metadata[metadata[,"Group"] %in% celltypes,]

days<- unique(metadata1[,"Day"])
#Exclude Day 6
days<- days[days!= "D6"]

fpkm_exp_cutoff<- 10;


##############################
get_gene_dotplot<- function(my_gene){

fpkm_gene<- fpkm_data1[fpkm_data1[,"Gene.Name"]==my_gene,]
stat_data<- data.frame(matrix(ncol=4, nrow=0));
colnames(stat_data)<- c("celltype", "day", "percent", "average");

for(i in 1:length(days)){

    sub_data<- data.frame(matrix(nrow=length(celltypes), ncol=4))
    colnames(sub_data)<- c("celltype", "day", "percent", "average");
    sub_data[,"day"]<- days[i]
    sub_data[,"celltype"]<- celltypes
   
    metadata_day<- metadata1[metadata1[,"Day"]==days[i],]
    for(j in 1:length(celltypes)){

        fpkm_sub<- fpkm_gene[,colnames(fpkm_gene) %in% rownames(metadata_day[metadata_day[,"Group"]==celltypes[j], ])] 
        fpkm_list<- as.numeric(fpkm_sub)
#Log2FPKM scale
        fpkm_list<- log2(fpkm_list+ 1);

        sub_data[j, "average"]<- mean(fpkm_list)
#cutoff log2(FPKM+1)
        sub_data[j, "percent"]<-round(length(which(fpkm_list>log2(fpkm_exp_cutoff+1)))/length(fpkm_list),3)
       }
    stat_data<- rbind(stat_data, sub_data)
   }

S1<- ggplot(stat_data, aes(x=factor(day, levels=days[order(as.numeric(gsub("D","",days)) )]), y=celltype, size=percent, color=average, group=day  )) + geom_point(alpha=0.8)+theme_classic() + ggtitle(paste0(my_gene) ) + xlab("d.p.f.") + ylab("Cell types");

S1<- S1+ scale_color_gradient(expression(paste("Mean of ",log[2],"RPKM", sep="")), low="mediumblue", high="red2", space="Lab", limit=c(0,11) )
S1<- S1+ labs(size=paste0("Cell percent of FPKM>", fpkm_exp_cutoff) )

S1<- S1+ theme(legend.direction="vertical", legend.box="horizontal", axis.text.x=element_text(angle=90, vjust=0.5, hjust=1) )
S1<- S1+ scale_x_discrete(labels=gsub("D", "", days))
S1<- S1+ scale_size(limits=c(0,1))

S1<- S1+ theme(plot.title=element_text(hjust=0.5, size=21, face="bold"), axis.text=element_text(size=12), axis.title=element_text(size=18) )
return(S1);
}
##############################

EG_p_list<- vector('list', length(EG_list))

for(i in 1:length(EG_list)){
  
    EG_p_list[[i]]<- get_gene_dotplot(EG_list[i])
   }

pdf(paste0("EG_dotplots.pdf"), width=17.2, height=5);

mat<- matrix(1:9, nrow=3, byrow=TRUE)

multiplot(plotlist=EG_p_list, cols=3, layout=mat)

dev.off();


GRG_p_list<- vector('list', length(GRG_list))

for(i in 1:length(GRG_list)){

    GRG_p_list[[i]]<- get_gene_dotplot(GRG_list[i])
   }


pdf(paste0("GRG_dotplots.pdf"), width=17.2, height=3.3);

mat<- matrix(c(1,2,3,4,5,0), nrow=2, byrow=TRUE)

multiplot(plotlist=GRG_p_list, cols=3, layout=mat)

dev.off();

###################### For specific genes (not EG or GRG)

gene_list<- c("AXIN2","E2F8", "MAPKAPK3", "NODAL", "SAV1", "NF2", "TBX3", "TGFBI", "PTPN14");

gene_p_list<- vector('list', length(gene_list))

for(i in 1:length(gene_list)){

    gene_p_list[[i]]<- get_gene_dotplot(gene_list[i]);
   }

print(paste0("gene list len: ", length(gene_list)) );
pdf("genes_dotplots.pdf", width=17.2, height=5);
mat<- matrix(c(1,2,3,4,5,6,7,8,9), nrow=3, byrow=TRUE);

multiplot(plotlist=gene_p_list, cols=3, layout=mat);
dev.off();



