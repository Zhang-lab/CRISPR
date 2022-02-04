library(ggplot2)
library(ggrepel)

data<- read.table("hTSC_mt_vs_hTSC_wt_DEG_full_list.csv", sep=",", header=T, row.names=1)

CPM_data<- read.table("CPM.csv", sep=",", header=T, row.names=1)

deg_data<- read.table("hTSC_mt_vs_hTSC_wt_DEG_significant_only.csv", sep=",", header=T, row.names=1)

index<- rep("Regularly expressed", nrow(data))
index[rownames(data) %in% rownames(deg_data[deg_data$log2FoldChange>0,])]<- "Highly expressed in KO"
index[rownames(data) %in% rownames(deg_data[deg_data$log2FoldChange<0,])]<- "Highly expressed in WT"

data_for_plot<- cbind.data.frame(data$log2FoldChange, rowMeans(CPM_data[rownames(data),]), index )

colnames(data_for_plot)<- c("LFC", "mean_CPM_WT_KO", "Index")
data_for_plot$Index<-  factor(data_for_plot$Index, level=c("Regularly expressed", "Highly expressed in KO", "Highly expressed in WT") )

data_for_plot$mean_CPM_WT_KO_log2<- log2(1+ data_for_plot$mean_CPM_WT_KO)

gene_list<- c("hla-g", "efna1", "fn1", "zfp57", "amotl2", "ceacam1", "cited2", "lepr", "loxl2", "mcam", "p2ry2", "psg11", "thbs1", "tle2", "ddit3", "fos", "met", "tet2", "unc5b", "tns1", "gli3", "bcl2", "cgb2", "cgb7", "fas", "egr1", "itgb8", "itga5", "itga1");

gene_list<- toupper(gene_list)

lfc_cutoff<- 0.5

p<- ggplot(data_for_plot, aes(x=mean_CPM_WT_KO_log2, y=LFC,col=Index, label= ifelse(rownames(data_for_plot) %in% gene_list,as.character(rownames(data_for_plot)), '') ))  + geom_point(size=0.4, alpha=0.5)

p<- p+ ggtitle("hTSC KO vs WT") + labs(x="log2(mean CPM of KO and WT)", ylab="LFC");
p<- p+ geom_hline(yintercept=as.numeric(lfc_cutoff),colour="black",linetype=2)+
       geom_hline(yintercept=-as.numeric(lfc_cutoff),colour="black",linetype=2)
p<- p+ theme_bw()+theme_classic()+ scale_colour_manual(values=c("#999999","#E69F00","#56B4E9"))

p<- p+ geom_text_repel(max.overlaps=Inf, size=5, nudge_y=0.1, nudge_x=0.1, box.padding=1, xlim=c(NA, Inf), ylim=c(-Inf, Inf) )

p<- p+ theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.5),
                text=element_text(size=12,family="Helvetica"),
                axis.title=element_text(face="bold"),
                axis.text.x=element_text(size=10,face="bold"),
                axis.text.y=element_text(size=10,face="bold"),
                legend.text=element_text(size=10,face="bold"),
                legend.title=element_text(size=10,face="bold"))

pdf("fig2i.pdf", width=10);
p;
dev.off();

