library(ggplot2)
library(tidyr)
library(ggsignif)

se <- function(x) sqrt(var(x)/length(x))

fc_data<- read.table("GSE138688_counts_table.txt", sep="\t", row.names=1,  header=T)

geneLength<- fc_data[,1]
readcountsRaw<- fc_data[,-1]
rpkm=apply(readcountsRaw,2,function(x) x/sum(x))*10^9/geneLength

rpkm_TEAD1<- rpkm["TEAD1",]

rpkm_data_TEAD1<- data.frame(matrix(nrow=2, ncol=5))
colnames(rpkm_data_TEAD1)<- c("primed", "naive", "TSC", "EVT", "STB");
rownames(rpkm_data_TEAD1)<- c("AN", "H9");

rpkm_data_TEAD1["AN", "primed"]<- rpkm_TEAD1["AN_primed"]
rpkm_data_TEAD1["AN", "naive"]<- rpkm_TEAD1["AN_naive"]
rpkm_data_TEAD1["AN", "TSC"]<- rpkm_TEAD1["AN_naive_TSC"]
rpkm_data_TEAD1["AN", "EVT"]<- rpkm_TEAD1["AN_EVT"]
rpkm_data_TEAD1["AN", "STB"]<- rpkm_TEAD1["AN_STB"]

rpkm_data_TEAD1["H9", "primed"]<- rpkm_TEAD1["H9_primed"]
rpkm_data_TEAD1["H9", "naive"]<- rpkm_TEAD1["H9_naive"]
rpkm_data_TEAD1["H9", "TSC"]<- rpkm_TEAD1["H9_naive_TSC"]
rpkm_data_TEAD1["H9", "EVT"]<- rpkm_TEAD1["H9_EVT"]
rpkm_data_TEAD1["H9", "STB"]<- rpkm_TEAD1["H9_STB"]

#Only for H9
t_rpkm_data_TEAD1<- data.frame( t(rpkm_data_TEAD1) )

p<- ggplot(t_rpkm_data_TEAD1, aes(x=factor(rownames(t_rpkm_data_TEAD1), levels=c("primed", "naive", "TSC", "EVT", "STB")) , y=H9, fill=rownames(t_rpkm_data_TEAD1) ) )
p<- p+ geom_bar(stat="identity")
p<- p+ labs(x="Condition", y="RPKM", title="H9 TEAD1 RNA-seq expression")

p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=15), axis.title=element_text(size=18) )
p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p<- p+ scale_x_discrete(labels=c("Primed hPSC","Naive hPSC", "Naive hTSC","EVT", "STB"));

p<- p+ theme(axis.text.x=element_text(angle=25, hjust=0.4, vjust=0.4))

pdf("fig3a_barplot_H9_TEAD1_mean_rpkm.pdf")
p
dev.off();

