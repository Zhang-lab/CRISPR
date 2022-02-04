library(ggplot2)
library(gridExtra)


EG_data<- read.table("EG_list", header=F)
GRG_data<- read.table("GRG_list", header=F)

EG_list<- toupper(EG_data[,1])
GRG_list<- toupper(GRG_data[,1])


rc_data<- read.table("merged_gene_name_expression.txt", sep="\t", row.names =1, header=T)

readcountsRaw<- rc_data[,-1]
geneLength<- rc_data[,1]

cpm=apply(readcountsRaw,2,function(x) x/sum(x)*10^6)
rpkm=apply(readcountsRaw,2,function(x) x/sum(x))*10^9/geneLength

BT5_mean<- rowMeans(rpkm[,c("BT5_hTSC_1", "BT5_hTSC_2")])

rpkm1<- cbind(rpkm, BT5_mean)

ec<- 1e-6;
log2rpkm<- data.frame(log2(rpkm1+ec) )

log2rpkm_EGs<- log2rpkm[rownames(log2rpkm) %in% EG_list,]
log2rpkm_GRGs<- log2rpkm[rownames(log2rpkm) %in% GRG_list,]


#################
get_plot<- function(my_data, my_col, my_title){

p<- ggplot(data = my_data, mapping = aes(x=my_col)) +
  geom_histogram(aes(y=..density..),fill="bisque",color="white",alpha=0.7, bins=50) +
  geom_density() +
  geom_rug() +
  labs(y='Density', x=expression(paste(Log[2],"(RPKM+1e-6)", sep="")), title=my_title ) 

p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=18), axis.title=element_text(size=20) )

return(p);
}
#################

log2rpkm_EGs_AN_naive_hTSC<- log2rpkm_EGs[,"AN_naive_hTSC"]
log2rpkm_EGs_H9_naive_hTSC<- log2rpkm_EGs[,"H9_naive_hTSC"]
log2rpkm_EGs_BT5_mean<- log2rpkm_EGs[,"BT5_mean"]

p1<- get_plot(log2rpkm_EGs, log2rpkm_EGs_BT5_mean, "BT5 hTSCs")
p2<- get_plot(log2rpkm_EGs, log2rpkm_EGs_H9_naive_hTSC, "H9 naive hTSCs")
p3<- get_plot(log2rpkm_EGs, log2rpkm_EGs_AN_naive_hTSC, "AN naive hTSCs")


pdf("figS1e_EG_density_BT5H9AN.pdf", width=13, height=5);

grid.arrange(p1, p2, p3, ncol=3)

dev.off();



log2rpkm_GRGs_AN_naive_hTSC<- log2rpkm_GRGs[,"AN_naive_hTSC"]
log2rpkm_GRGs_H9_naive_hTSC<- log2rpkm_GRGs[,"H9_naive_hTSC"]
log2rpkm_GRGs_BT5_mean<- log2rpkm_GRGs[,"BT5_mean"]

q1<- get_plot(log2rpkm_GRGs, log2rpkm_GRGs_BT5_mean, "BT5 hTSCs")
q2<- get_plot(log2rpkm_GRGs, log2rpkm_GRGs_H9_naive_hTSC, "H9 naive hTSCs")
q3<- get_plot(log2rpkm_GRGs, log2rpkm_GRGs_AN_naive_hTSC, "AN naive hTSCs")


pdf("figS8a_GRG_density_BT5H9AN.pdf", width=13, height=5);

grid.arrange(q1, q2, q3, ncol=3)

dev.off();

