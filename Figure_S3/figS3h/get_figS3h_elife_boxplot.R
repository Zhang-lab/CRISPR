#2021-04-09 Polish boxplots
library(tidyr)
library(ggplot2)
library(ggsignif)


file_name<- "EG_list"

EG_data<- read.table(file_name, header=F);

EGs<- toupper(EG_data[,1])

fc_data<- read.table("GSE138688_counts_table.txt", sep="\t", row.names=1,  header=T)

geneLength<- fc_data[,1]
readcountsRaw<- fc_data[,-1]
rpkm=apply(readcountsRaw,2,function(x) x/sum(x))*10^9/geneLength

#ec as 1 to add rpkm
log2rpkm<- log2(rpkm+1)


log2rpkm_EGs<- log2rpkm[rownames(log2rpkm) %in% EGs,]
log2rpkm_EGs_subset<- log2rpkm_EGs[,grep("naive_TSC|EVT|STB", colnames(log2rpkm_EGs))]

log2rpkm_EGs_subset<- log2rpkm_EGs_subset[, -grep("WIBR3_naive_TSC", colnames(log2rpkm_EGs_subset))]
print(dim(log2rpkm_EGs_subset))

rpkm_EGs_subset.long<-pivot_longer(data.frame(log2rpkm_EGs_subset), cols=1:ncol(log2rpkm_EGs_subset),names_to="Condition",values_to="Log2RPKM")

conditions<- colnames(log2rpkm_EGs_subset)
comparison<- list();
type_dat<- t(combn(conditions, 2))
for(k in 1:nrow(type_dat)) {
    comparison[[k]]<- c(as.character(type_dat[k,1]), as.character(type_dat[k,2]));
   }


sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}


Comparison<- list();
Comparison[[1]]<- c("AN_EVT", "AN_naive_TSC");
Comparison[[3]]<- c("AN_STB", "AN_naive_TSC");
Comparison[[2]]<- c("H9_EVT", "H9_naive_TSC");
Comparison[[4]]<- c("H9_STB", "H9_naive_TSC");

p<- ggplot(rpkm_EGs_subset.long, aes(x=Condition, y=Log2RPKM, fill=Condition)) + geom_boxplot(outlier.shape=NA)
p<- p+ geom_signif(comparisons= Comparison, map_signif_level= sigFunc, y_position=seq(9.5, 10.5, length.out=length(Comparison)))
p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=15), axis.title=element_text(size=20) )
p<- p+ lims(y=c(0, 11.2))
p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p<- p+ labs(y=expression(paste(Log[2], "RPKM", sep="")) )

p<- p+ scale_x_discrete(label=c("AN EVT", "AN naive TSC", "AN STB", "H9 EVT", "H9 naive TSC", "H9 STB"));
p<- p+ theme(axis.text.x = element_text(angle = 25, hjust=1))


pdf(paste0("figS3h_GSE138688_boxplot.pdf"), width=10);
p
dev.off()




