library(ggplot2)
library(tidyr)
library(ggsignif)
library(gridExtra)

EVT_data<- read.table("EVT", header=F, sep="\t");


data<- read.table("norm_result_k1", header=T, row.names=1, sep="\t")
readcountsRaw<- data[,-1]
cpm=apply(readcountsRaw,2,function(x) x/sum(x)*10^6)

cpm_mean<- cbind.data.frame(rowMeans(cpm[,1:6]), rowMeans(cpm[,7:8]), rowMeans(cpm[,9:14]), rowMeans(cpm[,15:16]), rowMeans(cpm[,17:22]), rowMeans(cpm[,23:24]))

 colnames(cpm_mean)<- c("mean_STB_KO", "mean_STB_WT", "mean_EVT_KO", "mean_EVT_WT", "mean_hTSC_KO", "mean_hTSC_WT")
cpm_mean<- cpm_mean[,c("mean_EVT_KO", "mean_EVT_WT")]

CPM_mean_log<- log2(1+cpm_mean)

sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}

sub_data<- CPM_mean_log[rownames(CPM_mean_log) %in% EVT_data[,1],]
CPM_mean_longer<- pivot_longer(sub_data, cols=1:ncol(sub_data), names_to="Type", values_to="CPM")

p<- ggplot(CPM_mean_longer, aes(x=Type, y=CPM)) + geom_boxplot()

p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p<- p+ labs(x="Condition", y=expression(paste(Log[2], "CPM",sep="")), title=paste0("EVT n=", nrow(sub_data)) )
p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, vjust=3, size=13, face="bold"), axis.text=element_text(size=10), axis.title=element_text(size=16) )

comparison<- list();
comparison[[1]]<- c(colnames(CPM_mean_log)[1:2]);

p<- p+ geom_signif(comparisons=comparison)
p<- p+ theme(axis.text.x = element_text(angle = 45, hjust=1))


pdf("fig4g_barplot_EVT.pdf");
p;

dev.off();











