#2020-12-30 use celltypes: CTB, EVT, STB 
#2020-12-18 log2(FPKM+1), change scale to 0-10 for log2(FPKM+1), only CTB, PrE, EPI
library(tidyr)
library(ggplot2)
library(ggsignif)
library(gridExtra)

file_name<- "EG_list"

fpkm_data<- read.table("GSE136447_555-samples-fpkm.txt",sep="\t", header=T, row.names=1)
fpkm_data1<- fpkm_data[grep("chr",fpkm_data[,"Reference"]),]

metadata<- read.table("41586_2019_1875_MOESM10_ESM.txt",sep="\t", row.names=1, header=T)

Celltypes<- c("CTB", "EVT", "STB");

metadata<- metadata[metadata[,"Group"] %in% Celltypes,]

EGinvivo<- read.table(file_name, header=F)
EGinvivo[,1]<- toupper(EGinvivo[,1])

sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}

Days<- unique(metadata[,"Day"])

##################
get_plot<- function(my_day, my_ymax, my_title, my_annot_min, my_annot_max, my_color){

celltypes<- unique(metadata[metadata[,"Day"]== my_day, "Group"])
fpkm_data_EGinvivo<- fpkm_data1[fpkm_data1[,"Gene.Name"] %in% EGinvivo[,1],]
EGinvivo_dat<- data.frame(matrix(nrow=nrow(fpkm_data_EGinvivo), ncol=length(celltypes) ))
colnames(EGinvivo_dat)<- celltypes;
for(j in 1:length(celltypes)){
    cells<- rownames(metadata[metadata[,"Day"]== my_day & metadata[,"Group"]==celltypes[j], ])
    fpkm_data_EGinvivo_sub<- fpkm_data_EGinvivo[,c(1,2, which(colnames(fpkm_data_EGinvivo) %in% cells))]
    fpkm_EGinvivo<- fpkm_data_EGinvivo_sub[,-c(1:2)]
    if(ncol(fpkm_data_EGinvivo_sub) >3){
      fpkm_EGinvivo<- rowMeans(fpkm_data_EGinvivo_sub[,-c(1:2)])
      }
#FPKM+1(pseudocount)
       log2fpkm_EGinvivo<- log2(fpkm_EGinvivo+1);
       EGinvivo_dat[,j]<- log2fpkm_EGinvivo;
    }
EGinvivo_dat.long<-pivot_longer(EGinvivo_dat, cols=1:length(celltypes),names_to="Celltypes",values_to="Log2FPKM")

comparison<- list();
type_dat<- t(combn(celltypes, 2))
for(k in 1:nrow(type_dat)) {
    comparison[[k]]<- c(as.character(type_dat[k,1]), as.character(type_dat[k,2]));
   }
p<- ggplot(EGinvivo_dat.long, aes(x=Celltypes, y=Log2FPKM, fill=Celltypes) ) + geom_boxplot(outlier.shape=NA) + geom_signif(comparisons= comparison, map_signif_level= sigFunc, y_position=seq( max(EGinvivo_dat.long[,2])+my_annot_min, max(EGinvivo_dat.long[,2])+my_annot_max, length.out=length(comparison))  ) +ggtitle(my_title)
p<- p+ ylab(expression(paste(Log[2],"(FPKM+1)") ) )
p<- p+ ylim(0, my_ymax)
p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=15), axis.title=element_text(size=18) )
p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p<- p+ scale_fill_manual(values=my_color)

return(p)
}
##################

p1<- get_plot("D10", 11, "10 d.p.f", -2, 1, c("#E69F00", "#009E73") )
p2<- get_plot("D12", 11.8, "12 d.p.f.", -3.2, -2, c("#E69F00", "#56B4E9", "#009E73"))

pdf("fig_S3f_S3g.pdf", width=10, height=5);
grid.arrange(p1, p2, nrow=1)
dev.off();













