library(ggplot2)
library(data.table)
library(tidyr)
library(rje)
library(gridExtra)
library(ggsignif)

group1_data<- read.table("filtered_peaks_3of3_qval5_min2_overlap_TSC_DAR_GREAT_out", sep="\t", header=F)
group2_data<- read.table("filtered_peaks_3of3_qval5_min2_overlap_shared_peaks_GREAT_out", sep="\t", header=F)

##################
get_gene_stat<- function(peak_data){
gene_stat<- data.frame(matrix(nrow=0, ncol=2), stringsAsFactors=F)

for(i in 1:nrow(peak_data)){

   str<- as.character(peak_data[i,4])
   str<- gsub(" ", "", str)
   if(str== "NONE"){
      next;
     }
   if(grepl(',', str)){
      str_vec<- strsplit(str, ",")[[1]]
      for(j in 1:length(str_vec)){
          gene_name<- strsplit(str_vec[j], "[(]")[[1]][1]
          dist<- abs(as.numeric(gsub(")", "", strsplit(str_vec[j], "[(]")[[1]][2])))
          if(!(gene_name %in% gene_stat[,1]) ){
             gene_stat<- rbind(gene_stat, data.frame(cbind(gene_name, dist), stringsAsFactors=F))
            }else{
              if(as.numeric(gene_stat[gene_stat[,1]==gene_name,2]) > dist){
                 gene_stat[gene_stat[,1]==gene_name,2]<- as.character(dist)
                }
             }
         }
     }else{
        gene_name<- strsplit(str, "[(]")[[1]][1]
        dist<- abs(as.numeric(gsub(")", "",strsplit(str, "[(]")[[1]][2] )))
         if(!(gene_name %in% gene_stat[,1]) ){
             gene_stat<- rbind(gene_stat, data.frame(cbind(gene_name, dist), stringsAsFactors=F))
            }else{
              if(as.numeric(gene_stat[gene_stat[,1]==gene_name,2]) > dist){
print(paste(gene_stat[gene_stat[,1]==gene_name,2], dist))
                 gene_stat[gene_stat[,1]==gene_name,2]<- as.character(dist)
                }
             }
       }
   }
return(gene_stat)
}
#################

group1_gene_stat<- get_gene_stat(group1_data);
group2_gene_stat<- get_gene_stat(group2_data);

elife_data<- read.table("GSE138688_counts_table.txt", header=T, sep="\t", row.names=1)

elife_readcountsRaw<- elife_data[,-1]
geneLength<- elife_data[,1]
elife_rpkm=apply(elife_readcountsRaw,2,function(x) x/sum(x))*10^9/geneLength
#get log2(rpkm+1)
elife_rpkm<- log2(elife_rpkm+1)
elife_rpkm_mean<- cbind(rowMeans(elife_rpkm[,c("AN_naive", "H9_naive")]), rowMeans(elife_rpkm[,c("AN_naive_TSC", "H9_naive_TSC")]), rowMeans(elife_rpkm[,c("AN_STB", "H9_STB")]), rowMeans(elife_rpkm[,c("AN_EVT", "H9_EVT")]));
colnames(elife_rpkm_mean)<- c("naive", "TSC", "STB", "EVT");

dist_cutoff<- 50000

elife_rpkm_mean_group1<- elife_rpkm_mean[rownames(elife_rpkm_mean) %in% group1_gene_stat[group1_gene_stat[,"dist"]<= dist_cutoff,"gene_name"],]
elife_rpkm_mean_group2<- elife_rpkm_mean[rownames(elife_rpkm_mean) %in% group1_gene_stat[group2_gene_stat[,"dist"]<= dist_cutoff,"gene_name"],]
#Exclude genes in group1
elife_rpkm_mean_group2<-elife_rpkm_mean_group2[!(rownames(elife_rpkm_mean_group2) %in% rownames(elife_rpkm_mean_group1)),]

sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}


################
get_boxplot<- function(my_data, my_title, my_sig_min, my_sig_max, my_ylim_max){

comparison<- list();
type_dat<- t(combn( colnames(my_data), 2))
for(k in 1:nrow(type_dat)) {
    comparison[[k]]<- c(as.character(type_dat[k,1]), as.character(type_dat[k,2]));
   }

my_data_long<- pivot_longer(data.frame(my_data), cols=1:ncol(my_data), names_to="type", values_to="rpkm");

p<- ggplot(my_data_long, aes(x=factor(type, levels=c("naive", "TSC", "EVT", "STB")), y=rpkm, fill=type) ) + geom_boxplot(outlier.shape=NA)
p<- p+ geom_signif(comparisons=comparison, map_signif_level=sigFunc, y_position=seq(my_sig_min,my_sig_max, length.out=length(comparison)) ) 
p<- p+ labs(x="Type", y=expression(paste(Log[2], "(RPKM+1)",sep="")), title=paste0(my_title, " n=", nrow(my_data)))
p<- p+ ylim(0, my_ylim_max)
p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=15), axis.title=element_text(size=18) )
return(p)
}
################

p1<- get_boxplot(elife_rpkm_mean_group1, "Group 1", 11, 13.8, 14)
p2<- get_boxplot(elife_rpkm_mean_group2, "Group 2", 9.5, 12, 12.6)

pdf("figS6e_boxplot_elife_rpkm_celltypes.pdf", width=14);
grid.arrange(p1, p2, nrow=1)
dev.off();

