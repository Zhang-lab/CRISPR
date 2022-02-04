#2021-03-26 CTB vs other cell types, no Group2 but group 2 exlucde group 1
#Consider genes even multiple genes corespond to a peak
library(ggplot2)
library(data.table)
library(tidyr)
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
#print(paste(gene_stat[gene_stat[,1]==gene_name,2], dist))
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

dist_cutoff<- 50000

group1_genes<- group1_gene_stat[group1_gene_stat[,"dist"]< dist_cutoff,"gene_name"]
group2_genes<- group2_gene_stat[group2_gene_stat[,"dist"]< dist_cutoff,"gene_name"]
#Exclude group1 genes
group2_genes_noGroup1<- group2_genes[!(group2_genes %in% group1_genes)]

CnT_genes<- unique(c(group1_genes, group2_genes))

TSC_vs_naive_deg_data<- read.table("TSC_vs_naive_DEG_significant_only.csv", sep=",", header=T, row.names=1)
TSC_vs_STB_deg_data<- read.table("TSC_vs_STB_DEG_significant_only.csv", sep=",", header=T, row.names=1)
TSC_vs_EVT_deg_data<- read.table("TSC_vs_EVT_DEG_significant_only.csv", sep=",", header=T, row.names=1)

TSC_vs_naive_up<- rownames(TSC_vs_naive_deg_data[TSC_vs_naive_deg_data[,"log2FoldChange"]>0,])
TSC_vs_naive_dn<- rownames(TSC_vs_naive_deg_data[TSC_vs_naive_deg_data[,"log2FoldChange"]<0,])

TSC_vs_STB_up<- rownames(TSC_vs_STB_deg_data[TSC_vs_STB_deg_data[,"log2FoldChange"]>0,])
TSC_vs_STB_dn<- rownames(TSC_vs_STB_deg_data[TSC_vs_STB_deg_data[,"log2FoldChange"]<0,])

TSC_vs_EVT_up<- rownames(TSC_vs_EVT_deg_data[TSC_vs_EVT_deg_data[,"log2FoldChange"]>0,])
TSC_vs_EVT_dn<- rownames(TSC_vs_EVT_deg_data[TSC_vs_EVT_deg_data[,"log2FoldChange"]<0,])

D12_CTB_vs_EPI_deg_data<- read.table("DEG_D12_CTBvsEPI_wName.txt", sep="\t", header=T, row.names=1)
D12_EVT_vs_CTB_deg_data<- read.table("DEG_D12_EVTvsCTB_wName.txt", sep="\t", header=T, row.names=1)
D12_STB_vs_CTB_deg_data<- read.table("DEG_D12_STBvsCTB_wName.txt", sep="\t", header=T, row.names=1)
CTB_D12_vs_D7_deg_data<- read.table("DEG_CTB_D12vsD7_wName.txt", sep="\t", header=T, row.names=1)

D12CTB_vs_EPI_up<-unique(D12_CTB_vs_EPI_deg_data[D12_CTB_vs_EPI_deg_data[,"CTBvsEPI"]>0,"geneName"])
D12CTB_vs_EPI_dn<-unique(D12_CTB_vs_EPI_deg_data[D12_CTB_vs_EPI_deg_data[,"CTBvsEPI"]<0,"geneName"])
#CTB vs others
D12CTB_vs_EVT_up<-unique(D12_EVT_vs_CTB_deg_data[D12_EVT_vs_CTB_deg_data[,"EVTvsCTB"]<0,"geneName"])
D12CTB_vs_EVT_dn<-unique(D12_EVT_vs_CTB_deg_data[D12_EVT_vs_CTB_deg_data[,"EVTvsCTB"]>0,"geneName"])

D12CTB_vs_STB_up<-unique(D12_STB_vs_CTB_deg_data[D12_STB_vs_CTB_deg_data[,"STBvsCTB"]<0,"geneName"])
D12CTB_vs_STB_dn<-unique(D12_STB_vs_CTB_deg_data[D12_STB_vs_CTB_deg_data[,"STBvsCTB"]>0,"geneName"])

CTBD12_vs_D7_up<-unique(CTB_D12_vs_D7_deg_data[CTB_D12_vs_D7_deg_data[,"D12vsD7"]>0,"geneName"])
CTBD12_vs_D7_dn<-unique(CTB_D12_vs_D7_deg_data[CTB_D12_vs_D7_deg_data[,"D12vsD7"]<0,"geneName"])

############
get_stacked_barplot<- function(my_deg_up, my_deg_dn, my_title, total_gene_num){

total_up_percent<- 100* length(my_deg_up)/total_gene_num
total_dn_percent<- 100* length(my_deg_dn)/total_gene_num

CnT_up_percent<- 100*length(which(my_deg_up %in% CnT_genes))/length(CnT_genes)
CnT_dn_percent<- 100*length(which(my_deg_dn %in% CnT_genes))/length(CnT_genes)

group1_up_percent<- 100*length(which(my_deg_up %in% group1_genes))/length(group1_genes)
group1_dn_percent<- 100*length(which(my_deg_dn %in% group1_genes))/length(group1_genes)

group2_up_percent<- 100*length(which(my_deg_up %in% group2_genes))/length(group2_genes)
group2_dn_percent<- 100*length(which(my_deg_dn %in% group2_genes))/length(group2_genes)

group2_noGroup1_up_percent<- 100*length(which(my_deg_up %in% group2_genes_noGroup1))/length(group2_genes_noGroup1)
group2_noGroup1_dn_percent<- 100*length(which(my_deg_dn %in% group2_genes_noGroup1))/length(group2_genes_noGroup1)

percents<- c(total_up_percent, CnT_up_percent, group1_up_percent, group2_noGroup1_up_percent, total_dn_percent, CnT_dn_percent, group1_dn_percent, group2_noGroup1_dn_percent)

condition<- c(rep("Up", 4), rep("Down", 4));

groups<- rep(c("total", "CnT", "group1", "group2_noGroup1"), 2)

my_data<- data.frame(condition, groups, percents)
my_data$groups<- factor(my_data$groups, levels=c("total", "CnT", "group1", "group2_noGroup1"))
my_data$condition<- factor(my_data$condition, levels=c("Up", "Down"))


p<- ggplot(my_data, aes(fill=condition, y=percents, x=groups)) + geom_bar(position="stack", stat="identity")
p<- p+ labs(x="Group", y="Percentage") + ggtitle(paste0(my_title  ) );
p<- p+ theme(axis.text.x = element_text(angle = 0))

p<- p+ theme(plot.title=element_text(hjust=0.5, size=18, face="bold"), axis.text=element_text(size=14), axis.title=element_text(size=18) )

p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p<- p+ scale_x_discrete(label=c("Total", "CnT", "Group 1", "Group 2"))
p<- p+ scale_fill_manual(values=c("navajowhite1", "navajowhite3"))

return(p)
}
############

p1<- get_stacked_barplot(TSC_vs_naive_up, TSC_vs_naive_dn, "hTSC vs naive hPSC", 56609)
p2<- get_stacked_barplot(TSC_vs_STB_up, TSC_vs_STB_dn, "hTSC vs STB", 56609)
p3<- get_stacked_barplot(TSC_vs_EVT_up, TSC_vs_EVT_dn, "hTSC vs EVT", 56609)

p4<- get_stacked_barplot(D12CTB_vs_EPI_up, D12CTB_vs_EPI_dn, "12 d.p.f. CTB vs 12 d.p.f. EPI", 64843)
p5<- get_stacked_barplot(D12CTB_vs_EVT_up, D12CTB_vs_EVT_dn, "12 d.p.f CTB vs 12 d.p.f. EVT", 64843)
p6<- get_stacked_barplot(D12CTB_vs_STB_up, D12CTB_vs_STB_dn, "12 d.p.f CTV vs 12. d.p.f STB", 64843)
p7<- get_stacked_barplot(CTBD12_vs_D7_up, CTBD12_vs_D7_dn, "CTBD12_vs_D7", 64843)


pdf("fig3jkm_elife_TSC_vs_DEG_stacked_barplot.pdf", width=15, height=5);
grid.arrange(p1, p2, p3, nrow=1);
dev.off();


pdf("figS6fgh_scRNA_DEG_stacked_barplot.pdf", width=15, height=5);
grid.arrange(p4, p6, p5, nrow=2);
dev.off();

