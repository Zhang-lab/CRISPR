library(ggplot2)
library(RColorBrewer)

Data<- read.table("filtered_peaks_3of3_qvalGt5_min2_hg38-all-gene.txt", sep="\t")

readcount_data<- read.table("GSE138688_counts_table.txt", sep="\t", header=T, row.names=1)

readcountsRaw<- readcount_data[,-1]
geneLength<- readcount_data[,1]

rpkm=apply(readcountsRaw,2,function(x) x/sum(x))*10^9/geneLength

gene_ls1k<- c();
gene_1kto2k<- c();
gene_2kto5k<- c();
gene_5kto10k<- c();
gene_10kto25k<- c();
gene_25kto50k<- c();
gene_50kto100k<- c();
gene_100kto200k<- c();
gene_200kto500k<- c();
gene_500kto1m<- c();

######################
group_gene<- function(my_gene, my_dist){

if(my_dist< 1000){
   gene_ls1k<<- c(gene_ls1k, my_gene);
  }
if(my_dist >= 1000 & my_dist< 2000){
   gene_1kto2k<<- c(gene_1kto2k, my_gene);
  }
if(my_dist >= 2000 & my_dist< 5000){
   gene_2kto5k<<- c(gene_2kto5k, my_gene);
  }
if(my_dist >= 5000 & my_dist< 10000){
   gene_5kto10k<<- c(gene_5kto10k, my_gene);
  }
if(my_dist >= 10000 & my_dist< 25000){
   gene_10kto25k<<- c(gene_10kto25k, my_gene);
  }
if(my_dist >= 25000 & my_dist< 50000){
   gene_25kto50k<<- c(gene_25kto50k, my_gene);
  }
if(my_dist >= 50000 & my_dist< 100000){
   gene_50kto100k<<- c(gene_50kto100k, my_gene);
  }
if(my_dist >= 100000 & my_dist< 200000){
   gene_100kto200k<<- c(gene_100kto200k, my_gene);
  }
if(my_dist >= 200000 & my_dist< 500000){
   gene_200kto500k<<- c(gene_200kto500k, my_gene);
  }
if(my_dist >= 500000 & my_dist< 1000000){
   gene_500kto1m<<- c(gene_500kto1m, my_gene);
  }
}
######################

max_dist<- 0;
for(i in 1:nrow(Data)){

   gene<- as.character(Data[i, 1])
   str<- as.character(Data[i, 2])

   if(grepl(',', str)){

      str_vec<- strsplit(str, ",")[[1]]

      for(j in 1:length(str_vec)){

          dist<- abs(as.numeric(strsplit(strsplit(str_vec[j], "[(]")[[1]][2], "[)]")[[1]][1]))
          if(dist > max_dist){
             max_dist<- dist;
            }
          group_gene(gene, dist);
         }
     }
   }

count_data<- cbind(length(gene_ls1k), length(gene_1kto2k), length(gene_2kto5k), length(gene_5kto10k), length(gene_10kto25k), length(gene_25kto50k), length(gene_50kto100k),length(gene_100kto200k), length(gene_200kto500k), length(gene_500kto1m) )
colnames(count_data)<- c("ls1k", "d1to2k", "d2to5k", "d5to10k", "d10to25k", "d25to50k", "d50to100k", "d100to200k", "d200to500k", "d500kto1m");

all_genes<- as.vector(Data[,1])

ec<- 1;

log2rpkm<- log2(rpkm+ec);

####################
get_gene_rpkm<- function(my_gene_bin){ 

gene_rpkms<- as.vector(log2rpkm[rownames(log2rpkm) %in% my_gene_bin, "H9_naive_TSC"])

return(gene_rpkms);
}
####################

rpkm_ls1k<- get_gene_rpkm(gene_ls1k);
rpkm_1kto2k<- get_gene_rpkm(gene_1kto2k);
rpkm_2kto5k<- get_gene_rpkm(gene_2kto5k);
rpkm_5kto10k<- get_gene_rpkm(gene_5kto10k);
rpkm_10kto25k<- get_gene_rpkm(gene_10kto25k);
rpkm_25kto50k<- get_gene_rpkm(gene_25kto50k);
rpkm_50kto100k<- get_gene_rpkm(gene_50kto100k);
rpkm_100kto200k<- get_gene_rpkm(gene_100kto200k);
rpkm_200kto500k<- get_gene_rpkm(gene_200kto500k);
rpkm_500kto1m<- get_gene_rpkm(gene_500kto1m);

rpkm_others<- as.vector(log2rpkm[!(rownames(log2rpkm) %in% all_genes), "H9_naive_TSC"])

exp_data<- cbind(rpkm_ls1k, "<1k");
exp_data<- rbind(exp_data,cbind(rpkm_1kto2k, "1-2k") );
exp_data<- rbind(exp_data,cbind(rpkm_2kto5k, "2-5k") );
exp_data<- rbind(exp_data,cbind(rpkm_5kto10k, "5-10k") );
exp_data<- rbind(exp_data,cbind(rpkm_10kto25k, "10-25k") );
exp_data<- rbind(exp_data,cbind(rpkm_25kto50k, "25-50k") );
exp_data<- rbind(exp_data,cbind(rpkm_50kto100k, "50-100k") );
exp_data<- rbind(exp_data,cbind(rpkm_100kto200k, "100-200k") );
exp_data<- rbind(exp_data,cbind(rpkm_200kto500k, "200-500k") );
exp_data<- rbind(exp_data,cbind(rpkm_500kto1m, "500k-1m") );

exp_data<- rbind(exp_data,cbind(rpkm_others, "Unbound genes") );

exp_data<- data.frame(exp_data);
colnames(exp_data)<- c("log2rpkm", "distTSS");
exp_data[,"log2rpkm"]<- as.numeric(as.character(exp_data[,"log2rpkm"]) )
exp_data[,"distTSS"]<- factor(exp_data[,"distTSS"], levels=c("<1k", "1-2k", "2-5k", "5-10k", "10-25k", "25-50k", "50-100k", "100-200k", "200-500k", "500k-1m", "Unbound genes"))


s<- ggplot(exp_data, aes(x=distTSS, y=log2rpkm, fill=distTSS))+ geom_boxplot( )
s<- s + labs(title="TEAD1 target gene expression", x="CnT peak distance to TSS", y=expression(paste(Log[2], "(RPKM+1)")) );

s<- s+ theme(axis.text.x=element_text(angle=45, hjust=1))
s<- s+ ylim(0,9)

s<- s+  scale_fill_manual(values=colorspace::lighten("slateblue4", seq(0.3, 1, length.out=11)) )

s<- s+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=14), axis.title=element_text(size=20) )

s<- s+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))



pdf("fig3i_filtered_CnT_qval5_min2_H9naiveTSC_gene_exp_bin.pdf");
s;
dev.off();

