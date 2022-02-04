#2020-12-11 Use log2FC rather than BF
#Include up/downstream genes
library(plyr)
suppressMessages(library(rtracklayer) )
#Include BT5=0 genes
deseq2_d6<- read.table("Day_6_vs_Day0_DEG_full_list.csv", sep=",", header=T, row.names=1)
deseq2_d12<- read.table("Day_12_vs_Day0_DEG_full_list.csv", sep=",", header=T, row.names=1)
deseq2_d18<- read.table("Day18_vs_Day0_DEG_full_list.csv", sep=",", header=T, row.names=1)

rownames(deseq2_d6)<- toupper(rownames(deseq2_d6))
rownames(deseq2_d12)<- toupper(rownames(deseq2_d12))
rownames(deseq2_d18)<- toupper(rownames(deseq2_d18))

z<- import("gencode.v27.annotation.gtf")
dat_z<- data.frame(z)

dat_gene<- dat_z[dat_z$type=="gene",]
dat_gene<- dat_gene[!grepl("_PAR_Y", dat_gene$gene_id),]

dat_gene<- dat_gene[order(dat_gene$seqnames, dat_gene$start),]

se <- function(x) sqrt(var(x)/length(x))

#Get a list of Growth-restrcting genes

GRG_data<- read.table("GRG_list", sep="\t", header=F)
bagel_genes<- toupper(GRG_data[,1])

count_data<- read.table("bagel_norm_count_rep", sep="\t", row.names=1, header=T)

genes<- sapply(strsplit(gsub("Br_", "",rownames(count_data)),"_"), head, 1)
count_data1<- cbind(count_data, genes)

count_data_mean<- aggregate(count_data1[,1:7], list(count_data1$genes), mean)
rownames(count_data_mean)<- count_data_mean[,"Group.1"]

rc_data<- cbind(count_data_mean[,"day0"], rowMeans(count_data_mean[,grep("6", colnames(count_data_mean))]),rowMeans(count_data_mean[,grep("12", colnames(count_data_mean))]), rowMeans(count_data_mean[,grep("18",colnames(count_data_mean))]))
colnames(rc_data)<- c("day0", "day6", "day12", "day18");

bagel_d18<- read.table("pr_d18_vs_d0", header=T, row.names=1, sep="\t")
bagel_d12<- read.table("pr_d12_vs_d0", header=T, row.names=1, sep="\t")
bagel_d6<- read.table("pr_d6_vs_d0", header=T, row.names=1, sep="\t")

#To calculate SE
sd_data<- cbind(apply(count_data_mean[,grep("6", colnames(count_data_mean))], 1, function(x) se(x)), apply(count_data_mean[,grep("12", colnames(count_data_mean))], 1, function(x) se(x)), apply(count_data_mean[,grep("18", colnames(count_data_mean))], 1, function(x) se(x)) )
colnames(sd_data)<- c("day6", "day12", "day18");
#https://stackoverflow.com/questions/15063287/add-error-bars-to-show-standard-deviation-on-a-plot-in-r

####begin study_gene
study_gene<- function(){

EGs<- rownames(rc_data)
EGs<- EGs[EGs %in% bagel_genes]

up_EGs<- c();
down_EGs<- c();

for(i in 1:length(EGs)){
if( nrow(dat_gene[dat_gene$gene_name==EGs[i],] ) >1){
print(paste0(i, " is ", dim(dat_gene[dat_gene$gene_name==EGs[i],]) ) );
  }
    for(j in 1:200){
        up_EG<- dat_gene[which(dat_gene$gene_name==EGs[i])-j,"gene_name"][1]     
         if(up_EG %in% rownames(deseq2_d6)){   
            up_EGs<- c(up_EGs, up_EG);
            break;
           }
       }
    for(j in 1:200){
        down_EG<- dat_gene[which(dat_gene$gene_name==EGs[i])+j,"gene_name"][1]
         if(down_EG %in% rownames(deseq2_d6)){
            down_EGs<- c(down_EGs, down_EG);
            break;
           }
       }
    }
EG_data_subset<- rc_data[rownames(rc_data) %in% bagel_genes,]
EG_sd_subset<- sd_data[rownames(sd_data) %in% bagel_genes,]

## Up & Down 
EG_data_subset<- rc_data[rownames(rc_data) %in% EGs,]
EG_up_data_subset<- rc_data[rownames(rc_data) %in% up_EGs,]
EG_down_data_subset<- rc_data[rownames(rc_data) %in% down_EGs,]

sd_sgrna<- colMeans(sd_data[rownames(sd_data)%in% EGs,])
sd_up_sgrna<- colMeans(sd_data[rownames(sd_data)%in% up_EGs,])
sd_down_sgrna<- colMeans(sd_data[rownames(sd_data)%in% down_EGs,])

EG_sgrna<- colMeans(EG_data_subset)
EG_up_sgrna<- colMeans(EG_up_data_subset)
EG_down_sgrna<- colMeans(EG_down_data_subset)

############# Use log2FC instead of BF

EG_BF<-c(0, mean(unlist(deseq2_d6[bagel_genes,"log2FoldChange"])), mean(unlist(deseq2_d12[bagel_genes,"log2FoldChange"])), mean(unlist(deseq2_d18[bagel_genes,"log2FoldChange"])) )
EG_up_BF<-c(0, mean(unlist(deseq2_d6[up_EGs,"log2FoldChange"])), mean(unlist(deseq2_d12[up_EGs,"log2FoldChange"])), mean(unlist(deseq2_d18[up_EGs,"log2FoldChange"])) )
EG_down_BF<-c(0, mean(unlist(deseq2_d6[down_EGs,"log2FoldChange"])), mean(unlist(deseq2_d12[down_EGs,"log2FoldChange"])), mean(unlist(deseq2_d18[down_EGs,"log2FoldChange"])) )



epsilon<- 0.02

Dat<- rbind(c(EG_up_sgrna, sd_up_sgrna, EG_up_BF), c(EG_sgrna, sd_sgrna, EG_BF), c(EG_down_sgrna, sd_down_sgrna, EG_down_BF));
colnames(Dat)<- c("sgrna_d0", "sgrna_d6", "sgrna_d12", "sgrna_d18",  "sd_d6", "sd_d12", "sd_d18", "BF_d0","BF_d6","BF_d12","BF_d18");

gene_list<- c("GCM1", "TGFBI", "PTPN14", "TET2", "NR6A1");

my_gene_list<- c();
for(i in 1:length(gene_list)){
    for(j in 1:200){
        up_gene<- dat_gene[which(dat_gene$gene_name==gene_list[i])-j, "gene_name"][1]
        if(up_gene %in% rownames(deseq2_d6)){
           my_gene_list<-c(my_gene_list, up_gene);
           break;
          }
       }
    my_gene_list<- c(my_gene_list, gene_list[i])
    for(j in 1:200){
        down_gene<- dat_gene[which(dat_gene$gene_name==gene_list[i])+j, "gene_name"][1]
        if(down_gene %in% rownames(deseq2_d6)){
           my_gene_list<-c(my_gene_list, down_gene);
           break;
          }
        }
   }

pdf("fig5b_S8c_GRG_exp_across_time.pdf", width=11.8, height=15.6);

par(mfrow=c(6,3))
par(mar=c(4.1, 8.1, 5.1, 6), xpd=T)

for(i in 1:nrow(Dat)){
#Set max of up/down/EG
ymax<- max(Dat[,grep("^sgrna_", colnames(Dat))]) + max(Dat[, grep("sd_", colnames(Dat))])+10;
    ymax<- round_any(ymax, 100, f=ceiling)
     EG_sgrna<- Dat[i, grep("^sgrna_", colnames(Dat))]
     sd_sgrna<- Dat[i, grep("^sd_", colnames(Dat))]
    EG_BF<- Dat[i, grep("BF_", colnames(Dat))]
    if(i==1){
      plot_main="Upstream genes"
      }
    if(i==2){
      plot_main="GRGs"
      }
    if(i==3){
      plot_main="Downstream genes"
      }
    plot(1:4, EG_sgrna, pch=21, xlab="", ylab="", type="l", col="brown", axes=F, main=plot_main, ylim=c(0,ymax), cex.main=3, lwd=7)
#    segments(2:4, EG_sgrna[-1] -sd_sgrna, 2:4, EG_sgrna[-1] +sd_sgrna, col="brown", lwd=7);
#    segments(2:4 -epsilon,EG_sgrna[-1]-sd_sgrna, 2:4 +epsilon, EG_sgrna[-1] - sd_sgrna,col="brown", lwd=7);
#    segments(2:4 -epsilon,EG_sgrna[-1]+sd_sgrna, 2:4 +epsilon, EG_sgrna[-1] + sd_sgrna,col="brown", lwd=7);
    axis(2, ylim=c(0, ymax), col="black", las=1, at=seq(0, ymax, by=100), cex.axis=1.5);
mtext(paste0("Normalized\nread counts"), side=2, line=3.5, cex=1.8)
    par(new=T);
    plot(1:4, EG_BF, pch=23, xlab="", ylab="", type="l", col="grey40", axes=F, ylim=c(-4,4), lwd=7)
    axis(4, ylim=c(-4, 4), col="grey40", las=1, col.axis="grey40", at=seq(-4, 4, by=2), cex.axis=1.5);
    mtext(bquote(paste('Log'['2'],'FC')), side=4, line=3.6, col="grey40", cex=1.8)            
    axis(1, at=1:4, labels=c(0,6,12,18), 10, cex.axis=1.5)
    mtext("Days", side=1, col="black", line=2.9, cex=1.8)
   }
#########################################
for(i in 1:length(my_gene_list)){
    my_gene<- my_gene_list[i]
    rc_data_subset<- rc_data[my_gene,]
############# Use log2FC instead of BF
     BF<- c(0,deseq2_d6[my_gene,"log2FoldChange"], deseq2_d12[my_gene,"log2FoldChange"],deseq2_d18[my_gene,"log2FoldChange"])
    my_data<- t(rbind(rc_data_subset, BF))
    colnames(my_data)<- c("gRNA#1", "log2FC")
    Days<- c(0, 6,12,18)
    my_data<- cbind(Days, my_data)
    my_data<- data.frame(my_data)
    rc_vec<- c(my_data[,2], my_data[,3] )
    day_vec<- rep(c(0,6,12,18),4)
    type_vec<- rep(c("gRNA_1"), each=4)
    my_data1<- data.frame(cbind(type_vec, day_vec, rc_vec) )
    my_data1[,"rc_vec"]<- as.numeric(as.character(my_data1[,"rc_vec"]))
    my_data_BF<- cbind(rep("log2FC",4), c(0,6,12,18), my_data[,3])
    colnames(my_data_BF)<- c("type_vec", "day_vec", "rc_vec")
    my_data2<- rbind(my_data1,my_data_BF)
#https://stackoverflow.com/questions/6142944/how-can-i-plot-with-2-different-y-axes 
    sd_gene_sgrna<- sd_data[my_gene,]
    gene_sgrna<- rc_data[my_gene,]
#Set ymax as max of up/down/gene
if(i%%3==0){
   ymax<- max(rc_data[my_gene_list[c(i-2,i-1,i)],])+max(sd_data[my_gene_list[c(i-2,i-1,i)],])+10
  }
if(i%%3==1){
   ymax<- max(rc_data[my_gene_list[c(i,i+1,i+2)],])+max(sd_data[my_gene_list[c(i,i+1,i+2)],])+10
  }
if(i%%3==2){
   ymax<- max(rc_data[my_gene_list[c(i-1,i,i+1)],])+max(sd_data[my_gene_list[c(i-1,i,i+1)],])+10
  }
    ymax<- round_any(ymax, 100, f=ceiling)
    interval<- 100;
    if(ymax > 500){interval<- 200};
    if(ymax > 1000){interval<- 500};
    plot(1:4, gene_sgrna, pch=21, xlab="", ylab="", type="l", col="brown", axes=F, main=my_gene, ylim=c(0,ymax), cex.main=3, lwd=7)    
#   segments(2:4, gene_sgrna[-1] -sd_gene_sgrna,2:4,gene_sgrna[-1]+sd_gene_sgrna,col="brown",lwd=7);
#    segments(2:4 -epsilon, gene_sgrna[-1] -sd_gene_sgrna, 2:4 +epsilon, gene_sgrna[-1] - sd_gene_sgrna, col="brown", lwd=7);
#    segments(2:4 -epsilon, gene_sgrna[-1] +sd_gene_sgrna, 2:4 +epsilon, gene_sgrna[-1] + sd_gene_sgrna, col="brown", lwd=7);
    axis(2, ylim=c(0, ymax), col="black", las=1, at=seq(0, ymax, by=interval), cex.axis=1.5);
    mtext(paste0("Normalized\nread counts"), side=2, line=3.5, cex=1.8)
    par(new=T)
    plot(1:4, as.numeric(as.character(my_data_BF[,"rc_vec"])), pch=23, xlab="", ylab="", type="l", col="grey40", axes=F, ylim=c(-4,4), lwd=7)
    axis(4, ylim=c(-4, 4), col="grey40", las=1, col.axis="grey40", at=seq(-4, 4, by=2), cex.axis=1.5 );
    mtext(bquote(paste('Log'['2'],'FC')), side=4, line=3.6, col="grey40", cex=1.8)             
    axis(1, at=1:4, labels=c(0,6,12,18), 10, cex.axis=1.5)
    mtext("Days", side=1, col="black", line=2.9, cex=1.8)
    }
dev.off();


}
####end study_gene

study_gene();


