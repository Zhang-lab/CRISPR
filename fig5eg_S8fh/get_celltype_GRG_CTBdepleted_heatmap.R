library(pheatmap)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(gplots)

#library(scater)

fpkm_data<- read.table("GSE136447_555-samples-fpkm.txt",sep="\t", header=T, row.names=1)
fpkm_data1<- fpkm_data[grep("chr",fpkm_data[,"Reference"]),]

metadata<- read.table("41586_2019_1875_MOESM10_ESM.txt",sep="\t", row.names=1, header=T)

hTSC_data<- read.table("GRG_list", header=F)
hTSC_gene<- toupper(hTSC_data[,1])

fpkm_data2<- fpkm_data1[fpkm_data1[,"Gene.Name"] %in% hTSC_gene,]
rownames(fpkm_data2)<- paste0(fpkm_data2[,"Gene.ID"], "_",fpkm_data2[, "Gene.Name"])

Celltypes<- c("CTB", "PrE", "EPI");

times_cutoff<- 2;
mean_fpkm_cutoff<- 2;

metadata<- metadata[metadata[,"Group"] %in% Celltypes,]

Days<- unique(metadata[,"Day"])
#Day6 only has 1 celltype
Days<- Days[-grep("D6", Days)]
#Only D10 & D12
Days<-c("D10", "D12")

for(i in 1:length(Days)){
print(i) 
    metadata_day<- metadata[metadata[,"Day"]==Days[i],]

    CTB_cols<-which(colnames(fpkm_data2)%in%rownames(metadata_day[metadata_day["Group"]=="CTB",]) )

    rowMeans_CTB<- rowMeans(fpkm_data2[,CTB_cols])

    rowMeans_data<- data.frame(matrix(ncol=1, nrow=nrow(fpkm_data2)))
    rowMeans_data[,1]<- rowMeans_CTB
    colnames(rowMeans_data)[1]<- "CTB"

#    cov_rows<- as.vector(which(rowMeans_CTB > mean_fpkm_cutoff))
#    selected_rows<- cov_rows;

    selected_rows<- 1:length(rowMeans_CTB)

    if("EPI" %in% metadata_day[,"Group"]){

       EPI_cols<-which(colnames(fpkm_data2)%in%rownames(metadata_day[metadata_day["Group"]=="EPI",]))  
       rowMeans_EPI<- rowMeans(fpkm_data2[,EPI_cols]) 
       diff_rows<- as.vector(which(rowMeans_CTB * times_cutoff < rowMeans_EPI) )
#Select EPI > mean_fpkm_cutoff  
       cov_rows<- as.vector(which(rowMeans_EPI > mean_fpkm_cutoff))

       selected_rows<- intersect(intersect(selected_rows, diff_rows), cov_rows);
      
       rowMeans_data<- cbind(rowMeans_data, rowMeans_EPI);
       colnames(rowMeans_data)[ncol(rowMeans_data)]<- "EPI"
      }
    if("PrE" %in% metadata_day[,"Group"]){

       Pre_cols<-which(colnames(fpkm_data2)%in%rownames(metadata_day[metadata_day["Group"]=="PrE",]))
       if(length(Pre_cols) >1){
           rowMeans_Pre<- rowMeans(fpkm_data2[,Pre_cols]) 
         }else{
           rowMeans_Pre<- fpkm_data2[,Pre_cols]
           }
       diff_rows<- as.vector(which(rowMeans_CTB * times_cutoff < rowMeans_Pre) )
#Select PrE > mean_fpkm_cutoff
       cov_rows<- as.vector(which(rowMeans_Pre > mean_fpkm_cutoff)) 

       selected_rows<- intersect(intersect(selected_rows, diff_rows), cov_rows);

       rowMeans_data<- cbind(rowMeans_data, rowMeans_Pre);
       colnames(rowMeans_data)[ncol(rowMeans_data)]<- "PrE"
      }
    rowMeans_data1<- rowMeans_data[selected_rows,]
    rowMeans_data2<- apply(rowMeans_data1, 2, function(x) as.numeric(x))
    rownames(rowMeans_data2)<- rownames(rowMeans_data1);

    log2_rowMeans_data2<- log2(rowMeans_data2 +1e-6)

    write.table(data.frame("name"=rownames(rowMeans_data2), rowMeans_data2), file=paste0(Days[i], "_CTBdepleted_meanFpkm.txt"), quote=F, sep="\t", row.names=F, col.names=T); 

    rownames(log2_rowMeans_data2)<- sapply(strsplit(rownames(log2_rowMeans_data2), "_"), tail,1)
 
    pdf(paste0("heatmap_", Days[i], "_CTBdepleted.pdf"), width=7, height=8);

    par(mar=c(3, 7,3, 19))
    colfunc<- colorRampPalette(c("blue", "white", "red"))

    heatmap.2(log2_rowMeans_data2, scale="none", trace="none", col=colfunc(15), dendrogram="row", cexCol=1, main=paste0(Days[i], "\nMean log2(FPKM+1e-6)\n n=", nrow(log2_rowMeans_data2)), margins=c(4,14), cexRow=0.9)

#    pheatmap(log2_rowMeans_data2, cluster_cols=F, angle_col=45, main=paste0(gsub("D","",Days[i]), " d.p.f (n=", nrow(log2_rowMeans_data2),")"))

    dev.off();
   }



