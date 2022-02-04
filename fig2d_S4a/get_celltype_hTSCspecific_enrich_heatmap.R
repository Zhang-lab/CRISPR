library(pheatmap)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(gplots)

fpkm_data<- read.table("GSE136447_555-samples-fpkm.txt",sep="\t", header=T, row.names=1)
fpkm_data1<- fpkm_data[grep("chr",fpkm_data[,"Reference"]),]

metadata<- read.table("41586_2019_1875_MOESM10_ESM.txt",sep="\t", row.names=1, header=T)

hTSC_data<- read.table("hTSC_specific_genes", header=F)
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
Days<- c("D10", "D12");

for(i in 1:length(Days)){
print(i) 
    metadata_day<- metadata[metadata[,"Day"]==Days[i],]

    CTB_cols<-which(colnames(fpkm_data2)%in%rownames(metadata_day[metadata_day["Group"]=="CTB",]) )

    rowMeans_CTB<- rowMeans(fpkm_data2[,CTB_cols])

    rowMeans_data<- data.frame(matrix(ncol=1, nrow=nrow(fpkm_data2)))
    rowMeans_data[,1]<- rowMeans_CTB
    colnames(rowMeans_data)[1]<- "CTB"

    cov_rows<- as.vector(which(rowMeans_CTB > mean_fpkm_cutoff))
    selected_rows<- cov_rows;

    if("EPI" %in% metadata_day[,"Group"]){

       EPI_cols<-which(colnames(fpkm_data2)%in%rownames(metadata_day[metadata_day["Group"]=="EPI",]))  
       rowMeans_EPI<- rowMeans(fpkm_data2[,EPI_cols]) 
       diff_rows<- as.vector(which(rowMeans_CTB > times_cutoff * rowMeans_EPI) )
       selected_rows<- intersect(selected_rows, diff_rows);
      
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
       diff_rows<- as.vector(which(rowMeans_CTB >times_cutoff * rowMeans_Pre) )
       selected_rows<- intersect(selected_rows, diff_rows);

       rowMeans_data<- cbind(rowMeans_data, rowMeans_Pre);
       colnames(rowMeans_data)[ncol(rowMeans_data)]<- "PrE"
      }
    rowMeans_data1<- rowMeans_data[selected_rows,]
    rowMeans_data2<- apply(rowMeans_data1, 2, function(x) as.numeric(x))
    rownames(rowMeans_data2)<- rownames(rowMeans_data1);

    log2_rowMeans_data2<- log2(rowMeans_data2 +1e-6)

    rownames(log2_rowMeans_data2)<- sapply(strsplit(rownames(log2_rowMeans_data2), "_"), tail,1)

    pdf(paste0("heatmap_", Days[i], "_CTBenrich_hTSCspecific.pdf"), height=8, width=2.5);

    colfunc<- colorRampPalette(c("blue", "white", "red"))

     pheatmap(log2_rowMeans_data2, cluster_cols=F, angle_col=45, main=paste0(gsub("D","",Days[i]), " d.p.f (n=", nrow(log2_rowMeans_data2),")"))

    dev.off();
    write.table(data.frame("name"=rownames(rowMeans_data2), rowMeans_data2), file=paste0(Days[i], "_CTBenriched_meanFpkm.txt"), quote=F, sep="\t", row.names=F, col.names=T); 
   }

