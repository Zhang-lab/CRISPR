#Plot PCA using DESeq2 package
library(ggfortify)
library(gridExtra)

data<- read.table("GSE138688_counts_table.txt", sep="\t", row.names=1, header=T)

hTSC_data<- read.table("EG_list", sep="\t", header=F)
hTSC_TF_data<- read.table("EG_TF", sep="\t", header=F)

hTSC_genes<- toupper(hTSC_data[,1])
hTSC_TF<- toupper(hTSC_TF_data[,1])

design_data<- read.table("design_table", sep="\t", header=T, row.names=1);

data1<- data[, rownames(design_data)]

cpm=apply(data1,2,function(x) x/sum(x)*10^6)
#Use genes w/ max CPM >1
cpm_gt1<- cpm[apply(cpm, 1, max) > 1,]
data2<- data1[rownames(cpm_gt1),]

################ Begin make_pca
make_pca<- function(my_data, my_title){

t_data2<- t(my_data);

pca_res <- prcomp(t_data2, scale. = TRUE)

cellType<- sapply(strsplit(colnames(my_data), "_"), tail,1)
cellLine<- sapply(strsplit(colnames(my_data), "_"), head,1)

df_out<- data.frame(pca_res$x)
df_out$cellType<- cellType
df_out$cellLine<- cellLine

percentage <- round(pca_res$sdev / sum(pca_res$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

p<- ggplot(df_out, aes(x=PC1, y=PC2, color=cellType, shape=cellLine)) + geom_point(size=5)
p<- p+ xlab(percentage[1]) + ylab(percentage[2])
p<- p+ labs(color="Cell type", shape="Cell line", title=my_title)

return(p);
}
################ end make_pca

p1<- make_pca(data2, paste0("Genes w/ max CPM>1, n=", nrow(data2)));

hTSC_data<- data2[rownames(data2) %in% hTSC_genes,]
hTSC_TF_data<- data2[rownames(data2) %in% hTSC_TF,]

p2<- make_pca(hTSC_data, paste0("EGs w/ max CPM>1, n=", nrow(hTSC_data)) );
p3<- make_pca(hTSC_TF_data, paste0("EG TFs w/ max CPM>1, n=", nrow(hTSC_TF_data)) );

pdf("fig1g_S3de_pca_EG_wDiffShape.pdf", width=9, height=7);

grid.arrange(p1, p2, p3, nrow=2)

dev.off();

############# Random sampling to get PCA

set.seed(42)

gene_random<- sample(1:nrow(cpm_gt1))
random_hTSC<- gene_random[1:nrow(hTSC_data)]
random_hTSC[1:3]
r1<- make_pca(cpm_gt1[random_hTSC,], paste0("Random to EG size n=", length(random_hTSC)) ); 

gene_random<- sample(1:nrow(cpm_gt1))
random_hTSC<- gene_random[1:nrow(hTSC_data)]
random_hTSC[1:3]
r2<- make_pca(cpm_gt1[random_hTSC,], paste0("Random to EG size n=", length(random_hTSC)) ); 

gene_random<- sample(1:nrow(cpm_gt1))
random_hTSC<- gene_random[1:nrow(hTSC_data)]
random_hTSC[1:3]
r3<- make_pca(cpm_gt1[random_hTSC,], paste0("Random to EG size n=", length(random_hTSC)) ); 


gene_random<- sample(1:nrow(cpm_gt1));
random_TF<- gene_random[1:nrow(hTSC_TF_data)]
random_TF[1:3]
s1<- make_pca(cpm_gt1[random_TF,], paste0("Random to EG TF size n=", length(random_TF)));

gene_random<- sample(1:nrow(cpm_gt1));
random_TF<- gene_random[1:nrow(hTSC_TF_data)]
random_TF[1:3]
s2<- make_pca(cpm_gt1[random_TF,], paste0("Random to EG TF size n=", length(random_TF)));

gene_random<- sample(1:nrow(cpm_gt1));
random_TF<- gene_random[1:nrow(hTSC_TF_data)]
random_TF[1:3]
s3<- make_pca(cpm_gt1[random_TF,], paste0("Random to EG TF size n=", length(random_TF)));


pdf("random_PCA_EGnTF.pdf", width=16, height=9);

grid.arrange(r1,r2,r3,s1, s2, s3, nrow=2)

dev.off()











