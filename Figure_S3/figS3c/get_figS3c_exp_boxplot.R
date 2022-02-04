library(ggplot2)
library(gridExtra)


rc_data<- read.table("merged_gene_name_expression.txt", sep="\t", row.names =1, header=T)

data<- read.table("bagel_d18vsd0_wClass.txt", sep="\t", header=T)
EG_data<- data[data[,"mean_RPKM_BT5_hTSC"]>0,]

CEG_data<- read.table("Hart_et_al_2017_CEGs.txt", sep="\t", header=T)
NEG_data<- read.table("Hart_et_al_2014_nonessential_genes.txt", sep="\t", header=T)



readcountsRaw<- rc_data[,-1]
geneLength<- rc_data[,1]

cpm=apply(readcountsRaw,2,function(x) x/sum(x)*10^6)
rpkm=apply(readcountsRaw,2,function(x) x/sum(x))*10^9/geneLength

BT5_mean<- rowMeans(rpkm[,c("BT5_hTSC_1", "BT5_hTSC_2")])

rpkm1<- cbind(rpkm, BT5_mean)

# rpkm+1
log2rpkm<- log2(rpkm1+1)

log2rpkm_unannotated<- log2rpkm[!(rownames(log2rpkm) %in% union(union(EG_data[,1], CEG_data[,1]), NEG_data[,1])),]

####################






####################

log2rpkm_EGs<- log2rpkm[rownames(log2rpkm)%in% EG_data[,1],]
log2rpkm_CEGs<- log2rpkm[rownames(log2rpkm)%in% CEG_data[,1],]
log2rpkm_NEGs<- log2rpkm[rownames(log2rpkm)%in% NEG_data[,1],]

rownames(log2rpkm_EGs)<- paste0("EG_",rownames(log2rpkm_EGs) )
rownames(log2rpkm_CEGs)<- paste0("CEG_",rownames(log2rpkm_CEGs) )

BT5_data<- data.frame(rbind(cbind("EG", log2rpkm_EGs[,"BT5_mean"]), cbind("CEG", log2rpkm_CEGs[,"BT5_mean"]), cbind("NEG", log2rpkm_NEGs[,"BT5_mean"]), cbind("Others", log2rpkm_unannotated[,"BT5_mean"])))
AN_naive_data<- data.frame(rbind(cbind("EG", log2rpkm_EGs[,"AN_naive_hTSC"]), cbind("CEG", log2rpkm_CEGs[,"AN_naive_hTSC"]), cbind("NEG", log2rpkm_NEGs[,"AN_naive_hTSC"]), cbind("Others", log2rpkm_unannotated[,"AN_naive_hTSC"])))
H9_naive_data<- data.frame(rbind(cbind("EG", log2rpkm_EGs[,"H9_naive_hTSC"]), cbind("CEG", log2rpkm_CEGs[,"H9_naive_hTSC"]), cbind("NEG", log2rpkm_NEGs[,"H9_naive_hTSC"]), cbind("Others", log2rpkm_unannotated[,"H9_naive_hTSC"])))

colnames(BT5_data)<- c("Type", "log2rpkm");
colnames(AN_naive_data)<- c("Type", "log2rpkm");
colnames(H9_naive_data)<- c("Type", "log2rpkm");
BT5_data[,"log2rpkm"]<- as.numeric(as.character(BT5_data[,"log2rpkm"]))
AN_naive_data[,"log2rpkm"]<- as.numeric(as.character(AN_naive_data[,"log2rpkm"]))
H9_naive_data[,"log2rpkm"]<- as.numeric(as.character(H9_naive_data[,"log2rpkm"]))

BT5_data[,"Type"]<- factor(BT5_data[,"Type"], levels=c("Others", "NEG", "CEG", "EG"));
AN_naive_data[,"Type"]<- factor(AN_naive_data[,"Type"], levels=c("Others", "NEG", "CEG", "EG"));
H9_naive_data[,"Type"]<- factor(H9_naive_data[,"Type"], levels=c("Others", "NEG", "CEG", "EG"));

p1<- ggplot(BT5_data, aes(x=Type, y=log2rpkm, color=Type)) + geom_boxplot();
p1<- p1+ coord_flip()
p1<- p1+ theme(legend.position="none", plot.title=element_text(hjust=0.5, vjust=-1, size=18, face="bold"), axis.text=element_text(size=13), axis.title=element_text(size=16) )
p1<- p1+ labs(title="BT5 hTSCs", x="", y=expression(paste(Log[2], "RPKM",sep="")))
p1<- p1+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1<- p1+ scale_color_manual(values=c("#999999", "#E69F00", "#8b1a1a", "#56B4E9"))
p1<- p1+ theme(plot.margin = unit(c(0,0.5,0,0), "lines"))
p1<- p1+ scale_x_discrete(labels=c("Others"="Other", "NEG"="Nonessential\ngenes", "CEG"="Core EG", "EG"= "hTSC EG"));


p2<- ggplot(AN_naive_data, aes(x=Type, y=log2rpkm, color=Type)) + geom_boxplot();
p2<- p2+ coord_flip()
p2<- p2+ theme(legend.position="none", plot.title=element_text(hjust=0.5, vjust=-1, size=18, face="bold"), axis.text=element_text(size=13), axis.title=element_text(size=16) )
p2<- p2+ labs(title="AN naive hTSCs", x="", y=expression(paste(Log[2], "RPKM",sep="")))
p2<- p2+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2<- p2+ scale_color_manual(values=c("#999999", "#E69F00", "#8b1a1a", "#56B4E9"))
p2<- p2+ theme(plot.margin = unit(c(0,0.5,0,0), "lines"))
p2<- p2+ scale_x_discrete(labels=c("Others"="Other", "NEG"="Nonessential\ngenes", "CEG"="Core EG", "EG"= "hTSC EG"));


p3<- ggplot(H9_naive_data, aes(x=Type, y=log2rpkm, color=Type)) + geom_boxplot();
p3<- p3+ coord_flip()
p3<- p3+ theme(legend.position="none", plot.title=element_text(hjust=0.5, vjust=-1,size=18, face="bold"), axis.text=element_text(size=13), axis.title=element_text(size=16) )
p3<- p3+ labs(title="H9 naive hTSCs", x="", y=expression(paste(Log[2], "RPKM",sep="")))
p3<- p3+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3<- p3+ scale_color_manual(values=c("#999999", "#E69F00", "#8b1a1a", "#56B4E9"))
p3<- p3+ theme(plot.margin = unit(c(0,0.5,0,0), "lines"))
p3<- p3+ scale_x_discrete(labels=c("Others"="Other", "NEG"="Nonessential\ngenes", "CEG"="Core EG", "EG"= "hTSC EG"));


pdf("figS3c_exp_boxplot.pdf", width=9.5, height=3.2);
grid.arrange(p1, p3, p2, nrow=1)


dev.off();



