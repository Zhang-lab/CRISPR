#2020-10-28 mean of counts to get FC
library(ggplot2)
library(plyr)
library(gridExtra)

nonessential_list<- read.table("Hart_et_al_2014_nonessential_genes.txt", header=T, sep="\t")

essential_list<- read.table("Hart_et_al_2017_CEGs.txt", header=T, sep="\t")

count_data<- read.table("mageck_bruello_count_mean.txt", header=T, sep="\t", row.names=1)

mean_table<- aggregate(count_data[,2:5], list(count_data$Gene), mean)
colnames(mean_table)[2:5]<- c("Day_0", "Day_6", "Day_12", "Day_18");

mean_table$Group.1<- gsub("1-Dec", "DEC1",mean_table$Group.1)
mean_table$Group.1<- gsub("10-Mar", "MARCH10",mean_table$Group.1)
mean_table$Group.1<- gsub("10-Sep", "SEPT10",mean_table$Group.1)
mean_table$Group.1<- gsub("11-Mar", "MARCH11",mean_table$Group.1)
mean_table$Group.1<- gsub("11-Sep", "SEPT11",mean_table$Group.1)
mean_table$Group.1<- gsub("12-Sep", "SEPT12",mean_table$Group.1)
mean_table$Group.1<- gsub("14-Sep", "SEPT14",mean_table$Group.1)
mean_table$Group.1<- gsub("15-Sep", "SEPT15",mean_table$Group.1)
mean_table$Group.1<- gsub("2-Mar", "MARC2",mean_table$Group.1)
mean_table$Group.1<- gsub("2-Sep", "SEPT2",mean_table$Group.1)
mean_table$Group.1<- gsub("3-Mar", "MARCH3",mean_table$Group.1)
mean_table$Group.1<- gsub("3-Sep", "SEPT3",mean_table$Group.1)
mean_table$Group.1<- gsub("4-Mar", "MARCH4",mean_table$Group.1)
mean_table$Group.1<- gsub("4-Sep", "SEPT4",mean_table$Group.1)
mean_table$Group.1<- gsub("5-Mar", "MARCH5",mean_table$Group.1)
mean_table$Group.1<- gsub("5-Sep", "SEPT5",mean_table$Group.1)
mean_table$Group.1<- gsub("6-Mar", "MARCH6",mean_table$Group.1)
mean_table$Group.1<- gsub("6-Sep", "SEPT6",mean_table$Group.1)
mean_table$Group.1<- gsub("7-Mar", "MARCH7",mean_table$Group.1)
mean_table$Group.1<- gsub("7-Sep", "SEPT7",mean_table$Group.1)
mean_table$Group.1<- gsub("8-Mar", "MARCH8",mean_table$Group.1)
mean_table$Group.1<- gsub("8-Sep", "SEPT8",mean_table$Group.1)
mean_table$Group.1<- gsub("9-Mar", "MARCH9",mean_table$Group.1)
mean_table$Group.1<- gsub("9-Sep", "SEPT9",mean_table$Group.1)



table_essential<- mean_table[mean_table$Group.1 %in% essential_list[,1],]
table_nonessential<- mean_table[mean_table$Group.1 %in% nonessential_list[,1],]

table_nontarget<- mean_table[grep("NO_CURENT", mean_table$Group.1),]

log2fc_d18_essential<- log2(table_essential$Day_18/table_essential$Day_0)
log2fc_d18_nonessential<- log2(table_nonessential$Day_18/table_nonessential$Day_0)

log2fc_d12_essential<- log2(table_essential$Day_12/table_essential$Day_0)
log2fc_d12_nonessential<- log2(table_nonessential$Day_12/table_nonessential$Day_0)

log2fc_d6_essential<- log2(table_essential$Day_6/table_essential$Day_0)
log2fc_d6_nonessential<- log2(table_nonessential$Day_6/table_nonessential$Day_0)

data_d18<- data.frame(rbind(cbind(log2fc_d18_essential, rep("Essential", length(log2fc_d18_essential) ) ), cbind(log2fc_d18_nonessential, rep("Nonessential", length(log2fc_d18_nonessential) ) ) ) )
colnames(data_d18)<- c("Log2FC", "Type")
data_d18[,"Log2FC"]<- as.numeric(as.character(data_d18[,"Log2FC"]))
mu_d18 <- ddply(data_d18, "Type", summarise, grp.mean=mean(Log2FC))

data_d12<- data.frame(rbind(cbind(log2fc_d12_essential, rep("Essential", length(log2fc_d12_essential) ) ), cbind(log2fc_d12_nonessential, rep("Nonessential", length(log2fc_d12_nonessential) ) ) ) )
colnames(data_d12)<- c("Log2FC", "Type")
data_d12[,"Log2FC"]<- as.numeric(as.character(data_d12[,"Log2FC"]))
mu_d12 <- ddply(data_d12, "Type", summarise, grp.mean=mean(Log2FC))

data_d6<- data.frame(rbind(cbind(log2fc_d6_essential, rep("Essential", length(log2fc_d6_essential) ) ), cbind(log2fc_d6_nonessential, rep("Nonessential", length(log2fc_d6_nonessential) ) ) ) )
colnames(data_d6)<- c("Log2FC", "Type")
data_d6[,"Log2FC"]<- as.numeric(as.character(data_d6[,"Log2FC"]))
mu_d6 <- ddply(data_d6, "Type", summarise, grp.mean=mean(Log2FC))



p1<- ggplot(data_d18, aes(x=Log2FC, fill=Type, color=Type )) + geom_density(alpha=0.2, size=1)
p1<- p1+geom_vline(data=mu_d18,aes(xintercept=grp.mean, color=Type),linetype="dashed",show.legend=F)
p1<- p1+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1<- p1+ theme(plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=12), axis.title=element_text(size=18), legend.title=element_blank(), legend.text=element_text(size=18) )
p1<- p1+ scale_fill_manual(values=c("brown", "grey60"))
p1<- p1+ scale_color_manual(values=c("brown", "grey60"))
p1<- p1+ labs(x=expression(paste(Log[2], "FC",sep="")), y="Density", title="Day 18")
p1<- p1+ theme(legend.position=c(0.2,0.9))


p2<- ggplot(data_d12, aes(x=Log2FC, fill=Type, color=Type )) + geom_density(alpha=0.2, size=1)
p2<- p2+geom_vline(data=mu_d12,aes(xintercept=grp.mean, color=Type),linetype="dashed",show.legend=F)
p2<- p2+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2<- p2+ theme(plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=12), axis.title=element_text(size=18), legend.title=element_blank(), legend.text=element_text(size=18) )
p2<- p2+ scale_fill_manual(values=c("brown", "grey60"))
p2<- p2+ scale_color_manual(values=c("brown", "grey60"))
p2<- p2+ labs(x=expression(paste(Log[2], "FC",sep="")), y="Density", title="Day 12")
p2<- p2+ theme(legend.position=c(0.2,0.9))


p3<- ggplot(data_d6, aes(x=Log2FC, fill=Type, color=Type )) + geom_density(alpha=0.2, size=1)
p3<- p3+geom_vline(data=mu_d6,aes(xintercept=grp.mean, color=Type),linetype="dashed",show.legend=F)
p3<- p3+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3<- p3+ theme(plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=12), axis.title=element_text(size=18), legend.title=element_blank(), legend.text=element_text(size=18) )
p3<- p3+ scale_fill_manual(values=c("brown", "grey60"))
p3<- p3+ scale_color_manual(values=c("brown", "grey60"))
p3<- p3+ labs(x=expression(paste(Log[2], "FC",sep="")), y="Density", title="Day 6")
p3<- p3+ theme(legend.position=c(0.2,0.9))


pdf("fig1c_S1c_FC_density.pdf", width=18, height=6);
grid.arrange(p1, p2, p3, nrow=1)
dev.off();



