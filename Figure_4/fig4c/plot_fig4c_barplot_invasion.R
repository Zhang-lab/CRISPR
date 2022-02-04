library(ggplot2)
library(tidyr)
library(ggbreak)

data<- read.table("invasion_assay.txt", sep="\t", row.names=1, header=T)
ctrl_mean<- mean(data[,1])
relative_data<- data/ctrl_mean
relative_data1<- pivot_longer(relative_data, cols=1:ncol(relative_data), names_to="Condition", values_to="Value")
relative_data2<- data.frame(relative_data1)
relative_data2[,2]<- as.numeric(as.character(relative_data2[,2]))

se <- function(x) sqrt(var(x)/length(x))

t_data<- t(relative_data)

mean_vec<- as.numeric(apply(t_data, 1, function(x) mean(x)) )
se_vec<- apply(t_data, 1, function(x) se(x))

t_data1<- data.frame(cbind(t_data, mean_vec, se_vec))
colnames(t_data1)[6:7]<- c("MEAN", "SE")

color_vec<- c("grey60", "coral3", "goldenrod3", "darkolivegreen3");

pdf("fig4c_invasion.pdf", width=6, height=7, onefile=F);

p<- ggplot(t_data1, aes(x=factor(rownames(t_data1), levels=c("WT", "KO1", "KO2", "KO3")) , y=MEAN, fill=color_vec, ymin=mean_vec-se_vec, ymax=mean_vec+se_vec) )
p<- p+ geom_bar(stat="identity", fill=color_vec)
p<- p+ geom_errorbar(width=0.2, position=position_dodge(0.05))
p<- p+ labs(y="Relative number of invading cells", x="", title="Invasion Assay")
p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=12), axis.title=element_text(size=18) )

p<- p+ scale_x_discrete(labels=c("WT", "TEAD1 KO#1", "TEAD1 KO#2", "TEAD1 KO#3"));

p<- p + geom_point(aes(y=R1),position= position_jitterdodge(jitter.width=0.8))
p<- p + geom_point(aes(y=R2),position= position_jitterdodge(jitter.width=0.8))
p<- p + geom_point(aes(y=R3),position= position_jitterdodge(jitter.width=0.8))
p<- p + geom_point(aes(y=R4),position= position_jitterdodge(jitter.width=0.8))
p<- p + geom_point(aes(y=R5),position= position_jitterdodge(jitter.width=0.8))

p<- p+ annotate("text", size=10, x=2, y=max(t_data1[2,1:3])+0.05, label="*")
p<- p+ annotate("text", size=10, x=3, y=max(t_data1[3,1:3])+0.05, label="*")
p<- p+ annotate("text", size=10, x=4, y=max(t_data1[4,1:3])+0.05, label="*")
p<- p+ ylim(0, 2.52)
p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p<- p+ scale_y_break(c(1.51, 2.49), scale=0.15)

p;
dev.off()









