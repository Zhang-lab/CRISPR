library(ggplot2)
library(tidyr)

data<- read.table("fig2g_data", sep="\t", row.names=1, header=T)

ctrl_mean<- mean(data[,1])

relative_data<- data/ctrl_mean

relative_data1<- pivot_longer(relative_data, cols=1:ncol(relative_data), names_to="Condition", values_to="Value")

relative_data2<- data.frame(relative_data1)
relative_data2[,2]<- as.numeric(as.character(relative_data2[,2]))

data_summary <- function(data, varname, groupnames){
require(plyr)
summary_func <- function(x, col){
c(mean = mean(x[[col]], na.rm=TRUE),
  sd = sd(x[[col]], na.rm=TRUE))
  }
data_sum<-ddply(data, groupnames, .fun=summary_func, varname)
ata_sum <- rename(data_sum, c("mean" = varname))
return(data_sum)
}

se <- function(x) sqrt(var(x)/length(x))

t_data<- t(relative_data)

mean_vec<- as.numeric(apply(t_data, 1, function(x) mean(x)) )
se_vec<- apply(t_data, 1, function(x) se(x))

t_data1<- data.frame(cbind(t_data, mean_vec, se_vec))
colnames(t_data1)[4:5]<- c("MEAN", "SE")

color_vec<- c("grey60", "coral3", "coral3", "goldenrod3", "goldenrod3", "darkolivegreen3", "darkolivegreen3")


pdf("fig2g.pdf", width=9, height=6);
p<- ggplot(t_data1, aes(x=factor(rownames(t_data1), levels=c("Ctrl_sgRNA", "SKP2_sgRNA_1", "SKP2_sgRNA_2", "TEAD1_sgRNA_1", "TEAD1_sgRNA_2", "TCAF1_sgRNA_1", "TCAF1_sgRNA_2")) , y=MEAN, fill=color_vec, ymin=mean_vec-se_vec, ymax=mean_vec+se_vec) )
p<- p+ geom_bar(stat="identity", fill=color_vec)
p<- p+ geom_errorbar(width=0.2, position=position_dodge(0.05))
p<- p+ labs(y="Relative cell count", x="", title="Day 6")

p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=12), axis.title=element_text(size=18) )

p<- p+ scale_x_discrete(labels=c("Control\nsgRNA", "SKP2\nsgRNA 1", "SKP2\nsgRNA 2", "TEAD1\nsgRNA 1", "TEAD1\nsgRNA 2", "TCAF1\nsgRNA 1", "TCAF1\nsgRNA 2") );

p<- p + geom_point(aes(y=Biological_Replicate_1),position= position_jitterdodge(jitter.width=0.8))
p<- p + geom_point(aes(y=Biological_Replicate_2),position= position_jitterdodge(jitter.width=0.8))
p<- p + geom_point(aes(y=Biological_Replicate_3),position= position_jitterdodge(jitter.width=0.8))

p<- p+ annotate("text", size=10, x=2, y=max(t_data1[2,1:3])+0.05, label="***")
p<- p+ annotate("text", size=10, x=3, y=max(t_data1[3,1:3])+0.05, label="***")
p<- p+ annotate("text", size=10, x=4, y=max(t_data1[4,1:3])+0.05, label="**")
p<- p+ annotate("text", size=10, x=5, y=max(t_data1[5,1:3])+0.05, label="***")
p<- p+ annotate("text", size=10, x=6, y=max(t_data1[6,1:3])+0.05, label="**")
p<- p+ annotate("text", size=10, x=7, y=max(t_data1[7,1:3])+0.05, label="*")
p<- p+ ylim(0, 1.2)
p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

dev.off()









