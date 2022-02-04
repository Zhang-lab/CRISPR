library(ggplot2)
library(gridExtra)

se <- function(x) sqrt(var(x)/length(x))

data<- read.table("qPCR.txt", sep="\t", header=T, row.names=1)

data1<- data[-c(1:2),]


############
get_barplot_gene<- function(my_gene, my_title, tag1, tag2, tag3, my_color){

data_gene<- data1[,grep(my_gene, colnames(data1))]

mean_vec<- as.numeric(apply(data_gene, 2, function(x) mean(x)) )
se_vec<- apply(data_gene, 2, function(x) se(x))

t_data_gene<- data.frame(t(rbind(data_gene,mean_vec, se_vec)))
colnames(t_data_gene)[4:5]<- c("MEAN", "SE");

color_vec<- c("grey60", rep(my_color,3))

p<- ggplot(t_data_gene, aes(x=factor(rownames(t_data_gene), levels=rownames(t_data_gene)), y=MEAN, fill=color_vec, ymin=mean_vec-se_vec, ymax=mean_vec+se_vec))
p<- p+ geom_bar(stat="identity", fill=color_vec)
p<- p+ geom_errorbar(width=0.2, position=position_dodge(0.05))

p<- p+ labs(y="Relative expression", x="", title=my_title)
p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=13), axis.title=element_text(size=18) )

p<- p+ scale_x_discrete(labels=c("WT", paste0("TEAD1\nKO #1"), paste0("TEAD1\nKO #2"), paste0("TEAD1\nKO #3")));

p<- p + geom_point(aes(y=Rep_1),position= position_jitterdodge(jitter.width=0.8))
p<- p + geom_point(aes(y=Rep_2),position= position_jitterdodge(jitter.width=0.8))
p<- p + geom_point(aes(y=Rep_3),position= position_jitterdodge(jitter.width=0.8))

p<- p+ annotate("text", size=10, x=2, y=max(t_data_gene[2,1:3])+0.05, label=tag1);
p<- p+ annotate("text", size=10, x=3, y=max(t_data_gene[3,1:3])+0.05, label=tag2);
p<- p+ annotate("text", size=10, x=4, y=max(t_data_gene[4,1:3])+0.05, label=tag3);

p<- p+ ylim(0, max(data_gene)+0.14)
p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

return(p);
}
############

p1<- get_barplot_gene("TEAD4", "TEAD4", "**", "**", "**", "coral3");
p2<- get_barplot_gene("SKP2", "SKP2", "***", "*", "***", "goldenrod3");
p3<- get_barplot_gene("KRT8", "KRT8", "", "**", "**", "darkolivegreen3");
p4<- get_barplot_gene("SDC1", "SDC1", "**", "***", "*", "rosybrown3");
p5<- get_barplot_gene("CGB", "CGB", "***", "***", "", "steelblue3");
p6<- get_barplot_gene("HLAG", "HLA-G", "", "", "", "darkorchid3");

pdf("figS5d_qPCR_barplot.pdf", width=11, height=6.1);
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3)
dev.off();


