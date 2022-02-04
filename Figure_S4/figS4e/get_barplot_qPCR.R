library(ggplot2)
library(gridExtra)

se <- function(x) sqrt(var(x)/length(x))


data<- read.table("qPCR_data", sep="\t", header=T, row.names=1)


######
get_label<- function(my_pval){
print(my_pval);
if(my_pval < 0.001){
   return("***");
  }
if(my_pval > 0.001 & my_pval < 0.01){
   return("**");
  }
if(my_pval > 0.01 & my_pval < 0.05){
   return("*");
  }
}
######



#####################
get_barplot<- function(my_gene, my_pval_1, my_pval_2, my_color){

my_data<- data[,grep(my_gene, colnames(data))]

ctrl_mean<- mean(my_data[,1])

relative_data<- my_data/ctrl_mean

t_data<- t(relative_data)

mean_vec<- as.numeric(apply(t_data, 1, function(x) mean(x)) )
se_vec<- apply(t_data, 1, function(x) se(x))

t_data1<- data.frame(cbind(t_data, mean_vec, se_vec))
colnames(t_data1)[4:5]<- c("MEAN", "SE")

color_vec<- c("grey60", rep(my_color,2))
p<- ggplot(t_data1, aes(x=factor(rownames(t_data1), levels=rownames(t_data1)) , y=MEAN, fill=color_vec, ymin=mean_vec-se_vec, ymax=mean_vec+se_vec) )
p<- p+ geom_bar(stat="identity", fill=color_vec)
p<- p+ geom_errorbar(width=0.2, position=position_dodge(0.05))

p<- p+ labs(y="Relative expression", x="", title=my_gene)

p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=12), axis.title=element_text(size=18) )

p<- p+ scale_x_discrete(labels=c("Control\nsgRNA", rep(paste0(my_gene, "\nsgRNA"),2)));

p<- p + geom_point(aes(y=Rep_1),position= position_jitterdodge(jitter.width=0.8))
p<- p + geom_point(aes(y=Rep_2),position= position_jitterdodge(jitter.width=0.8))
p<- p + geom_point(aes(y=Rep_3),position= position_jitterdodge(jitter.width=0.8))

p<- p+ annotate("text", size=10, x=2, y=max(t_data1[2,1:3])+0.05, label=get_label(my_pval_1));
p<- p+ annotate("text", size=10, x=3, y=max(t_data1[3,1:3])+0.05, label=get_label(my_pval_2));

p<- p+ ylim(0, 1.2)
p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

return(p);
}
##########################


p1<- get_barplot("TEAD1", 0.000174819, 0.000138908, "coral3");
p2<- get_barplot("SKP2", 0.020187053,0.003199302, "goldenrod3");
p3<- get_barplot("TCAF1", 0.005227982,0.013576238, "darkolivegreen3");
p4<- get_barplot("PTPN14", 0.000926346,0.001503466, "rosybrown3");
p5<- get_barplot("TET2", 0.041293211,0.011697001, "steelblue3");


pdf("qPCR_barplot.pdf", width=10.5);
grid.arrange(p1, p2, p3, p4, p5, ncol=3)

dev.off();













