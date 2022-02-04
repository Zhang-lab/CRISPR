#https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
library(ggplot2)
library(viridis)
library(gridExtra)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data_d6<- read.table("pr_d6_vs_d0", sep="\t", header=T, row.names=1)
data_d12<- read.table("pr_d12_vs_d0", sep="\t", header=T, row.names=1)
data_d18<- read.table("pr_d18_vs_d0", sep="\t", header=T, row.names=1)

##################### get scatter density plot in ggplot
get_density_plot<- function(my_data, x1, x2, x1_lab, x2_lab, my_title, text_dist){

#Add R-square
mod1<- lm(x1 ~ x2)
modsum<- summary(mod1);
r2<- modsum$adj.r.squared
r<- sqrt(r2)
mylabel = bquote(italic(R) == .(format(r, digits=3)))
R<- round(r, 3)
scale_min<- min(x1,x2); 
scale_max<- max(x1,x2); 
my_data$Density<- get_density(x1, x2, n=200)

p<- ggplot(my_data) + geom_point(aes(x=x1, y=x2, color=Density)) + scale_color_viridis()
p<- p+ xlim(scale_min, scale_max) + ylim(scale_min, scale_max)
p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

p<- p+ labs(y=x2_lab, x=x1_lab, title=my_title)
p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=28, face="bold"), axis.text=element_text(size=15.5), axis.title=element_text(size=21) )
mylabel

p<- p+ annotate("text", size=8, x=scale_max-text_dist, y=scale_min, label=bquote(italic(R) == .(format(r, digits=3))) )
return(p)
}
#####################

get_density_plot_fixed<- function(my_data, x1, x2, x1_lab, x2_lab, my_title){

#Add R-square
mod1<- lm(x1 ~ x2)
modsum<- summary(mod1);
r2<- modsum$adj.r.squared
r<- sqrt(r2)
mylabel = bquote(italic(R) == .(format(r, digits=3)))
R<- round(r, 3)
scale_min<- 0
scale_max<- 1
my_data$Density<- get_density(x1, x2, n=200)

p<- ggplot(my_data) + geom_point(aes(x=x1, y=x2, color=Density)) + scale_color_viridis()
p<- p+ xlim(scale_min, scale_max) + ylim(scale_min, scale_max)
p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

p<- p+ labs(y=x2_lab, x=x1_lab, title=my_title)
p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=28, face="bold"), axis.text=element_text(size=15.5), axis.title=element_text(size=21) )
mylabel

p<- p+ annotate("text", size=8, x=scale_max-0.1, y=scale_min, label=bquote(italic(R) == .(format(r, digits=3))) )

return(p)
}
#################


### For BF
BF_data<- data.frame(cbind(data_d6[,"BF"], data_d12[rownames(data_d6),"BF"], data_d18[rownames(data_d6),"BF"]) )
colnames(BF_data)<- c("d6", "d12", "d18") 

p1<- get_density_plot(BF_data, BF_data[,"d6"], BF_data[,"d12"], "Day 6", "Day 12", "BF",4);
p2<- get_density_plot(BF_data, BF_data[,"d6"], BF_data[,"d18"], "Day 6", "Day 18", "BF",4);
p3<- get_density_plot(BF_data, BF_data[,"d12"], BF_data[,"d18"], "Day 12", "Day 18", "BF",4);

pdf("figS2a_bagel_BF_density_plots.pdf", width=22.6, height=7);
grid.arrange(p1, p2,p3, nrow=1)
dev.off();


######## For FDR
FDR_data<- data.frame(cbind(data_d6[,"FDR"], data_d12[rownames(data_d6),"FDR"], data_d18[rownames(data_d6),"FDR"]) )
colnames(FDR_data)<- c("d6", "d12", "d18")

s1<- get_density_plot(FDR_data, FDR_data[,"d6"], FDR_data[,"d12"], "Day 6", "Day 12", "FDR", 0.05);
s2<- get_density_plot(FDR_data, FDR_data[,"d6"], FDR_data[,"d18"], "Day 6", "Day 18", "FDR",0.05);
s3<- get_density_plot(FDR_data, FDR_data[,"d12"], FDR_data[,"d18"], "Day 12", "Day 18", "FDR",0.05);

pdf("figS2a_bagel_FDR_density_plots.pdf", width=21.5);
grid.arrange(s1, s2,s3, nrow=1)
dev.off();

############## For Recall
Recall_data<- data.frame(cbind(data_d6[,"Recall"], data_d12[rownames(data_d6),"Recall"], data_d18[rownames(data_d6),"Recall"]) )
colnames(Recall_data)<- c("d6", "d12", "d18")

t1<- get_density_plot_fixed(Recall_data, Recall_data[,"d6"], Recall_data[,"d12"], "Day 6", "Day 12", "Recall");
t2<- get_density_plot_fixed(Recall_data, Recall_data[,"d6"], Recall_data[,"d18"], "Day 6", "Day 18", "Recall");
t3<- get_density_plot_fixed(Recall_data, Recall_data[,"d12"], Recall_data[,"d18"], "Day 12", "Day 18", "Recall");

pdf("figS2a_bagel_Recall_density_plots.pdf", width=21.5);
grid.arrange(t1, t2,t3, nrow=1)
dev.off();


############## For Precision
Precision_data<- data.frame(cbind(data_d6[,"Precision"], data_d12[rownames(data_d6),"Precision"], data_d18[rownames(data_d6),"Precision"]) )
colnames(Precision_data)<- c("d6", "d12", "d18")

q1<- get_density_plot(Precision_data, Precision_data[,"d6"], Precision_data[,"d12"], "Day 6", "Day 12", "Precision",0.1);
q2<- get_density_plot(Precision_data, Precision_data[,"d6"], Precision_data[,"d18"], "Day 6", "Day 18", "Precision", 0.1);
q3<- get_density_plot(Precision_data, Precision_data[,"d12"], Precision_data[,"d18"], "Day 12", "Day 18", "Precision", 0.1);

pdf("figS2a_bagel_Precision_density.pdf", width=21.5);
grid.arrange(q1, q2,q3, nrow=1)
dev.off();






