#2020-12-15 With replicate check between Day0 replicates
#2021-02-20 Polish scatter density plot w/ ggplot

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

data<- read.table("bagel_norm_count_rep_wDay0", sep="\t", row.names=1, header=T)

Gene<- gsub("_1|_2|_3|_4|Br_", "", rownames(data) )
data1<- cbind(Gene, data)

mean_table<- aggregate(data1[,2:11], list(data1$Gene), mean)

day0_mean<- mean_table[,"day0"]


##################### get scatter density plot in ggplot
get_density_plot<- function(my_data, x1, x2, x1_lab, x2_lab, my_title, text_dist){

#Add R-square
mod1<- lm(x1 ~ x2)
modsum<- summary(mod1);
r2<- modsum$adj.r.squared
r<- sqrt(r2)
#mylabel = bquote(italic(R) == .(format(r, digits=3)))
mylabel = bquote(italic(R) == .(format(r, digits=3, format="f")))
R<- round(r, 3)
scale_min<- min(x1,x2);
scale_max<- ceiling(max(x1,x2) );
my_data$Density<- get_density(x1, x2, n=200)

p<- ggplot(my_data) + geom_point(aes(x=x1, y=x2, color=Density)) + scale_color_viridis()
p<- p+ xlim(scale_min, scale_max) + ylim(scale_min, scale_max)
p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

p<- p+ labs(y=x2_lab, x=x1_lab, title=my_title)
p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=28, face="bold"), axis.text=element_text(size=18), axis.title=element_text(size=21) )
print(mylabel)

p<- p+ annotate("text", size=8, x=scale_max-text_dist, y=scale_min, label=bquote(italic(R) == .(formatC(r, digits=3, format="f"))) )
p;

return(p)
}
#####################

log2_mean_table<- log2(1+ mean_table[,-1]);

p1<- get_density_plot(log2_mean_table, log2_mean_table[,"Day_6_Rep_1"], log2_mean_table[,"Day_6_Rep_2"], "Replicate 1", "Replicate 2", "Day 6", 1)
p2<- get_density_plot(log2_mean_table, log2_mean_table[,"Day_12_Rep_1"], log2_mean_table[,"Day_12_Rep_2"], "Replicate 1", "Replicate 2", "Day 12", 1)
p3<- get_density_plot(log2_mean_table, log2_mean_table[,"Day18Rep1"], log2_mean_table[,"Day18Rep2"], "Replicate 1", "Replicate 2", "Day 18", 1)

p4<- get_density_plot(log2_mean_table, log2_mean_table[,"Day0_Rep_1"], log2_mean_table[,"Day0Rep1"], "Replicate 1", "Replicate 2", "Day 0", 0.5)
p5<- get_density_plot(log2_mean_table, log2_mean_table[,"Day0_Rep_1"], log2_mean_table[,"Day0Rep2"], "Replicate 1", "Replicate 3", "Day 0", 0.5)
p6<- get_density_plot(log2_mean_table, log2_mean_table[,"Day0Rep1"], log2_mean_table[,"Day0Rep2"], "Replicate 2", "Replicate 3", "Day 0", 0.5)



pdf("figS1a_replicate_norm_log2readcount_day6_12_18.pdf", width=21)
grid.arrange(p1, p2,p3, nrow=1)

dev.off();

pdf("figS1a_replicate_norm_log2readcount_day0.pdf", width=21)
grid.arrange(p4, p5,p6, nrow=1)

dev.off();


