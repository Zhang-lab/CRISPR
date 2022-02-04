
library(tidyr)
library(ggplot2)
library(gridExtra)

library(plotly)


data<- read.table("figS5f", header=T) 

data_to_plot<- data.frame(pivot_longer(data, cols=2:5, names_to="Sample", values_to="Value"), stringsAsFactors=F)


data_to_plot$Value<- data_to_plot$Value/100
data_to_plot$Type<- as.character(data_to_plot$Type)
data_to_plot$Type[data_to_plot$Type=="G01"]<- "G0/G1";
data_to_plot$Type[data_to_plot$Type=="G2M"]<- "G2/M";
data_to_plot$Type<- factor(data_to_plot$Type, levels=c("G0/G1", "S", "G2/M"));
data_to_plot$Sample<- factor(data_to_plot$Sample, levels=c("WT", "R1", "R2", "R3"))

p<- ggplot(data_to_plot, aes(x=Sample, y=Value, fill=Sample)) + geom_bar(stat="identity")

p<- p+ facet_wrap(~Type);

p<- p+ theme(plot.title=element_text(hjust=0.5, size=18, face="bold"), axis.text=element_text(size=14), axis.title=element_text(size=18), legend.text=element_text(size=14), legend.title=element_text(size=16)  )

p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

p<- p+ scale_fill_manual(values=c("gray84", "wheat", "wheat", "wheat"))

p<- p+ labs(x="TEAD KO", y="Percentage");

p<- p+ scale_y_continuous(labels = scales::percent)

p<- p+  theme(strip.text.x = element_text(size = 20))


pdf("barplot_figS5f.pdf");
p;
dev.off();

