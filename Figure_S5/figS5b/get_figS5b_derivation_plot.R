
library(ggplot2)
library(tidyr)

library(gridExtra)


data<- read.table("count_data", header=T, row.names=1)

data1<- cbind.data.frame(data, rownames(data))

colnames(data1)[ncol(data1)]<- "type"

input_data<- data.frame(pivot_longer(data1[1:4,], cols=1:7, names_to="Day", values_to="Count"))

input_data1<- cbind.data.frame(input_data, data.frame(pivot_longer(data1[5:8,], cols=1:7, names_to="Day", values_to="SD"))[,"SD"])

colnames(input_data1)[4]<- "SD"

input_data1$Day<- gsub("Day", "",input_data1$Day)


input_data1$Day<- factor(input_data1$Day, levels=c(0, 7, 9, 13, 16, 19, 25))

input_data1$Day<- as.numeric(as.character(input_data1$Day))

input_data1$type<- factor(input_data1$type, levels=c("WT", "TEAD1_KO1", "TEAD1_KO2", "TEAD1_KO3"))

input_data2<- na.omit(input_data1)

p<- ggplot(input_data2, aes(x=Day, y=Count, group=type, color=type)) + geom_line() + geom_point() 

p<- p+ geom_errorbar(aes(ymin=Count -SD, ymax=Count +SD), width=2.5, position=position_dodge(0.04))

p<- p+ theme_classic()

p<- p+ labs(y="Cumulative cell count")

y_max<- max(input_data1$Count)+ max(input_data1$SD)+100

p<- p+ scale_x_continuous(expand = c(0, 0), breaks=c(0, 7, 9,13, 16, 19, 25), limits=c(0, 26))

p <- p + labs(color="")

p<- p+ scale_colour_discrete(labels=c("WT", "TEAD1 KO #1", "TEAD1 KO #2", "TEAD1 KO #3"))



p1<- p+ scale_y_continuous(expand = c(0, 0), limits = c(0, y_max))  


p2<- p+ scale_y_continuous(trans='log10', expand = c(0, 0) )


pdf("figS5b_hTSC_derivation_cell_counts.pdf", height=5, width=5);
grid.arrange( p2, nrow=1)

dev.off();



