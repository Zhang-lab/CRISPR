
library(tidyr)
library(ggplot2)
library(gridExtra)

##########
get_plot<- function(my_file_name){

data<- read.table(my_file_name, header=T, sep="\t")
title<- my_file_name
data_to_plot<- data.frame(pivot_longer(data, cols=2:5, names_to="Day", values_to="Percent") )

data_to_plot$Day<- factor(data_to_plot$Day, levels=c("Day0", "Day6", "Day12", "Day18"))
data_to_plot$Type<- factor(data_to_plot$Type, levels=c("WT", "In-frame", "Frameshift"));

p<- ggplot(data_to_plot, aes(fill=Type, x=Day, y=Percent)) + geom_bar(position="fill", stat="identity") 

p<- p+ labs(x="Day", y="Percentage") + ggtitle(paste0(title  ) );

p<- p+ theme(plot.title=element_text(hjust=0.5, size=18, face="bold"), axis.text=element_text(size=14), axis.title=element_text(size=18), legend.text=element_text(size=14), legend.title=element_text(size=16)  )

p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

p<- p+ scale_x_discrete(label=c("0", "6", "12", "18"))

p<- p+ geom_text(aes(label = paste0(format(Percent*100, digits=5),"%")), 
            position = position_stack(vjust = 0.5), size = 6)

p<- p+ scale_fill_manual(values=c("gray84", "wheat", "lightblue"))
p<- p+ scale_y_continuous(labels = scales::percent)
return(p);
}
##########


p1<- get_plot("PTPN14_R1");
p2<- get_plot("PTPN14_R2");


q1<- get_plot("TET2_R1");
q2<- get_plot("TET2_R2");


pdf("figS8e_percent_stacked_barchart.pdf", width=11.5, height=9.5);

grid.arrange(p1,q1, p2, q2, nrow=2);

dev.off();


