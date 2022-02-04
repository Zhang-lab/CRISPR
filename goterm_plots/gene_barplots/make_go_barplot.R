library(ggplot2)

library(RColorBrewer)
library(RPMG)

library(Hmisc)


args=commandArgs(T)

tag<- args[1];

data<- read.table(tag, header=T, sep="\t")

data1<- data.frame(cbind(1:nrow(data),as.character(data[,"Name",]),-log10(as.numeric(as.character(data[,"q.value.FDR.B.H"])))) ) 

colnames(data1)<- c("pos", "name", "log10Q");

data1[,"pos"]<- factor(as.numeric(as.character(data1[,"pos"])), levels=nrow(data):1)
data1[,"log10Q"]<- as.numeric(as.character(data1[,"log10Q"]))
data1[,"name"]<- capitalize(as.character(data1[,"name"]) )

p<- ggplot(data1, aes(x=pos, y=log10Q) ) + geom_bar(stat="identity", fill=alpha("steelblue",0.3) )

p<- p+ coord_flip()

p<- p+ geom_text(aes(label=name), y=0.01, hjust=0, size=6)

p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

p<- p+ labs(x="Name", y=expression(paste(-Log[10], "q-value",sep="")), title=gsub("_", " ",tag) )

p<- p+  scale_x_discrete(labels = NULL, breaks = NULL)

p<- p+ theme(legend.position="none", plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=14), axis.title=element_text(size=20) )

p<- p+ scale_y_continuous(limits=c(0,NA), expand=c(0,0))

pdf_filename<- paste0(tag, ".pdf");

pdf(pdf_filename, width=9) 
p

dev.off();







