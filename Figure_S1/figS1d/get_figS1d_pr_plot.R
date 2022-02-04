#https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
#library(PRROC)


data_d6<- read.table("pr_d6_vs_d0", sep="\t", header=T, row.names=1)

data_d12<- read.table("pr_d12_vs_d0", sep="\t", header=T, row.names=1)

data_d18<- read.table("pr_d18_vs_d0", sep="\t", header=T, row.names=1)




pdf("figS1d_bagel_pr_plot.pdf")

par(mfrow=c(2,2))

plot(data_d6[,"Recall"], data_d6[,"Precision"], type='l', lwd=1, col='brown', xlab="Recall", ylab="Precision", main="day 6 vs day 0")

plot(data_d12[,"Recall"], data_d12[,"Precision"], type='l', lwd=1, col='brown', xlab="Recall", ylab="Precision", main="day 12 vs day 0")

plot(data_d18[,"Recall"], data_d18[,"Precision"], type='l', lwd=1, col='brown', xlab="Recall", ylab="Precision", main="day 18 vs day 0")

dev.off();

