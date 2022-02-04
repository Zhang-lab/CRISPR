library(gridExtra)
library(ggplot2)

cnt_data<- read.table("merged_cnt_rup_on_all_peak.bed", header=T)

atac_data<- read.table("ATAC_counts_table.txt", header=T, row.names=1)

atac_cpm<- apply(atac_data[,-1], 2,function(x) x/sum(x)*10^6)

atac_cpm_mean<- cbind(rowMeans(atac_cpm[,c("naive_hESC_1", "naive_hESC_2", "naive_hESC_3", "naive_hESC_1")]), rowMeans(atac_cpm[,c("AN_hTSC_1", "AN_hTSC_2", "H9_hTSC_1", "H9_hTSC_2")]) )
 
colnames(atac_cpm_mean)<- c("naive", "TSC");

rownames(cnt_data)<- paste0(cnt_data[,1], ":", cnt_data[,2], "-", cnt_data[,3])
cnt_cpm<- rowMeans(apply(cnt_data[,4:6], 2, function(x) x/sum(x)*10^6))

cnt_atac_intersect<- read.table("CnT_atac_intersect", header=F)

naive_DAR_data<- read.table("naive_uniq_data_naive.bg", header=F)
TSC_DAR_data<- read.table("TSC_uniq_data_TSC.bg", header=F)
shared_atac_data<- read.table("shared_data_naive.bg", header=F)

rownames(naive_DAR_data)<-paste0(naive_DAR_data[,1],":",naive_DAR_data[,2], "-", naive_DAR_data[,3])
rownames(TSC_DAR_data)<- paste0(TSC_DAR_data[,1], ":", TSC_DAR_data[,2], "-", TSC_DAR_data[,3])
rownames(shared_atac_data)<- paste0(shared_atac_data[,1], ":", shared_atac_data[,2], "-", shared_atac_data[,3])

type<- rep("", nrow(cnt_atac_intersect)) 

for(i in 1:nrow(cnt_atac_intersect)){

    tag<-paste0(cnt_atac_intersect[i,4], ":", cnt_atac_intersect[i,5],"-",cnt_atac_intersect[i,6]);
 
    if(tag %in% rownames(naive_DAR_data)){
       type[i]<- "Group 3"
      }
   if(tag %in% rownames(TSC_DAR_data)){
       type[i]<- "Group 1"
      }
   if(tag %in% rownames(shared_atac_data)){
       type[i]<- "Group 2"
      }
   }

cnt_tags<- paste0(cnt_atac_intersect[,1], ":", cnt_atac_intersect[,2], "-", cnt_atac_intersect[,3])
atac_tags<- paste0(cnt_atac_intersect[,4], "_", cnt_atac_intersect[,5], "_", cnt_atac_intersect[,6])

cnt_atac_intersect1<- cbind(cnt_atac_intersect, type, cnt_cpm[cnt_tags], atac_cpm_mean[atac_tags,1:2])

colnames(cnt_atac_intersect1)[8:10]<- c("cnt_cpm", "atac_naive_cpm", "atac_TSC_cpm");

##################### get scatter density plot in ggplot
get_scatter_plot<- function(my_data, X1, X2, x1_lab, x2_lab, my_title, text_dist){

x1<- log2(1+X1)
x2<- log2(1+X2)

X1_group1<- X1[which(my_data[,"type"]=="Group 1")]
X1_group2<- X1[which(my_data[,"type"]=="Group 2")]
X2_group1<- X2[which(my_data[,"type"]=="Group 1")]
X2_group2<- X2[which(my_data[,"type"]=="Group 2")]

R_group1<- round(sqrt(summary(lm(X1_group1 ~ X2_group1))$adj.r.squared), 3)
R_group2<- round(sqrt(summary(lm(X1_group2 ~ X2_group2))$adj.r.squared), 3)

#Add R-square
mod1<- lm(x1 ~ x2)
modsum<- summary(mod1);
r2<- modsum$adj.r.squared
r<- sqrt(r2)
mylabel = bquote(italic(R) == .(format(r, digits=3, format="f")))
R<- round(r, 3)
scale_min<- min(x1,x2);
scale_max<- max(x1,x2);

p<- ggplot(my_data) + geom_point(aes(x=x1, y=x2, color=type)) 
p<- p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

p<- p+ labs(y=x2_lab, x=x1_lab, title=my_title)
p<- p+ theme(plot.title=element_text(hjust=0.5, size=25, face="bold"), axis.text=element_text(size=12), axis.title=element_text(size=18) )
p<- p+ guides(color=guide_legend(title="Group"))

print(mylabel)

return(p)
}
#####################

s1<- get_scatter_plot(cnt_atac_intersect1, cnt_atac_intersect1[,"cnt_cpm"], cnt_atac_intersect1[,"atac_naive_cpm"], "TEAD1 CnT signal", "ATAC naive data", "log2(CPM+1)", 1.2)

s2<- get_scatter_plot(cnt_atac_intersect1, cnt_atac_intersect1[,"cnt_cpm"], cnt_atac_intersect1[,"atac_TSC_cpm"], "TEAD1 CnT signal", "hTSC ATAC-seq signal", "log2(CPM+1)", 1.2)



pdf("fig3h_cnt_atac_log2cpm.pdf", width=8);
s2;
dev.off();

