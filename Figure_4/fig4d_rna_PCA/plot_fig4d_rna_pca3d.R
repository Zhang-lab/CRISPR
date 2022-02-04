
library(factoextra)
library(ggfortify)
library(gridExtra)
library(dplyr)
library(magrittr)
library(plotly)
library(htmlwidgets)
library(data.table);


data<- read.table("rna_norm_result_k1", header=T, row.names=1)

d<- apply(data[,-1]+1,2,function(x) x/sum(x)*10^6)
t_data1<-t(log2(d))

condition<- rep("KO", nrow(t_data1))

condition[grepl("wt", rownames(t_data1))]<- "WT"

celltype<- sapply(strsplit(rownames(t_data1), "_"), head, 1)


#######generate_3Dpca
generate_3Dpca<- function(my_pca_plot_data, my_field, field_title, my_pca){

title<- paste0("By ", field_title);


sdev_df <- data.frame(components = paste0("PC", 1:length(my_pca$sdev)), value = my_pca$sdev ** 2)
sdev_df$value <- sdev_df$value / sum(sdev_df$value) * 100
pc1_perc <- sdev_df$value[1] %>% round(., digits = 2)
pc2_perc <- sdev_df$value[2] %>% round(., digits = 2)
pc3_perc <- sdev_df$value[3] %>% round(., digits = 2)



my_pca_plot_data$field<- as.factor(my_field)
my_pca_plot_data$condition<- as.factor(my_pca_plot_data$condition)
my_pca_plot_data$sample_id<- rownames(my_pca_plot_data)

x<- plot_ly(my_pca_plot_data, x=~PC1, y=~PC2, z=~PC3, color=~field, symbol=~condition, symbols=c("x", "circle"),size=I(150) ) %>%
 add_markers(text= ~sample_id) %>%
 layout(title= title,
 scene=list(xaxis= list(title= paste0("PC1 (", unique(pc1_perc), "%)") ),
 yaxis= list(title=paste0("PC2 (", unique(pc2_perc), "%)") ),
 zaxis= list(title=paste0("PC3 (", unique(pc3_perc), "%)") ) )
)

htmlwidgets::saveWidget(as_widget(x), selfcontained = F, paste0("/home/shuhua/Public_html/theunissen_CRISPR/hTSC_EVT_STB_atac_rna_pca/rna_pca_k1_by_",field_title,".html") )
}
#############===end generate_3Dpca

pca_res<- prcomp(t_data1, scale. =T)

pca_plot_data<- data.frame(pca_res$x[,1:3]);

pca_plot_data$condition<- condition
pca_plot_data$celltype<- celltype

generate_3Dpca(pca_plot_data, as.vector(pca_plot_data[,"celltype"]), "celltype", pca_res);



