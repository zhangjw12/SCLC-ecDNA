library(imcRtools);library(tidyverse);library(stringr);library(readr)
library(cytomapper);library(CATALYST);library(scater);library(batchelor)
library(dittoSeq);library(patchwork);library(pheatmap);library(scran)
library(bluster);library(dittoSeq);library(viridis);library(BiocParallel)
mytheme = theme(plot.title = element_text(hjust = 0.5),
                text=element_text(size=10),
                axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
                axis.text.y = element_text(size = 8,color = "black"),
                legend.text = element_text(size=6,color = "black"),
                legend.title = element_text(size=8),
                legend.key.size = unit(11, "pt"))

#data
spe_fastMNN <- readRDS("/home/ylab/IMC/3_clean_data/spe_fastMNN.RDS")

######### spe_fastMNN dim: 40 103341
##Figure 5A
setwd("/home/ylab/SCLC/code_ocean/data")
pdf("../results/Fig5A.pdf",width = 220/25.4, height = 160/25.4)
dittoDimPlot(spe_fastMNN, var = "celltype",
 reduction.use = "UMAP_mnnCorrected", size = 0.2,do.label = TRUE,labels.repel = TRUE ,labels.size = 3,labels.highlight = F) +
 ggtitle("Cluster cell type on UMAP after correction")+
 scale_color_manual(values= c(  "T cell" = "#E69F00",
  "Ki67+ Tumor" = "#f87171",
  "HLA-DR+ Macrophage" = "#009E73",
  "PD-L1+ Macrophage" = "#F0E442",
  "Macrophage" = "#0072B2",
  "CD163+ Macrophage" = "#D55E00",
  "HIF1a+ Macrophage" = "#CC79A7",
  "B cell" = "#0b1e29",
  "Neutrophil" = "#AD7700",
  "Undefined" = "#565353",
  "Ecad+ Tumor" = "#FFBE2D",
  "VEGFA+ Tumor" = "#929292",
  "Endothelial" = "#850016",
  "Fibroblast" = "#007756",
  "COL1A1+ Fibroblast" = "#D5C711",
  "FAP+ Fibroblast" = "#A04700",
  "Myofibroblast" = "#B14380",
  "Tumor" = "#80C7EF"
))
dev.off()

##Figure 5B
cluster_mean <- aggregateAcrossCells(as(spe_fastMNN, "SingleCellExperiment"),
 ids = spe_fastMNN$celltype,
 statistics = "mean",
 use.assay.type = "normalized",
 subset_row = rowData(spe_fastMNN)$use_channel)
custom_colors <- c(
  "T cell" = "#E69F00",
  "Ki67+ Tumor" = "#f87171",
  "HLA-DR+ Macrophage" = "#009E73",
  "PD-L1+ Macrophage" = "#F0E442",
  "Macrophage" = "#0072B2",
  "CD163+ Macrophage" = "#D55E00",
  "HIF1a+ Macrophage" = "#CC79A7",
  "B cell" = "#0b1e29",
  "Neutrophil" = "#AD7700",
  "Undefined" = "#565353",
  "Ecad+ Tumor" = "#FFBE2D",
  "VEGFA+ Tumor" = "#929292",
  "Endothelial" = "#850016",
  "Fibroblast" = "#007756", 
  "COL1A1+ Fibroblast" = "#D5C711",
  "FAP+ Fibroblast" = "#A04700",
  "Myofibroblast" = "#B14380",
  "Tumor" = "#80C7EF"
)
set.seed((240440))
genes<-c("Pan-cytokeratin","VEGFA","Ki-67","E-cadherin","CD31","aSMA","FAP","CollagenI","CD15",
"CD68", "HIF1a","HLA-DR","PD-L1","CD163","CD3","CD20","CD45")
pdf("../results/Fig5B.pdf")
dittoHeatmap(cluster_mean,genes,
    heatmap.colors = colorRampPalette(c("darkblue", "white", "firebrick3"))(300),
 assay = "normalized", cluster_cols = T,cluster_rows = F,
 scale = "row",
 annot.by = c("celltype", "ncells"),
annotation_colors = list(celltype = custom_colors))
dev.off()

#landscape 5C
pdf(file.path("/mnt/c/Users/ylab/Desktop/New", paste0( "5_proportion.pdf")),
    width = 120/25.4, height = 100/25.4)
ggplot(df_counts_all, aes(x = proportion, y = sample_id, fill = cluster_celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cluster Cell Type Proportions",
       x = "Cluster Cell Type",
       y = "Proportion") +
  scale_fill_manual(values=c("#E31A1C","#ffffb3","#30d724","#b8b1e5","#207aae","#d070d0","#D3D3D3")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Calculate the mean of the arsinh-transformed counts per cell phenotype


celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),
 ids = spe$celltype, statistics = "mean",
 use.assay.type = "exprs",
 subset_row = rowData(spe)$use_channel)

dittoHeatmap(celltype_mean,
 assay = "exprs", cluster_cols = TRUE,cluster_rows = F,
 scale = "row", 
 heatmap.colors = viridis(100),
 annot.by = c("celltype", "ncells"))
# Calculate the mean of the arsinh-transformed counts per annotated 


#################landscape S5A
#boxplot
library(ggpubr)
matrix<-assay(spe_fastMNN, "normalized")%>%as.data.frame()%>%t()
df_long <- as.data.frame(matrix)
df_long$sample_id <- rownames(df_long)
df_long <- tidyr::pivot_longer(df_long, -sample_id, names_to = "marker", values_to = "expression")
head(df_long)
#S5A expression
df_long$sample_id<-sub("_[^_]*$","",df_long$sample_id)
mean<-df_long
for (marker1 in unique(df_long$marker)) {
  tmp<-mean%>%filter(marker==marker1)%>%group_by(sample_id)%>%summarize(mean_expression=median(expression))
  tmp<-tmp%>%mutate(ecMYC = ifelse(substr(tmp$sample_id,1,8) %in% c("Patient1", "Patient2", "Patient3"), "yes", "no"))
  p1<-ggboxplot(tmp, x="ecMYC", y="mean_expression",color="ecMYC",
    add = "jitter",palette  =  "nejm",add.params  = list (fill  =  "white"))+  
    stat_compare_means(aes(group = ecMYC), label = "p.format",size=2,method="wilcox.test")+labs(y =marker1)
    theme(axis.title.y=element_text(size=6),axis.text.y=element_text(size=6),axis.text.x=element_text(size=6))+mytheme
  ggsave(paste0("/home/ylab/IMC/4_output/T_cell/mean_2/", marker1, "_2_group.pdf"), p1,  width = 40/25.4, height = 70/25.4, dpi = 300)
}
#S5A Tcell
library(ggpubr)
cell<-colData(spe_fastMNN)%>%as.data.frame()
df_pro<-data.frame(sample_id = character(),
                     cell_type = character(),
                     proportion = numeric())
for (i in unique(cell$sample_id)) {
  for (j in unique(cell$celltype)) {
  sample_data <- subset(cell, sample_id == i)
 total_counts<-table(sample_data$sample_id==i)%>%sum()
  sample_counts <- sum(sample_data$celltype==j)
  sample_proportions <- sample_counts / total_counts
  df_pro<-rbind(df_pro,data.frame(sample_id = i,
                                     cell_type = j,
                                     proportion = sample_proportions))
  }
}
for (marker1 in unique(df_pro$cell_type)) {
  tmp<-df_pro%>%filter(df_pro$cell_type==marker1)%>%group_by(sample_id)%>%summarize(mean_expression=median(proportion))
  tmp<-tmp%>%mutate(ecMYC = ifelse(substr(tmp$sample_id,1,8) %in% c("Patient1", "Patient2", "Patient3"), "yes", "no"))
  p1<-ggboxplot(tmp, x="ecMYC", y="mean_expression",color="ecMYC",
    add = "jitter",palette  =  "nejm",add.params  = list (fill  =  "white"))+  
    stat_compare_means(aes(group = ecMYC), label = "p.format",size=2, method = "wilcox.test")+labs(y =marker1)
    theme(axis.title.y=element_text(size=6),axis.text.y=element_text(size=6),axis.text.x=element_text(size=6))
  ggsave(paste0("./4_output/T_cell/", marker1, "_2_group.pdf"), p1,  width = 40/25.4, height = 70/25.4, dpi = 300)
}



for (marker1 in unique(df3$celltype)) {
  tmp <-   df_counts_all %>% 
    filter(cluster_celltype == marker1) %>% 
    group_by(sample_id) %>% 
    summarize(mean_expression = mean(proportion))
  tmp<-tmp%>%mutate(ecMYC = ifelse(substr(tmp$sample_id,1,8) %in% c("Patient1", "Patient2", "Patient3"), "yes", "no"))
  p1 <- ggboxplot(tmp, x = "ecMYC", y = "mean_expression", color = "ecMYC",
                  add = "jitter", palette = "nejm",
                  add.params = list(fill = "white")) +  
    stat_compare_means(aes(group = ecMYC), label = "p.format", size = 2, method = "wilcox.test") +
    labs(y = marker1) +
    theme(axis.title.y = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6)) +
    mytheme
  ggsave(paste0("./4_output/landscape/", marker1, "_2_group.png"), p1, width = 40/25.4, height = 60/25.4, dpi = 400)
}
pdf(file.path("/mnt/c/Users/ylab/Desktop/New", paste0( "S5B_CN.pdf")),
    width = 120/25.4, height = 100/25.4)
pheatmap(mat,
 color = colorRampPalette(c("dark blue", "white", "dark red"))(100),
 scale = "column")
dev.off()