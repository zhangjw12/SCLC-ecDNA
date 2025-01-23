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