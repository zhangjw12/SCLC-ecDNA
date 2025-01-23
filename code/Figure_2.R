library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggrepel)
library(pheatmap)
library(ggExtra)
library(tidyverse)

##Figure 2A
load("../data/CCLE.Rdata")
annotation_colors<-c(list(EC_type=c(ecMYC="#66C2A5",ecMYCN="#FC8D62",ecMYCL1="#8DA0CB",ecNA="grey")),
    list(ecMYC=c(yes="#1B9E77",no='grey')))
pdf("../results/Fig2A.pdf",width=30/2.54,height=15/2.54)
pheatmap(heatmap,cluster_cols=T,cluster_rows=T,
    scale="row",
    color = colorRampPalette(c("darkblue", "white", "firebrick3"))(300),
    clustering_distance_cols = 'euclidean',
    annotation_col=annotation_col,
    annotation_colors=annotation_colors,
    breaks=c(seq(-1.5,1.5,by=0.01)),
    cutree_cols=2,
    legend=T,
    show_colnames=F,
    fontsize_row=8,
    label_row=T
)
dev.off()

#Figure 2B
EMT_p <- ggplot(EMT,mapping = aes(x=GOCC_MHC_CLASS_I_PEPTIDE_LOADING_COMPLEX,y=GOBP_CELL_CYCLE_DNA_REPLICATION_INITIATION,color=ecMYC))+
    geom_point(size=1.5,alpha=0.75)+
    geom_smooth (method = "lm",se=FALSE,colour="#5f5f5f",alpha=0.5)+
    theme_bw()+
    guides(fill = guide_legend(label=F))+
    theme(axis.title.x = element_text(size = 4,color = "black"))+
     theme(axis.title.y = element_text(size = 4,color = "black"))+
     theme(axis.text.x = element_text(size = 5,color = "black"))+
     theme(axis.text.y = element_text(size = 5,color = "black"))+
    labs(title = NULL,
         x = "GOCC_MHC_CLASS_I_PEPTIDE_LOADING_COMPLEX",  
         y = "GOBP_CELL_CYCLE_DNA_REPLICATION_INITIATION")+
         theme(legend.position="none")+
    scale_discrete_manual(values=c("grey","#E41A1C"),
                          aesthetics = 'colour',
                          labels = c("ecMYC-","ecMYC+")) 
cor.test(EMT$GOCC_MHC_CLASS_I_PEPTIDE_LOADING_COMPLEX,EMT$GOBP_CELL_CYCLE_DNA_REPLICATION_INITIATION)
pdf("../results/Fig2B.pdf",width=7/2.54,height=7/2.54)
ggMarginal(EMT_p, type = "density", groupColour = TRUE, groupFill = TRUE)
dev.off()

#Figure 2E
#PCA
load("../data/cell_line_sequencing.Rdata")
pca_data <- vst(dds_524, blind = FALSE) %>%
  plotPCA(intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
PCA <- vst(dds_524, blind = FALSE) %>%
  plotPCA(intgroup = "group", returnData = TRUE) %>%
  ggplot(aes(PC1, PC2, color = group)) +
  geom_point(size = 1) +  
  theme_bw() +  
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c("C" = "darkblue", 
                               "G" = "lightgreen", 
                               "H" = "firebrick3")) + 
  labs(title = NULL, 
       x = paste("PC1 (", percent_var[1], "%)", sep = ""), 
       y = paste("PC2 (", percent_var[2], "%)", sep = "")) 
pdf("../results/Fig3E.pdf",width=7.2/2.54,height=6/2.54)
PCA
dev.off()

##Figure 2F
res_524_hu <-
  dds_524 %>%
  results(name = "group_H_vs_C") %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble()
cut_off_pvalue = 1.301
cut_off_logFC = 1
a <- res_524_hu %>%
  mutate(
    A = log2(baseMean + 1),  
    M = log2FoldChange,
    p = -log10(padj),  
    change = case_when(
      p > cut_off_pvalue & abs(log2FoldChange) >= cut_off_logFC ~ 
        case_when(log2FoldChange > cut_off_logFC ~ 'Up', log2FoldChange < -cut_off_logFC ~ 'Down', TRUE ~ 'No difference'),
      TRUE ~ 'No difference'
    )
  ) %>%drop_na()%>%
  left_join(exp_mat_524 %>% dplyr::select(gene_id, gene_name), by = "gene_id") 
#MA_plot
p_ma <- a %>%
  mutate(ensg_id = str_extract(gene_id, "[^|]+")) %>%
  filter(ensg_id %in% protein_coding$V1) %>%
  ggplot(aes(x = A, y = M, colour = change),position = "jitter", ) +
  geom_point(data = a[a$change == "No difference", ],alpha = 1, size = 0.5) +
  geom_point(data = a[a$change == "Down", ],alpha = 1, size = 0.5 ) +
  geom_point(data = a[a$change == "Up", ],alpha = 1, size = 0.5 ) +
  scale_color_manual(values = c("darkblue", "grey", "firebrick3")) +
  geom_hline(yintercept = 0, lty = 4, col = "black", lwd = 0.8) +  # 0 line for fold change
  labs(x = "log2(baseMean+1)", y = "log2(FoldChange)", title = "MA Plot") +
  theme(axis.title = element_text(size = 20, face = "bold")) +
  theme_bw() +
  coord_cartesian(ylim = c(-10, 10))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank()) +
  geom_text_repel(data = a %>% filter(gene_name %in% c("HLA-A", "HLA-B", "HLA-C")),
                  aes(label = gene_name), 
                  fontface = "bold", 
                  color = "black",
                  box.padding = 0.5,   
                  point.padding = 0,
                  segment.size = 0.5, 
                  nudge_x = 3,
                  nudge_y = 3,
                  max.overlaps = 5,
                  force = 10, 
                  min.segment.length = 1   )
pdf("../results/Fig2F.pdf",width=7/2.54,height=7/2.54)
p_ma
dev.off()