library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(stringr)
library(ggpubr)

##PKUPH_cohort
exp_mat <- read_csv("../data/PKUPH_cohort/sample40_count_matrix.csv") %>%
  mutate(
    ensg_id = str_extract(gene_id, "[^.]+"),
    gene_name = str_extract(gene_id, "[!|].+") %>% substr(2, stop = nchar(.))
  )
col_data<-read.csv("../data/PKUPH_cohort/col_data.csv",row.names=1)
col_data$group <- factor(col_data$group, levels = c("N", "P"))
dds <- DESeq2::DESeqDataSetFromMatrix(
  exp_mat %>%
    dplyr::select(-ensg_id, -gene_name) %>%
    column_to_rownames("gene_id"),
  colData = col_data %>%
    column_to_rownames("sample"),
  ~ group
) %>% DESeq2::DESeq()
res<-dds %>%
  results() %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")%>%
  as_tibble() %>%
  drop_na()

#GSEA
gene_list_for_GSEA <- res %>%
  mutate(entrezgene_id = mapIds(org.Hs.eg.db, 
                                keys = sub("\\..*", "", gene_id), 
                                column = "ENTREZID", 
                                keytype = "ENSEMBL", 
                                multiVals = "first")) %>%
  filter(!is.na(entrezgene_id)) %>%
  distinct(entrezgene_id, .keep_all = TRUE) %>%  
  filter(!is.na(log2FoldChange)) %>%  
  dplyr::select(entrezgene_id, log2FoldChange) %>%  
  arrange(desc(log2FoldChange)) %>%  
  deframe() 
gsea_results <- gseGO(
  geneList = gene_list_for_GSEA, 
  OrgDb = org.Hs.eg.db, 
  ont = "ALL",
  minGSSize = 10, 
  maxGSSize = 500, 
  pvalueCutoff = 0.05, 
  verbose = TRUE
)
##Figure 4A
pdf("../results/Fig4A.pdf",width=22/2.54,height=14/2.54)
dotplot(gsea_results, showCategory = 10,split=".sign")+facet_wrap(~.sign,scales="free")
dev.off()
##Figure 4B
pdf("../results/Fig4B.pdf",,width=18/2.54,height=14/2.54)
gseaplot2(gsea_results, geneSetID =c("GO:0019814","GO:0003823","GO:0050853"), color = c("#003B5C","#D3C19B","#F5A623"),)
dev.off()

##Figure 4C
load("../data/George_GSEA.Rdata")
ggplot(df, aes(x = NES, y = reorder(Description, NES), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(high= "blue",  low= "red") +
  labs(title = "GSEA (George et al.)", x = "NES", y = "Description")+theme_bw()
ggsave("../results/Fig4C.pdf",height=100,width=160,units="mm",dpi=400)
dev.off()

##Figure 4D
load("../data/MCP_results.Rdata")
#T cell
ggboxplot(data_combined, x="group", y="T cells",color="group",
    add = "jitter",palette  =  "nejm",
    add.params  = list (fill  =  "white"))+ 
    theme(axis.title.y=element_text(size=12),axis.text.y=element_text(size=12),
    legend.text = element_text(size = 8), legend.title = element_text(size = 8),)+
    stat_compare_means(aes(group = group),method = "wilcox.test", label = "p.format",size=2)+
    theme(axis.title.x =element_blank())
ggsave("../results/Fig4D_MCP_T_cell.pdf",height=120,width=60,units="mm",dpi=400)
#NK cells
ggboxplot(data_combined, x="group", y="NK cells",color="group",
    add = "jitter",palette  =  "nejm",
    add.params  = list (fill  =  "white"))+ 
    theme(axis.title.y=element_text(size=12),axis.text.y=element_text(size=12),
    legend.text = element_text(size = 8), legend.title = element_text(size = 8),)+
    stat_compare_means(aes(group = group),method = "wilcox.test", label = "p.format",size=2)+
    theme(axis.title.x =element_blank())
ggsave("../results/Fig4D_MCP_NK cells.pdf",height=120,width=60,units="mm",dpi=400)
#Monocytic lineage
ggboxplot(data_combined, x="group", y="Monocytic lineage",color="group",
    add = "jitter",palette  =  "nejm",
    add.params  = list (fill  =  "white"))+ 
    theme(axis.title.y=element_text(size=12),axis.text.y=element_text(size=12),
    legend.text = element_text(size = 8), legend.title = element_text(size = 8),)+
    stat_compare_means(aes(group = group),method = "wilcox.test", label = "p.format",size=2)+
    theme(axis.title.x =element_blank())
ggsave("../results/Fig4D_MCP_Mono.pdf",height=120,width=60,units="mm",dpi=400)