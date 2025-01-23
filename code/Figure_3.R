library(pheatmap)
library(tidyverse)

sample_info_SCLC<-read.csv("../data/PKUPH_cohort/table_s3.csv",check.names = FALSE,row.names=1)
clinical_data <- sample_info_SCLC[, c("Ki67","TTF-1", "cMYC", "NMYC", "LMYC", "MYC", "MYCN", "MYCL")]

pdf("../results/Fig3A.pdf",width=10/2.54,height=15/2.54)
pheatmap(clinical_data,
         annotation_row = sample_info_SCLC[, c("stage","sex","age")],
         annotation_col = NULL,  
         cluster_rows = F,    # 对行进行聚类
         cluster_cols = FALSE,   # 不对列进行聚类
         color = colorRampPalette(c("#f7fdff","#1a1ad3"))(10),  # 自定义颜色
         main = NA,
         fontsize = 8 ,
         border_color = NA # 字体大小
)
dev.off()
