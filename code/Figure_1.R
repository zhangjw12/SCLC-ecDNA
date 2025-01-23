library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)

load("ecDNA_cell_line_results.Rdata")
mytheme = theme(plot.title = element_text(hjust = 0.5),
                text=element_text(size=10),
                axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
                axis.text.y = element_text(size = 8,color = "black"),
                legend.text = element_text(size=6,color = "black"),
                legend.title = element_text(size=8),
                legend.key.size = unit(11, "pt"))

##figure 1A
ecDNA_distribution<-ggplot(ecDNA,mapping = aes(x=TotalIntervalSize,y=AverageAmplifiedCopyCount,color=subtype))+
    geom_point(size=1.5,alpha=0.75)+
    theme_bw()+
    labs(title =  bquote("ecDNA distribution"),
         x = 'Log10[size(Mb)+1]',
         y = "Mean copy number") +
          scale_discrete_manual(values=c("grey","#E41A1C"),
                          aesthetics = 'colour',
                          labels = c("Others","SCLC"))+mytheme
pdf("../results/Fig1A.pdf",width=4,height=4)
ecDNA_distribution
dev.off()
#violin plot
pdf("../results/Fig1A2.pdf",width=2,height=6)
ggviolin(ecDNA, x="subtype", y="AverageAmplifiedCopyCount",
     add = "boxplot",fill  =  "subtype",palette  =  c("grey", "#E41A1C"),alpha=0.75,
     add.params  = list (fill  =  "white"),)+mytheme+
     stat_compare_means(aes(group = subtype), label = "p.signif",size=6)
dev.off()

##figure 1B
figure_B_bar_data$cancer_type <- factor(figure_B_bar_data$cancer_type,levels = c("SCLC","Others"))
ggplot(figure_B_bar_data,aes(reorder(DepMap_ID,copyCount,decreasing = F),copyCount,fill =cancer_type ))+
  geom_bar(stat = "identity")+theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  annotate("text",size = 8/2.45, x = 30 , y = 27,label = "CNV>20",colour="firebrick3") +
  geom_hline(yintercept=20,linetype = "dashed",linewidth = 0.3,color = "firebrick3")+
  scale_fill_manual(values = c("firebrick3","#aaaaaa"))+
  theme_classic()+mytheme
ggsave("../results/Fig1B.pdf",height=95,width=130,units="mm",dpi=400)
#pie plot
pdf("../results/Fig1B2.pdf.pdf",width= 10/2.54,height=10/2.54)
pie(figure_B_pie_data$frequency,labels = figure_B_pie_data$DNA,
    col = c("#ba6124","#de752d","#f09971","#aaaaaa"),
    border = "white",edges = 200,init.angle = 90)
dev.off()

##figure 1C
ggplot(figure_data,aes(x=ecDNA,y=sample_percent,fill=cancer_type))+
  geom_bar(stat = "identity",width = 0.5)+
  scale_fill_manual(values = c("#bfbfbf","firebrick3"))+
  guides(fill=guide_legend(reverse = T))+
  scale_y_continuous(breaks = NULL)+
  theme_void()+mytheme
ggsave("../results/Fig1C.pdf",height=95,width=130,units="mm",dpi=400)

##figure 1D
data_melt$lineage<-factor(data_melt$lineage,levels = c("peripheral_nervous_system","plasma_cell","soft_tissue","blood",
     "upper_aerodigestive","unirary_tract","esophagous","liver","pancreas",
     "uterus","central_neurvous_system","skin","gastric","colorectal","ovary","breast","lung","SCLC","NSCLC"))
data_melt$oncogene<-factor(data_melt$oncogene,levels = c("MYC","MYCN","MYCL1","ERBB2","CCND1","CCND2","CCNE1","CDK4","CDK12","FGFR1","EGFR","KRAS","NRAS","AKT2"))
ecDNA_pan <- ggplot(data_melt, aes(x = oncogene, y =lineage, size = count,color=mean))+
     scale_colour_gradient(low = "blue", high = "red")+
     theme_bw()+
     theme(axis.text.x=element_text(angle=45, face="italic",hjust=1))+mytheme+
     geom_point()
pdf("../results/Fig1D.pdf",width=7.5,height=5)
ecDNA_pan
dev.off()