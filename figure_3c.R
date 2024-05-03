setwd("/scbio7/mengwei/NNI_brain/data/")
library(pbmcapply)
library(Seurat)
library(scDblFinder)
library(ggplot2)
library(stringr)
library(dittoSeq)
library(reshape2)
library(ggpubr)
library(dplyr)
library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler)
library("gridExtra")
library(ggpubr)
library(patchwork)
options(bitmapType = "cairo")
library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler)
`%!in%` <- Negate('%in%')


rna.data = readRDS("MDE.rds")
DimPlot(rna.data, group.by = "cell.type")

rna.data = subset(rna.data, cell.type %in% c("Microglia"))

activatedMicroglia = c("CT010467.1","EEF1A1","EEF1B2","EEF2","FAU","GM10076","GM10269","GM10443","GM10736","GM14303","GM2000","GM23935","GM3511","GM5735","GM7618","LARS2","MIR6236","MT-CO1","MT-CYTB","MT-RNR2","RPL10","RPL13","RPL14","RPL17","RPL18A","RPL21","RPL27A","RPL28","RPL31","RPL32","RPL35","RPL36","RPL37","RPL37A","RPL38","RPL41","RPL5","RPL9","RPLP0","RPLP1","RPLP2","RPS10-PS2","RPS11","RPS14","RPS16","RPS17","RPS18","RPS19","RPS23","RPS25","RPS27","RPS27A","RPS28","RPS29","RPS3","RPS6","RPS9","RPSA","TMSB4X")

rna.data = AddModuleScore(
  object = rna.data,
  ctrl = 30,
  features = list(activatedMicroglia),
  name = "activatedMicroglia",
  seed = 1
)

my_comparisons <- list( c("CTRL", "ADD"), c("CTRL", "DLB"), c("CTRL", "PDD") )
plot.data = data.frame(score = rna.data$activatedMicroglia1, group = rna.data$group)
p = ggboxplot(plot.data, x = "group", y = "score", palette = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"),
          color = "group")+ 
  stat_compare_means(comparisons = my_comparisons) + NoLegend()
ggsave("../figure/activatedMicroglia.pdf", width = 4, height = 4)



agedMicroglia = c("ALDOA","ARF3","ARL2","ARL8A","ATP5B","ATP5D","ATP5E","ATP5G1","ATP5G2","ATP5G3","ATP5H","ATP5J2","ATP5O","CALM1","CDC37","CHCHD2","COX4I1","COX5B","COX6B1","COX7A2","COX7B","CYBA","EIF1","ENO1","EWSR1","FCER1G","FKBP8","GABARAP","GAPDH","GAS5","GLRX3","GM10123","GPI1","GPX1","GPX4","HCFC1R1","HINT1","KARS","MBP","MIF","MPC1","MRFAP1","MYL12B","MYL6","NDUFB7","NDUFS5","NDUFS7","NPC1","NUTF2","OAZ2","PARK7","PDHB","PEBP1","PKM","PPIA","PRDX1","PRDX5","PSMB2","PSMB4","PTMA","RAN","RHOBTB1","RNH1","RPS27L","SDC4","SH3BGRL3","SLC25A3","SLC25A4","SNX3","SOD1","SRP14","SRPK1","STUB1","TECPR1","TMEM256","TOMM7","TPI1","TUBB4B","TWF1","TYROBP","UBB","UBL5","UQCR10","UQCR11","UQCRQ","YBX1")

rna.data = AddModuleScore(
  object = rna.data,
  ctrl = 30,
  features = list(agedMicroglia),
  name = "agedMicroglia",
  seed = 1
)
my_comparisons <- list( c("CTRL", "ADD"), c("CTRL", "DLB"), c("CTRL", "PDD") )
plot.data = data.frame(score = rna.data$agedMicroglia1, group = rna.data$group)
p = ggboxplot(plot.data, x = "group", y = "score", palette = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"),
              color = "group")+ 
  stat_compare_means(comparisons = my_comparisons) + NoLegend()

ggsave("../figure/agedMicroglia.pdf", width = 4, height = 4)



DAM = c("ANXA2","ANXA5","APOE","ATP6V0C","B2M","C1QB","CAPG","CD52","CD63","CD63-PS","CD74","CLEC7A","CRIP1","CST7","CTSB","CTSS","CTSZ","CYBB","FAM20C","FTH1","FTL1","FTL1-PS1","FTL1-PS2","GM12164","GM7541","H2-D1","H2-K1","IFITM3","LGALS1","LGALS3","LYZ2","MIR692-1","SPP1","TSPO","VIM")

rna.data = AddModuleScore(
  object = rna.data,
  ctrl = 30,
  features = list(DAM),
  name = "DAM",
  seed = 1
)
my_comparisons <- list( c("CTRL", "ADD"), c("CTRL", "DLB"), c("CTRL", "PDD") )
plot.data = data.frame(score = rna.data$DAM1, group = rna.data$group)
p = ggboxplot(plot.data, x = "group", y = "score", palette = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"),
              color = "group")+ 
  stat_compare_means(comparisons = my_comparisons) + NoLegend()
p
ggsave("../figure/DAM.pdf", width = 4, height = 4)



Homeostatic = c("1700017B05RIK","CCR5","CD164","CMTM6","CSF1R","CX3CR1","ECSCR","FCRLS","GPR34","HEXB","IFNGR1","LAPTM5","LDHB","LGMN","LPCAT2","LRRC3","MGAT4A","OLFML3","P2RY12","P2RY13","PTGS1","SELPLG","SIGLECH","SIRPA","SLC2A5","SRGAP2","TMEM119","UBQLN1","UPK1B","VDAC1","VSIR")

rna.data = AddModuleScore(
  object = rna.data,
  ctrl = 30,
  features = list(Homeostatic),
  name = "Homeostatic",
  seed = 1
)
my_comparisons <- list( c("CTRL", "ADD"), c("CTRL", "DLB"), c("CTRL", "PDD") )
plot.data = data.frame(score = rna.data$Homeostatic1, group = rna.data$group)
p = ggboxplot(plot.data, x = "group", y = "score", palette = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"),
              color = "group")+ 
  stat_compare_means(comparisons = my_comparisons) + NoLegend()
p
ggsave("../figure/Homeostatic.pdf", width = 4, height = 4)




MyelinIngMicroglia = c("TCF12","SHTN1","DOCK10","ZEB2","PIP4K2A","MBP","TMTC2","SLC44A1","CTNNA3","CLDN11","PLP1","ST18","RNF220","TMEM144","PDE4B","MAN2A1","ERBB2IP","FRYL","QKI","NCKAP5","PLD1")

rna.data = AddModuleScore(
  object = rna.data,
  ctrl = 30,
  features = list(MyelinIngMicroglia),
  name = "MyelinIngMicroglia",
  seed = 1
)
my_comparisons <- list( c("CTRL", "ADD"), c("CTRL", "DLB"), c("CTRL", "PDD") )
plot.data = data.frame(score = rna.data$MyelinIngMicroglia1, group = rna.data$group)
p = ggboxplot(plot.data, x = "group", y = "score", palette = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"),
              color = "group")+ 
  stat_compare_means(comparisons = my_comparisons) + NoLegend()
p
ggsave("../figure/MyelinIngMicroglia.pdf", width = 4, height = 4)

#######################################################################

rna.data = readRDS("MDE.rds")
DimPlot(rna.data, group.by = "cell.type")

p = VlnPlot(rna.data, c("PLP1", "MBP", "ST18", "PTPRC", "APOE", "CD74"), group.by = "cell.type", pt.size = 0, cols = c("#3f51b5", "#e32365", "#ff7d00")) + NoLegend()
ggsave("../figure/MIM_gene.pdf", width = 8, height = 6)

