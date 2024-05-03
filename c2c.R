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
library(CellChat)
library(patchwork)
library(ggpubr)
library(patchwork)
options(bitmapType = "cairo")
library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler)
`%!in%` <- Negate('%in%')


rna.data = readRDS("MDE.rds")
rna.data = subset(rna.data, cell.type != "ME")


unique(rna.data$group)
rna = subset(rna.data, group == "ADD")

cellchat = createCellChat(object = rna@assays$RNA@data, meta = rna@meta.data, group.by = "cell.type")
CellChatDB = CellChatDB.human
CellChatDB.use = CellChatDB
cellchat@DB = CellChatDB.use

cellchat = subsetData(cellchat)
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)
cellchat = computeCommunProb(cellchat, type = "triMean")
cellchat = filterCommunication(cellchat, min.cells = 10)
cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)

saveRDS(cellchat, "../data/cell_chat/ADD.rds", compress = F)


control = readRDS("cell_chat/CTRL.rds")
dlb = readRDS("cell_chat/DLB.rds")
pdd = readRDS("cell_chat/PDD.rds")
add = readRDS("cell_chat/ADD.rds")


plot.data = data.frame(
  group = c("CTRL","CTRL","ADD","ADD","DLB","DLB","PDD","PDD"),
  type = c( "weight","count", "weight","count", "weight","count", "weight", "count"),
  val = c(control@net$weight[1,2] + control@net$weight[2,1], control@net$count[1,2] + control@net$count[2,1],
          add@net$weight[1,2] + add@net$weight[2,1], add@net$count[1,2] + add@net$count[2,1],
          dlb@net$weight[1,2] + dlb@net$weight[2,1], dlb@net$count[1,2] + dlb@net$count[2,1],
          pdd@net$weight[1,2] + pdd@net$weight[2,1], pdd@net$count[1,2] + pdd@net$count[2,1])
)


plot.data.sub = plot.data[which(plot.data$type == "count"),]
p = ggbarplot(plot.data.sub, x = "group", y = "val",
          fill = "group",
          palette = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"),
          x.text.angle = 90   
) + NoLegend() + xlab("") + ylab("Count of interaction")

ggsave("../figure/interaction_count.pdf",p,width = 5, height = 5)


plot.data.sub = plot.data[which(plot.data$type == "weight"),]
p = ggbarplot(plot.data.sub, x = "group", y = "val",
          fill = "group",
          palette = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"),
          x.text.angle = 90   
) + NoLegend() + xlab("") + ylab("Weight of interaction")
ggsave("../figure/interaction_weight.pdf",p,width = 5, height = 5)
