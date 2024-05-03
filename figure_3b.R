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
rna.data = subset(rna.data, cell.type == "Microglia")
exp.count = rowSums(x = GetAssayData(object = rna.data, slot = "counts") > 0) 
exp.count = exp.count/ncol(rna.data)
genes = names(which(exp.count > 0.05))
genes = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")



mm = read.csv("../figure/DLB_vs_control_microglia.txt", sep = "\t")
go_enrich <- enrichGO(gene = mm$gene[which(mm$avg_log2FC > 0)],
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

a = go_enrich@result
a = a[which(a$p.adjust < 0.01),]
a$type = "positive DEG"
pdf("../figure/ME_pathway_pos.pdf", height = 8, width = 9)
treeplot(pairwise_termsim(go_enrich))
dev.off()


go_enrich <- enrichGO(gene = mm$gene[which(mm$avg_log2FC < 0)],
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
b = go_enrich@result
b = b[which(b$p.adjust < 0.01),]
b$type = "negative DEG"
pdf("../figure/ME_pathway_neg.pdf", height = 8, width = 9)
treeplot(pairwise_termsim(go_enrich))
dev.off()

a = rbind(a,b)
write.table(a, "../figure/MIM_GO.txt", sep = "\t", quote = F, col.names = T, row.names = F)

neg = bitr(mm$gene[which(mm$avg_log2FC < 0)], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
pos = bitr(mm$gene[which(mm$avg_log2FC > 0)], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

res = list(pos, neg, genes)
saveRDS(res, "tmp.rds")
