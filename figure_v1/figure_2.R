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
options(bitmapType = "cairo")
`%!in%` <- Negate('%in%')

rna.data = readRDS("rna_annotated.rds")
rna.data$cell.type[which(rna.data$cell.type == "microglia MDE")] = "microglia"

unique(rna.data$group)
rna.sub = subset(rna.data, group %in% c("CTRL", "PDD"))
rna.sub = SetIdent(rna.sub, value = "cell.type")
a = data.frame(table(rna.sub$cell.type, rna.sub$group))
a = a[which(a$Freq > 30),]
a = a[which(a$Var1 %in% names(which(table(a$Var1) > 1))),]

a = as.character(unique(a$Var1))


mm = pbmclapply(
  a, function(ct) {
    m = FindMarkers(rna.sub, ident.1 = "PDD", ident.2 = "CTRL", group.by = "group", logfc.threshold = 0.5, subset.ident = ct)
    m$gene = rownames(m)
    m$disease = "PDD"
    m$cell_type = ct
    return(m)
  }, mc.cores = 20
)
mm = do.call(rbind, mm)

write.table(mm, "../figure/DEG_PDD.txt", sep = "\t", quote = F, row.names = F, col.names = T)



rna.sub = subset(rna.data, group %in% c("CTRL", "ADD"))
rna.sub = SetIdent(rna.sub, value = "cell.type")
a = data.frame(table(rna.sub$cell.type, rna.sub$group))
a = a[which(a$Freq > 30),]
a = a[which(a$Var1 %in% names(which(table(a$Var1) > 1))),]

a = as.character(unique(a$Var1))


mm = pbmclapply(
  a, function(ct) {
    m = FindMarkers(rna.sub, ident.1 = "ADD", ident.2 = "CTRL", group.by = "group", logfc.threshold = 0.5, subset.ident = ct)
    m$gene = rownames(m)
    m$disease = "ADD"
    m$cell_type = ct
    return(m)
  }, mc.cores = 20
)
mm = do.call(rbind, mm)

write.table(mm, "../figure/DEG_ADD.txt", sep = "\t", quote = F, row.names = F, col.names = T)




rna.sub = subset(rna.data, group %in% c("CTRL", "DLB"))
rna.sub = SetIdent(rna.sub, value = "cell.type")
a = data.frame(table(rna.sub$cell.type, rna.sub$group))
a = a[which(a$Freq > 30),]
a = a[which(a$Var1 %in% names(which(table(a$Var1) > 1))),]

a = as.character(unique(a$Var1))


mm = pbmclapply(
  a, function(ct) {
    m = FindMarkers(rna.sub, ident.1 = "DLB", ident.2 = "CTRL", group.by = "group", logfc.threshold = 0.5, subset.ident = ct)
    m$gene = rownames(m)
    m$disease = "DLB"
    m$cell_type = ct
    return(m)
  }, mc.cores = 20
)
mm = do.call(rbind, mm)

write.table(mm, "../figure/DEG_DLB.txt", sep = "\t", quote = F, row.names = F, col.names = T)



mm1 = read.csv("../figure/DEG_ADD.txt", sep = "\t")
mm2 = read.csv("../figure/DEG_DLB.txt", sep = "\t")
mm3 = read.csv("../figure/DEG_PDD.txt", sep = "\t")

mm = rbind(mm1, mm2, mm3)
mm$type = "positive"
mm$type[which(mm$avg_log2FC < 0)] = "negative"

mm$group = paste0(mm$disease, "-",mm$cell_type)
mm = mm[which(mm$p_val_adj < 0.01),]
plot.data = data.frame(table(mm$group, mm$type))
a = data.frame(table(mm$group))
a = a[order(a$Freq),]

plot.data$Var1 = factor(as.character(plot.data$Var1), levels = a$Var1)
p = ggplot(plot.data, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")

ggsave("../figure/DEG_pos_neg.pdf", p, width = 12, height = 7)

mm.sub = mm[which(mm$gene == "APOE"),]
rna.sub = subset(rna.data, cell.type == "microglia")

p = VlnPlot(rna.sub, "APOE", group.by = "group") + xlab("") + NoLegend()
ggsave("../figure/apoe.pdf", p, width = 7, height = 6)



mm.sub = mm[which(mm$gene == "SNCA"),]
rna.sub = subset(rna.data, cell.type == "L5b EN")
p = VlnPlot(rna.sub, "SNCA", group.by = "group") + xlab("") + NoLegend()
ggsave("../figure/snca.pdf", p, width = 7, height = 6)


mm.sub = mm[which(mm$gene == "BIN1"),]
rna.sub = subset(rna.data, cell.type == "astrocyte")
p = VlnPlot(rna.sub, "BIN1", group.by = "group") + xlab("") + NoLegend()
ggsave("../figure/bin1.pdf", p, width = 7, height = 6)



mm %>%
  group_by(cell_type, disease) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10

length(unique(top10$gene))

plot.data = dcast(group ~ gene, data = top10, value.var = "avg_log2FC", fill = 0)
rownames(plot.data) = plot.data[,1]
plot.data = plot.data[,-1]
pheatmap::pheatmap(plot.data, scale = "column")


mm.sub = mm[which(mm$gene %in% names(sort(-table(mm$gene))[1:50])),]
plot.data = dcast(group ~ gene, data = mm.sub, value.var = "avg_log2FC", fill = 0)
rownames(plot.data) = plot.data[,1]
plot.data = plot.data[,-1]

pheatmap::pheatmap(plot.data, scale = "none", breaks = seq(-4, 4, by = 0.05),filename = "../figure/common_deg.pdf",
                   color = colorRampPalette(c("blue", "white", "red"))(161), width = 10, height = 8)



library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler)


a = unique(mm$group)
pathways = pbmclapply(
  a, function(g){
    mm.sub = mm[which(mm$group == g),]
    go_enrich <- enrichGO(gene = mm.sub$gene,
                          OrgDb = org.Hs.eg.db, 
                          keyType = 'SYMBOL',
                          readable = T,
                          ont = "BP",
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)
    go_enrich = go_enrich@result
    go_enrich$group = g
    go_enrich = go_enrich[which(go_enrich$p.adjust <0.05),]
    
    return(go_enrich)
  }, mc.cores = 45
)
pathways = do.call(rbind, pathways)

pathways.sub = pathways[which(pathways$ID %in% names(sort(-table(pathways$ID))[1:20])),]
pathways.sub$val = -log(pathways.sub$pvalue)

plot.data = dcast(Description ~ group, data = pathways.sub, value.var = "val", fill = 0)
rownames(plot.data) = plot.data[,1]
plot.data = plot.data[,-1]

pheatmap::pheatmap(plot.data, scale = "none",filename = "../figure/common_pathway.pdf",
                   width = 10, height = 7,
                   color = colorRampPalette(c("blue", "white", "red"))(161))
####################################################################

go_enrich <- enrichGO(gene = top10$gene,
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

pdf("../figure/common_gene_pathway.pdf", height = 8, width = 12)
treeplot(pairwise_termsim(go_enrich))
dev.off()

treeplot(pairwise_termsim(go_enrich))
    

#####################################################################
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

mm1 = read.csv("../figure/DEG_ADD.txt", sep = "\t")
mm2 = read.csv("../figure/DEG_DLB.txt", sep = "\t")
mm3 = read.csv("../figure/DEG_PDD.txt", sep = "\t")

mm = rbind(mm1, mm2, mm3)
mm$type = "positive"
mm$type[which(mm$avg_log2FC < 0)] = "negative"

mm = mm[which(mm$p_val_adj < 0.01),]


ct = unique(mm$cell_type)

pbmclapply(
  ct, function(cc) {
    mm.sub = mm[which(mm$cell_type == cc),]
    
    cc = sub("\\/", "_", cc)
    venn.diagram(
      x = list(mm.sub$gene[which(mm.sub$disease == "PDD")], 
               mm.sub$gene[which(mm.sub$disease == "ADD")], 
               mm.sub$gene[which(mm.sub$disease == "DLB")]),
      category.names = c("PDD", "ADD" , "DLB"), 
      fill = myCol,
      filename = paste0(cc, ".png"),imagetype = "png"
    )
  }
)


########################gene disease specificity##############################
gene = unique(mm$gene[which(mm$disease == "ADD")])

res = pbmclapply(
  gene, function(g) {
    g = "PLP1"
    a = max(sum(mm$gene == g & mm$disease == "ADD"),1)
    b = max(sum(mm$gene == g) - a, 1)
    c = sum(mm$disease == "ADD") - a
    d = 10113 - a -b - c
    fisher.test(matrix(c(a,b,c,d), nrow = 2))$estimate
  }
)
res = unlist(res)
res = data.frame(specifity = res, row.names = gene)
res.add = res
res.add$disease = "ADD"
res.add$gene = rownames(res.add)

gene = unique(mm$gene[which(mm$disease == "PDD")])

res = pbmclapply(
  gene, function(g) {
    a = max(sum(mm$gene == g & mm$disease == "PDD"),1)
    b = max(sum(mm$gene == g) - a, 1)
    c = sum(mm$disease == "PDD") - a
    d = 10113 - a -b - c
    fisher.test(matrix(c(a,b,c,d), nrow = 2))$estimate
  }
)
res = unlist(res)

res = data.frame(specifity = res, row.names = gene)
res.pdd = res
res.pdd$disease = "PDD"
res.pdd$gene = rownames(res.pdd)

gene = unique(mm$gene[which(mm$disease == "DLB")])

res = pbmclapply(
  gene, function(g) {
    a = max(sum(mm$gene == g & mm$disease == "DLB"),1)
    b = max(sum(mm$gene == g) - a, 1)
    c = sum(mm$disease == "DLB") - a
    d = 10113 - a -b - c
    fisher.test(matrix(c(a,b,c,d), nrow = 2))$estimate
  }
)
res = unlist(res)
res = data.frame(specifity = res, row.names = gene)
res.dlb = res
res.dlb$disease = "DLB"
res.dlb$gene = rownames(res.dlb)


res = rbind(res.add, res.dlb, res.pdd)
write.table(res, "../figure/gene_disease_specifity.txt", sep = "\t", row.names = F)

rna.data = readRDS("rna_annotated.rds")
VlnPlot(rna.data, "GMDS-DT", group.by = "group", pt.size = 0) + xlab("") + NoLegend()


res %>%
  group_by(disease) %>%
  slice_max(specifity, n = 10) %>%
  ungroup() -> top10

rna.data$type = paste0(rna.data$group, " - ", rna.data$cell.type)


rna.average = AverageExpression(rna.data, group.by = "type")
rna.average = rna.average$raw
rna.average = rna.average[top10$gene,]

rna.average = sweep(rna.average, 1, FUN="/", rowMeans(rna.average))



p = ggbarplot(top10, x = "gene", y = "specifity",
           fill = "disease",
           color = "white",  
           palette = "jco",
           sort.val = "asc",  
           sort.by.groups = TRUE,
           x.text.angle = 90,
           ylab = "Specificity",
           xlab = F,
)
p
ggsave("../figure/gene_disease_specifity.pdf", p , width = 9, height = 6)

p = VlnPlot(rna.data, "PLP1", group.by = "group", pt.size = 0) + xlab("") 
ggsave("../figure/PLP1.pdf", p , width = 6, height = 5)


