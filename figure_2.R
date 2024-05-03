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
rna.data$sex = NA
rna.data$sex[which(rna.data$sample =="CTRL_4")] = "F"
rna.data$sex[which(rna.data$sample =="CTRL_5")] = "M"
rna.data$sex[which(rna.data$sample =="CTRL_6")] = "M"
rna.data$sex[which(rna.data$sample =="CTRL_7")] = "F"
rna.data$sex[which(rna.data$sample =="ADD_4")] = "F"
rna.data$sex[which(rna.data$sample =="ADD_5")] = "M"
rna.data$sex[which(rna.data$sample =="ADD_6")] = "F"
rna.data$sex[which(rna.data$sample =="ADD_7")] = "F"
rna.data$sex[which(rna.data$sample =="DLB_2")] = "M"
rna.data$sex[which(rna.data$sample =="DLB_4")] = "F"
rna.data$sex[which(rna.data$sample =="DLB_5")] = "M"
rna.data$sex[which(rna.data$sample =="DLB_6")] = "F"
rna.data$sex[which(rna.data$sample =="DLB_7")] = "M"
rna.data$sex[which(rna.data$sample =="DLB_8")] = "F"
rna.data$sex[which(rna.data$sample =="PDD_1")] = "F"
rna.data$sex[which(rna.data$sample =="PDD_2")] = "M"
rna.data$sex[which(rna.data$sample =="PDD_4")] = "F"
rna.data$sex[which(rna.data$sample =="PDD_5")] = "M"
rna.data$sex[which(rna.data$sample =="PDD_7")] = "F"
rna.data$sex[which(rna.data$sample =="PDD_9")] = "M"

unique(rna.data$group)
rna.sub = subset(rna.data, group %in% c("CTRL", "PDD"))
rna.sub = SetIdent(rna.sub, value = "cell.type")
a = data.frame(table(rna.sub$cell.type, rna.sub$group))
a = a[which(a$Freq > 30),]
a = a[which(a$Var1 %in% names(which(table(a$Var1) > 1))),]

a = as.character(unique(a$Var1))


mm = pbmclapply(
  a, function(ct) {
    m = FindMarkers(rna.sub, ident.1 = "PDD", ident.2 = "CTRL", group.by = "group", logfc.threshold = 0.5, subset.ident = ct, test.use = "MAST", latent.vars = "sex")
    m$gene = rownames(m)
    m$disease = "PDD"
    m$cell_type = ct
    return(m)
  }, mc.cores = 20
)
mm = do.call(rbind, mm)
mm = mm[which(mm$p_val_adj < 0.01),]
write.table(mm, "../figure/DEG_PDD.txt", sep = "\t", quote = F, row.names = F, col.names = T)



rna.sub = subset(rna.data, group %in% c("CTRL", "ADD"))
rna.sub = SetIdent(rna.sub, value = "cell.type")
a = data.frame(table(rna.sub$cell.type, rna.sub$group))
a = a[which(a$Freq > 30),]
a = a[which(a$Var1 %in% names(which(table(a$Var1) > 1))),]

a = as.character(unique(a$Var1))


mm = pbmclapply(
  a, function(ct) {
    m = FindMarkers(rna.sub, ident.1 = "ADD", ident.2 = "CTRL", group.by = "group", logfc.threshold = 0.5, subset.ident = ct, test.use = "MAST", latent.vars = "sex")
    m$gene = rownames(m)
    m$disease = "ADD"
    m$cell_type = ct
    return(m)
  }, mc.cores = 20
)
mm = do.call(rbind, mm)
mm = mm[which(mm$p_val_adj < 0.01),]

write.table(mm, "../figure/DEG_ADD.txt", sep = "\t", quote = F, row.names = F, col.names = T)

dim(rna.data)


rna.sub = subset(rna.data, group %in% c("CTRL", "DLB"))
rna.sub = SetIdent(rna.sub, value = "cell.type")
a = data.frame(table(rna.sub$cell.type, rna.sub$group))
a = a[which(a$Freq > 30),]
a = a[which(a$Var1 %in% names(which(table(a$Var1) > 1))),]

a = as.character(unique(a$Var1))


mm = pbmclapply(
  a, function(ct) {
    m = FindMarkers(rna.sub, ident.1 = "DLB", ident.2 = "CTRL", group.by = "group", logfc.threshold = 0.5, subset.ident = ct, test.use = "MAST", latent.vars = "sex")
    m$gene = rownames(m)
    m$disease = "DLB"
    m$cell_type = ct
    return(m)
  }, mc.cores = 20
)
mm = do.call(rbind, mm)
mm = mm[which(mm$p_val_adj < 0.01),]
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

plot.data$disease = str_match(plot.data$Var1, "^(.*)-")[,2]
plot.data$cell.type = str_match(plot.data$Var1, "-(.*)$")[,2]

plot.data.pos = plot.data[which(plot.data$Var2 == "positive"),]
plot.data.pos = dcast(plot.data.pos, cell.type ~ disease, value.var = "Freq", fill = 0)
rownames(plot.data.pos) = plot.data.pos[,1]
plot.data.pos = plot.data.pos[,-1]

plot.data.neg = plot.data[which(plot.data$Var2 == "negative"),]
plot.data.neg = dcast(plot.data.neg, cell.type ~ disease, value.var = "Freq", fill = 0)
rownames(plot.data.neg) = plot.data.neg[,1]
plot.data.neg = plot.data.neg[,-1]


library(ComplexHeatmap)
library(circlize)

ha1 = HeatmapAnnotation(bar1 = anno_barplot(as.numeric(colSums(plot.data.neg)), ylim = c(0,1500)))
ha2 = HeatmapAnnotation(bar2 = anno_barplot(as.numeric(colSums(plot.data.pos)), ylim = c(0,1500)))
rowAnnotation = rowAnnotation(bar3 = anno_barplot(as.numeric(rowSums(plot.data.neg+plot.data.pos)), ylim = c(0,2000)))



col_fun = colorRamp2(c(0,600), c("#eeeeee", "#f44336"))

ht1 = Heatmap(as.matrix(plot.data.neg), name = "DEG(pos)", cluster_rows = F, cluster_columns = F, col = col_fun, top_annotation = ha1,
              left_annotation = rowAnnotation)
ht2 = Heatmap(as.matrix(plot.data.pos), name = "DEG(neg)", cluster_rows = F, cluster_columns = F, col = col_fun, top_annotation = ha2)
ht_list = ht1 + ht2 
pdf("../figure/DEG_pos_neg.pdf", width = 8, height = 7)
draw(ht_list)
dev.off()


############################
mm.sub = mm[which(mm$gene %in% names(sort(-table(mm$gene))[1:50])),]
plot.data = dcast(group ~ gene, data = mm.sub, value.var = "avg_log2FC", fill = 0)
rownames(plot.data) = plot.data[,1]
plot.data = plot.data[,-1]

ht = Heatmap(t(plot.data))
draw(ht)
genes = colnames(plot.data)[row_order(ht)]


pos = as.numeric(apply(t(plot.data)[genes,], 1, function(i){
  length(which(i > 0))
}))
neg = as.numeric(apply(t(plot.data)[genes,], 1, function(i){
  length(which(i < 0))
}))
rowAnno = rowAnnotation(bar3 = anno_barplot(cbind(pos, neg), 
                                            gp = gpar(fill = c("#f44336","#3f51b5"), col = c("white","white"))))

col_fun = colorRamp2(c(-2, 0, 2), c("#3f51b5", "#eeeeee", "#f44336"))
plot.data.ADD = plot.data[grep("ADD", rownames(plot.data)),]
plot.data.ADD = t(plot.data.ADD)
colnames(plot.data.ADD) = str_match(colnames(plot.data.ADD), "-(.*$)")[,2]
ht1 = Heatmap(plot.data.ADD[genes,], name = "ADD", col = col_fun, cluster_rows = F, show_row_names=F, column_names_gp = gpar(fontsize = 8))

plot.data.DLB = plot.data[grep("DLB", rownames(plot.data)),]
plot.data.DLB = t(plot.data.DLB)
colnames(plot.data.DLB) = str_match(colnames(plot.data.DLB), "-(.*$)")[,2]
ht2 = Heatmap(plot.data.DLB[genes,], name = "DLB", col = col_fun, cluster_rows = F, column_names_gp = gpar(fontsize = 8))


plot.data.PDD = plot.data[grep("PDD", rownames(plot.data)),]
plot.data.PDD = t(plot.data.PDD)
colnames(plot.data.PDD) = str_match(colnames(plot.data.PDD), "-(.*$)")[,2]
ht3 = Heatmap(plot.data.PDD[genes,], name = "PDD", col = col_fun, cluster_rows = F, right_annotation = rowAnno,
              row_names_side = "right", row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))

pdf("../figure/common_deg.pdf", width = 7, height = 7)
draw(ht1 + ht2 + ht3)
dev.off()
#########################################


####################################################################
mm %>%
  group_by(cell_type, disease) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

go_enrich <- enrichGO(gene = top10$gene,
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
go_enrich@result = pathways
tree = treeplot(pairwise_termsim(go_enrich))
pdf("../figure/common_gene_pathway.pdf", width = 10, height = 6)
tree
dev.off()

tree.data = data.frame(item = tree[["data"]][["label"]], cluster = tree[["data"]][["group"]])
tree.data = tree.data[1:30,]
tree.data = tree.data[order(tree.data$cluster),]
rownames(tree.data) = tree.data$item
tree.data = tree.data[rev(c("regulation of neuron apoptotic process","regulation of postsynapse organization","postsynapse organization","positive regulation of neuron apoptotic process","regulation of synapse structure or activity","regulation of synapse organization","negative regulation of neuron projection development","negative regulation of cell projection organization","ligand-gated ion channel signaling pathway","ionotropic glutamate receptor signaling pathway","glutamate receptor signaling pathway","regulation of ion transmembrane transport","cerebellum morphogenesis","cerebellar cortex morphogenesis","hindbrain morphogenesis","cerebellar cortex development","regulation of cell junction assembly","regulation of synapse assembly","cell junction assembly","synapse assembly","regulation of nervous system development","regulation of neuron projection development","synapse organization","homophilic cell adhesion via plasma membrane adhesion molecules","cell-cell adhesion via plasma-membrane adhesion molecules","retinal ganglion cell axon guidance","negative regulation of synapse organization","negative regulation of cell junction assembly","negative chemotaxis","axon development")),]

##########################################

library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler)

mm$group = paste0(mm$group,"-",mm$type)
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
  }, mc.cores = 1
)
pathways = do.call(rbind, pathways)
pathways$val = -log(pathways$pvalue)



pathways.sub = pathways[which(pathways$Description %in% tree.data$item),]
pathways.sub$val = -log(pathways.sub$pvalue)

pathways.sub$Description = factor(pathways.sub$Description, levels = tree.data$item)
pathways.sub$group = paste0(str_match(pathways.sub$group, "-.*-(.*$)")[,2], "-", str_match(pathways.sub$group, "(^.*-.*)-.*$")[,2])

for (i in 1:nrow(pathways.sub)) {
  pathways.sub$GeneRatio[i] = eval(parse(text = pathways.sub$GeneRatio[i]))
}

pathways.sub$GeneRatio = as.numeric(pathways.sub$GeneRatio)

library(viridis)

p = ggplot(pathways.sub, aes(x = group, y = Description)) + 
  geom_point(aes(size = GeneRatio, fill = val), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 0.5), range = c(1,10), breaks = c(0.1, 0.2, 0.5)) + 
  scale_fill_gradient(low="blue", high="red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
ggsave("../figure/common_pathway.pdf", p, width = 11, height = 8)


############################################################
library(VennDiagram)

mm1 = read.csv("../figure/DEG_ADD.txt", sep = "\t")
mm2 = read.csv("../figure/DEG_DLB.txt", sep = "\t")
mm3 = read.csv("../figure/DEG_PDD.txt", sep = "\t")

mm = rbind(mm1, mm2, mm3)
mm$type = "positive"
mm$type[which(mm$avg_log2FC < 0)] = "negative"

mm$group = paste0(mm$gene, "-",mm$cell_type, "-",mm$type)


set1 = mm$group[which(mm$disease == "ADD")]
set2 = mm$group[which(mm$disease == "DLB")]
set3 = mm$group[which(mm$disease == "PDD")]

venn.diagram(
  x = list(set1, set2, set3), imagetype = "png", fill = c("#8CDCDA", "#B1D877", "#F16A70"),
  category.names = c("ADD" , "DLB " , "PDD"), width = 1000, height = 1000, resolution = 300,
  filename = '../figure/deg_all_venn.png',
  output = T
)


ct = unique(mm$cell_type)

pbmclapply(
  ct, function(cc) {
    mm.sub = mm[which(mm$cell_type == cc),]
    
    cc = sub("\\/", "_", cc)
    venn.diagram(
      x = list(mm.sub$group[which(mm.sub$disease == "ADD")], 
               mm.sub$group[which(mm.sub$disease == "DLB")], 
               mm.sub$group[which(mm.sub$disease == "PDD")]),
      category.names = c("ADD", "DLB" , "PDD"), 
      fill = c("#8CDCDA", "#B1D877", "#F16A70"),  width = 1000, height = 1000, resolution = 300,
      filename = paste0("../figure/deg_overlap/", cc, ".png"),imagetype = "png"
    )
  }
)

##########################################
mm1 = read.csv("../figure/DEG_ADD.txt", sep = "\t")
mm2 = read.csv("../figure/DEG_DLB.txt", sep = "\t")
mm3 = read.csv("../figure/DEG_PDD.txt", sep = "\t")

mm = rbind(mm1, mm2, mm3)
mm$type = "positive"
mm$type[which(mm$avg_log2FC < 0)] = "negative"

mm$group = paste0(mm$gene, "-",mm$cell_type, "-",mm$type)

rna.sub = subset(rna.data, cell.type == "Inhibitory neuron")
rna.sub$group = factor(rna.sub$group, levels = c("CTRL", "ADD", "DLB", "PDD"))


p = VlnPlot(rna.sub, "SNCA", group.by = "group", cols = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70")) + NoLegend() +
  xlab("")
ggsave("../figure/snca.pdf", width = 3, height = 4)


p = VlnPlot(rna.data, "SNCA", group.by = "cell.type", split.by = "group", pt.size = 0, 
        cols = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70")) + xlab("") +
  theme(axis.text=element_text(size=10),
         axis.title=element_text(size=10,face="bold"))
ggsave("../figure/snca_all.pdf", width = 9, height = 4)


rna.sub = subset(rna.data, cell.type == "Microglia")
rna.sub$group = factor(rna.sub$group, levels = c("CTRL", "ADD", "DLB", "PDD"))


p = VlnPlot(rna.sub, "APOE", group.by = "group", cols = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70")) + NoLegend() +
  xlab("")
ggsave("../figure/apoe.pdf", width = 3, height = 4)

p = VlnPlot(rna.data, "APOE", group.by = "cell.type", split.by = "group", pt.size = 0, 
            cols = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70")) + xlab("") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"))
p
ggsave("../figure/apoe_all.pdf", width = 9, height = 4)

rna.sub = subset(rna.data, cell.type == "Astrocyte")
rna.sub$group = factor(rna.sub$group, levels = c("CTRL", "ADD", "DLB", "PDD"))


p = VlnPlot(rna.sub, "BIN1", group.by = "group", cols = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70")) + NoLegend() +
  xlab("")
ggsave("../figure/bin1.pdf", width = 3, height = 4)

p = VlnPlot(rna.data, "BIN1", group.by = "cell.type", split.by = "group", pt.size = 0, 
            cols = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70")) + xlab("") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"))
ggsave("../figure/bin1_all.pdf", width = 9, height = 4)

##########################################

mm1 = read.csv("../figure/DEG_ADD.txt", sep = "\t")
mm2 = read.csv("../figure/DEG_DLB.txt", sep = "\t")
mm3 = read.csv("../figure/DEG_PDD.txt", sep = "\t")

mm = rbind(mm1, mm2, mm3)
mm$type = "positive"
mm$type[which(mm$avg_log2FC < 0)] = "negative"

mm$group = paste0(mm$disease, "-",mm$cell_type)


combinations = expand.grid(1:length(unique(mm$group)), 1:length(unique(mm$group)))
# combinations = combinations[combinations$Var1 < combinations$Var2, , drop = FALSE]
combinations$Var1 = unique(mm$group)[combinations$Var1]
combinations$Var2 = unique(mm$group)[combinations$Var2]


res = pbmclapply(
  1:nrow(combinations), function(i){
    gene1 = mm$gene[which(mm$group == combinations$Var1[i])]
    gene2 = mm$gene[which(mm$group == combinations$Var2[i])]
    
    a = length(intersect(gene1, gene2)) + 1
    b = length(setdiff(gene1, gene2)) + 1
    c = length(setdiff(gene2, gene1)) + 1
    d = length(setdiff(setdiff(unique(mm$gene), gene1), gene2))
    return(a*d/(b*c))
  }
)
combinations$or = unlist(res)
combinations$or = log2(combinations$or)
combinations$or[which(combinations$Var1 == combinations$Var2)] = NA

plot.data = dcast(combinations, Var1 ~ Var2, value.var = "or")
rownames(plot.data) = plot.data[,1]
plot.data = plot.data[,-1]



disease.col = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70")
names(disease.col) = c("CTRL", "ADD", "DLB", "PDD")

cell.type.col = c("#c0b1a5", "#f3de2c", "#5fa8d3", "#ff7d00", "#57cc99", "#e32365", 
                  "#ecaecf", "#2ec4b6")
names(cell.type.col) = levels(rna.data$cell.type)

row.anno = data.frame(disease = str_match(rownames(plot.data), "^(.*)-")[,2],
                      cell.type = str_match(rownames(plot.data), "-(.*$)")[,2],
                      row.names = rownames(plot.data))

row.anno.col = list(
  disease = disease.col, cell.type = cell.type.col
)

pdf("../figure/deg_shared_or.pdf", width = 7, height = 5)
pheatmap::pheatmap(plot.data, annotation_row = row.anno, annotation_colors = row.anno.col)
dev.off()



########################gene disease specificity##############################
rna.data = readRDS("rna_annotated.rds")
rna.data = SetIdent(rna.data, value = "group")

rna.data$sex = NA
rna.data$sex[which(rna.data$sample =="CTRL_4")] = "F"
rna.data$sex[which(rna.data$sample =="CTRL_5")] = "M"
rna.data$sex[which(rna.data$sample =="CTRL_6")] = "M"
rna.data$sex[which(rna.data$sample =="CTRL_7")] = "F"
rna.data$sex[which(rna.data$sample =="ADD_4")] = "F"
rna.data$sex[which(rna.data$sample =="ADD_5")] = "M"
rna.data$sex[which(rna.data$sample =="ADD_6")] = "F"
rna.data$sex[which(rna.data$sample =="ADD_7")] = "F"
rna.data$sex[which(rna.data$sample =="DLB_2")] = "M"
rna.data$sex[which(rna.data$sample =="DLB_4")] = "F"
rna.data$sex[which(rna.data$sample =="DLB_5")] = "M"
rna.data$sex[which(rna.data$sample =="DLB_6")] = "F"
rna.data$sex[which(rna.data$sample =="DLB_7")] = "M"
rna.data$sex[which(rna.data$sample =="DLB_8")] = "F"
rna.data$sex[which(rna.data$sample =="PDD_1")] = "F"
rna.data$sex[which(rna.data$sample =="PDD_2")] = "M"
rna.data$sex[which(rna.data$sample =="PDD_4")] = "F"
rna.data$sex[which(rna.data$sample =="PDD_5")] = "M"
rna.data$sex[which(rna.data$sample =="PDD_7")] = "F"
rna.data$sex[which(rna.data$sample =="PDD_9")] = "M"

disease.specific.mm = pbmclapply(
  unique(rna.data$cell.type), function(ct) {
    rna.sub = subset(rna.data, cell.type == ct)
    mm = FindAllMarkers(rna.sub, logfc.threshold = 1, only.pos = T, latent.vars = "sex", test.use = "MAST")
    mm$cell.type = ct
    mm = mm[which(mm$gene %in% names(which(table(mm$gene) == 1))),]
    # mm = mm[which(mm$pct.1 > 0.5 & mm$pct.2 < 0.1),]
    return(mm)
  }, mc.cores = 20
)
disease.specific.mm = do.call(rbind, disease.specific.mm)
disease.specific.mm = disease.specific.mm[which(disease.specific.mm$pct.1 -disease.specific.mm$pct.2 > 0.5),]


FeaturePlot(rna.data, "GRID2", split.by = "group")


########################gene disease specificity##############################
rna.data = readRDS("rna_annotated.rds")
rna.data = SetIdent(rna.data, value = "group")
mm = FindAllMarkers(rna.sub, logfc.threshold = 1, only.pos = T)
mm = mm[which(mm$p_val_adj < 0.01),]

mm.cell.type = read.csv("../figure/deg_cell_type.txt", sep = "\t")

mm = mm[which(mm$gene %!in% mm.cell.type$gene),]


rna.sub = subset(rna.data, sample %in% c("CTRL_4","CTRL_7","ADD_4","ADD_6","ADD_7","DLB_4","DLB_6","DLB_8","PDD_1","PDD_4","PDD_7"))
mm = mm[which(mm$p_val_adj < 0.01),]

mm2 = FindAllMarkers(rna.sub, logfc.threshold = 1, only.pos = T)

p = VlnPlot(rna.sub, "XIST", group.by = "group", cols = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"), pt.size = 0) + NoLegend() +
  xlab("")
ggsave("../figure/XIST_female.pdf", width = 3, height = 4)


rna.sub = subset(rna.data, sample %in% c("CTRL_5","CTRL_6","ADD_5","DLB_2","DLB_5","DLB_7","PDD_2","PDD_5","PDD_9"))
p = VlnPlot(rna.sub, "XIST", group.by = "group", cols = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"), pt.size = 0) + NoLegend() +
  xlab("")
ggsave("../figure/XIST_male.pdf", width = 3, height = 4)


p = VlnPlot(rna.data, "XIST", group.by = "group", cols = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"), pt.size = 0) + NoLegend() +
  xlab("")
ggsave("../figure/XIST_all.pdf", width = 3, height = 4)


########################gene disease specificity##############################
