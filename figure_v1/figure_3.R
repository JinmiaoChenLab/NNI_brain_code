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

options(bitmapType = "cairo")
`%!in%` <- Negate('%in%')


rna.data = readRDS("rna_annotated.rds")
rna.data = subset(rna.data, cell.type %in% c("OPC", "oligodendrocyte", "astrocyte", "microglia", "endothelial cell/pericyte", "microglia MDE"))
DefaultAssay(rna.data) = "RNA"
rna.data[["RNA"]] = rna.data[["raw"]] 
rna.data[["raw"]] = NULL

rna.list = SplitObject(rna.data, split.by = "sample")
for (i in 1:20) {
  rna.list[[i]] = FindVariableFeatures(rna.list[[i]])
}

features = SelectIntegrationFeatures(rna.list, nfeatures = 3000)
rna.bind = pbmclapply(
  1:20, function(i) {
    rna = readRDS(paste0("../FastIntegrationTmp/inte/inte_", i, ".rds"))
    rna = rna[intersect(rownames(rna), features),]
    return(rna)
  }, mc.cores = 20
)
rna.bind = do.call(rbind, rna.bind)

rna.data[["inte"]] = CreateAssayObject(rna.bind[,Cells(rna.data)])
DefaultAssay(rna.data) = "inte"

rna.data = ScaleData(rna.data, features = features)
rna.data = RunPCA(rna.data, features = features)
rna.data = RunUMAP(rna.data, dims = 1:30)
rna.data = FindNeighbors(rna.data, dims = 1:30)
rna.data = FindClusters(rna.data, algorithm = 2)


p = DimPlot(rna.data, label = T, group.by = "cell.type", raster = T, raster.dpi = c(1000,1000),
            repel = T, label.box = T) + NoLegend() + ggtitle("Cell type")
ggsave("../figure/glia_umap.pdf", p, width = 7, height = 7)

DefaultAssay(rna.data) = "RNA"

FeaturePlot(rna.data, "APOE")

DimPlot(rna.data, label = T)
rna.data = subset(rna.data, seurat_clusters %!in% c(21,13,15,8))
rna.data$cell.type[which(rna.data$seurat_clusters == 12)] = "microglia MDE"
rna.data$cell.type[which(rna.data$seurat_clusters == 5)] = "microglia"

saveRDS(rna.data, "glia_annotated.rds")


plot.data = data.frame(table(rna.data$cell.type, rna.data$sample))
plot.data$ratio = NA

for (i in 1:nrow(plot.data)) {
  plot.data$ratio[i] = plot.data$Freq[i]/sum(plot.data$Freq[which(plot.data$Var2 == plot.data$Var2[i])])
} 

plot.data$disease = str_match(plot.data$Var2, "^(.*)_")[,2]
plot.data$disease = factor(plot.data$disease, levels = c("CTRL", "ADD", "DLB", "PDD"))

plot.data.sub  = plot.data[which(plot.data$Var1 == "microglia MDE"),]

my_comparisons <- list( c("ADD", "CTRL"), c("ADD", "CTRL"), c("DLB", "CTRL"), c("PDD", "CTRL"))

p = ggboxplot(plot.data.sub, "disease", "ratio",
                           color = "disease", 
                           add = "jitter")  + NoLegend() 
ggsave("../figure/ratio_MDE.pdf", p, width = 5, height = 6)

FeaturePlot(rna.data, "MSR1")
DefaultAssay(rna.data) = "RNA"
rna.data = SetIdent(rna.data, value = "cell.type")

mm = FindMarkers(rna.data, ident.1 = "microglia MDE", ident.2 = "microglia", logfc.threshold = 0.5)
write.table(mm, "../figure/deg_MDE.txt", sep = "\t", quote = F, row.names = F, col.names = T)



########################################
rna.data = readRDS("glia_annotated.rds")
rna.data = subset(rna.data, cell.type == "microglia")
rna.data = SetIdent(rna.data, value = "group")
mm = FindMarkers(rna.data, ident.1 = "DLB", ident.2 = "CTRL", logfc.threshold = 0.5, only.pos = T)
mm = mm[which(mm$p_val_adj < 0.01),]
mm = mm[order(-mm$avg_log2FC),]
mm = mm[1:50,]
VlnPlot(rna.data, "MSR1", group.by = "group", sort = "increasing")

go_enrich <- enrichGO(gene = rownames(mm),
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

pdf("../figure/DLB_marco_pathway.pdf", height = 8, width = 12)
treeplot(pairwise_termsim(go_enrich))
dev.off()

rna.data = ScaleData(rna.data, features = rownames(mm))
DoHeatmap(rna.data, features = rownames(mm)[1:50], raster = T)

rna.average = AverageExpression(rna.data, group.by = "group")
rna.average = rna.average$RNA[rownames(mm)[1:50],]

pdf("../figure/DLB_marco_marker.pdf")
pheatmap::pheatmap(rna.average, scale = "row")
dev.off()
VlnPlot(rna.data, "APOE", group.by = "group", sort = "increasing")


table(rna.data$group)
