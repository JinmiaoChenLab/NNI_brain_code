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
`%!in%` <- Negate('%in%')


rna.data = readRDS("rna_annotated.rds")
unique(rna.data$cell.type)
rna.data = subset(rna.data, cell.type %in% c("Oligodendrocyte", "Microglia"))


DefaultAssay(rna.data) = "RNA"
rna.data[["RNA"]] = rna.data[["raw"]] 
rna.data[["raw"]] = NULL

rna.list = SplitObject(rna.data, split.by = "sample")
for (i in 1:20) {
  rna.list[[i]] = FindVariableFeatures(rna.list[[i]])
}

features = SelectIntegrationFeatures(rna.list, nfeatures = 1000)
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
rna.data = RunUMAP(rna.data, dims = 1:10)
rna.data = FindNeighbors(rna.data, dims = 1:10)
rna.data = FindClusters(rna.data, algorithm = 2)

DefaultAssay(rna.data) = "RNA"
FeaturePlot(rna.data, c("ST18", "ARHGAP24"), min.cutoff =  "q10", max.cutoff = "q90")

DimPlot(rna.data, label = T, group.by = "seurat_clusters")

mm = FindMarkers(rna.data, ident.1 = 12, logfc.threshold = 1, only.pos = T)

rna.all = readRDS("rna_annotated.rds")
FeaturePlot(rna.all, "SKAP1")
rna.all$cell.type = factor(rna.all$cell.type, levels = unique(rna.all$cell.type))
colors = data.frame(
  color = c("#c0b1a5", "#f3de2c", "#5fa8d3", "#ff7d00", "#57cc99", "#e32365", 
            "#ecaecf", "#2ec4b6"),
  cell.type = levels(rna.all$cell.type)
)

rna.data = subset(rna.data, seurat_clusters != 12)


rna.data = readRDS("glia_annotated.rds")
mm = read.csv("../figure/deg_cell_type.txt", sep = "\t")
mm = mm[which(mm$p_val_adj < 0.01),]
mm.macro = mm$gene[which(mm$cluster == "Microglia" & mm$avg_log2FC > 0)]
mm.oligo = mm$gene[which(mm$cluster == "Oligodendrocyte" & mm$avg_log2FC > 0)]

DimPlot(rna.data, label = T)
rna.data = SetIdent(rna.data, value = "seurat_clusters")
rna.data$cell.type = "Oligodendrocyte"
rna.data$cell.type[which(rna.data$seurat_clusters %in% c(9))] = "Microglia"
rna.data$cell.type[which(rna.data$seurat_clusters %in% c(27))] = "ME"

rna.data = AddModuleScore(
  object = rna.data,
  features = list(mm.macro),
  # features = list(phagocytosis.gene),
  name = "macro",
  seed = 1
)

rna.data = AddModuleScore(
  object = rna.data,
  features = list(mm.oligo),
  # features = list(myelination.gene),
  name = "oligo",
  seed = 1
)

DimPlot(rna.data, label = T, group.by = "cell.type")
rna.data$cell.type = "Oligodendrocyte"
rna.data$cell.type[which(rna.data$seurat_clusters %in% c(5,9,11))] = "Microglia"

plot.data = data.frame(macro = rna.data$macro1, oligo = rna.data$oligo1, cell.type = rna.data$cell.type)


macro = plot.data$macro[which(plot.data$cell.type == "Microglia")]
# macro.model = Mclust(plot.data$macro)
macro.cut = qnorm(0.005, mean(macro), sd(macro))

oligo = plot.data$oligo[which(plot.data$cell.type == "Oligodendrocyte")]
# oligo.model = Mclust(plot.data$oligo)
oligo.cut = qnorm(0.005, mean(oligo), sd(oligo))

# plot.data$cell.type[which(plot.data$macro > 0.1 & plot.data$oligo > 0.25)] = "ME"
plot.data$cell.type[which(plot.data$macro > macro.cut & plot.data$oligo > oligo.cut)] = "ME"
# plot.data$cell.type[which(plot.data$macro < macro.cut & plot.data$oligo < oligo.cut)] = "remove"

xdensity = ggplot(plot.data, aes(macro, fill=cell.type)) + 
  geom_density(alpha=.5) + xlab("") +
  theme(legend.position = "none") + scale_fill_manual(values = c("#3f51b5", "#e32365", "#ff7d00"))

ydensity = ggplot(plot.data, aes(oligo, fill=cell.type)) + 
  geom_density(alpha=.5) + xlab("") +
  theme(legend.position = "none")+coord_flip() + scale_fill_manual(values = c("#3f51b5", "#e32365", "#ff7d00"))




scatterPlot = ggplot(plot.data, aes(x=macro, y=oligo, color=cell.type)) +
  geom_point() + NoLegend() + scale_color_manual(values = c("#3f51b5", "#e32365", "#ff7d00"))

p = xdensity + plot_spacer() + scatterPlot + ydensity + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

p
ggsave("../figure/MDE_identification.pdf", p, width = 5, height = 5)

table(plot.data$cell.type)
########################################################################
rna.data$cell.type = plot.data$cell.type

rna.data = readRDS("glia_annotated.rds")

plot.data = data.frame(table(rna.data$cell.type, rna.data$sample))
plot.data$ratio = NA

for (i in 1:nrow(plot.data)) {
  plot.data$ratio[i] = plot.data$Freq[i]/(plot.data$Freq[which(plot.data$Var2 == plot.data$Var2[i] & plot.data$Var1 == "Microglia")]+plot.data$Freq[which(plot.data$Var2 == plot.data$Var2[i] & plot.data$Var1 == "Oligodendrocyte")])
} 



plot.data$disease = str_match(plot.data$Var2, "^(.*)_")[,2]
plot.data$disease = factor(plot.data$disease, levels = c("CTRL", "ADD", "DLB", "PDD"))

plot.data.sub  = plot.data[which(plot.data$Var1 == "ME"),]
# write.table(plot.data.sub, "../figure/MDE_ratio.txt", sep = "\t", quote = F, row.names = F, col.names = T)

my_comparisons <- list( c("ADD", "CTRL"), c("ADD", "CTRL"), c("DLB", "CTRL"), c("PDD", "CTRL"))

colors.disease = data.frame(
  color = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"),
  cell.type = c("CTRL", "ADD", "DLB", "PDD")
)

rna.all = readRDS("rna_annotated.rds")
plot.data.sub$ratio = plot.data.sub$Freq/table(rna.all$sample)[plot.data.sub$Var2]
plot.data.sub = plot.data.sub[which(plot.data.sub$Var2 %!in% c("PDD_1","PDD_2","PDD_4","DLB_4")),]

p = ggboxplot(plot.data.sub, "disease", "ratio",
              color = "disease", palette = colors.disease$color,
              add = "jitter")  + NoLegend() + xlab("") + ylab("MDE ratio")

p
ggsave("../figure/ratio_MDE.pdf", p, width = 5, height = 5)

########################################################################

p =DimPlot(rna.data, group.by = "cell.type", cols = c("#3f51b5", "#e32365", "#ff7d00")) 
p
ggsave("../figure/umap_MDE.pdf", p, width = 6.5, height = 5)


saveRDS(rna.data, "glia_annotated.rds", compress = F)


########################################################################
rna.data = readRDS("glia_annotated.rds")
DefaultAssay(rna.data)

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


rna.data = subset(rna.data, cell.type == "Microglia")
rna.data = SetIdent(rna.data, value = "group")
mm = FindMarkers(rna.data, ident.1 = "DLB", test.use = "MAST", latent.vars = "sex")
mm = mm[which(mm$p_val_adj < 0.01),]
mm$gene = rownames(mm)

write.table(mm, "../figure/DLB_vs_control_microglia.txt", sep = "\t", quote = F, row.names = F)

rna.data = subset(rna.data, group %in% c("CTRL", "DLB"))
rna.data$group = factor(as.character(rna.data$group), levels = c("CTRL", "DLB"))
p = VlnPlot(rna.data, "MSR1", group.by = "group", cols = c( "#4D4D4D","#B1D877"))+ NoLegend() +
  xlab("")
p
ggsave("../figure/MSR1.pdf", width = 3, height = 4)


rna.data = readRDS("glia_annotated.rds")
DefaultAssay(rna.data)
rna.sub = subset(rna.data, cell.type == "Oligodendrocyte")


p = FeaturePlot(rna.data, features = "HSPA1A", min.cutoff = "q10", max.cutoff = "q90")
ggsave("../figure/HSPA1A_feature_plot.pdf", width = 5, height = 4.5)

DimPlot(rna.data, label = T)

rna.data$cell.type[which(rna.data$seurat_clusters == 8)] = "HSPA1A+ Oligodendrocyte"

rna.data$cell.type = factor(rna.data$cell.type, levels = unique(rna.data$cell.type))
p = DimPlot(rna.data, group.by = "cell.type", cols = c("#ff7d00", "#3f51b5", "#009688" , "#e32365"))
ggsave("../figure/HSPA1A_cell_type.pdf", width = 6.4, height = 4.5)



plot.data = data.frame(table(rna.data$cell.type, rna.data$sample))
plot.data$ratio = NA
for (i in 1:nrow(plot.data)) {
  plot.data$ratio[i] = plot.data$Freq[i]/(plot.data$Freq[which(plot.data$Var2 == plot.data$Var2[i] & plot.data$Var1 == "Oligodendrocyte")])
} 
plot.data$disease = str_match(plot.data$Var2, "^(.*)_")[,2]
plot.data$disease = factor(plot.data$disease, levels = c("CTRL", "ADD", "DLB", "PDD"))


plot.data.sub  = plot.data[which(plot.data$Var1 == "HSPA1A+ Oligodendrocyte"),]
# write.table(plot.data.sub, "../figure/MDE_ratio.txt", sep = "\t", quote = F, row.names = F, col.names = T)

my_comparisons <- list( c("ADD", "CTRL"), c("ADD", "CTRL"), c("DLB", "CTRL"), c("PDD", "CTRL"))

colors.disease = data.frame(
  color = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"),
  cell.type = c("CTRL", "ADD", "DLB", "PDD")
)

rna.all = readRDS("rna_annotated.rds")
plot.data.sub$ratio = plot.data.sub$Freq/table(rna.all$sample)[plot.data.sub$Var2]
plot.data.sub = plot.data.sub[which(plot.data.sub$Var2 %!in% c("PDD_1","PDD_2","PDD_4","DLB_4")),]

p = ggboxplot(plot.data.sub, "disease", "ratio",
              color = "disease", palette = colors.disease$color,
              add = "jitter")  + NoLegend() + xlab("") + ylab("MDE ratio")

p
table(rna.data$sample[which(rna.data$seurat_clusters == 8)])

rna.sub = subset(rna.data, cell.type == "Oligodendrocyte")
rna.sub = subset(rna.sub, group %in% c("CTRL", "DLB"))
p = VlnPlot(rna.sub, "HSPA1A", group.by = "group", cols = c("#4D4D4D","#B1D877"))+ NoLegend() +
  xlab("")
ggsave("../figure/HSPA1A.pdf", width = 3, height = 4)

