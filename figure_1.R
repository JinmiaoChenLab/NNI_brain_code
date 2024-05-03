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


rna.data$cell.type[grep("IN", rna.data$cell.type)] = "Inhibitory neuron"
rna.data$cell.type[grep("EN", rna.data$cell.type)] = "Excitatory neuron"
rna.data$cell.type[grep("oligodendrocyte", rna.data$cell.type)] = "Oligodendrocyte"
rna.data$cell.type[grep("astrocyte", rna.data$cell.type)] = "Astrocyte"
rna.data$cell.type[grep("microglia", rna.data$cell.type)] = "Microglia"
rna.data$cell.type[grep("endothelial cell/pericyte", rna.data$cell.type)] = "Endothelial cell"
rna.data$cell.type[grep("ependymal cell", rna.data$cell.type)] = "Ependymal cell"

DimPlot(rna.data, group.by = "cell.type", label = T)

rna.data$cell.type = factor(rna.data$cell.type, levels = unique(rna.data$cell.type))
colors = data.frame(
  color = c("#c0b1a5", "#f3de2c", "#5fa8d3", "#ff7d00", "#57cc99", "#e32365", 
           "#ecaecf", "#2ec4b6"),
  cell.type = levels(rna.data$cell.type)
)


p = DimPlot(rna.data, label = T, group.by = "cell.type", raster = F, cols = colors$color,
            repel = T) + NoLegend() + ggtitle("Cell type")
ggsave("../figure/global_umap.pdf", p, width = 7, height = 7)



slices <- as.numeric(table(rna.data$cell.type))
lbls <- names(table(rna.data$cell.type))
pct <- round(slices/sum(slices)*100)
lbls <- paste0(lbls, "(",pct)
# add percents to labels
lbls <- paste0(lbls,"%",")") # ad % to labels

pdf("../figure/cell_pct_pie.pdf", width = 8, height = 6)
pie(slices, labels = lbls, main="Percentage of Cell Type", col = colors$color)
dev.off()


#####################
table(rna.data$group)

colors.disease = data.frame(
  color = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"),
  cell.type = c("CTRL", "ADD", "DLB", "PDD")
)


plot.data = data.frame(table(rna.data$cell.type, rna.data$sample))
plot.data$ratio = NA

for (i in 1:nrow(plot.data)) {
  plot.data$ratio[i] = plot.data$Freq[i]/sum(plot.data$Freq[which(plot.data$Var2 == plot.data$Var2[i])])
} 

plot.data$disease = str_match(plot.data$Var2, "^(.*)_")[,2]
plot.data$disease = factor(plot.data$disease, levels = c("CTRL", "ADD", "DLB", "PDD"))
plot.data = plot.data[which(plot.data$Var2 %!in% c("PDD_1","PDD_2","PDD_4","DLB_4")),]


plot.list = list()

for (i in 1:length(unique(plot.data$Var1))) {
  
  plot.data.sub  = plot.data[which(plot.data$Var1 == unique(plot.data$Var1)[i]),]
  
  plot.list[[i]] = ggboxplot(plot.data.sub, "disease", "ratio",
                             color = "disease", palette = colors.disease$color,
                             add = "jitter")  + NoLegend() + ggtitle(unique(plot.data$Var1)[i])
}
plot.list[[8]] =  plot.list[[8]] + ylim(c(0,0.03))
p = cowplot::plot_grid(plotlist = plot.list, ncol = 4)
p
ggsave("../figure/pct_disease.pdf", p, width = 12, height = 7)



rna.data = SetIdent(rna.data, value = "cell.type")
rna.sub = subset(rna.data, downsample = 1000)
levels(rna.sub@active.ident)
DefaultAssay(rna.sub) = "raw"
mm = FindAllMarkers(rna.sub, logfc.threshold = 1, only.pos = T)
mm %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 2) %>%
  ungroup() -> top10


p = DotPlot(object = rna.sub, features = c("SYT1", top10$gene)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab("") 
write.table(mm, "../figure/deg_cell_type.txt", sep = "\t", quote = F, row.names = F, col.names = T)
ggsave("../figure/global_marker.pdf", p, width = 8, height = 6)


rna.data$group = factor(rna.data$group, levels = c("CTRL", "ADD", "DLB", "PDD"))
p = DimPlot(rna.data, label = T, split.by = "group", group.by = "cell.type", ncol = 2, repel = T, cols = colors$color, label.size = 3) + NoLegend() + ggtitle("Cell Type")

p
ggsave("../figure/global_umap_disease.pdf", p, width = 12, height = 12)


saveRDS(rna.data, "rna_annotated.rds", compress = F)


####################################
colors.disease = data.frame(
  color = c("#4D4D4D", "#8CDCDA", "#B1D877", "#F16A70"),
  cell.type = c("CTRL", "ADD", "DLB", "PDD")
)


plot.data = data.frame(table(rna.data$cell.type, rna.data$sample))
plot.data$ratio = NA

for (i in 1:nrow(plot.data)) {
  plot.data$ratio[i] = plot.data$Freq[i]/sum(plot.data$Freq[which(plot.data$Var2 == plot.data$Var2[i])])
} 

plot.data$disease = str_match(plot.data$Var2, "^(.*)_")[,2]
plot.data$disease = factor(plot.data$disease, levels = c("CTRL", "ADD", "DLB", "PDD"))
plot.data = plot.data[which(plot.data$Var2 %!in% c("PDD_1","PDD_2","PDD_4","DLB_4")),]
plot.data = dcast(plot.data, Var1 ~ Var2, value.var = "ratio")
rownames(plot.data) = plot.data[,1]
plot.data = plot.data[,-1]
library(ggbiplot)

pc <- prcomp(t(plot.data),
             center = T,
             scale. = F)
plot.data = data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2], disease = str_match(colnames(plot.data), "^(.*)_")[,2])
plot.data$disease = factor(plot.data$disease, levels = c("CTRL", "ADD", "DLB", "PDD"))

p = ggplot(plot.data, aes(x=PC1, y=PC2, color=disease)) + 
  geom_point(size=6) + scale_color_manual(values=colors.disease$color)

ggsave("../figure/pca_pct.pdf", width = 5, height = 4.2)

