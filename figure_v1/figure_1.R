setwd("/scbio7/mengwei/NNI_brain/data/")
library(pbmcapply)
library(Seurat)
library(scDblFinder)
library(ggplot2)
library(stringr)
library(dittoSeq)
library(reshape2)
library(ggpubr)
options(bitmapType = "cairo")
`%!in%` <- Negate('%in%')


rna.data = readRDS("rna_annotated.rds")
rna.data$cell.type[which(rna.data$cell.type == "PVALN IN")] = "PVALB IN"
saveRDS(rna.data, "rna_annotated.rds", compress = F)


rna.data$cell.type[which(rna.data$cell.type == "microglia MDE")] = "microglia"


p = DimPlot(rna.data, label = T, group.by = "cell.type", raster = T, raster.dpi = c(1000,1000),
            repel = T, label.box = T) + NoLegend() + ggtitle("Cell type")
ggsave("../figure/global_umap.pdf", p, width = 7, height = 7)


p = DimPlot(rna.data, label = T, split.by = "group", raster = T, group.by = "cell.type",
        raster.dpi = c(800,800), ncol = 2, repel = T, label.size = 3) + NoLegend()
ggsave("../figure/global_umap_disease.pdf", p, width = 12, height = 12)




plot.data = data.frame(table(rna.data$cell.type, rna.data$sample))
plot.data$ratio = NA

for (i in 1:nrow(plot.data)) {
  plot.data$ratio[i] = plot.data$Freq[i]/sum(plot.data$Freq[which(plot.data$Var2 == plot.data$Var2[i])])
} 

plot.data$disease = str_match(plot.data$Var2, "^(.*)_")[,2]
plot.data$disease = factor(plot.data$disease, levels = c("CTRL", "ADD", "DLB", "PDD"))


plot.list = list()

for (i in 1:length(unique(plot.data$Var1))) {
  
  plot.data.sub  = plot.data[which(plot.data$Var1 == unique(plot.data$Var1)[i]),]
  
  plot.list[[i]] = ggboxplot(plot.data.sub, "disease", "ratio",
                             color = "disease", 
                             add = "jitter")  + NoLegend() + ggtitle(unique(plot.data$Var1)[i])
}
p = cowplot::plot_grid(plotlist = plot.list, ncol = 5)

ggsave("../figure/pct_disease.pdf", p, width = 15, height = 12)


rna.data = SetIdent(rna.data, value = "cell.type")
rna.sub = subset(rna.data, downsample = 1000)
mm = FindAllMarkers(rna.sub, logfc.threshold = 0.5, only.pos = T)

write.table(mm, "../figure/deg_cell_type.txt", sep = "\t", quote = F, row.names = F, col.names = T)

unique(rna.data$cell.type)
rna.sub@active.ident = factor(rna.sub$cell.type, levels=c("SST IN", "CXCL14 IN", "SV2C IN", "LUZP2 IN", "PVALN IN", 
                                                            "MEIS2 IN", "VIP IN", "L2/3 EN", "L4 EN", "L4/5 EN", "L5b EN",
                                                          "L5/6 EN", "L6 EN", "OPC", "oligodendrocyte", "astrocyte", "microglia",
                                                          "endothelial cell/pericyte", "ependymal cell"))
rna.sub[["RNA"]] = rna.sub[["raw"]]
DefaultAssay(rna.sub) = "RNA"
rna.sub[["raw"]] = NULL
p = DotPlot(object = rna.sub, features = c("SYT1", "GAD2", "SST", "CXCL14", "SV2C", "LUZP2", "PVALB",
                                            "MEIS2", "VIP", "SLC17A7", "ENC1", "RORB", "PLCH1", "TSHZ2", 
                                           "SEMA3E", "ZNF804B", "VCAN", "PLP1", "SLC14A1", "CD74",
                                           "RGS5", "PECAM1", "CFAP299"
                                            )) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab("") 
p



ggsave("../figure/global_marker.pdf", p, width = 12, height = 8)


