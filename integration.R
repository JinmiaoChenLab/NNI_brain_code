setwd("/scbio7/mengwei/NNI_brain/")
library(pbmcapply)
library(Seurat)
library(scDblFinder)
library(ggplot2)
library(stringr)
options(bitmapType = "cairo")
`%!in%` <- Negate('%in%')
fs = list.files("raw")

rna.list = pbmclapply(
  fs, function(f) {
    rna = Read10X(paste0("raw/",f,"/",f,"/outs/filtered_feature_bc_matrix/")) 
    rna = CreateSeuratObject(rna)
    rna[["percent.mt"]] = PercentageFeatureSet(rna, pattern = "^MT-")
    rna = subset(rna, subset = nFeature_RNA > 500 & percent.mt < 5)
    sce = scDblFinder(rna@assays$RNA@counts)
    rna$doublet = sce@colData@listData[["scDblFinder.class"]]
    rna = subset(rna, doublet == "singlet")
    rna = NormalizeData(rna)
    rna = FindVariableFeatures(rna)
    rna$sample = f
    rna = RenameCells(rna, new.names = paste0(Cells(rna), "--", rna$sample))
    return(rna)
  }, mc.cores = 20
)

saveRDS(rna.list,  "data/rna_list.rds", compress = F)

rna.list = readRDS("data/rna_list.rds")

library(FastIntegration)
BuildIntegrationFile(rna.list = rna.list, tmp.dir = "./", nCores = 30, nfeatures = 5000)
FastFindAnchors(tmp.dir = "./", nCores = 100)

genes = readRDS("./FastIntegrationTmp/raw/1.rds")
genes = rownames(genes)
idx = split(1:length(genes), cut(1:length(genes), 20, labels = FALSE))

pbmclapply(
  1:20, function(i) {
    rna.integrated = FastIntegration(tmp.dir = "./", npcs = 1:30, slot = "data",cut.low = -0.1,
                                     features.to.integrate = genes[idx[[i]]])
    saveRDS(rna.integrated, paste0("./FastIntegrationTmp/inte/inte_", i, ".rds"), 
            compress = F)
  }, mc.cores = 20
)



features = readRDS("./FastIntegrationTmp/others/features.rds")
rna.bind = pbmclapply(
  1:20, function(i) {
    rna = readRDS(paste0("./FastIntegrationTmp/inte/inte_", i, ".rds"))
    rna = rna[intersect(rownames(rna), features),]
    return(rna)
  }, mc.cores = 20
)
rna.bind = do.call(rbind, rna.bind)


rna.bind = CreateSeuratObject(rna.bind)
rna.bind = ScaleData(rna.bind, features = features)
rna.bind = RunPCA(rna.bind, features = features)
rna.bind = RunUMAP(rna.bind, dims = 1:30)
rna.bind = FindNeighbors(rna.bind, dims = 1:30)
rna.bind = FindClusters(rna.bind, algorithm = 2)


rna.list = readRDS("data/rna_list.rds")
rna.list = merge(rna.list[[1]], rna.list[2:20])

rna.bind = AddMetaData(rna.bind, rna.list@meta.data)
rna.bind[["raw"]] = CreateAssayObject(rna.list@assays$RNA@counts)
DefaultAssay(rna.bind) = "raw"
rna.bind = NormalizeData(rna.bind)
rna.bind$group = str_match(rna.bind$sample, "^(.*)_")[,2]

library(DISCOtoolkit)
rna.average = AverageExpression(rna.bind)
ref.data = readRDS("/disco_500t/mengwei/disco_subatlas/idiopathic_pulmonary_fibrosis/DISCOtmp/ref_data.rds")
ref.deg = readRDS("/disco_500t/mengwei/disco_subatlas/idiopathic_pulmonary_fibrosis/DISCOtmp/ref_deg.rds")
predict.ct = CELLiDCluster(rna = rna.average$raw, ref.data = ref.data, ref.deg = ref.deg, n.predict = 3)
rna.bind$predict.ct = predict.ct$predict_cell_type_1[as.numeric(rna.bind$seurat_clusters)]

DimPlot(rna.bind, label = T)

rna.bind$nFeature_raw
FeaturePlot(rna.bind, c("ST18", "PLP1", "PTPRC"), min.cutoff = "q10", max.cutoff = "q90")

# OPC
FeaturePlot(rna.bind, "VCAN", min.cutoff = "q10", max.cutoff = "q90")
# OLI
FeaturePlot(rna.bind, "PLP1", min.cutoff = "q10", max.cutoff = "q90")
# MIC
FeaturePlot(rna.bind, "CD74", min.cutoff = "q10", max.cutoff = "q90")
# AST
FeaturePlot(rna.bind, "SLC14A1", min.cutoff = "q10", max.cutoff = "q90")
# END/PER
FeaturePlot(rna.bind, c("RGS5", "PECAM1"), min.cutoff = "q10", max.cutoff = "q90")
# Neuron
FeaturePlot(rna.bind, c("SYT1"), min.cutoff = "q10", max.cutoff = "q90")
# IN/EN
FeaturePlot(rna.bind, c("GAD2", "SLC17A7"), min.cutoff = "q10", max.cutoff = "q90")

# SST IN
FeaturePlot(rna.bind, c("SST"), min.cutoff = "q10", max.cutoff = "q90")

# L5b EN
FeaturePlot(rna.bind, c("TSHZ2"), min.cutoff = "q10", max.cutoff = "q90")
# L6 EN
FeaturePlot(rna.bind, c("NTNG2"), min.cutoff = "q10", max.cutoff = "q90")
# L2/3 EN
FeaturePlot(rna.bind, c("ENC1"), min.cutoff = "q10", max.cutoff = "q90")
# L4/5 EN
FeaturePlot(rna.bind, c("TSHZ2", "PLCH1"), min.cutoff = "q10", max.cutoff = "q90")
# L4 EN
FeaturePlot(rna.bind, c("RORB"), min.cutoff = "q10", max.cutoff = "q90")
# L5/6 EN
FeaturePlot(rna.bind, c("SEMA3E"), min.cutoff = "q10", max.cutoff = "q90")

FeaturePlot(rna.bind, c("TSHZ1"), min.cutoff = "q10", max.cutoff = "q90")

DimPlot(rna.bind, group.by = "group")


saveRDS(rna.bind, "data/rna_integrated.rds", compress = F)
mm = FindMarkers(rna.bind, ident.1 = c(22),
                 logfc.threshold = 2)

rna.bind$cell.type = NA
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(7,36))] = "OPC"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(0,1,3,20))] = "oligodendrocyte"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(9))] = "microglia"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(2,12,13))] = "astrocyte"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(23,35))] = "Endo/Pericyte"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(27))] = "microglia MDE"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(22))] = "ependymal cell"

#myelin debris engulfment
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(16))] = "SST IN"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(14,33))] = "PVALN IN"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(19))] = "SV2C IN"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(10))] = "VIP IN"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(24))] = "CXCL14 IN"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(17))] = "LUZP2 IN"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(4,8,34))] = "MEIS2 IN"


rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(32))] = "L5b EN"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(30))] = "L6 EN"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(5,6,21))] = "L2/3 EN"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(31,18,26,15))] = "L4/5 EN"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(11))] = "L4 EN"
rna.bind$cell.type[which(rna.bind$seurat_clusters %in% c(25))] = "L5/6 EN"


rna.bind = subset(rna.bind, seurat_clusters %!in% c(28,29))

DimPlot(rna.bind, group.by = "cell.type", label = T)
DimPlot(rna.bind, group.by = "predict.ct", label = T) + NoLegend()


table(rna.bind$predict.ct[which(rna.bind$seurat_clusters == 22)])
saveRDS(rna.bind, "data/rna_annotated.rds", compress = F)
