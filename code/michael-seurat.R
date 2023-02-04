## Single cell analysis on H3K36 methylation data: Seurat analysis
## Author: Simona Cristea
## Date: November 2022
## Publication: H3K36 methylation maintains cell identity by regulating opposing lineage programs

source(here::here("libraries.R"))


###########################################################################
###########################################################################
###                                                                     ###
###              READ IN MONOCLE OBJECTS AND CREATE SEURAT              ###
###                                                                     ###
###########################################################################
###########################################################################
monocle.list <- readRDS(file = "../data/rdata/michael-monocle-matrices-raw.RDS")

michael <- CreateSeuratObject(monocle.list$m, min.cells = 0, min.features = 0)
michael <- AddMetaData(michael, as.data.frame(monocle.list$ann.data))



## normalize
michael <- SCTransform(michael, verbose = TRUE)
michael <- RunPCA(michael, features = VariableFeatures(object = michael, assay = "SCT"), assay = "SCT")
michael <- RunUMAP(michael, dims = 1:30, verbose = TRUE)
michael <- FindNeighbors(michael, dims = 1:30, verbose = TRUE)
michael <- FindClusters(michael, verbose = TRUE)
michael <- CellCycleScoring(michael, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)



############################################################################
############################################################################
###                                                                      ###
###                               PLOTTING                               ###
###                                                                      ###
############################################################################
############################################################################
pdf("../plots/seurat/umap-sample.pdf")
DimPlot(michael, group.by = "sample", label = TRUE)
dev.off()

pdf("../plots/seurat/umap-clusters.pdf")
DimPlot(michael, label = TRUE)
dev.off()

pdf("../plots/seurat/umap-day.pdf")
DimPlot(michael, group.by = "day", label = TRUE)
dev.off()

pdf("../plots/seurat/umap-histone.pdf")
DimPlot(michael, group.by = "histone", label = TRUE)
dev.off()

pdf("../plots/seurat/umap-qc.pdf", width = 12)
FeaturePlot(michael, features = c("percent.mt", "nFeature_SCT"))
dev.off()

pdf("../plots/seurat/umap-cellcycle.pdf")
DimPlot(michael, group.by = "Phase")
dev.off()

pdf("../plots/seurat/ridge-percentmt.pdf")
RidgePlot(michael, features = c("percent.mt"))
dev.off()

pdf("../plots/seurat/ridge-n_FeatureSCT.pdf")
RidgePlot(michael, features = c("nFeature_SCT"))
dev.off()


#### plot markers
genes_michael <- read_excel("../../tables/Genes_for_scRNAseq.xlsx")
genes_michael <- as.data.frame(genes_michael)
genes_michael <- toupper(genes_michael$gene)
for (g in genes_michael){
  pdf(paste0("../plots/seurat/markers/genes-", g, ".pdf"), width = 8)
  print(FeaturePlot(michael, features = g))
  dev.off()
}

pdf("../plots/seurat/markers/genes-CD24A.pdf", width = 8)
FeaturePlot(michael, features = "CD24A")
dev.off()

#saveRDS(michael, file = "../../rdata/michael-seurat.RDS")


library(viridis)
DefaultAssay(michael) <- "RNA"
michael <- NormalizeData(michael)
Idents(michael) <- michael$sample
genes.dot.plot <- rev(toupper(c("Vim", "Fn1", "Prrx1", "Col1a2", "Twist1", "Zeb1", "Mef2c", "Nanog",
                                "Sall4", "Pou5f1", "Cldn4", "Cldn6", "Cdh1", "Krt8", "Krt18", "Epcam")))
pdf("../../plots/seurat/dot-plot-short.pdf", width = 10)
DotPlot(michael, features = unique(genes.dot.plot), idents = c("H3-d2", "K36M-d2", "H3-d4", "K36M-d4"),
        assay = "RNA") + scale_color_viridis() + coord_flip()
dev.off()



###########################################################################
###########################################################################
###                                                                     ###
###                            DPT DIFFUSION                            ###
###                                                                     ###
###########################################################################
###########################################################################
library(destiny)

## read in dpt object
michael.dpt <- readRDS("../../rdata/michael-dpt.RDS")
plot(michael.dpt$dpt.histones@dm$DC1,michael.dpt$dpt.histones@dm$DC2)

write.table(cbind(names(michael.dpt$dpt.histones@dm$DC1), michael.dpt$dpt.histones@dm$DC1, michael.dpt$dpt.histones@dm$DC2), file = "../tables/michael-DPT-histones.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)




