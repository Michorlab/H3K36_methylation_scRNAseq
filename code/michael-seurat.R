#### done from the monocle CDS 

## load libraries
library(Seurat)
library(tidyverse)

## read in Monocle object
monocle.list <- readRDS(file = "../../rdata/michael-monocle-matrices-raw.RDS")

## create Seurat from it
michael <- CreateSeuratObject(monocle.list$m, min.cells = 0, min.features = 0)
michael <- AddMetaData(michael, as.data.frame(monocle.list$ann.data))

## normalize
michael <- SCTransform(michael, verbose = TRUE)
michael <- RunPCA(michael, features = VariableFeatures(object = michael, assay = "SCT"), assay = "SCT")
michael <- RunUMAP(michael, dims = 1:30, verbose = TRUE)
michael <- FindNeighbors(michael, dims = 1:30, verbose = TRUE)
michael <- FindClusters(michael, verbose = TRUE)
michael <- CellCycleScoring(michael, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)


## plotting
pdf("../../michael-oncohistones/plots/seurat/umap-sample.pdf")
DimPlot(michael, group.by = "sample", label = TRUE)
dev.off()

pdf("../../michael-oncohistones/plots/seurat/umap-clusters.pdf")
DimPlot(michael, label = TRUE)
dev.off()

pdf("../../michael-oncohistones/plots/seurat/umap-day.pdf")
DimPlot(michael, group.by = "day", label = TRUE)
dev.off()

pdf("../../michael-oncohistones/plots/seurat/umap-histone.pdf")
DimPlot(michael, group.by = "histone", label = TRUE)
dev.off()

pdf("../../michael-oncohistones/plots/seurat/umap-qc.pdf", width = 12)
FeaturePlot(michael, features = c("percent.mt", "nFeature_SCT"))
dev.off()

pdf("../../michael-oncohistones/plots/seurat/umap-cellcycle.pdf")
DimPlot(michael, group.by = "Phase")
dev.off()

pdf("../../michael-oncohistones/plots/seurat/ridge-percentmt.pdf")
RidgePlot(michael, features = c("percent.mt"))
dev.off()

pdf("../../michael-oncohistones/plots/seurat/ridge-n_FeatureSCT.pdf")
RidgePlot(michael, features = c("nFeature_SCT"))
dev.off()


#### plot markers
genes_michael <- read_excel("../../tables/Genes_for_scRNAseq.xlsx")
genes_michael <- as.data.frame(genes_michael)
genes_michael <- toupper(genes_michael$gene)
for (g in genes_michael){
  pdf(paste0("../../michael-oncohistones/plots/seurat/markers/genes-", g, ".pdf"), width = 8)
  print(FeaturePlot(michael, features = g))
  dev.off()
}

pdf("../../michael-oncohistones/plots/seurat/markers/genes-CD24A.pdf", width = 8)
FeaturePlot(michael, features = "CD24A")
dev.off()

#saveRDS(michael, file = "../../rdata/michael-seurat.RDS")
library(tidyverse)
library(Seurat)
michael <- readRDS("../../rdata/michael-seurat.RDS")


## differential expression
michael <- NormalizeData(michael, assay = "RNA")
Idents(michael) <- michael$sample
p.val.thresh <- 0.1
de <- list()
de$d2 <- FindMarkers(michael, ident.1 = "K36M-d2", ident.2 = "H3-d2", assay = "RNA")
de$d4 <- FindMarkers(michael, ident.1 = "K36M-d4", ident.2 = "H3-d4", assay = "RNA")
de$d6 <- FindMarkers(michael, ident.1 = "K36M-d6", ident.2 = "H3-d6", assay = "RNA")
de$d8 <- FindMarkers(michael, ident.1 = "K36M-d8", ident.2 = "H3-d8", assay = "RNA")
for (i in 1:length(de)){
  de[[i]] <- cbind("gene" = rownames(de[[i]]), de[[i]])
  de[[i]] <- as.data.frame(de[[i]])
  de[[i]] <- de[[i]][which(de[[i]]$p_val_adj < p.val.thresh),]
  de[[i]] <- de[[i]][order(abs(de[[i]]$avg_log2FC), decreasing = TRUE),]
  write.table(de[[i]], paste0("../../michael-oncohistones/tables/differential-expression/seurat/K36MvsH3-", names(de)[i], ".txt"),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
}

no.clusters.integrated <- length(unique(Idents(integrated)))
p.val.thresh <- 0.1
de.markers.integrated <- list()
for (c in 0:(no.clusters.integrated-1)){
  print(paste("performing differential expression for cluster", c))
  de.markers.integrated[[c+1]] <- FindMarkers(integrated, ident.1 = c, ident.2 = NULL, assay = "RNA")
  names(de.markers.integrated)[c+1] <- paste0("cluster", c)
  de.markers.integrated[[c+1]] <- cbind("gene" = rownames(de.markers.integrated[[c+1]]), de.markers.integrated[[c+1]])
  de.markers.integrated[[c+1]] <- as.data.frame(de.markers.integrated[[c+1]])
  de.markers.integrated[[c+1]] <- de.markers.integrated[[c+1]][which(de.markers.integrated[[c+1]]$p_val_adj < p.val.thresh),]
  de.markers.integrated[[c+1]] <- de.markers.integrated[[c+1]][order(de.markers.integrated[[c+1]]$avg_logFC, decreasing = TRUE),]
  write.table(de.markers.integrated[[c+1]], quote = FALSE, sep = "\t", row.names = FALSE,
              file = paste0("../../TGFBRi study/Simona scRNAseq/new-tables/", sample, "/differential-expression/integrated/RNA/res0.41/", names(de.markers.integrated)[c+1], ".txt"))
}

