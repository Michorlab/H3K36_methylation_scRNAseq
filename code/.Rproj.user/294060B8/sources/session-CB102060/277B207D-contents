####### not used in the end

## load libraries
library(Seurat)
library(tidyverse)
library(bannerCommenter)
source("../funcs-masaki.R")
sample <- "new"



###########################################################################
###########################################################################
###                                                                     ###
###                         LOAD AND CLEAN DATA                         ###
###                                                                     ###
###########################################################################
###########################################################################
dataset_loc <- "~/Dropbox (Michor)/Harvard/projects/michael-michael/sequencing/cellranger-data/single-samples-michael/"
samples <- c("MHb10-1", "MHb10-2", "MHb10-3", "MHb10-4", "MHb10-5", "MHb10-6", "MHb10-7", "MHb10-8", "MHb10-9", "MHb10-10")

michael.data <- sapply(samples, function(i){
  d10x <- Read10X(file.path(dataset_loc,i,"outs/"), )
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
michael.data <- do.call("cbind", michael.data)
rownames(michael.data) <- toupper(rownames(michael.data))
dim(michael.data)

michael <- CreateSeuratObject(
  michael.data,
  min.cells = 10,
  min.features = 300,
  names.field = 2,
  names.delim = "\\-")
dim(michael)
rm(michael.data)

michael <- AddMetaData(object = michael, metadata = sapply(rownames(michael@meta.data), function(x){paste(strsplit(x, "-")[[1]][2:3], collapse= "-")}), col.name = "sample")
michael@meta.data$condition <- recode(michael@meta.data$sample, "MHb10-1" = "d0",
                                     "MHb10-2" = "iPSC", "MHb10-3" = "H3_d2", "MHb10-4" = "K36M_d2",
                                     "MHb10-5" = "H3_d4", "MHb10-6" = "K36M_d4", "MHb10-7" = "H3_d6",
                                     "MHb10-8" = "K36M_d6", "MHb10-9" = "H3_d8", "MHb10-10" = "K36M_d8")
michael@meta.data$day <- recode(michael@meta.data$sample, "MHb10-1" = "d0","MHb10-2" = "iPSC", "MHb10-3" = "d2", "MHb10-4" = "d2",
                                "MHb10-5" = "d4", "MHb10-6" = "d4", "MHb10-7" = "d6", "MHb10-8" = "d6", 
                                "MHb10-9" = "d8", "MHb10-10" = "d8")
michael@meta.data$histone <- recode(michael@meta.data$sample, "MHb10-1" = "d0", "MHb10-2" = "iPSC", "MHb10-3" = "H3", "MHb10-4" = "K36M",
                                    "MHb10-5" = "H3", "MHb10-6" = "K36M", "MHb10-7" = "H3", "MHb10-8" = "K36M", 
                                    "MHb10-9" = "H3", "MHb10-10" = "K36M")
Idents(michael) <- michael@meta.data$condition


## calculate proportion of mitochondrial genes with Seurat
michael[["percent.mt"]] <- PercentageFeatureSet(michael, pattern = "^MT-")

pdf(paste0("../../michael-oncohistones/plots/", sample, "/histmt.pdf"))
hist(michael$percent.mt)
dev.off()


## remove all cells that have mitochondria percentage higher than a given threshold
percent.mito.thresh <- 20
michael <- michael[,-which(michael$percent.mt > percent.mito.thresh)]
dim(michael)
table(michael$condition)



###########################################################################
###########################################################################
###                                                                     ###
###                              NORMALIZE                              ###
###                                                                     ###
###########################################################################
michael <- SCTransform(michael, verbose = TRUE)
michael <- RunPCA(michael, features = VariableFeatures(object = michael, assay = "SCT"), assay = "SCT")
michael <- FindNeighbors(michael, dims = 1:30, verbose = TRUE, assay = "SCT", reduction = "pca", graph.name = "SNN.SCT")
michael <- FindClusters(michael, verbose = TRUE, graph.name = "SNN.SCT")
michael <- RunUMAP(michael, dims = 1:30, verbose = TRUE, reduction.name = "UMAP.SCT", assay = "SCT")
michael <- CellCycleScoring(michael, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)



pdf(paste0("../../michael-oncohistones/plots/", sample, "/umap-sample.pdf"))
DimPlot(michael, group.by = "sample", label = TRUE)
dev.off()

pdf(paste0("../../michael-oncohistones/plots/", sample, "/umap-condition.pdf"))
DimPlot(michael, group.by = "condition", label = TRUE)
dev.off()

pdf(paste0("../../michael-oncohistones/plots/", sample, "/umap-clusters.pdf"))
DimPlot(michael, label = TRUE)
dev.off()

pdf(paste0("../../michael-oncohistones/plots/", sample, "/umap-day.pdf"))
DimPlot(michael, group.by = "day", label = TRUE)
dev.off()

pdf(paste0("../../michael-oncohistones/plots/", sample, "/umap-histone.pdf"))
DimPlot(michael, group.by = "histone", label = TRUE)
dev.off()

pdf(paste0("../../michael-oncohistones/plots/", sample, "/umap-qc.pdf"), width = 12)
FeaturePlot(michael, features = c("percent.mt", "nFeature_SCT"))
dev.off()

pdf(paste0("../../michael-oncohistones/plots/", sample, "/ridge-percentmt.pdf"))
RidgePlot(michael, features = "percent.mt")
dev.off()

pdf(paste0("../../michael-oncohistones/plots/", sample, "/ridge-nFeature_SCT.pdf"))
RidgePlot(michael, features = "nFeature_SCT")
dev.off()

pdf(paste0("../../michael-oncohistones/plots/", sample, "/umap-cellcycle.pdf"))
DimPlot(michael, group.by = "Phase")
dev.off()



###########################################################################
###########################################################################
###                                                                     ###
###                               MARKERS                               ###
###                                                                     ###
###########################################################################
###########################################################################
genes_michael <- read_excel("../../tables/Genes_for_scRNAseq.xlsx")
genes_michael <- as.data.frame(genes_michael)
genes_michael <- toupper(genes_michael$gene)
for (g in genes_michael){
  pdf(paste0("../../michael-oncohistones/plots/", sample, "/genes-", g, ".pdf"), width = 8)
  print(FeaturePlot(michael, features = g))
  dev.off()
}



###########################################################################
###########################################################################
###                                                                     ###
###                       DIFFERENTIAL EXPRESSION                       ###
###                                                                     ###
###########################################################################
###########################################################################
## SCT
p.val.thresh <- 0.1
Idents(michael) <- michael$seurat_clusters
no.clusters.here <- length(unique(Idents(michael)))
de.markers.scale <- list()
for (c in 0:(no.clusters.here-1)){
  print(paste("performing differential expression for cluster", c))
  de.markers.scale[[c+1]] <- get_diff_expr_seurat(michael, ident.1 = c, ident.2 = NULL)
  de.markers.scale[[c+1]] <- de.markers.scale[[c+1]][order(de.markers.scale[[c+1]]$avg_diff, decreasing = TRUE),]
  names(de.markers.scale)[c+1] <- paste0("cluster", c)
  write.table(de.markers.scale[[c+1]], quote = FALSE, sep = "\t", row.names = FALSE,
              file = paste0("../../michael-oncohistones/tables/differential-expression/", sample, "/SCT/", names(de.markers.scale)[c+1], ".txt"))
}


