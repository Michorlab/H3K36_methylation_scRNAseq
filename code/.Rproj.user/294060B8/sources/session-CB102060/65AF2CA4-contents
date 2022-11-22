source(here::here("libraries.R"))



############################################################################
############################################################################
###                                                                      ###
###                         READ IN CLEANED DATA                         ###
###                                                                      ###
############################################################################
############################################################################
seuratX <- readRDS(file = "../data/rdata/michael-initial.RDS")


## remove low genes clusters and contamination (clusters 9,10,11): 372 cells
idx.to.remove <- which(cds@clusters$UMAP$clusters %in% c(9,10,11))
#saveRDS(idx.to.remove, file = "../../rdata/michael-idx.to.remove.RDS")
cds <- cds[, -idx.to.remove]

## remove low genes clusters (new cluster 10): 903 cells
idx.to.remove.2 <- which(cds@clusters$UMAP$clusters %in% c(10))
#saveRDS(idx.to.remove.2, file = "../../rdata/michael-idx.to.remove.2.RDS")
cds <- cds[, -idx.to.remove.2]


## extract RNA count data from Seurat
expression_matrix <- as(as.matrix(seuratX@assays$RNA@counts), 'sparseMatrix')
cell_metadata <- seuratX@meta.data
cell_metadata <- cbind("barcode" =rownames(cell_metadata), cell_metadata)
gene_metadata <- data.frame(gene_short_name = row.names(expression_matrix), row.names = row.names(expression_matrix))


## construct monocle cds
cds <- new_cell_data_set(expression_data = expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
rm(expression_matrix)

## correct label
cds@colData$sample[which(cds@colData$sample == "d9")] <- "iPS"


## view data
pData(cds)
fData(cds)


## preprocess data
cds <- preprocess_cds(cds, num_dim = 100)
pdf("../../michael-oncohistones/plots/monocle/variance.explained.pdf")
plot_pc_variance_explained(cds)
dev.off()


## reduce dimension and plot
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
pdf("../../michael-oncohistones/plots/monocle/umap-sample.pdf")
plot_cells(cds, color_cells_by = "sample", label_cell_groups = FALSE, group_label_size = 5)
dev.off()
pdf("../../michael-oncohistones/plots/monocle/umap-histone.pdf")
plot_cells(cds, color_cells_by = "histone", label_cell_groups = TRUE, group_label_size = 5)
dev.off()
pdf("../../michael-oncohistones/plots/monocle/umap-percent.mt.pdf")
plot_cells(cds, color_cells_by = "percent.mt", label_cell_groups = FALSE)
dev.off()
pdf("../../michael-oncohistones/plots/monocle/umap-detected.pdf")
plot_cells(cds, color_cells_by = "detected", label_cell_groups = FALSE)
dev.off()
pdf("../../michael-oncohistones/plots/monocle/umap-sum.pdf")
plot_cells(cds, color_cells_by = "sum", label_cell_groups = FALSE)
dev.off()


## cluster cells and plot
cds <- cluster_cells(cds, resolution=1e-5)

#pdf("../../michael-oncohistones/plots/monocle/umap-cluster.pdf")
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = TRUE)
#dev.off()
#pdf("../../michael-oncohistones/plots/monocle/umap-partition.pdf")
plot_cells(cds, color_cells_by = "partition", label_cell_groups = TRUE)
#dev.off()

head(cds@clusters$UMAP$partitions)
head(cds@clusters$UMAP$clusters)


pd <- pData(cds.regressed)
pd$partitions <- cds.regressed@clusters$UMAP$partitions
pd <- as_tibble(pd)

pd <- pData(cds)
pd$partitions <- cds@clusters$UMAP$partitions
pd$clusters <- cds@clusters$UMAP$clusters
pd <- as_tibble(pd)

pdf("../../michael-oncohistones/plots/monocle/density_percentmt-partitions.pdf")
ggplot(pd, aes (x = percent.mt)) + geom_density() + facet_wrap(~partitions)
dev.off()
pdf("../../michael-oncohistones/plots/monocle/density_detected-partitions.pdf")
ggplot(pd, aes (x = detected)) + geom_density() + facet_wrap(~partitions)
dev.off()

pdf("../../michael-oncohistones/plots/monocle/density_percentmt-clusters.pdf")
ggplot(pd, aes (x = percent.mt)) + geom_density() + facet_wrap(~clusters)
dev.off()
pdf("../../michael-oncohistones/plots/monocle/density_detected-clusters.pdf")
ggplot(pd, aes (x = detected)) + geom_density() + facet_wrap(~clusters)
dev.off()


## plot markers
genes_michael <- read_excel("../../tables/Genes_for_scRNAseq.xlsx")
genes_michael <- as.data.frame(genes_michael)
genes_michael <- toupper(genes_michael$gene)
for (g in genes_michael){
  pdf(paste0("../../michael-oncohistones/plots/monocle/markers/genes-", g, ".pdf"), width = 8)
  print(plot_cells(cds, genes=g, color_cells_by = "sample"))
  dev.off()
}


## learn graph and plot (each partition becomes its own trajectory)
cds <- learn_graph(cds)
pdf("../../michael-oncohistones/plots/monocle/umap-clusters-trajectory.pdf")
plot_cells(cds, group_label_size = 4)
dev.off()


## detailed trajectory plots
# leaves (or gray circles with numbers in them) correspond to different outcomes (i.e. cell fate) of the trajectory
# branch nodes (black circles) are points in which cells can travel to one of several outcomes
pdf("../../michael-oncohistones/plots/monocle/umap-trajectory.pdf")
plot_cells(cds, color_cells_by = "sample", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, group_label_size = 4)
dev.off()

pdf("../../michael-oncohistones/plots/monocle/umap-trajectory-detailed.pdf")
plot_cells(cds, color_cells_by = "sample", label_groups_by_cluster=FALSE, label_leaves = TRUE, label_branch_points=TRUE, group_label_size = 4)
dev.off()


## loupe browser
write.table(file = "~/Downloads/michael-umap-monocle.txt", sep = "\t", quote = FALSE, row.names = TRUE,
            reducedDims(cds)$UMAP)
#cbind("barcode" = colnames(cds), cds@clustersreductions$UMAP.SCT@cell.embeddings[,c(1:2)]))
write.table(file = "~/Downloads/michael-metadata-monocle.txt", sep = "\t", quote = FALSE, row.names = FALSE,
            cbind("barcode" = colnames(cds), cds@colData))


#### add cell cycle inference 
# michael <- readRDS(file = "../../rdata/michael-seurat.RDS")
# all.equal(unname(michael$barcode), colData(cds)$barcode)
# all.equal(dim(michael)[2], dim(cds)[2])
# colData(cds)$G2M.Score <- michael$G2M.Score
# colData(cds)$S.Score <- michael$S.Score
# colData(cds)$Phase <- michael$Phase

pdf("../../michael-oncohistones/plots/monocle/umap-Phase.pdf")
plot_cells(cds, color_cells_by = "Phase", show_trajectory_graph = FALSE, label_cell_groups = FALSE)
dev.off()

pdf("../../michael-oncohistones/plots/monocle/umap-S.Score.pdf")
plot_cells(cds, color_cells_by = "S.Score", show_trajectory_graph = FALSE, label_cell_groups = FALSE)
dev.off()

pdf("../../michael-oncohistones/plots/monocle/umap-G2M.Score.pdf")
plot_cells(cds, color_cells_by = "G2M.Score", show_trajectory_graph = FALSE, label_cell_groups = FALSE)
dev.off()





#### load data
library(monocle3)
library(tidyverse)
library(viridis)
library(scales)
cds <- readRDS("../../rdata/michael-monocle.RDS")


#### compare with atlas
atlas <- read.table("../../michael-oncohistones/papers/gkx054_Supplementary_Data/final_table_for_supp.tsv", 
                    header = TRUE, sep = "\t")
atlas$name <- toupper(atlas$name)
atlas <- atlas[which(atlas$name != ""),]
pseudo <- cbind("gene" = atlas$name, "ESCsV6.5" = atlas$ESCs..V6.5., "MEF" = atlas$MEF, "Keratinocyte" = atlas$Keratinocyte)
pseudo <- as.data.frame(pseudo)
common <- intersect(rownames(cds), pseudo$gene)
pseudo <- pseudo[match(common, pseudo$gene),]
cds <- cds[match(common, rownames(cds)),]
all.equal(pseudo$gene, rownames(cds))

cds.normalized.counts <- normalized_counts(cds)
md <- pData(cds)


## compute correlations on all genes of atlas
cds.correlations <- apply(cds.normalized.counts, 2, function(x){
  cor.ESCsV6.5 <- cor(x, as.numeric(pseudo$ESCsV6.5))
  cor.MEFs <- cor(x, as.numeric(pseudo$MEF))
  return(c(cor.ESCsV6.5, cor.MEFs))
  #cor.Keratinocyte <- cor(x, as.numeric(pseudo$Keratinocyte))
  #return(c(cor.ESCsV6.5, cor.MEFs, cor.Keratinocyte))
})
#rownames(cds.correlations) <- c("ESCsV6.5", "MEFs", "Keratinocyte")
rownames(cds.correlations) <- c("ESCsV6.5", "MEFs")
pData(cds)$ESCsV6.5 <- cds.correlations["ESCsV6.5",]
pData(cds)$MEFs <- cds.correlations["MEFs",]
#pData(cds)$Keratinocyte <- cds.correlations["Keratinocyte",]


cds.correlations.zscore <- apply(cds.correlations, 1, function(x){
  return((x - mean(x))/sd(x))
})
cds.correlations.zscore <- t(cds.correlations.zscore)


no.breaks.col <- 1000
gene <- "POU5F1"
vir <- viridis(no.breaks.col)
unique.samples <- unique(md$sample)
vals.gene <- cds.normalized.counts[which(rownames(cds.normalized.counts) == gene),]
cols <- vir[as.numeric(cut(vals.gene, breaks = no.breaks.col))]
#rbPal <- colorRampPalette(c('red','blue'))
pdf(paste0("../../michael-oncohistones/plots/projections/cor.atlas.allgenes.zscore-enh-", gene, ".pdf"), width = 12, height = 9.5)
par(mfrow=c(3,4))
unique.samples <- c("d0", "H3-d2", "H3-d4", "H3-d6", "H3-d8", "K36M-d2", "K36M-d4", "K36M-d6", "K36M-d8", "iPS")
for (u in unique.samples){
  print(plot(cds.correlations.zscore["ESCsV6.5",], cds.correlations.zscore["MEFs",], col = alpha("gray", alpha = 0.1), 
             #xlab = "z-score of correlation to ESCsV6.5", ylab = "z-score of correlation to MEFs", 
             xlab = "", ylab = "",
             main = u, cex = 0.7, pch = 20))
  #cols <- rbPal(10)[as.numeric(cut(vals.gene, breaks = 10))]
  cols.here <- cols[which(md$sample == u)]
  points(cds.correlations.zscore["ESCsV6.5",which(md$sample == u)], 
         cds.correlations.zscore["MEFs",which(md$sample == u)], col = cols.here, cex = 0.7, pch = 20)
}
dev.off()



pdf("../../michael-oncohistones/plots/projections/cor.allgenes.ESCsV6.5.pdf")
plot_cells(cds, color_cells_by = "ESCsV6.5", show_trajectory_graph = FALSE)
dev.off()
pdf("../../michael-oncohistones/plots/projections/cor.allgenes.MEFs.pdf")
plot_cells(cds, color_cells_by = "MEFs", show_trajectory_graph = FALSE)
dev.off()

unique.samples <- unique(md$sample)
#pdf("../../michael-oncohistones/plots/projections/cor.atlas.allgenes.pdf")
par(mfrow=c(3,4))
unique.samples <- c("d0", "H3-d2", "H3-d4", "H3-d6", "H3-d8", "K36M-d2", "K36M-d4", "K36M-d6", "K36M-d8", "iPS")
for (u in unique.samples){
  md.partial <- md[-(which(md$sample == u)),]
  cds.correlations.partial <- cds.correlations[,-which(md$sample == u)]
  print(plot(cds.correlations.partial["ESCsV6.5",], cds.correlations.partial["MEFs",], col = "gray", 
             xlab = "correlation to ESCsV6.5", ylab = "correlation to MEFs", main = u, xlim = c(0,1), ylim = c(0,1)))
  points(cds.correlations["ESCsV6.5",which(md$sample == u)], 
         cds.correlations["MEFs",which(md$sample == u)], col = "blue")
}
#dev.off()


## compute correlations on genes of interest only on atlas
genes.of.interest <- rownames(cds)[order(apply(cds.normalized.counts,1,sd), decreasing = TRUE)[1:3000]]
cds.correlations.short <- apply(cds.normalized.counts[match(genes.of.interest, rownames(cds.normalized.counts)),], 2, function(x){
  cor.ESCsV6.5 <- cor(x, as.numeric(pseudo$ESCsV6.5[match(genes.of.interest, pseudo$gene)]))
  cor.MEFs <- cor(x, as.numeric(pseudo$MEF[match(genes.of.interest, pseudo$gene)]))
  return(c(cor.ESCsV6.5, cor.MEFs))
})
rownames(cds.correlations.short) <- c("ESCsV6.5", "MEFs")
pData(cds)$ESCsV6.5.short <- cds.correlations.short["ESCsV6.5",]
pData(cds)$MEFs.short <- cds.correlations.short["MEFs",]



#### do pseudo-bulk
md <- pData(cds)
cds.normalized.counts <- normalized_counts(cds)

pseudo.d0 <- rowSums(cds.normalized.counts[,which(md$sample == "d0")])
pseudo.iPS <- rowSums(cds.normalized.counts[,which(md$sample == "iPS")])


## compute correlations on all genes on pseudo-bulk
cds.correlations <- apply(cds.normalized.counts, 2, function(x){
  cor.d0 <- cor(x, as.numeric(pseudo.d0))
  cor.iPS <- cor(x, as.numeric(pseudo.iPS))
  return(c(cor.d0, cor.iPS))
})
rownames(cds.correlations) <- c("d0", "iPS")
pData(cds)$cord0 <- cds.correlations["d0",]
pData(cds)$coriPS <- cds.correlations["iPS",]

pdf("../../michael-oncohistones/plots/projections/cor.allgenes.d0.pdf")
plot_cells(cds, color_cells_by = "cord0", show_trajectory_graph = FALSE)
dev.off()
pdf("../../michael-oncohistones/plots/projections/cor.allgenes.iPS.pdf")
plot_cells(cds, color_cells_by = "coriPS", show_trajectory_graph = FALSE)
dev.off()

cds.correlations <- apply(cds.normalized.counts, 2, function(x){
  cor.d0 <- cor(x, as.numeric(pseudo.d0))
  cor.iPS <- cor(x, as.numeric(pseudo.iPS))
  return(c(cor.d0, cor.iPS))
})
rownames(cds.correlations) <- c("d0", "iPS")
pData(cds)$cord0 <- cds.correlations["d0",]
pData(cds)$coriPS <- cds.correlations["iPS",]

unique.samples <- unique(md$sample)
pdf("../../michael-oncohistones/plots/projections/cor.diff.expr.allgenes.pdf")
par(mfrow=c(3,4))
unique.samples <- c("d0", "H3-d2", "H3-d4", "H3-d6", "H3-d8", "K36M-d2", "K36M-d4", "K36M-d6", "K36M-d8", "iPS")
for (u in unique.samples){
  md.partial <- md[-(which(md$sample == u)),]
  cds.correlations.partial <- cds.correlations[,-which(md$sample == u)]
  print(plot(cds.correlations.partial["d0",], cds.correlations.partial["iPS",], col = "gray", 
             xlab = "correlation to pseudobulk d0", ylab = "correlation to pseudobulk iPS", main = u, xlim = c(0.3,1), ylim = c(0.3,1)))
  points(cds.correlations["d0",which(md$sample == u)], 
         cds.correlations["iPS",which(md$sample == u)], col = "blue")
}
dev.off()




## compute correlations on genes of interest only on pseudo-bulk
#genes.of.interest <- rownames(cds)[order(apply(cds.normalized.counts,1,sd), decreasing = TRUE)[1:3000]]
genes.up <- read.table(file = "../../michael-oncohistones/tables/bulk/H3_3_pluri_rpkm.tsv", sep = "\t", header = TRUE)
genes.up <- genes.up[order(genes.up$logFC, decreasing = TRUE),]
genes.down <- read.table(file = "../../michael-oncohistones/tables/bulk/H3_3_fibro_rpkm.tsv", sep = "\t", header = TRUE)
genes.down <- genes.down[order(genes.down$logFC),]

no.genes <- 500
genes.up <- genes.up$gene[1:no.genes]
genes.down <- genes.down$gene[1:no.genes]
genes.of.interest <- c(genes.up, genes.down)
#genes.of.interest <- genes.down
genes.of.interest <- toupper(genes.of.interest)
genes.of.interest <- intersect(rownames(cds.normalized.counts), genes.of.interest)


cds.correlations.short <- apply(cds.normalized.counts[match(genes.of.interest, rownames(cds.normalized.counts)),], 2, function(x){
  cor.d0 <- cor(x, as.numeric(pseudo.d0[match(genes.of.interest, names(pseudo.d0))]))
  cor.iPS <- cor(x, as.numeric(pseudo.iPS[match(genes.of.interest, names(pseudo.iPS))]))
  return(c(cor.d0, cor.iPS))
})
rownames(cds.correlations.short) <- c("d0", "iPS")
pData(cds)$cord0.short <- cds.correlations.short["d0",]
pData(cds)$coriPS.short <- cds.correlations.short["iPS",]




#### compute correlations on atlas using differentially expressed genes on the single cell data
atlas <- read.table("../../michael-oncohistones/papers/gkx054_Supplementary_Data/final_table_for_supp.tsv", 
                    header = TRUE, sep = "\t")
atlas$name <- toupper(atlas$name)
atlas <- atlas[which(atlas$name != ""),]
pseudo <- cbind("gene" = atlas$name, "ESCsV6.5" = atlas$ESCs..V6.5., "MEF" = atlas$MEF, "Keratinocyte" = atlas$Keratinocyte)
pseudo <- as.data.frame(pseudo)
common <- intersect(rownames(cds), pseudo$gene)
pseudo <- pseudo[match(common, pseudo$gene),]
cds <- cds[match(common, rownames(cds)),]
all.equal(pseudo$gene, rownames(cds))

cds.normalized.counts <- normalized_counts(cds)
md <- pData(cds)

diff.genes.single.cell <- read.table("../../michael-oncohistones/tables/differential-expression/seurat/d0-iPS.txt", sep = "\t", header = TRUE)
genes.of.interest <- diff.genes.single.cell$gene
genes.of.interest <- intersect(rownames(cds.normalized.counts), genes.of.interest)

cds.correlations.short <- apply(cds.normalized.counts[match(genes.of.interest, rownames(cds.normalized.counts)),], 2, function(x){
  cor.ESCsV6.5 <- cor(x, as.numeric(pseudo$ESCsV6.5[match(genes.of.interest, pseudo$gene)]))
  cor.MEFs <- cor(x, as.numeric(pseudo$MEF[match(genes.of.interest, pseudo$gene)]))
  return(c(cor.ESCsV6.5, cor.MEFs))
})
rownames(cds.correlations.short) <- c("ESCsV6.5", "MEFs")
pData(cds)$ESCsV6.5.short <- cds.correlations.short["ESCsV6.5",]
pData(cds)$MEFs.short <- cds.correlations.short["MEFs",]


cds.correlations.zscore <- apply(cds.correlations.short, 1, function(x){
  return((x - mean(x))/sd(x))
})
cds.correlations.zscore <- t(cds.correlations.zscore)


no.breaks.col <- 1000
gene <- "EPCAM"
vir <- viridis(no.breaks.col)
unique.samples <- unique(md$sample)
vals.gene <- cds.normalized.counts[which(rownames(cds.normalized.counts) == gene),]
cols <- vir[as.numeric(cut(vals.gene, breaks = no.breaks.col))]
#rbPal <- colorRampPalette(c('red','blue'))
pdf(paste0("../../michael-oncohistones/plots/projections/cor.atlas.single.cell.diff.expr-enh-", gene, ".pdf"), width = 12, height = 9.5)
par(mfrow=c(3,4))
unique.samples <- c("d0", "H3-d2", "H3-d4", "H3-d6", "H3-d8", "K36M-d2", "K36M-d4", "K36M-d6", "K36M-d8", "iPS")
for (u in unique.samples){
  print(plot(cds.correlations.zscore["ESCsV6.5",], cds.correlations.zscore["MEFs",], col = alpha("gray", alpha = 0.1), 
             #xlab = "z-score of correlation to ESCsV6.5", ylab = "z-score of correlation to MEFs", 
             xlab = "", ylab = "",
             main = u, cex = 0.7, pch = 20))
  #cols <- rbPal(10)[as.numeric(cut(vals.gene, breaks = 10))]
  cols.here <- cols[which(md$sample == u)]
  points(cds.correlations.zscore["ESCsV6.5",which(md$sample == u)], 
         cds.correlations.zscore["MEFs",which(md$sample == u)], col = cols.here, cex = 0.7, pch = 20)
}
dev.off()



## differential expression (takes too long to finish)
#c <- fit_models(cds, model_formula_str = "~histone + day")
#write.table(c[,c(1:3,6)], file = "~/Downloads/diff.expr.monocle.histone.day.txt")


## signature scores
md <- pData(cds)
cds.normalized.counts <- normalized_counts(cds)
somatic <- c("VIM", "CDH2", "MGP", "THBS1", "COL1A2", "THY1", "PRRX2", "ACTA2", "ZEB1", "TWIST1")
somatic <- toupper(unique(somatic))
epithelial <- c("EPCAM", "CHD1", "Cldn3", "Cldn4", "Cldn6", "Cldn7", "Ocln", "Krt8", "Krt18", "Crb3")
epithelial <- toupper(unique(epithelial))
pluripotency <- c("Pou5f1", "Nanog", "Sox2", "Sall4", "Zfp42", "Alpl", "Sall4", "Dppa3", "Tfap2c", "Lin28a")
pluripotency <- toupper(unique(pluripotency))

pData(cds)$somatic.score <- colSums(cds.normalized.counts[match(somatic, rownames(cds.normalized.counts)),])
pData(cds)$epithelial.score <- colSums(cds.normalized.counts[match(epithelial, rownames(cds.normalized.counts)),])
pData(cds)$pluripotency.score <- colSums(cds.normalized.counts[match(pluripotency, rownames(cds.normalized.counts)),])

pdf("../../michael-oncohistones/plots/monocle/somatic.score.pdf")
plot_cells(cds, color_cells_by = "somatic.score", show_trajectory_graph = FALSE) +  
  scale_color_gradient(low = "antiquewhite3", high = "brown4")
dev.off()

pdf("../../michael-oncohistones/plots/monocle/epithelial.score.pdf")
plot_cells(cds, color_cells_by = "epithelial.score", show_trajectory_graph = FALSE) +  
  scale_color_gradient(low = "antiquewhite3", high = "darkblue")
dev.off()

pdf("../../michael-oncohistones/plots/monocle/pluripotency.score.pdf")
plot_cells(cds, color_cells_by = "pluripotency.score", show_trajectory_graph = FALSE) +  
  scale_color_gradient(low = "antiquewhite3", high = "darkgreen")
dev.off()


#### load data
library(monocle3)
library(tidyverse)
library(viridis)
cds <- readRDS("../../rdata/michael-monocle.RDS")
#m <- t(t(exprs(cds)) /  pData(cds)[, 'Size_Factor'])
#saveRDS(list("ann.data" = ann.data, "m" = m), file = "../../rdata/michael-monocle-matrices.RDS")
ann.data <- colData(cds)
m <- cds@assays@data$counts
saveRDS(list("ann.data" = ann.data, "m" = m), file = "../../rdata/michael-monocle-matrices-raw.RDS")
