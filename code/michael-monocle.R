## Single cell analysis on H3K36 methylation data: monocle3 analysis
## Author: Simona Cristea
## Date: November 2022
## Publication: H3K36 methylation maintains cell identity by regulating opposing lineage programs

source(here::here("libraries.R"))



############################################################################
############################################################################
###                                                                      ###
###                         READ IN CLEANED DATA                         ###
###                                                                      ###
############################################################################
############################################################################
seuratX <- readRDS(file = "../data/rdata/michael-initial.RDS")



############################################################################
############################################################################
###                                                                      ###
###                  EXTRACT RNA COUNT DATA FROM SEURAT                  ###
###                                                                      ###
############################################################################
############################################################################
expression_matrix <- as(as.matrix(seuratX@assays$RNA@counts), 'sparseMatrix')
cell_metadata <- seuratX@meta.data
cell_metadata <- cbind("barcode" =rownames(cell_metadata), cell_metadata)
gene_metadata <- data.frame(gene_short_name = row.names(expression_matrix), row.names = row.names(expression_matrix))



###########################################################################
###########################################################################
###                                                                     ###
###                        CONSTRUCT MONOCLE CDS                        ###
###                                                                     ###
###########################################################################
###########################################################################
cds <- new_cell_data_set(expression_data = expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
rm(expression_matrix)

## correct label
cds@colData$sample[which(cds@colData$sample == "d9")] <- "iPS"


## view data
pData(cds)
fData(cds)



############################################################################
############################################################################
###                                                                      ###
###                       PREPROCESS USING MONOCLE                       ###
###                                                                      ###
############################################################################
############################################################################
cds <- preprocess_cds(cds, num_dim = 100)
pdf("../plots/monocle/variance.explained.pdf")
plot_pc_variance_explained(cds)
dev.off()


## reduce dimension and plot
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
pdf("../plots/monocle/umap-sample.pdf")
plot_cells(cds, color_cells_by = "sample", label_cell_groups = FALSE, group_label_size = 5)
dev.off()
pdf("../plots/monocle/umap-histone.pdf")
plot_cells(cds, color_cells_by = "histone", label_cell_groups = TRUE, group_label_size = 5)
dev.off()
pdf("../plots/monocle/umap-percent.mt.pdf")
plot_cells(cds, color_cells_by = "percent.mt", label_cell_groups = FALSE)
dev.off()
pdf("../plots/monocle/umap-detected.pdf")
plot_cells(cds, color_cells_by = "detected", label_cell_groups = FALSE)
dev.off()
pdf("../plots/monocle/umap-sum.pdf")
plot_cells(cds, color_cells_by = "sum", label_cell_groups = FALSE)
dev.off()



############################################################################################
############################################################################################
###                                                                                      ###
###  ROUND 1: REMOVE LOW GENES AND CONTAMINATION (CLUSTERS 9,10,11: 372 CELLS) AND REDO  ###
###                                                                                      ###
############################################################################################
############################################################################################
# uncomment next 3 lines only for the first normalization & preprocessing run
#idx.to.remove <- which(cds@clusters$UMAP$clusters %in% c(9,10,11))
#saveRDS(idx.to.remove, file = "../data/rdata/michael-idx.to.remove.RDS")
#cds <- cds[, -idx.to.remove]



######################################################################################
######################################################################################
###                                                                                ###
###  ROUND 2: REMOVE LOW GENES AND CONTAMINATION (CLUSTER 10: 903 CELLS) AND REDO  ###
###                                                                                ###
######################################################################################
######################################################################################
# uncomment next 3 lines only for the second normalization & preprocessing run
## remove low genes clusters (new cluster 10): 903 cells
#idx.to.remove.2 <- which(cds@clusters$UMAP$clusters %in% c(10))
#saveRDS(idx.to.remove.2, file = "../data/rdata/michael-idx.to.remove.2.RDS")
#cds <- cds[, -idx.to.remove.2]



############################################################################
############################################################################
###                                                                      ###
###                        CLUSTER CELLS AND PLOT                        ###
###                                                                      ###
############################################################################
############################################################################
cds <- cluster_cells(cds, resolution=1e-5)

pdf("../plots/monocle/umap-cluster.pdf")
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = TRUE)
dev.off()
pdf("../plots/monocle/umap-partition.pdf")
plot_cells(cds, color_cells_by = "partition", label_cell_groups = TRUE)
dev.off()

head(cds@clusters$UMAP$partitions)
head(cds@clusters$UMAP$clusters)


pd <- pData(cds.regressed)
pd$partitions <- cds.regressed@clusters$UMAP$partitions
pd <- as_tibble(pd)

pd <- pData(cds)
pd$partitions <- cds@clusters$UMAP$partitions
pd$clusters <- cds@clusters$UMAP$clusters
pd <- as_tibble(pd)

pdf("../plots/monocle/density_percentmt-partitions.pdf")
ggplot(pd, aes (x = percent.mt)) + geom_density() + facet_wrap(~partitions)
dev.off()
pdf("../plots/monocle/density_detected-partitions.pdf")
ggplot(pd, aes (x = detected)) + geom_density() + facet_wrap(~partitions)
dev.off()

pdf("../plots/monocle/density_percentmt-clusters.pdf")
ggplot(pd, aes (x = percent.mt)) + geom_density() + facet_wrap(~clusters)
dev.off()
pdf("../plots/monocle/density_detected-clusters.pdf")
ggplot(pd, aes (x = detected)) + geom_density() + facet_wrap(~clusters)
dev.off()


## plot markers
genes_michael <- readxl::read_excel("../tables/Genes_for_scRNAseq.xlsx")
genes_michael <- as.data.frame(genes_michael)
genes_michael <- toupper(genes_michael$gene)
for (g in genes_michael){
  pdf(paste0("../plots/monocle/markers/genes-", g, ".pdf"), width = 8)
  print(plot_cells(cds, genes=g, color_cells_by = "sample"))
  dev.off()
}



############################################################################
############################################################################
###                                                                      ###
###                         LEARN GRAPH AND PLOT                         ###
###                                                                      ###
############################################################################
############################################################################
# each partition is its own trajectory
cds <- learn_graph(cds)
pdf("../plots/monocle/umap-clusters-trajectory.pdf")
plot_cells(cds, group_label_size = 4)
dev.off()


## detailed trajectory plots
# leaves (or gray circles with numbers in them) correspond to different outcomes (i.e. cell fate) of the trajectory
# branch nodes (black circles) are points in which cells can travel to one of several outcomes
pdf("../plots/monocle/umap-trajectory.pdf")
plot_cells(cds, color_cells_by = "sample", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, group_label_size = 4)
dev.off()

pdf("../plots/monocle/umap-trajectory-detailed.pdf")
plot_cells(cds, color_cells_by = "sample", label_groups_by_cluster=FALSE, label_leaves = TRUE, label_branch_points=TRUE, group_label_size = 4)
dev.off()



###########################################################################
###########################################################################
###                                                                     ###
###                        PROJECT ON BULK ATLAS                        ###
###                                                                     ###
###########################################################################
###########################################################################
#### read in atlas data
atlas <- read.table("../tables/final_table_for_supp.tsv", header = TRUE, sep = "\t")
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
})
rownames(cds.correlations) <- c("ESCsV6.5", "MEFs")
pData(cds)$ESCsV6.5 <- cds.correlations["ESCsV6.5",]
pData(cds)$MEFs <- cds.correlations["MEFs",]


## do z-score
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
pdf(paste0("../plots/projections/cor.atlas.allgenes.zscore-enh-", gene, ".pdf"), width = 12, height = 9.5)
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


unique.samples <- unique(md$sample)
pdf("../plots/projections/cor.atlas.allgenes.pdf")
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


###########################################################################
###########################################################################
###                                                                     ###
###                        HETEROGENEITY BOXPLOT                        ###
###                                                                     ###
###########################################################################
###########################################################################
um <- reducedDims(cds)$UMA
md <- pData(cds)

## separate per samples
unique.samples <- c("d0", "H3-d2", "H3-d4", "H3-d6", "H3-d8", "K36M-d2", "K36M-d4", "K36M-d6", "K36M-d8", "iPS")
distance.mats <- list()
ms <- list()
for (u in unique.samples){
  print(u)
  idx <- which(md$sample == u)
  um.now <- um[idx,]
  distance.mats[[length(distance.mats)+1]] <- dist(um.now)
  distance.mats[[length(distance.mats)]] <- as.matrix(distance.mats[[length(distance.mats)]])
  ms[[length(distance.mats)]] <- rowMeans(distance.mats[[length(distance.mats)]])
  ms[[length(distance.mats)]] <- cbind(ms[[length(distance.mats)]], rep(u, length(ms[[length(distance.mats)]])))
}

ms.all.points <- do.call(rbind, ms)
ms.all.points <- as_data_frame(ms.all.points)
colnames(ms.all.points) <- c("average_distance", "sample")
ms.all.points$average_distance <- as.numeric(ms.all.points$average_distance)
ms.all.points$histone <- ms.all.points$sample
ms.all.points$histone <- dplyr::recode(ms.all.points$histone,
                                       "H3-d2" = "H3", "H3-d4" = "H3", "H3-d6" = "H3", "H3-d8" = "H3",
                                       "K36M-d2" = "K36M", "K36M-d4" = "K36M", "K36M-d6" = "K36M", 
                                       "K36M-d8" = "K36M")
ms.all.points$sample <- factor(ms.all.points$sample, levels = c("d0", "H3-d2", "K36M-d2", "H3-d4", "K36M-d4", "H3-d6", "K36M-d6", "H3-d8", "K36M-d8", "iPS"))


ms.all.points.short.short <- ms.all.points %>% filter(average_distance <= 7.5)
pdf("../plots/monocle/avg_distance_no_outliers.2.pdf")
ggplot(ms.all.points.short.short, aes (x = sample, y = average_distance)) + geom_boxplot(na.rm = TRUE, notch = TRUE, outlier.size = 0.2, 
                                                                                         aes(fill = histone)) + theme_minimal()
dev.off()

saveRDS(cds, file = "../data/rdata/michael-monocle.RDS")
