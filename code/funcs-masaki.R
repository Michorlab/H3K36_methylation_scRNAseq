## add log fold change from the SCT data slot for the genes resulting from differential expression of SCT scale.data slot
get_diff_expr_seurat <- function(seurat.obj, ident.1, ident.2, logfc.thresh = 0.25, p.val.thresh = 0.1){
  
  genes.scaled <- FindMarkers(seurat.obj, ident.1 = ident.1, ident.2 = ident.2, assay = "SCT", slot = "scale.data", verbose = TRUE)
  genes.scaled <- cbind("gene" = rownames(genes.scaled), genes.scaled)
  
  genes.data <- FindMarkers(seurat.obj, ident.1 = ident.1, ident.2 = ident.2, assay = "SCT", slot = "data", verbose = TRUE, logfc.threshold = logfc.thresh)
  genes.data <- cbind("gene" = rownames(genes.data), genes.data)
  
  genes.scaled$avg_logFC <- genes.data$avg_log2FC[match(genes.scaled$gene, genes.data$gene)]
  genes.scaled <- genes.scaled %>% arrange(desc(abs(avg_logFC)))
  genes.scaled <- genes.scaled[which(genes.scaled$p_val_adj < p.val.thresh),]
  genes.scaled$avg_logFC[which(is.na(genes.scaled$avg_logFC))] <- 0
  
  #return(list("genes.scaled" = genes.scaled, "genes.data" = genes.data))
  return(genes.scaled)
}