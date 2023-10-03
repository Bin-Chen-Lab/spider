get_gene_coexp <- function(seurat_data, save_path){
  
  DefaultAssay(seurat_data) = "RNA"
  seurat_data <- NormalizeData(seurat_data, normalization.method = "LogNormalize", scale.factor = 10000)
  gene_coexp <- cor(t(as.matrix(seurat_data@assays$RNA@data)), method = 'pearson')
  write.csv(gene_coexp, paste0(save_path, 'query_gene_coexpression.csv'))
  
  return(gene_coexp)
}