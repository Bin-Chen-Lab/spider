#' @export
get_gene_coexp <- function(seurat_data, save_path, use_python_path){
  
  library(reticulate)
  
  if(!is.null(use_python_path)){
    use_python(use_python_path, required = T)
  }
  
  pd <- import('pandas', convert = FALSE)
  
  DefaultAssay(seurat_data) = "RNA"
  seurat_data <- NormalizeData(seurat_data, normalization.method = "LogNormalize", scale.factor = 10000)
  gene_coexp <- pd$DataFrame(t(as.matrix(seurat_data@assays$RNA@data)))$corr(method='pearson')
  gene_coexp$index = rownames(seurat_data)
  gene_coexp$columns = rownames(seurat_data)
  gene_coexp$to_csv(paste0(save_path, 'query_gene_coexpression.csv'))
  
  return(gene_coexp)
}

