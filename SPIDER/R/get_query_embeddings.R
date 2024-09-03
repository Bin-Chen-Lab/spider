#' Generating query transcriptome embeddings
#'
#' Generating query transcriptome embeddings via scArches-SCANVI.
#' @param seurat_data Query transcriptomes after prepocessing. Including Seurat log normalization, clustering and umap reductions. You must also use seurat_data[["study"]] = ... to specify the batch IDs for all cells. 
#' @param save_path The path to the directory where you want to save your generated query transcriptome embeddings. This should be the same directory as where you later save your protein abundance prediction results.
#' @param use_pretrain Whether to use pretrained SPIDER model weights or not. If yes, set use_pretrain = 'T', otherwise set use_pretrain = 'F'.
#' @param SPIDER_model_file_path The path to the saved pretrained SPIDER model weights.
#' @param scarches_path The path to the directory where the scArches package is downloaded.
#' @param use_python_path The path to the specific version of python for R reticulate. This parameter is only needed if you use a separate version of python for R reticulate from your default python configuration for reticulate. It will automatically pass this parameter to reticulate's "use_python" function. Otherwise just set this parameter to NULL.

#' @return Generated query transcriptome embeddings will be saved to your indicated directory save_path.

#' @export
get_query_embeddings <- function(seurat_data, save_path, use_pretrain, SPIDER_model_file_path, scarches_path, use_python_path){
  
  library(Seurat)
  library(dplyr)
  library(anndata)
  
  library(reticulate)
  if(!is.null(use_python_path)){
  use_python(use_python_path, required = T)
  }
  sc <- import('scanpy', convert = FALSE)
  scvi <- import('scvi', convert = FALSE)
  if(!is.null(scarches_path)){
    setwd(scarches_path)
  }
  sca <- import('scarches', convert = FALSE)
  torch <- import('torch', convert = FALSE)
  remove_sparsity <- import('scarches.dataset.trvae.data_handling', convert = FALSE)
  plt <- import('matplotlib.pyplot', convert = FALSE)
  np <- import('numpy', convert = FALSE)
  os <- import('os', convert = FALSE)
  
  DefaultAssay(seurat_data) = "RNA"
  RNA <- seurat_data #The query transcriptomes

  setwd("../")
  if(use_pretrain == 'T'){
    save_top_1000_HVGs = read.csv(paste0(SPIDER_model_file_path, 'training_combined_6_datasets_RNA_SCANVI_latent_128dim_20230115/top_1000_HVGs_training_combined_6_datasets_RNA_SCANVI_latent_128dim_20230115.csv'), stringsAsFactors = F, row.names = 1)$x
  }
  if(use_pretrain == 'F'){
    save_top_1000_HVGs = read.csv(paste0(SPIDER_model_file_path, 'reference_SCANVI_top_HVGs.csv'), stringsAsFactors = F, row.names = 1)$x
  }
  setwd("scarches-0.4.0/")

  if(length(setdiff(save_top_1000_HVGs, rownames(RNA))) != 0){
    non_exist_genes <- setdiff(save_top_1000_HVGs, rownames(RNA)) 
    non_exist_gene_exp <- as.data.frame(matrix(0, nrow = length(non_exist_genes), ncol = ncol(RNA)))
    rownames(non_exist_gene_exp) = non_exist_genes
    colnames(non_exist_gene_exp) = colnames(RNA)
    RNA <- as.sparse(rbind(as.data.frame(RNA@assays$RNA@counts), non_exist_gene_exp))
    RNA <- CollapseSpeciesExpressionMatrix(RNA)
    RNA <- CreateSeuratObject(counts = RNA)
  }
  
  RNA[['study']] = seurat_data[['study']]
  RNA[['cell_type']] = 'Unknown'
  
  RNA <- NormalizeData(RNA, normalization.method = "LogNormalize", scale.factor = 10000)
  
  RNA <- RNA[save_top_1000_HVGs] #If there are HVGs that are not existed in the query data set's genes, use zeros to impute.
  print(RNA)
  
  scvi$settings$progress_bar_style = 'tqdm'
  adata <- sc$AnnData(
    X   = t(as.matrix(GetAssayData(RNA,slot='counts',assay = "RNA"))), #scVI requires raw counts
    obs = RNA[[]],
    var = GetAssay(RNA)[[]]
  )
  print(adata)
  
  #SCANVI model setting:
  sc$settings$set_figure_params(dpi=as.integer(200), frameon=FALSE)
  sc$set_figure_params(dpi=as.integer(200))
  sc$set_figure_params(figsize=c(4, 4))
  torch$set_printoptions(precision=as.integer(3), sci_mode=FALSE, edgeitems=as.integer(7))
  
  condition_key = 'study'
  cell_type_key = 'cell_type'

  setwd("../")
  if(use_pretrain == 'T'){
  ref_model_path = paste0(SPIDER_model_file_path, 'training_combined_6_datasets_RNA_SCANVI_latent_128dim_20230115/')
  }
  if(use_pretrain == 'F'){
    ref_model_path = paste0(SPIDER_model_file_path)
  }
  setwd("scarches-0.4.0/")
  
  #Perform surgery on reference model and train on query dataset without cell type labels
  model = sca$models$SCANVI$load_query_data(
    adata,
    ref_model_path,
    freeze_dropout = TRUE
  )
  
  model$`_unlabeled_indices` = np$arange(adata$n_obs)
  model$`_labeled_indices` = list()
  
  print(paste0("Labelled Indices: ", length(model$`_labeled_indices`)))
  print(paste0("Unlabelled Indices: ", length(model$`_unlabeled_indices`)))
  
  model$train(
    max_epochs=as.integer(100),
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=as.integer(10)
  )
  
  query_latent = sc$AnnData(model$get_latent_representation())
  query_latent$obs$predictions = model$predict()

  query_latent = model$get_latent_representation()
  query_latent = as.matrix(query_latent)
  rownames(query_latent) = colnames(RNA)
  setwd("../")
  write.csv(query_latent, paste0(save_path, 'query_embeddings.csv'))
  
  return(query_latent)
  
}
