#generate reference embeddings using scArches-SCANVI:
get_reference_embeddings <- function(reference, 
                                     SPIDER_model_file_path, 
                                     scarches_path,
                                     use_cell_type = 'SingleR',
                                     use_python_path,
                                     reference_celltype_save_path,
                                     reference_cell_type = NULL){
  
  #reference_cell_type: A meta table with the same format as "query_celltype_SingleR.csv". The row names are cell IDs and there should be one column named 'final_celltype' containing all self-defined cell types for every cell. Otherwise if users set use_cell_type = 'SingleR', they do not need to provide any matrix for the query_cell_type parameter. Default is NULL.
  #reference_celltype_save_path: The path where you previously save your reference SingleR cell types. Use NULL if you didn't set use_cell_type='SingleR'
  #-------------------------------------------------------------------------
  cat('Start reference data embedding...\n')
  
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
  scvi$settings$progress_bar_style = 'tqdm'
  
  DefaultAssay(reference) = "RNA"
  
  if(use_cell_type == 'SingleR'){
    reference[['cell_type']] = read.csv(paste0(reference_celltype_save_path, 'reference_celltype_SingleR.csv'), stringsAsFactors = F, row.names = 1, check.names = F)$final_celltype
  }else{
    reference[['cell_type']] = reference_cell_type
  }
  #---------------------------------------------------------------
  RNA <- NormalizeData(reference, normalization.method = "LogNormalize", scale.factor = 10000)
  RNA <- FindVariableFeatures(RNA, selection.method = "vst", assay = "RNA", nfeatures = 1000)
  top1000 <- head(VariableFeatures(RNA, assay = "RNA"), 1000)
  RNA <- RNA[top1000]
  print(RNA)
  
  
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
  
  #Create SCANVI model and train it on fully labelled reference dataset:
  sca$dataset$setup_anndata(adata, batch_key=condition_key, labels_key=cell_type_key)
  
  vae = sca$models$SCVI(
    adata,
    n_layers=as.integer(2),
    n_latent=as.integer(128),
    encode_covariates=TRUE,
    deeply_inject_covariates=FALSE,
    use_layer_norm="both",
    use_batch_norm="none")
  
  vae$train()
  
  scanvae = sca$models$SCANVI$from_scvi_model(vae, "Unknown")
  print(paste0("Labelled Indices: ", length(scanvae$`_labeled_indices`)))
  print(paste0("Unlabelled Indices: ", length(scanvae$`_unlabeled_indices`)))
  
  scanvae$train(max_epochs=as.integer(20))
  
  #Getting the latent representation and visualization
  reference_latent = scanvae$get_latent_representation()
  
  reference_latent = as.matrix(reference_latent)
  rownames(reference_latent) = colnames(RNA)
  
  setwd(SPIDER_model_file_path)
  
  scanvae$save(SPIDER_model_file_path, overwrite=TRUE)
  
  write.csv(reference_latent, 'reference_embeddings.csv')
  
  write.csv(as.character(py_to_r(adata$var$index)), 'reference_SCANVI_top1000_HVGs.csv')
  
  return(reference_latent)
  
}


