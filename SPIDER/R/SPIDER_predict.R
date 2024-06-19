#' Surface protein abundance prediction
#'
#' Predict surface protein abundance on provided query transcriptomes.
#' @param seurat_data Query transcriptomes after prepocessing. Including Seurat log normalization, clustering and umap reductions. You must also use seurat_data[["study"]] = ... to specify the batch IDs for all cells. 
#' @param tissue The tissue source of the query transcriptomes. Default options include: 'bone marrow', 'brain', 'blood', 'pleura', 'peritoneum'.
#' @param disease The disease state of the query transcriptomes. Default options include: 'healthy', 'mesothelioma', 'glioblastoma', 'leukemia'.
#' @param SPIDER_model_file_path The path to the saved pretrained SPIDER model weights.
#' @param use_cell_type The method name of cell type annotations. Default option is 'SingleR'. If you intend to use the pre-trained SPIDER (use_pretrain = 'T'), you should set use_cell_type to default.
#' @param query_cell_type A data frame specifying cell types for the query dataset. The row names should be cell IDs and there should be one column named 'final_celltype' specifying your defined cell type for every cell. This parameter is only needed if users do not set use_cell_type = 'SingleR'. Otherwise if users set use_cell_type = 'SingleR', they do not need to provide any data frame for the query_cell_type parameter. Default is NULL.
#' @param protein A vector specifying corresponding gene names of proteins to be predicted. Default option is 'All', where all surface proteins listed in UniProt will be predicted.
#' @param use_pretrain Whether to use pretrained SPIDER model weights or not. If yes, set use_pretrain = 'T', otherwise set use_pretrain = 'F'.
#' @param save_path The path to the directory where you want to save your prediction results.
#' @param use_python_path The path to the specific version of python for R reticulate. This parameter is only needed if you use a separate version of python for R reticulate from your default python configuration for reticulate. It will automatically pass this parameter to reticulate's "use_python" function. Otherwise just set this parameter to NULL.
#' @param scarches_path The path to the directory where the scArches package is downloaded.
#' @param all_trainable_proteins_gene_names A data frame matching all trainable proteins' (i.e., seen proteins in your reference set) names to their corresponding gene names. Columns must contain column names: 'consistent_protein_name' specifying the protein names, and 'gene_name' specifying the corresponding gene names. This parameter is only needed if users set use_pretrain = 'F', otherwise just set this parameter to NULL.
#' @param file_A_B_C_matching_table A data frame specifying the path to the reference data split by their context factors. Must containing column names: 'file_B' specifying all file paths of embedded reference transcriptomes, 'file_C' specifying all corresponding file paths of reference CLR normalized ADT counts, 'tissue_class' specifying corresponding reference tissue classes using numeric numbers (0, 1, 2, ...), 'disease_class' specifying corresponding reference disease classes using numeric numbers (0, 1, 2, ...), 'cell_type_class' specifying corresponding reference cell type classes using numeric numbers (0, 1, 2, ...). Make sure that each row only contain one cell type, one disease and one tissue correctly describing the states of cells in the corresponding transcriptome and ADT files. This parameter is only needed if users set use_pretrain = 'F', otherwise just set this parameter to NULL. 


#' @return Predicted surface protein abundance and confidence scores will be saved to your indicated directory save_path.


SPIDER_predict <- function(seurat_data, 
                           tissue = 'pancreas', 
                           disease = 'healthy', 
                           SPIDER_model_file_path,
                   use_cell_type = 'SingleR', 
                   query_cell_type = NULL,
                   protein = 'All', 
                   use_pretrain = 'T',
                   save_path,
                   use_python_path,
                   scarches_path,
                   all_trainable_proteins_gene_names = NULL,
                   file_A_B_C_matching_table = NULL,
                          n_ensemble_members = 8){ 

  #seurat_data: Query transcriptomes after prepocessing. Including Seurat normalization, clustering and umap reductions. You should use seurat_data[["study"]] = ... to specify the batch IDs. 
  #SPIDER_model_file_path: The path to the saved SPIDER model weights and one-hot encoded features files.
  #tissue: 'bone marrow', 'brain', 'blood', 'pleura', 'peritoneum'
  #disease: 'healthy', 'mesothelioma', 'glioblastoma', 'leukemia'
  #use_cell_type: 'SingleR', self-defined_cell_type
  #query_cell_type: A meta table with the same format as "query_celltype_SingleR.csv". The row names are cell IDs and there should be one column named 'final_celltype' containing all self-defined cell types for every cell. Otherwise if users set use_cell_type = 'SingleR', they do not need to provide any matrix for the query_cell_type parameter. Default is NULL.
  #protein: Corresponding gene names of proteins to be predicted. 'All', self-defined_protein_list.
  #use_python_path: Defined path to the reticulate python. Default is NULL.
  #all_trainable_proteins_gene_names: A table matching all seen proteins' names to corresponding gene names. columns: 'gene_name', 'consistent_protein_name'
  #file_A_B_C_matching_table: A table containing column names: 'file_B', 'file_C', 'tissue_class', 'disease_class', 'cell_type_class'.
  #------------------------------------------
  library(Seurat)
  library(dplyr)
  
  all_trainable_proteins_gene_names_6_training_sets = all_trainable_proteins_gene_names 
  #------------------------------------------
  #Assign cell types for the query dataset for later SPIDER prediction.
  cat('Check query cell type annotation...\n')
  
  if(use_cell_type == 'SingleR'){
    query_celltype <- get_SingleR_celltype(seurat_data, save_path, specify = 'query')
  }
  
  if(use_pretrain == 'T' & use_cell_type != 'SingleR'){
    stop('Error: Please set use_cell_type to "SingleR" when using pretrained model')
  }
  
  if(use_pretrain == 'F' & use_cell_type != 'SingleR'){
    if(sum(colnames(seurat_data) != rownames(query_cell_type)) == 0){ #examine if the provided self-defined cell type file has the same cell IDs as the query transcriptomes.
      query_celltype <- query_cell_type
    }else{
      stop('Error: Provided query cell type file contains different cell IDs from the query transcriptomes')
    }
  }
  
  #------------------------------------------
  #Use scArches-SCANVI to embed the reference (Only when use_pretrain = 'F')/query transcriptomes into 128 dimensions.
  
  cat('Start query data embedding...\n')
  
  query_embeddings <- as.matrix(get_query_embeddings(seurat_data, save_path, use_pretrain, SPIDER_model_file_path, scarches_path, use_python_path))
  
  #------------------------------------------
  #Generate protein-protein similarity for the query dataset.
  
  cat('Prepare unseen...\n')
  
  query_gene_coexp <- get_gene_coexp(seurat_data, save_path)

  #------------------------------------------
  #Run SPIDER model to predict seen proteins on the query dataset, and save confidence scores.
  cat('Predict seen proteins...\n')
  
  if(use_pretrain == 'T'){
    #training_set_cell_type_class = read.csv(paste0(SPIDER_model_file_path, 'celltype_class_combined_training_6_datasets.csv'), stringsAsFactors = F, check.names = F, row.names = 'cell_type')
    training_epoch = read.csv(paste0(SPIDER_model_file_path, 'training_epoch_289_proteins_cell_features_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_DNN_SCANVI_128dim_20230115.csv'), stringsAsFactors = F, check.names = F, row.names = 'protein_name')
    all_trainable_proteins_gene_names_6_training_sets = read.csv(paste0(SPIDER_model_file_path, 'protein_gene_names_union_289_DNNs_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_0.6_20230115.csv'), stringsAsFactors = F, check.names = F, row.names = 1)
    file_A_B_C_matching_table = read.csv(paste0(SPIDER_model_file_path, 'file_B_C_multitask_DNN_celltype_tissue_disease_cell_features_seen_proteins_all_cell_matching_table_combined_training_6_datasets_SCANVI_128dim_20230115.csv'), stringsAsFactors = F, check.names = F, row.names = 1)
  }else{
    #training_set_cell_type_class = ...
    training_epoch = read.csv(paste0(SPIDER_model_file_path, 'retrain_record_epochs.csv'), stringsAsFactors = F, check.names = F, row.names = 'protein_name')
    all_trainable_proteins_gene_names_6_training_sets = all_trainable_proteins_gene_names_6_training_sets
    if(is.null(all_trainable_proteins_gene_names_6_training_sets)){
      stop('Error: Please input a table for all_trainable_proteins_gene_names when not using pretrained model')
    }
    file_A_B_C_matching_table = file_A_B_C_matching_table
    if(is.null(file_A_B_C_matching_table)){
      stop('Error: Please input a table for file_A_B_C_matching_table when not using pretrained model')
    }
  }
  
  if(protein == 'All'){
    if(use_pretrain == 'T'){
    all_protein_list = read.csv(paste0(SPIDER_model_file_path, 'union_trainable_protein_names_from_combined_6_training_sets_20230115.csv'), stringsAsFactors = F, row.names = 1, check.names = F)$x
    }
    if(use_pretrain == 'F'){
      all_protein_list = all_trainable_proteins_gene_names_6_training_sets$consistent_protein_name
    }
  }else{
    protein2 <- intersect(all_trainable_proteins_gene_names_6_training_sets$gene_name, protein)
    if(length(protein2) > 0){
      all_protein_list <- NULL
      for (p in protein2) {
        all_protein_list =  c(all_protein_list, filter(all_trainable_proteins_gene_names_6_training_sets, gene_name == p)$consistent_protein_name)
      }
    }else{
      all_protein_list <- NULL
    }
  }
  
  library(reticulate)
  if(!is.null(use_python_path)){
    use_python(use_python_path, required = T)
  }
  
  setwd(SPIDER_model_file_path)
  setwd('../python/')
  SPIDER <- reticulate::import("SPIDER", convert = F)
  
  if(length(all_protein_list) == 1){
    all_protein_list = list(all_protein_list)
  }
    
  if(length(all_protein_list) > 0){
  SPIDER$predict_seen$predict_seen(tissue, 
                                   disease, 
                                   save_path, 
                                   all_protein_list, 
                                   SPIDER_model_file_path, 
                                   training_epoch, 
                                   use_cell_type, 
                                   query_cell_type, 
                                   use_pretrain,
                                   file_A_B_C_matching_table)
  }
  #------------------------------------------
  #Run SPIDER model to predict unseen proteins on the query dataset, and save confidence scores.
  
  cat('Predict unseen proteins...\n')
  
  if(use_cell_type == 'SingleR'){
    cell_type_file = paste0(save_path, 'query_celltype_SingleR.csv')
  }else{
    cell_type_file = query_cell_type
  }
  
  if(use_pretrain == 'T'){
    all_trainable_proteins_gene_names_6_training_sets = read.csv(paste0(SPIDER_model_file_path, 'protein_gene_names_union_289_DNNs_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_0.6_20230115.csv'), stringsAsFactors = F, check.names = F, row.names = 1)
    training_epoch = read.csv(paste0(SPIDER_model_file_path, 'training_epoch_289_proteins_cell_features_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_DNN_SCANVI_128dim_20230115.csv'), stringsAsFactors = F, check.names = F, row.names = 'protein_name')
    file_A_B_C_matching_table = read.csv(paste0(SPIDER_model_file_path, 'file_B_C_multitask_DNN_celltype_tissue_disease_cell_features_seen_proteins_all_cell_matching_table_combined_training_6_datasets_SCANVI_128dim_20230115.csv'), stringsAsFactors = F, check.names = F, row.names = 1)
  }else{
    all_trainable_proteins_gene_names_6_training_sets = all_trainable_proteins_gene_names_6_training_sets
    training_epoch = read.csv(paste0(SPIDER_model_file_path, 'retrain_record_epochs.csv'), stringsAsFactors = F, check.names = F, row.names = 'protein_name')
    if(is.null(all_trainable_proteins_gene_names_6_training_sets)){
      stop('Error: Please input a table for all_trainable_proteins_gene_names when not using pretrained model')
    }
    file_A_B_C_matching_table = file_A_B_C_matching_table
    if(is.null(file_A_B_C_matching_table)){
      stop('Error: Please input a table for file_A_B_C_matching_table when not using pretrained model')
    }
  }
  
  if(protein == 'All'){
    all_test_proteins_gene_names = read.csv(paste0(SPIDER_model_file_path, 'unseen_file/uniprot_cell_membrane_protein_gene_list.csv'), stringsAsFactors = F, check.names = F, row.names = 1)$surface_protein_gene
    }else{
    all_test_proteins_gene_names = protein #gene names
  }
  
  test_gene_coexp_matrix <- read.csv(paste0(save_path, 'query_gene_coexpression.csv'), stringsAsFactors = F, row.names = 1, check.names = F)
  
  library(reticulate)
  if(!is.null(use_python_path)){
    use_python(use_python_path, required = T)
  }
  
  setwd(SPIDER_model_file_path)
  setwd('../python/')
  SPIDER <- reticulate::import("SPIDER", convert = F)
  
  if(length(all_test_proteins_gene_names) == 1){
    all_test_proteins_gene_names = list(all_test_proteins_gene_names)
  }
  
  if(length(intersect(all_test_proteins_gene_names, rownames(test_gene_coexp_matrix)) > 0)){
  SPIDER$predict_unseen$predict_unseen(tissue, 
                                       disease, 
                                       save_path, 
                                       SPIDER_model_file_path, 
                                       use_pretrain, 
                                       cell_type_file, 
                                       all_trainable_proteins_gene_names_6_training_sets, 
                                       training_epoch, 
                                       all_test_proteins_gene_names,
                                       file_A_B_C_matching_table,
                                      as.integer(n_ensemble_members))
  }
  
  #-------------------------------------------
  cat('Finished!\n')
}




