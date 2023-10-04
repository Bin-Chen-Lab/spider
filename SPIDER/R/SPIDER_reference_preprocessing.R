#------------------------------------------
#Reference set preprocessing: 1. Check cell type annotation. 2. Embed transcriptomes.
SPIDER_reference_preprocessing <- function(reference = NULL,
                         use_cell_type = 'SingleR',
                         reference_cell_type = NULL,
                         save_path,
                         SPIDER_model_file_path,
                         scarches_path){
  
  #reference: The Seurat object of a reference CITE-seq dataset provided by the user (If using multiple Merge all your datasets together). Use reference[["batch"]] = ... to specify batch IDs. And if the user choose to use self-defined cell types for the reference set, they must also contain a column "cell_type" in the meta table.
  #reference_cell_type: A meta table with the same format as "query_celltype_SingleR.csv". The row names are cell IDs and there should be one column named 'final_celltype' containing all self-defined cell types for every cell. Otherwise if users set use_cell_type = 'SingleR', they do not need to provide any matrix for the query_cell_type parameter. Default is NULL.
  #------------------------------------------
  #Assign cell types for the reference (Only when use_pretrain = 'F') for scArches-SCANVI embedding.
  cat('Check reference cell type annotation...\n')
  
  if(use_cell_type == 'SingleR'){
    reference_celltype <- get_SingleR_celltype(reference, save_path, specify = 'reference')
  }
  
  if(use_cell_type != 'SingleR'){
    if(sum(colnames(reference) != rownames(reference_cell_type)) == 0){ #examine if the provided self-defined cell type file has the same cell IDs as the query transcriptomes.
      reference_celltype <- reference_cell_type
    }else{
      stop('Error: Provided reference cell type file contains different cell IDs from the reference transcriptomes')
    }
  }
  
  #------------------------------------------

}






