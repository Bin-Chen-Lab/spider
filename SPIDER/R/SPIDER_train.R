#------------------------------------------
#Train on reference set.
SPIDER_train <- function(SPIDER_model_file_path,
                         file_A_B_C_matching_table,
                         all_protein_list,
                         use_python_path){
  
  #------------------------------------------
  #Train the SPIDER model for seen proteins (Only when use_pretrain = 'F'), save model weights and confidence scores.
  
  cat('Start training...\n')
  
  library(reticulate)
  if(!is.null(use_python_path)){
    use_python(use_python_path, required = T)
  }
  
  setwd(SPIDER_model_file_path)
  setwd('../python/')
  SPIDER <- reticulate::import("SPIDER", convert = F)
  
  SPIDER$train_model$train_model(file_A_B_C_matching_table, all_protein_list, SPIDER_model_file_path)

}





