# spider
#First load your query transcriptomes in R studio (Seurat object):
load('.../RNA_ADT_QC_SeuratObject_GSM5025059_20230115.RData')
RNA[['study']] = 'GSM5025059' 

setwd('/SPIDER/')
library(devtools)
load_all(".")

#Use pretrained SPIDER to directly predict on your query transcriptomes:
SPIDER_predict (           RNA, 
                           tissue = ..., 
                           disease = ..., 
                           SPIDER_model_file_path = '/home/ubuntu/single_cell/SPIDER/SPIDER_weight/',
                           use_cell_type = 'SingleR', 
                           query_cell_type = NULL,
                           protein = 'All', 
                           use_pretrain = 'T',
                           save_path = ...,
                           use_python_path = NULL,
                           scarches_path = NULL,
                           all_trainable_proteins_gene_names_6_training_sets = NULL,
                           file_A_B_C_matching_table = NULL)

#-------------------------------------------------------------------------------------------------
#For SPIDER retraining: 

#First run: reference <- SPIDER_reference_preprocessing(...)

#Then manually select cells for consisting your reference set (e.g., only using cell types with >=100 cells).

#Then run: get_reference_embeddings(...)

#Then manually prepare file for:
#file_A_B_C_matching_table: A table containing column names: 'file_B', 'file_C', 'tissue_class', 'disease_class', 'cell_type_class'.
#all_protein_list: Names of all proteins in your reference set.

#Then run: SPIDER_train(...)

#Then run: SPIDER_predict(...)

