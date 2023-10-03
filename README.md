# spider
#First load your query transcriptomes in R studio (Seurat object):<br />
load('.../RNA_ADT_QC_SeuratObject_GSM5025059_20230115.RData') <br />
RNA[['study']] = 'GSM5025059' <br />

setwd('/SPIDER/') <br />
library(devtools) <br />
load_all(".") <br />

#Use pretrained SPIDER to directly predict on your query transcriptomes: <br />
SPIDER_predict (           RNA,  <br />
                           tissue = ...,  <br />
                           disease = ..., <br />
                           SPIDER_model_file_path = '/home/ubuntu/single_cell/SPIDER/SPIDER_weight/',<br />
                           use_cell_type = 'SingleR', <br />
                           query_cell_type = NULL,<br />
                           protein = 'All', <br />
                           use_pretrain = 'T',<br />
                           save_path = ...,<br />
                           use_python_path = NULL,<br />
                           scarches_path = NULL,<br />
                           all_trainable_proteins_gene_names_6_training_sets = NULL,<br />
                           file_A_B_C_matching_table = NULL)<br />

#-------------------------------------------------------------------------------------------------<br />
#For SPIDER retraining: <br />

#First run: reference <- SPIDER_reference_preprocessing(...)<br />

#Then manually select cells for consisting your reference set (e.g., only using cell types with >=100 cells).<br />

#Then run: get_reference_embeddings(...)<br />

#Then manually prepare file for:<br />
#file_A_B_C_matching_table: A table containing column names: 'file_B', 'file_C', 'tissue_class', 'disease_class', 'cell_type_class'.<br />
#all_protein_list: Names of all proteins in your reference set.<br />

#Then run: SPIDER_train(...)<br />

#Then run: SPIDER_predict(...)<br />

