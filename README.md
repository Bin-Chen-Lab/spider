# SPIDER
#First load your query transcriptomes in R studio (Seurat object):<br />
load('.../RNA_ADT_QC_SeuratObject_GSM5025059_20230115.RData') <br />
RNA[['study']] = 'GSM5025059' <br />

library(devtools) <br />
devtools::install_github(repo = 'Bin-Chen-Lab/spider', <br />
                         subdir = '/SPIDER') <br />

#Use SPIDER to predict on your query transcriptomes: <br />
SPIDER_predict (           RNA,  <br />
                           tissue = ...,  <br />
                           disease = ..., <br />
                           SPIDER_model_file_path = ...,<br />
                           use_cell_type = 'SingleR', <br />
                           query_cell_type = NULL,<br />
                           protein = 'All', <br />
                           use_pretrain = 'T',<br />
                           save_path = ...,<br />
                           use_python_path = NULL,<br />
                           scarches_path = NULL,<br />
                           all_trainable_proteins_gene_names_6_training_sets = NULL,<br />
                           file_A_B_C_matching_table = NULL)<br />
