# SPIDER
SPIDER (surface protein prediction using deep ensembles from single-cell RNA-seq) is a context-agnostic zero-shot deep ensemble model, which enables the large-scale prediction of cell surface protein abundance. Users can install and implement SPIDER with R.

# SPIDER intallation
library(devtools) <br />
devtools::install_github(repo = 'Bin-Chen-Lab/spider', subdir = '/SPIDER') <br />

# SPIDER usage with sample data
#Make sure you have the scArches downloaded first for python, and reticulate downloaded for R. <br />
#First download the "sample_query.RData" from our github. <br />

load("sample_query.RData") <br />

library(SPIDER) <br />

#Use SPIDER to predict on your query transcriptomes: <br />
SPIDER_predict (           seurat_data = RNA  <br />
                           tissue = 'pancreas',  <br />
                           disease = 'healthy', <br />
                           SPIDER_model_file_path = paste0(.libPaths(), '/SPIDER/SPIDER_weight/'),<br />
                           use_cell_type = 'SingleR', <br />
                           query_cell_type = NULL,<br />
                           protein = 'All', <br />
                           use_pretrain = 'T',<br />
                           save_path = ..., #enter your path where you want to save prediction results <br />
                           use_python_path = NULL, #If you're using a specific python path for reticulate, enter this path. <br />
                           scarches_path = NULL, #If you're using a specific path for scArches, enter this path. <br />
                           all_trainable_proteins_gene_names_6_training_sets = NULL,<br />
                           file_A_B_C_matching_table = NULL)<br />
