# SPIDER
SPIDER (surface protein prediction using deep ensembles from single-cell RNA-seq) is a context-agnostic zero-shot deep ensemble model, which enables the large-scale prediction of cell surface protein abundance. 

# Installation
Before installing SPIDER, make sure you have installed the dependent R and python packages. To ensure this, first run the R lines as in: <br />
https://github.com/Bin-Chen-Lab/spider/blob/fbcd525b52eb66a72b6257946bb30d1dff46737e/dependencies/R%20package%20dependencies.R <br />
Then run the commands in terminal as in: <br />
https://github.com/Bin-Chen-Lab/spider/blob/fbcd525b52eb66a72b6257946bb30d1dff46737e/dependencies/python%20package%20dependencies <br />

Install the SPIDER package: <br />
Open R studio, type the following lines: <br />
```
devtools::install_github(repo = 'Bin-Chen-Lab/spider', subdir = '/SPIDER')
``` 

Then in the terminal, type the following commands:
```
cd '/user_path' #enter a filepath representing the directory where you want to store the later downloaded contents.
mkdir SPIDER_python
cd SPIDER_python
git init
git remote add -f origin https://github.com/Bin-Chen-Lab/spider.git
git config core.sparseCheckout true
echo "SPIDER_python/" >> .git/info/sparse-checkout
git pull origin main
```

# SPIDER usage with sample data
In R studio, load the sample query transcriptomes:
```
library(SPIDER)
data("sample_query")
```
Use SPIDER to predict on the sample query transcriptomes:
```
SPIDER_predict ( seurat_data = RNA,
                 tissue = 'pancreas',
                 disease = 'healthy',
                 SPIDER_model_file_path = '/user_path/SPIDER_python/SPIDER_python/SPIDER_weight/', #The user_path here should be the same as the user_path in "cd '/user_path'" in the previous installation part.
                 use_cell_type = 'SingleR',
                 query_cell_type = NULL,
                 protein = 'All',
                 use_pretrain = 'T', #Using pretrained SPIDER
                 save_path = ..., #enter a path to the directory where you want to save your prediction results.
                 use_python_path = ..., #If you use a specific version of python that separates from your default python configuration for reticulate, indicate the path to it. It will pass this parameter to reticulate's "use_python" function. Otherwise just set this parameter to NULL.
                 scarches_path = '/scarches_user_path/scarches-0.4.0/', #The scarches_user_path here should be the same as the scarches_user_path in "cd '/scarches_user_path'" in the previous scArches installation part.
                 all_trainable_proteins_gene_names_6_training_sets = NULL, #If you're using pretrained model, set this parameter to NULL
                 file_A_B_C_matching_table = NULL ) #If you're using pretrained model, set this parameter to NULL
```
