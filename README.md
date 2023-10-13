# SPIDER
SPIDER (surface protein prediction using deep ensembles from single-cell RNA-seq) is a context-agnostic zero-shot deep ensemble model, which enables the large-scale prediction of cell surface protein abundance. 

# Conda environment
Before installing SPIDER, make sure you have installed the dependent R and python packages. You can do this by using our environment file, which will create a conda environment named "SPIDER" with the required dependencies: <br /> <br />
In the terminal, first type the following commands. '/user_path' is a filepath representing the directory where you want to store the later downloaded contents: <br />
```
cd '/user_path' 
mkdir SPIDER
cd SPIDER
git init
git remote add -f origin https://github.com/Bin-Chen-Lab/spider.git
git config core.sparseCheckout true
echo "SPIDER_python/" >> .git/info/sparse-checkout
git pull origin main
cd '/user_path/SPIDER/SPIDER_python/SPIDER_env'
conda env create -f SPIDER_environment_test_basic_all.yaml
```
Also, in the terminal, type the following commands. '/scarches_user_path' is a filepath representing the directory where you want to store the later downloaded dependent scArches package:
```
cd '/scarches_user_path' 
wget https://github.com/theislab/scarches/archive/refs/tags/v0.4.0.zip
unzip 'v0.4.0.zip'
```

# Alternative way of installing dependencies
Alternatively, you can choose to install these R and python dependent packages by first run the R lines as in: <br />
https://github.com/Bin-Chen-Lab/spider/blob/fbcd525b52eb66a72b6257946bb30d1dff46737e/dependencies/R%20package%20dependencies.R <br />
Then run the commands in terminal as in: <br />
https://github.com/Bin-Chen-Lab/spider/blob/fbcd525b52eb66a72b6257946bb30d1dff46737e/dependencies/python%20package%20dependencies <br />

# Installation
Install the SPIDER package: <br />
Open R studio, type the following lines: <br />
```
devtools::install_github(repo = 'Bin-Chen-Lab/spider', subdir = '/SPIDER')
``` 

# SPIDER usage with sample data
If you have created the conda environment as previously described, first enter the conda environment by typing the following command in the terminal:
```
conda activate SPIDER
```

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
                 SPIDER_model_file_path = '/user_path/SPIDER/SPIDER_python/SPIDER_weight/', 
                 use_cell_type = 'SingleR',
                 query_cell_type = NULL,
                 protein = 'All',
                 use_pretrain = 'T', #Using pretrained SPIDER
                 save_path = ..., 
                 use_python_path = ..., 
                 scarches_path = '/scarches_user_path/scarches-0.4.0/',
                 all_trainable_proteins_gene_names = NULL, 
                 file_A_B_C_matching_table = NULL ) 
```
Note that:<br />
SPIDER_model_file_path: The "user_path" here should be the same as your "user_path" in "cd '/user_path'" in the previous # Conda environment part.
save_path: Enter a path to the directory where you want to save your prediction results. <br />
use_python_path: If you use a specific version of python that separates from your default python configuration for reticulate, indicate the path to it. It will pass this parameter to reticulate's "use_python" function. Otherwise just set this parameter to NULL. <br />
scarches_path: The scarches_user_path here should be the same as the scarches_user_path in "cd '/scarches_user_path'" in the previous scArches installation part. <br />
all_trainable_proteins_gene_names: If you're using pretrained model, set this parameter to NULL. <br />
file_A_B_C_matching_table: If you're using pretrained model, set this parameter to NULL <br />
You can also type the following line in R to access the help file and check more details: <br />
```
help(SPIDER_predict)
```
