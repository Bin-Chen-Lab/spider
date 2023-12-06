# SPIDER
SPIDER is a zero-shot model, which enables the large-scale prediction of cell surface protein abundance from single-cell transcriptomes. 

# Conda environment
Before installing SPIDER, make sure you have installed the dependent R and python packages. To do this, you can use our environment file, which will create a conda environment named "SPIDER" with the required dependencies: <br /> <br />
In the terminal, first type the following commands. '/user_path' is a filepath representing the directory where you want to store the later downloaded contents. The downloaded folder will also contain SPIDER's pretrained weights: <br />
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

Your query transcriptome dataset should only contain one tissue and one disease. If your dataset is the combination of multiple tissues and diseases, you should first split your data into multiple subsets and run SPIDER on each of them separately. We provide a sample query dataset For users to explore. In R studio, load the sample query transcriptomes:
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
                 save_path = '/.../', 
                 use_python_path = '/...', 
                 scarches_path = '/scarches_user_path/scarches-0.4.0/',
                 all_trainable_proteins_gene_names = NULL, 
                 file_A_B_C_matching_table = NULL ) 
```
Note that:<br /><br />
SPIDER_model_file_path: The "user_path" here should be the same as your "user_path" in "cd '/user_path'" in the previous # Conda environment part. <br /><br />
save_path: Enter a path to the directory where you want to save your prediction results. <br /><br />
use_python_path: If you use a specific version of python that separates from your default python configuration for reticulate, indicate the path to it. It will pass this parameter to reticulate's "use_python" function. Otherwise just set this parameter to NULL. <br /><br />
scarches_path: The scarches_user_path here should be the same as the scarches_user_path in "cd '/scarches_user_path'" in the previous scArches installation part. <br /><br />
all_trainable_proteins_gene_names: If you're using pretrained model, set this parameter to NULL. <br /><br />
file_A_B_C_matching_table: If you're using pretrained model, set this parameter to NULL <br /><br />
You can also type the following line in R to access the help file and check more details: <br />
```
help(SPIDER_predict)
```
# Reproducibility
To find code to reproduce the results we generated in the manuscript, please visit [this separate github repository](https://github.com/Bin-Chen-Lab/spider_analysis/), which provides all code necessary to reproduce our results.
