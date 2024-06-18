# SPIDER
SPIDER (surface protein prediction using deep ensembles from single-cell RNA-seq) is a context-agnostic zero-shot deep ensemble model, which enables the large-scale prediction of cell surface protein abundance. 

# Step 1: Installation of dependency packages
Before installing SPIDER, make sure you have installed the dependent R and python packages. You can do this by using our environment file, which will create a conda environment named "SPIDER" with the required dependencies: <br /> <br />

## 1.1
In the terminal, first type the following commands. The downloaded folder will also contain SPIDER's pretrained weights: <br />
```
mkdir SPIDER
cd SPIDER
git init
git remote add -f origin https://github.com/Bin-Chen-Lab/spider.git
git config core.sparseCheckout true
echo "SPIDER_python/" >> .git/info/sparse-checkout
git pull origin main
``` 

## 1.2 
Users can use our yaml file to conveniently install all the dependency packages. To do this, in your terminal, type the following:

```
cd 'SPIDER_python/SPIDER_env'
conda env create -f SPIDER_environment_test_basic_all.yaml
```

## 1.3 (optional)
For some users who fail to execute our yaml file in 1.2 due to system incompatibility, we also provide a manual way for these users to download dependency packages. These codes of manual installation have the same effect as the yaml file in 1.2. First, create a conda environment with specified R and python versions by typing the following commands in your terminal:

```
conda create -n SPIDER python=3.9.2 
conda install conda-forge::r-base=4.1.1
```

Then, in your terminal, run the following commands in the following link file to manually install all the python dependency packages: <br />
https://github.com/Bin-Chen-Lab/spider/blob/fbcd525b52eb66a72b6257946bb30d1dff46737e/dependencies/python%20package%20dependencies <br />

Then, in your R studio, run the R lines as in the following link file to manually install all the R dependency packages: <br />
https://github.com/Bin-Chen-Lab/spider/blob/fbcd525b52eb66a72b6257946bb30d1dff46737e/dependencies/R%20package%20dependencies.R <br />


## 1.4 
In your terminal, type the following commands to download the scArches package. It will create a folder "scarches-0.4.0":
```
wget https://github.com/theislab/scarches/archive/refs/tags/v0.4.0.zip
unzip 'v0.4.0.zip'
``` 
(Alternatively, you can also simply open the link in your browser and directly download the folder without wget.)

# Step 2: Installation of SPIDER
You should first complete step 1 before you do this step 2. After you have created the conda environment as previously described, first enter the conda environment by typing the following command in the terminal:
```
conda activate SPIDER
```
To install the SPIDER package: <br />
Open your R studio, type the following lines in R studio: <br />
```
devtools::install_github(repo = 'Bin-Chen-Lab/spider', subdir = '/SPIDER')
``` 

# Step 3: SPIDER usage with sample data

In R studio, load the sample query transcriptomes:
```
library(SPIDER)
data("sample_query")
```

In R studio, use SPIDER to predict on the sample query transcriptomes:
```
SPIDER_predict ( seurat_data = RNA,
                 tissue = 'pancreas',
                 disease = 'healthy',
                 SPIDER_model_file_path = '.../SPIDER/SPIDER_python/SPIDER_weight/', 
                 use_cell_type = 'SingleR',
                 query_cell_type = NULL,
                 protein = 'All',
                 use_pretrain = 'T', #Using pretrained SPIDER
                 save_path = '/.../', 
                 use_python_path = '/...', 
                 scarches_path = '.../scarches-0.4.0/',
                 all_trainable_proteins_gene_names = NULL, 
                 file_A_B_C_matching_table = NULL,
                 n_ensemble_members = 8 ) 
```
Note that:<br /><br />

SPIDER_model_file_path: This is the absolute path to the "SPIDER_weight" folder, a sub-folder stored in the "SPIDER" folder which you previously created in step1.1.1. <br /><br />

save_path: Enter a path to the directory where you want to save your prediction results. <br /><br />

use_python_path: If you use a specific version of python that separates from your default python configuration for reticulate, indicate the path to it. It will pass this parameter to reticulate's "use_python" function. Otherwise just set this parameter to NULL. If you don't know how to locate your python path, see Q1 in "frequently asked questions".<br /><br />

scarches_path: This is the absolute path to the "scarches-0.4.0" folder which you previously downloded in the scArches installation part. <br /><br />

all_trainable_proteins_gene_names: If you're using pretrained model, set this parameter to NULL. <br /><br />

file_A_B_C_matching_table: If you're using pretrained model, set this parameter to NULL. <br /><br />

n_ensemble_members: The number of ensemble members. The default setting is 8 as in our paper. <br /><br />

You can also type the following line in R to access the help file and check more details: <br />
```
help(SPIDER_predict)
``` 
# Frequently asked questions
#### Q1: 
I have already executed the codes from step 1 & 2 without enconutering error, but when I run step 3, why do I still enconter the following error? 
```
Error in py_module_import(module, convert = convert) : 
  ModuleNotFoundError: No module named 'scanpy'
```
#### A1: 
This is likely because you set the "use_python_path" parameter to NULL in step 3, however, you have multiple python installed on your computer, and the default python is not the same one as you installed in your SPIDER environment. To solve this problem, you should specify the correct python path using the "use_python_path" parameter. You can check the correct python path by doing the following:

Open python, and type the following codes in python:
```
import sys 
sys.path[1]
```
It should return a path in the format of '.../SPIDER/lib/python39.zip'. <br /> You should set your "use_python_path" parameter as '.../SPIDER/bin/python' <br /> (the "..." part keep the same).

#### Q2: 
When I run the commands in 1.2, why do I encounter the following error?
```
PackagesNotFoundError: The following packages are not available from current channels:

  - bioconductor-singler
```
#### A2: 
This is likely because you use osx-arm64 for your computer system, which is incompatible with the Bioconda approach of installing the dependency packages. You may need to manually install the dependency packages as in 1.3 instead of 1.2. 

# Reproducibility
To find code to reproduce the results we generated in the manuscript, please visit [this separate github repository](https://github.com/Bin-Chen-Lab/spider_analysis/), which provides all code necessary to reproduce our results.

# Citation
please cite this repo https://github.com/Bin-Chen-Lab/spider with authors (Ruoqiao Chen (chenruo4@msu.edu), Jiayu Zhou, Bin Chen (chenbi12@msu.edu), Michigan State University).
