# SPIDER
SPIDER (surface protein prediction using deep ensembles from single-cell RNA-seq) is a context-agnostic zero-shot deep ensemble model, which enables the large-scale prediction of cell surface protein abundance (e.g., can predict the abundance for >2,500 proteins) from single-cell transcriptomes. 

# Step 1: Installation of dependency packages
Before installing SPIDER, you will need to install all the dependent R and python packages ([dependent python and python packages requirements](https://github.com/Bin-Chen-Lab/spider/blob/main/dependencies/python%20package%20dependencies), and [dependent R and R packages requirements](https://github.com/Bin-Chen-Lab/spider/blob/fbcd525b52eb66a72b6257946bb30d1dff46737e/dependencies/R%20package%20dependencies.R)). We provide a convenient way for you to install these dependencies via our environment file, which will create a conda environment named "SPIDER" with the required dependencies (If your computer does not have conda, you should go the the [conda website](https://conda.io/projects/conda/en/latest/index.html) to install conda first): <br /> <br />

## 1.1
In your computer's terminal, first type the following commands. The downloaded folder will also contain SPIDER's pretrained weights: <br />
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
For users who are NOT using osx-arm64 for their computer system, they can use our yaml file to conveniently install all the dependency packages. To do this, in your terminal, type the following:

```
conda env create -f SPIDER_python/SPIDER_env/SPIDER_environment_test_basic_all.yaml
```

## 1.3 (optional) For users who use osx-arm64 for their computer system
osx-arm64 is incompatible with the Bioconda approach of installation, therefore, these users may fail to execute our yaml file in 1.2. Also see Q2 in the "frequently asked questions" section below. If you use osx-arm64 and directly run the yaml file as in 1.2, you are likely to encounter the error shown in Q2. Instead, you can replace the codes in 1.2 with the following codes:
```
CONDA_SUBDIR=osx-64 conda env create -f SPIDER_python/SPIDER_env/SPIDER_environment_test_basic_all_osx-arm64.yaml
```

## 1.4 
In your terminal, type the following commands to download the scArches package in the "SPIDER" folder. It will create another folder "scarches-0.4.0" there:
```
wget https://github.com/theislab/scarches/archive/refs/tags/v0.4.0.zip
unzip 'v0.4.0.zip'
``` 
(Alternatively, you can also simply open the link in your browser and directly download the zip file without wget. Munually unzip the file and put the "scarches-0.4.0" folder in your "SPIDER" folder.)

# Step 2: Installation of SPIDER
You should first complete step 1 before you do this step 2. After you have created the conda environment with all the dependency packages installed as previously described, first open R in the activated conda environment by typing the following commands in the terminal:
```
conda activate SPIDER
R
```
Then install the SPIDER package in R by typing the following lines in R: <br />
```
devtools::install_github(repo = 'Bin-Chen-Lab/spider', subdir = '/SPIDER')
```
Your system may ask you "Enter one or more numbers, or an empty line to skip updates", just enter an empty line to skip updates.

# Step 3: SPIDER usage with sample data
First, let's create an empty folder for saving your results. In your computer's terminal, use the following command to create another folder named SPIDER_results in your "SPIDER" folder, and then open R in the activated conda environment: <br />
```
mkdir SPIDER_results
conda activate SPIDER
R
```

Then in R (opened in the activated conda environment), load our sample query transcriptomes:
```
library(SPIDER)
data("sample_query")
```

In R, use SPIDER to predict on the sample query transcriptomes:
```
SPIDER_predict ( seurat_data = RNA,
                 tissue = 'pancreas',
                 disease = 'healthy',
                 SPIDER_model_file_path = paste0(getwd(), '/SPIDER_python/SPIDER_weight/'), 
                 use_cell_type = 'SingleR',
                 query_cell_type = NULL,
                 protein = 'All',
                 use_pretrain = 'T', #Using pretrained SPIDER
                 save_path = paste0(getwd(), '/SPIDER_results/'), 
                 use_python_path = NULL, 
                 scarches_path = paste0(getwd(), '/scarches-0.4.0/'),
                 all_trainable_proteins_gene_names = NULL, 
                 file_A_B_C_matching_table = NULL,
                 n_ensemble_members = 8 ) 
```
Note that you may need to modify "use_python_path = NULL" here according to your system setting. SPIDER will pass this parameter to [reticulate](https://rstudio.github.io/reticulate/)'s [use_python](https://rstudio.github.io/reticulate/reference/use_python.html) function. You can use "use_python_path = NULL" if your default python configuration for reticulate is the same as the python installed in your SPIDER conda environment. However, if your default python configuration for reticulate is different from the python installed in your SPIDER conda environment, you'll need to modify the "use_python_path" parameter to indicate the path to the python installed in your SPIDER conda environment. If you don't know how to locate your python path, see Q1 in "frequently asked questions" section below.<br /><br /> 

For other parameters:<br /><br />

tissue: The name of the source tissue of your transcriptome data (If your data contain multiple tissues, subset your data by tissue and run SPIDER separately on each subset). Use help(SPIDER_predict) to read more about this parameter. If your data's corresponding tissue is NOT among the 5 default tissues ('bone marrow', 'brain', 'blood', 'pleura', 'peritoneum'), use a new name that represents your data's corresponding tissue. <br /><br />

disease: The name of the disease state of your transcriptome data (If your data contain multiple diseases, subset your data by disease and run SPIDER separately on each subset). Use help(SPIDER_predict) to read more about this parameter. If your data's corresponding disease is NOT among the 4 default diseases ('healthy', 'mesothelioma', 'glioblastoma', 'leukemia'), use a new name that represents your data's corresponding disease.

SPIDER_model_file_path: This is the ABSOLUTE path to the "SPIDER_weight" folder, a sub-folder stored in your "SPIDER" folder. Avoid using the "~" symbol to locate your path. <br /><br />

save_path: This is the ABSOLUTE path to the folder where you want to save your prediction results. Avoid using the "~" symbol to locate your path.<br /><br />

scarches_path: This is the ABSOLUTE path to the "scarches-0.4.0" folder, a sub-folder stored in your "SPIDER" folder. Avoid using the "~" symbol to locate your path. <br /><br />

all_trainable_proteins_gene_names: If you're using pretrained model, set this parameter to NULL. <br /><br />

file_A_B_C_matching_table: If you're using pretrained model, set this parameter to NULL. <br /><br />

n_ensemble_members: The number of ensemble members. The default setting is 8 as in our paper. <br /><br />

You can also type the following line in R to access the help file and check more details: <br />
```
help(SPIDER_predict)
```

# Step 4: Downstream applications with SPIDER's output files

The output files from SPIDER will be stored in your specified directory. The file "all_seen_proteins_predicted.csv" contains the predicted surface abundance for all the seen proteins. The file "all_unseen_proteins_predicted.csv" contains the predicted surface abundance for all the unseen proteins. The file "confidence_score_all_unseen_proteins.csv" contains the estimated prediction confidence for all the unseen proteins.

# Frequently asked questions
#### Q1: 
I have already executed the codes from step 1 & 2 without enconutering error, but when I run step 3, why do I still enconter the following errors such as: 
```
the No non-system installation of Python could be found.
Would you like to download and install Miniconda? [Y/n]
```
and
```
Error in py_module_import(module, convert = convert) : 
  ModuleNotFoundError: No module named 'scanpy'
```
#### A1: 
These errors are because you set the "use_python_path" parameter to NULL in step 3, however, you have multiple pythons on your computer, and your default python configuration for [reticulate](https://rstudio.github.io/reticulate/) is different from the python installed in your SPIDER conda environment. To solve this problem, in the "use_python_path" parameter, you should specify the path to the python installed in your SPIDER conda environment. SPIDER will pass this parameter to [reticulate](https://rstudio.github.io/reticulate/)'s [use_python](https://rstudio.github.io/reticulate/reference/use_python.html) function. If you don't know how to locate your python path, you can locate it by doing the following:

First open python in your activated conda environment by typing the following command in your terminal:
```
conda activate SPIDER
python
```

Then type the following lines in python:
```
import sys 
sys.path[1]
```
It should return a path in the format of '.../SPIDER/lib/python39.zip'. <br /> You should set your "use_python_path" parameter as '.../SPIDER/bin/python' <br /> (the "..." parts are the same thing).

#### Q2: 
When I run the commands in 1.2, why do I encounter the following error?
```
PackagesNotFoundError: The following packages are not available from current channels:

  - bioconductor-singler
```
#### A2: 
This is likely because you use osx-arm64 for your computer system, which is incompatible with the Bioconda approach of installing the dependency packages. You should run the commands following step 1.3 instead of 1.2.

#### Q3:
Can I run SPIDER on mouse scRNA-seq data?

#### A3:
Yes, you can run SPIDER on other species' scRNA-seq data besides human data. But note that if you choose to use our pretrained SPIDER model (i.e., use_pretrain = 'T') to directly predict on other species' data, you will need to convert the gene names in your scRNA-seq data to human gene names (uppercase letters) first before you run SPIDER on your data.

# Reproducibility
To find code to reproduce the results we generated in the manuscript, please visit [this separate github repository](https://github.com/Bin-Chen-Lab/spider_analysis/), which provides all code necessary to reproduce our results.

# Citation
please cite this repo https://github.com/Bin-Chen-Lab/spider with authors (Ruoqiao Chen (chenruo4@msu.edu), Jiayu Zhou, Bin Chen (chenbi12@msu.edu), Michigan State University).


