R packages requirements:
SingleR
celldex
SingleCellExperiment
Seurat (version=4.3.0)
anndata
reticulate

#------------------------------------------------------------------------------------
#Demo code of installation:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("SingleCellExperiment") 

install.packages("Seurat", version = '4.3.0')
install.packages("anndata")
install.packages("reticulate")
