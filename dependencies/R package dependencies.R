if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("SingleCellExperiment") 

install.packages("Seurat", version = '4.4.0')
install.packages("anndata")
install.packages("reticulate")
