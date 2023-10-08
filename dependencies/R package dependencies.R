if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("SingleCellExperiment") 

install.packages("Seurat")
install.packages("anndata")