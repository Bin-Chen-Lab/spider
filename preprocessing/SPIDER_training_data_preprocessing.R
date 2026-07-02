library(Seurat)
library(dplyr)

import_RNA_count_file <- function(RNA_count_file, all_cell_type_file, dataset_id,
                                  cell_type, select_cell_type_file, tissue){
  for(t in tissue){
    
    if(t == 'BM'){
      load(RNA_count_file)
      WTA <- subset(WTA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
      features = rownames(WTA)[rowSums(WTA@assays$RNA@counts > 0) >= 5]
      WTA <- subset(WTA, features = features)
      
      RNA <- as.sparse(as.matrix(WTA@assays$RNA@counts))
      RNA <- CollapseSpeciesExpressionMatrix(RNA)
      RNA <- CreateSeuratObject(counts = RNA)
    }
    
    if(t == 'liver'|t == 'PBMC'){
      RNA <- as.sparse(read.csv(RNA_count_file, stringsAsFactors = F, row.names = 1, check.names = F))
      RNA <- CollapseSpeciesExpressionMatrix(RNA)
      RNA <- CreateSeuratObject(counts = RNA)
      RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = '^MT\\.') #Manually check
      RNA <- subset(RNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
      features = rownames(RNA)[rowSums(RNA@assays$RNA@counts > 0) >= 5]
      RNA <- subset(RNA, features = features)
    }
    
    if(t == 'peritoneum'|t == 'pleura'|t == 'BM_leukemia'|t == 'GBM'){
      RNA <- as.sparse(read.csv(RNA_count_file, stringsAsFactors = F, row.names = 1, check.names = F))
      RNA <- CollapseSpeciesExpressionMatrix(RNA)
      RNA <- CreateSeuratObject(counts = RNA)
      RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = '^MT-') #Manually check
      RNA <- subset(RNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
      features = rownames(RNA)[rowSums(RNA@assays$RNA@counts > 0) >= 5]
      RNA <- subset(RNA, features = features)
    }
    
    load(all_cell_type_file)
    cell_type_info <- get(cell_type)
    cell_type_info$cell_id = rownames(cell_type_info)
    RNA[['cell_type']] = cell_type_info$final_celltype
    RNA[['study']] = dataset_id
    
    #Select only cell types with >100 cells.
    select_cell_type <- read.csv(select_cell_type_file, stringsAsFactors = F, row.names = 1, check.names = F)
    select_cell_type = select_cell_type[, 1]
    select_cell_type_info = filter(cell_type_info, final_celltype == select_cell_type[1])
    for (i in select_cell_type[2:length(select_cell_type)]) {
      tmp = filter(cell_type_info, final_celltype == i)
      select_cell_type_info = rbind(select_cell_type_info, tmp)
    }
    cells = select_cell_type_info$cell_id
    RNA <- subset(RNA, cells = cells)
    return(RNA)
  }
}


#Run function to import RNA raw count files:
#BM (Triana et al.):
dataset_id = 'TS_NI_2021_BM'
assign(paste0(dataset_id, '_seuratobj'), import_RNA_count_file(RNA_count_file = 'WTA.rda',
                                                               all_cell_type_file = 'training_6_datasets_all_celltype_meta_SingleR_20230115.RData',
                                                               dataset_id = 'TS_NI_2021_BM',
                                                               cell_type = 'singler_13165_BM_meta',
                                                               select_cell_type_file = 'select_threshold_100_cell_type_training_6_datasets_SingleR_20230115.csv',
                                                               tissue = 'BM'))


#GBM:
dataset_id = 'GSM4972212'
assign(paste0(dataset_id, '_seuratobj'), import_RNA_count_file(RNA_count_file = 'GSM4972212_Citeseq_Human.GBM.R2_5_ND8.filtered.RNA.feature.bc.matrix.csv',
                                                               all_cell_type_file = 'training_6_datasets_all_celltype_meta_SingleR_20230115.RData',
                                                               dataset_id = 'GSM4972212',
                                                               cell_type = 'singler_GBM_meta',
                                                               select_cell_type_file = 'select_threshold_100_cell_type_training_6_datasets_SingleR_20230115.csv',
                                                               tissue = 'GBM'))


#PBMC (Hao et al.):
for(i in c('P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8')){
  for (j in c(0, 2, 7)) {
    
    dataset_id = paste(i, j, 'GSM5008738', sep = '_')
    
    assign(paste0(dataset_id, '_seuratobj'), import_RNA_count_file (RNA_count_file = paste('./PBMC (GSE164378)/GSM5008737 (RNA_3P)/RNA', i, j, '.csv', sep = '_'),
                                                                    all_cell_type_file = 'training_6_datasets_all_celltype_meta_SingleR_20230115.RData',
                                                                    dataset_id = paste(i, j, 'GSM5008738', sep = '_'),
                                                                    cell_type = paste('PBMC_meta', i, j, sep = '_'),
                                                                    select_cell_type_file = 'select_threshold_100_cell_type_training_6_datasets_SingleR_20230115.csv',
                                                                    tissue = 'PBMC'))
    
  }
}

#peritoneum:
dataset_id = 'GSM5242793'
assign(paste0(dataset_id, '_seuratobj'), import_RNA_count_file(RNA_count_file = './peritoneum_pleura (GSE172155)/peritoneum_GSE172155_RNA.csv',
                                                               all_cell_type_file = 'training_6_datasets_all_celltype_meta_SingleR_20230115.RData',
                                                               dataset_id = 'GSM5242793',
                                                               cell_type = 'singler_peritoneum_GSE172155_meta',
                                                               select_cell_type_file = 'select_threshold_100_cell_type_training_6_datasets_SingleR_20230115.csv',
                                                               tissue = 'peritoneum'))

#pleura:
dataset_id = 'GSM5242791'
assign(paste0(dataset_id, '_seuratobj'), import_RNA_count_file(RNA_count_file = './peritoneum_pleura (GSE172155)/pleura_GSE172155_RNA.csv',
                                                               all_cell_type_file = 'training_6_datasets_all_celltype_meta_SingleR_20230115.RData',
                                                               dataset_id = 'GSM5242791',
                                                               cell_type = 'singler_pleura_GSE172155_meta',
                                                               select_cell_type_file = 'select_threshold_100_cell_type_training_6_datasets_SingleR_20230115.csv',
                                                               tissue = 'pleura'))


#leukemia BM:
dataset_id = 'GSE143363'
assign(paste0(dataset_id, '_seuratobj'), import_RNA_count_file(RNA_count_file = './BM leukemia (GSE143363)/BM_leukemia_GSE143363_RNA.csv',
                                                               all_cell_type_file = 'training_6_datasets_all_celltype_meta_SingleR_20230115.RData',
                                                               dataset_id = 'GSE143363',
                                                               cell_type = 'singler_GSE143363_meta',
                                                               select_cell_type_file = 'select_threshold_100_cell_type_training_6_datasets_SingleR_20230115.csv',
                                                               tissue = 'BM_leukemia'))

