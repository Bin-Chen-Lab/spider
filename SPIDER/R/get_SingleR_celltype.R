#' @export
get_SingleR_celltype <- function(seurat_data, save_path, specify = 'query'){
  
  #specify: 'query', 'reference'
  #---------------------------------------------------------------------------------
  library(Seurat)
  library(celldex)
  library(SingleR)
  library(dplyr)
  
  DefaultAssay(seurat_data) = "RNA"
  #----------------------------------------------------------------------------------
  cells_classification = SingleR(test = as.SingleCellExperiment(DietSeurat(seurat_data)), 
                                 ref=celldex::HumanPrimaryCellAtlasData(),
                                 labels=celldex::HumanPrimaryCellAtlasData()$label.fine)
  cells_classification2 = SingleR(test = as.SingleCellExperiment(DietSeurat(seurat_data)), 
                                  ref=celldex::HumanPrimaryCellAtlasData(),
                                  labels=celldex::HumanPrimaryCellAtlasData()$label.main)
  seurat_data[['pruned.labels']]=cells_classification$pruned.labels
  seurat_data[['pruned.labels2']]=cells_classification2$pruned.labels
  RNA_QC_meta <- data.frame(cell_type_HCAfinelabels = seurat_data[['pruned.labels']]$pruned.labels, cell_type_HCAmainlabels = seurat_data[['pruned.labels2']]$pruned.labels, umap1 = seurat_data@reductions$umap@cell.embeddings[,1], umap2 = seurat_data@reductions$umap@cell.embeddings[,2], stringsAsFactors = FALSE, row.names = rownames(seurat_data@reductions$umap@cell.embeddings))

  #------------------------------------------------------------------------------------
  RNA_QC_meta$study = seurat_data@meta.data$study 
  
  #filtering cell types:
    tmp <- RNA_QC_meta
    tmp$final_celltype = tmp$cell_type_HCAfinelabels
    tmp$final_celltype[grep('B_cell', tmp$cell_type_HCAfinelabels)] = 'B cell'
    tmp$final_celltype[grep('CMP', tmp$cell_type_HCAfinelabels)] = 'CMP'
    tmp$final_celltype[grep('Chondrocytes:', tmp$cell_type_HCAfinelabels)] = 'Chondrocyte'
    tmp$final_celltype[grep('DC:', tmp$cell_type_HCAfinelabels)] = 'Dendritic cell'
    tmp$final_celltype[grep('Endothelial_cells', tmp$cell_type_HCAfinelabels)] = 'Endothelial cell'
    tmp$final_celltype[grep('Epithelial_cells', tmp$cell_type_HCAfinelabels)] = 'Epithelial cell'
    tmp$final_celltype[grep('GMP', tmp$cell_type_HCAfinelabels)] = 'GMP'
    tmp$final_celltype[grep('Hepatocytes', tmp$cell_type_HCAfinelabels)] = 'Hepatocyte'
    tmp$final_celltype[grep('HSC', tmp$cell_type_HCAfinelabels)] = 'HSC'
    tmp$final_celltype[grep('Macrophage:', tmp$cell_type_HCAfinelabels)] = 'Macrophage'
    tmp$final_celltype[grep('Monocyte:', tmp$cell_type_HCAfinelabels)] = 'Monocyte'
    tmp$final_celltype[grep('Neurons:', tmp$cell_type_HCAfinelabels)] = 'Neuron'
    tmp$final_celltype[grep('Neutrophil:', tmp$cell_type_HCAfinelabels)] = 'Neutrophil'
    tmp$final_celltype[grep('NK_cell', tmp$cell_type_HCAfinelabels)] = 'NK cell'
    tmp$final_celltype[grep('Pre-B_cell_CD34-', tmp$cell_type_HCAfinelabels)] = 'CD34- Pre-B cell'
    tmp$final_celltype[grep('Pro-B_cell_CD34+', tmp$cell_type_HCAfinelabels)] = 'CD34+ Pro-B cell'
    tmp$final_celltype[grep('T_cell:CD4+', tmp$cell_type_HCAfinelabels)] = 'CD4+ T cell'
    tmp$final_celltype[grep('T_cell:CD8+', tmp$cell_type_HCAfinelabels)] = 'CD8+ T cell'
    tmp$final_celltype[grep('T_cell:gamma-delta', tmp$cell_type_HCAfinelabels)] = 'gamma-delta T cell'
    tmp$final_celltype[grep('T_cell:Treg', tmp$cell_type_HCAfinelabels)] = 'CD4+ T cell'
    tmp$final_celltype[grep('Tissue_stem_cells:', tmp$cell_type_HCAfinelabels)] = 'Tissue stem cell'
    RNA_QC_meta = tmp
    if(specify == 'query'){
    write.csv(RNA_QC_meta, paste0(save_path, 'query_celltype_SingleR.csv'))
    }
    if(specify == 'reference'){
      write.csv(RNA_QC_meta, paste0(save_path, 'reference_celltype_SingleR.csv'))
    }
    return(RNA_QC_meta)
  
}
