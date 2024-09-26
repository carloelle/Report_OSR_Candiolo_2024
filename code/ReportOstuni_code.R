#preprocessing lee et al. 2020

library(Matrix)
library(data.table)

lee_k<-fread('/home/carlo/Desktop/transferdata/lee/lee_korean/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt',sep='\t',header = TRUE)

genes_k=as.character(unlist(lee_k[1:33694,1]))
cells_k=as.character(names(unlist(lee_k[1,2:6390])))
lee_k<-as.matrix(lee_k[1:33694,2:6390])

lee_k<-Matrix(lee_k,sparse = TRUE)

rownames(lee_k)=genes_k
colnames(lee_k)=cells_k

library(Seurat)

lee_k_metadata<-read.table('/home/carlo/Desktop/transferdata/lee/lee_korean/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt',sep='\t',header = TRUE)
rownames(lee_k_metadata)=lee_k_metadata$Index

lee_k<-CreateSeuratObject(counts = lee_k,assay = 'RNA',meta.data = lee_k_metadata)
rm(cells_k,genes_k,lee_k_metadata)

lee_b<-fread('/home/carlo/Desktop/transferdata/lee/lee_belgian/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt',sep='\t',header=TRUE)

genes_b=as.character(unlist(lee_b[1:33694,1]))
cells_b=as.character(names(unlist(lee_b[1,2:27415])))
lee_b<-as.matrix(lee_b[1:33694,2:27415])

lee_b<-Matrix(lee_b,sparse = TRUE)

rownames(lee_b)=genes_b
colnames(lee_b)=cells_b

lee_b_metadata<-read.table('/home/carlo/Desktop/transferdata/lee/lee_belgian/GSE144735_processed_KUL3_CRC_10X_annotation.txt',sep = '\t',header=TRUE)
rownames(lee_b_metadata)=lee_b_metadata$Index

lee_b<-CreateSeuratObject(counts = lee_b,assay = 'RNA',meta.data = lee_b_metadata)
rm(cells_b,genes_b,lee_b_metadata)

#BiocManager::install("harmony")

lee_k@meta.data$Cohort<-rep('Korean',6389)
lee_b@meta.data$Cohort<-rep('Belgian',27414)

#check
rownames(lee_b@assays$RNA@features@.Data)==rownames(lee_k@assays$RNA@features@.Data)
#align everything before merging!

rownames(lee_b@assays$RNA@layers$counts)=rownames(lee_b@assays$RNA@features@.Data)
rownames(lee_k@assays$RNA@layers$counts)=rownames(lee_k@assays$RNA@features@.Data)

lee_b@assays$RNA@layers$counts<-lee_b@assays$RNA@layers$counts[rownames(lee_k@assays$RNA@layers$counts),]
rownames(lee_b@assays$RNA@features@.Data)=rownames(lee_k@assays$RNA@features@.Data)

matrix=as.matrix(cbind(lee_b@assays$RNA@layers$counts,
                       lee_k@assays$RNA@layers$counts))

colnames(matrix)<-c(rownames(lee_b@meta.data),
                    rownames(lee_k@meta.data))
rownames(matrix)<-rownames(lee_b@assays$RNA@features@.Data)


lee2<-CreateSeuratObject(counts = matrix,meta.data = rbind(lee_b@meta.data,
                                                           lee_k@meta.data))
rm(matrix,lee_b,lee_k)


library(dplyr)
lee2<-NormalizeData(lee2)%>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE,features = VariableFeatures(lee2)) %>% 
  RunPCA(pc.genes = VariableFeatures(lee2), npcs = 20, verbose = FALSE)

lee2<-subset(lee2,features = VariableFeatures(lee2))


library(harmony)
lee2<-RunHarmony(lee2,'Cohort',lambda=NULL)


lee2 <- lee2 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()



library(reticulate)
use_python('/home/carlo/.local/share/r-miniconda/envs/giotto_env/bin/python', required = TRUE)
load("~/Desktop/transferdata/wanglong/signfinal.Robj")

library(circlize)
library(data.table)
library(Giotto)


library(tidyr)
process_PAGEanalysis_sc <- function(seurat_object, signatures_all,only_fibro=F) {
  if(only_fibro==T){signatures_all <- signatures_all[!names(signatures_all) %in% c("CCM", "EMRM")]}
  raw_exprs <- seurat_object@assays$RNA@layers$counts
  rownames(raw_exprs) <- rownames(seurat_object@assays$RNA@features@.Data)
  colnames(raw_exprs) <- rownames(seurat_object@meta.data)
  #spatial_locs <- as.data.table(GetTissueCoordinates(seurat_object)[,1:2])
  #colnames(spatial_locs) <- c("x", "y")
  
  myGiottoObj <- createGiottoObject(raw_exprs = raw_exprs)
  myGiottoObj <- normalizeGiotto(gobject = myGiottoObj)
  
  # Create signature matrix for initial PAGE analysis
  all_signatures <- names(signatures_all)
  signature_matrix_complete <- makeSignMatrixPAGE(
    sign_names = all_signatures,
    sign_list = lapply(all_signatures, function(sign) signatures_all[[sign]])
  )
  
  # Run initial PAGE enrichment analysis on all signatures
  myGiottoObj_initial <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 5,
    output_enrichment = c("original")
  )
  
   #Run PAGE enrichment analysis with p-values
  myGiottoObj_pval <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 5,
    p_value = TRUE,
    reverse_log_scale = FALSE,
    return_gobject = FALSE,
    n_times=500
  )
  
  # Extract p-value results
  reshaped_pval <- myGiottoObj_pval$DT %>%
    select(cell_ID, cell_type, pval) %>%  
   pivot_wider(names_from = cell_type, values_from = pval, values_fill = 0) 
  
  pzscore_df<-as.data.frame(reshaped_pval)
  colnames(pzscore_df)=paste0(colnames(pzscore_df), "_pvalPAGE")
  seurat_object@meta.data<-cbind(seurat_object@meta.data, pzscore_df)
  
  zscore_df<-as.data.frame(myGiottoObj_initial@spatial_enrichment$PAGE)
  zscore_df <- as.data.frame(zscore_df)
  colnames(zscore_df) <- paste0(colnames(zscore_df), "PAGE")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, zscore_df)
  
  return(seurat_object)
}

lee2<-process_PAGEanalysis_sc(seurat_object = lee2,signatures_all = sign)

lee2<-AddModuleScore(lee2,features = sign,name = names(sign),ctrl = 10)

lee2_md<-lee2@meta.data


library(ggplot2)

pvalnames<-colnames(lee2_md)[grepl('_pvalPAGE$',colnames(lee2_md)) & 
                  !grepl('cell_',colnames(lee2_md))]

pvaldf<-lee2_md[,pvalnames]<0.05
colnames(pvaldf)=paste(pvalnames,'_L',sep = '')
lee2_md<-cbind(lee2_md,pvaldf)

total_cells <- table(lee2_md$Cell_subtype)

gene_modules <- c("EMRPAGE", "MF1PAGE", "MF2PAGE", "cCAFPAGE", "iCAFPAGE",
                  "mCAFPAGE", "mrCAFPAGE", "myCAFPAGE", "vCAFPAGE")

gene_modules_module <- c("EMR1","MF13","MF24","cCAF5","iCAF6","mCAF7","mrCAF2","myCAF8","vCAF9")

pval_columns <- c("EMR_pvalPAGE_L", "MF1_pvalPAGE_L", "MF2_pvalPAGE_L", "cCAF_pvalPAGE_L",
                  "iCAF_pvalPAGE_L", "mCAF_pvalPAGE_L", "mrCAF_pvalPAGE_L",
                  "myCAF_pvalPAGE_L", "vCAF_pvalPAGE_L")

stopifnot(length(gene_modules_module) == length(pval_columns))

df <- lee2_md 
summary_list <- list()

for (i in seq_along(gene_modules_module)) {
  gene_module <- gene_modules_module[i]
  pval_column <- pval_columns[i]
  
  df_subset <- df %>%
    group_by(Cell_subtype) %>%
    summarise(
      Avg_Expression = mean(.data[[gene_module]], na.rm = TRUE),
      Percent_Positive = sum(.data[[pval_column]] == TRUE, na.rm = TRUE) / n() * 100
    ) %>%
    mutate(Gene_Module = gene_module)
  
  summary_list[[i]] <- df_subset
}


df_summary <- bind_rows(summary_list)



pdf('dotplot_cellsubtype_percentSTATtestsize_module1.pdf')
ggplot(df_summary, aes(x = Gene_Module, y = Cell_subtype)) +
  geom_point(aes(size = Percent_Positive, color = Avg_Expression)) +
  scale_color_gradient2(low ='white', mid='lightgray',high = "darkred",limits=c(-0.3,1),na.value = 'gray50',midpoint = 0) +
  theme_minimal() +
  labs(
    x = "Mean Module value",
    y = "Cell subtype",
    title = "DotPlot of Gene Module Expression by Cell Subtype",
    size = "Percent of Statistically \nPositive Cells (perm. on PAGE)",
    color = "Module Expression"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

