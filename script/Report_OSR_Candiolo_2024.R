#Supporting Code for Report, Collaboration between OSR and IRCC Candiolo, Part of PhD Thesis of Carlo Leonardi, PhD 

#Preprocessing lee et al. 2020 - single cell data atlas-level

library(Matrix)
library(data.table)

lee_k<-fread('~/Report_OSR_Candiolo_2024/data/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt',sep='\t',header = TRUE)

genes_k=as.character(unlist(lee_k[1:33694,1]))
cells_k=as.character(names(unlist(lee_k[1,2:6390])))
lee_k<-as.matrix(lee_k[1:33694,2:6390])

lee_k<-Matrix(lee_k,sparse = TRUE)

rownames(lee_k)=genes_k
colnames(lee_k)=cells_k

library(Seurat)

lee_k_metadata<-read.table('~/Report_OSR_Candiolo_2024/data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt',sep='\t',header = TRUE)
rownames(lee_k_metadata)=lee_k_metadata$Index

lee_k<-CreateSeuratObject(counts = lee_k,assay = 'RNA',meta.data = lee_k_metadata)
rm(cells_k,genes_k,lee_k_metadata)

lee_b<-fread('~/Report_OSR_Candiolo_2024/data/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt',sep='\t',header=TRUE)

genes_b=as.character(unlist(lee_b[1:33694,1]))
cells_b=as.character(names(unlist(lee_b[1,2:27415])))
lee_b<-as.matrix(lee_b[1:33694,2:27415])

lee_b<-Matrix(lee_b,sparse = TRUE)

rownames(lee_b)=genes_b
colnames(lee_b)=cells_b

lee_b_metadata<-read.table('~/Report_OSR_Candiolo_2024/data/GSE144735_processed_KUL3_CRC_10X_annotation.txt',sep = '\t',header=TRUE)
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


library(lisi)
lisi<-compute_lisi(t(lee2@assays$RNA@layers$data),meta_data = lee2@meta.data,label_colnames='Cohort')

lee2$LISI<-lisi$Cohort


data_to_plot <- data.frame(
  Cell_subtype = lee2@meta.data$Cell_subtype,
  Cohort = lee2@meta.data$Cohort,
  LISI = lee2@meta.data$LISI
)

# Ensure the cell types and cohorts are factors for plotting
data_to_plot$Cell_subtype <- factor(data_to_plot$Cell_subtype)
data_to_plot$Cohort <- factor(data_to_plot$Cohort)

# Summarize by mean LISI value for each Cell_subtype-Cohort combination
library(dplyr)
summary_data <- data_to_plot %>%
  group_by(Cell_subtype, Cohort) %>%
  summarise(mean_LISI = mean(LISI, na.rm = TRUE)) %>%
  ungroup()

# Reshape data for heatmap
library(tidyr)
heatmap_data <- summary_data %>%
  pivot_wider(names_from = Cohort, values_from = mean_LISI)



library(ggplot2)

# Melt the data into long format if needed for ggplot
library(reshape2)
long_heatmap_data <- melt(heatmap_data, id.vars = "Cell_subtype")

# Plot the heatmap

long_heatmap_data<-long_heatmap_data%>%
  filter(!is.na(value))%>%
  filter(value>1.000796)

setwd('~/Report_OSR_Candiolo_2024/script/')
p1<-ggplot(long_heatmap_data, aes(x = variable, y = Cell_subtype, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "darkred",mid = 'orange',na.value = 'white',midpoint = 1.2) +
  labs(x = "Cohort", y = "Cell Subtype", fill = "Mean LISI") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line())

save(p1,file = 'p1.Robj')



library(reticulate)
use_python('/home/carlo/.local/share/r-miniconda/envs/giotto_env/bin/python', required = TRUE)
load("~/Report_OSR_Candiolo_2024/data/signfinal.Robj")


CRIS<-read.csv('~/Report_OSR_Candiolo_2024/data/CRIS_sign.csv',header = T)
CRIS<-data.frame(Genes=CRIS$X.1[3:dim(CRIS)[1]],CRIS=CRIS$X.2[3:dim(CRIS)[1]])


sign$CRIS_A<-CRIS%>%filter(CRIS=='CRIS-A')%>%select(Genes)%>%unlist()%>%unname()
sign$CRIS_B<-CRIS%>%filter(CRIS=='CRIS-B')%>%select(Genes)%>%unlist()%>%unname()
sign$CRIS_C<-CRIS%>%filter(CRIS=='CRIS-C')%>%select(Genes)%>%unlist()%>%unname()
sign$CRIS_D<-CRIS%>%filter(CRIS=='CRIS-D')%>%select(Genes)%>%unlist()%>%unname()
sign$CRIS_E<-CRIS%>%filter(CRIS=='CRIS-E')%>%select(Genes)%>%unlist()%>%unname()
emr<-sign$EMR
sign$EMR=NULL
sign$iCAF_Cords=NULL
sign$mCAF_Cords=NULL
sign$dCAF_Cords=NULL
sign$apCAF_Cords=NULL
sign$pericyte_Cords=NULL
sign$vCAFrCAF_Cords=NULL
sign$tCAF_Cords=NULL
sign$MF3=NULL
sign$MF4=NULL
sign$dCAF=NULL
sign$EMR=emr

library(circlize)
library(data.table)
library(Giotto)


library(tidyr)
#we adapt PAGE Enrichment scores for single cell
process_PAGEanalysis_sc <- function(seurat_object, signatures_all) {
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

gene_modules <- c("CRIS_APAGE","CRIS_BPAGE","CRIS_CPAGE","CRIS_DPAGE","CRIS_EPAGE","EMRPAGE","MF1PAGE","MF2PAGE",
                  "cCAFPAGE","iCAFPAGE","mCAFPAGE","mrCAFPAGE","myCAFPAGE","vCAFPAGE")

gene_modules_module <- c("mrCAF1","MF12","MF23","cCAF4","iCAF5","mCAF6","myCAF7","vCAF8","CRIS_A9","CRIS_B10","CRIS_C11","CRIS_D12","CRIS_E13","EMR14")

pval_columns <- c("CRIS_A_pvalPAGE_L","CRIS_C_pvalPAGE_L","mCAF_pvalPAGE_L","myCAF_pvalPAGE_L","CRIS_B_pvalPAGE_L","vCAF_pvalPAGE_L",
                  "cCAF_pvalPAGE_L","MF1_pvalPAGE_L","CRIS_E_pvalPAGE_L","mrCAF_pvalPAGE_L","CRIS_D_pvalPAGE_L","MF2_pvalPAGE_L","iCAF_pvalPAGE_L","EMR_pvalPAGE_L")

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

df_summary$Gene_Module<-factor(df_summary$Gene_Module,levels = c("CRIS_A9","CRIS_B10","CRIS_C11","CRIS_D12","CRIS_E13","EMR14","mrCAF1","MF12","iCAF5","mCAF6","MF23","cCAF4","myCAF7","vCAF8"))


p2<-ggplot(df_summary, aes(x = Gene_Module, y = Cell_subtype)) +
  geom_point(aes(size = Percent_Positive, color = Avg_Expression)) +
  scale_color_gradient2(low ='white', mid='white',high = "darkred",limits=c(-0.3,1),na.value = 'gray50',midpoint = -0.1) +
  theme_minimal() +
  labs(
    x = "Mean Module value",
    y = "Cell subtype",
    title = "Real Gene Expression vs. Relative Expression \nin scRNAseq (atlas-level)",
    size = "Percent of Statistically \nPositive Cells (perm. on PAGE)",
    color = "Module Expression"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line())
save(p2,file = 'p2.Robj')




rm(lee2,lee2_md)


#### part on weighted correlation matrix and overlap analysis 
library(Seurat)
library(SeuratObject)

#load datasets related to Wang Long et al. 
load("~/Report_OSR_Candiolo_2024/data/GSMobjfinal.Robj")


# Custom subsetting function for Visium objects
custom_subset_visium <- function(seurat_object, cells) {
  # Subset the counts matrix
  subset_counts <- seurat_object@assays$Spatial@layers$counts[, cells, drop = FALSE]
  # Create a new Seurat object with the subsetted counts matrix and metadata
  new_seurat_object <- CreateSeuratObject(counts = subset_counts, meta.data = seurat_object@meta.data[cells, ])
  
  
  # Subset the image data if present
  if (!is.null(seurat_object@images)) {
    new_seurat_object@images <- seurat_object@images
    for (image_name in names(seurat_object@images)) {
      new_seurat_object@images[[image_name]]@boundaries$centroids@coords <- seurat_object@images[[image_name]]@boundaries$centroids@coords[cells, , drop = FALSE]
    }
  }
  
  # Subset the graph data if present
  if (!is.null(seurat_object@graphs)) {
    graph_object <- as.Graph(subset_counts)
    new_seurat_object@graphs <- list(graph = graph_object)
  }
  
  return(new_seurat_object)
}

library(Matrix)
allspatial<-CreateSeuratObject(counts = Matrix(cbind(GSM7058756_C1@assays$Spatial@layers$counts,
                                                     GSM7058757_C2@assays$Spatial@layers$counts,
                                                     GSM7058758_C3@assays$Spatial@layers$counts,
                                                     GSM7058759_C4@assays$Spatial@layers$counts),sparse = T),
                               meta.data = as.data.frame(rbind(GSM7058756_C1@meta.data,
                                                               GSM7058757_C2@meta.data,
                                                               GSM7058758_C3@meta.data,
                                                               GSM7058759_C4@meta.data)))

#allspatial1<-allspatial@meta.data
#pdf('allspatial_wl.pdf')
#ggplot(allspatial1,aes(x=log10(Fibroblasts_C2L),y=Fibroblasts_C2L_percentage,colour = orig.ident))+geom_point()
#dev.off()

#from here, we see that if we want to establish an inter-patient cutoff, 20% and 0.85 quantile of the joint distribution is acceptable --> it gives us dataset-specific enrichment 

spotstoconsider<-intersect(rownames(allspatial@meta.data)[allspatial$Fibroblasts_C2L>quantile(allspatial$Fibroblasts_C2L,.85)],
                           rownames(allspatial@meta.data)[allspatial$Fibroblasts_C2L_percentage>20])

GSM7058759_C4_fibroenriched<-custom_subset_visium(GSM7058759_C4,intersect(rownames(GSM7058759_C4@meta.data),spotstoconsider))
GSM7058756_C1_fibroenriched<-custom_subset_visium(GSM7058756_C1,intersect(rownames(GSM7058756_C1@meta.data),spotstoconsider))
GSM7058757_C2_fibroenriched<-custom_subset_visium(GSM7058757_C2,intersect(rownames(GSM7058757_C2@meta.data),spotstoconsider))
GSM7058758_C3_fibroenriched<-custom_subset_visium(GSM7058758_C3,intersect(rownames(GSM7058758_C3@meta.data),spotstoconsider))
rm(allspatial)

load('~/Report_OSR_Candiolo_2024/data/ValdolivasSTObj.Robj')

#These datasets were published in Valdolias et al., alongside with a deconvolution using unmatched data from Lee et al. 

deconvKorean<-read.csv('~/Report_OSR_Candiolo_2024/data/DeconvKorean_cell_density_q05.csv')
deconvBelgian<-read.csv('~/Report_OSR_Candiolo_2024/data/DeconvBelgian_cell_density_q05.csv')

colnames(deconvKorean)=gsub('q05_spot_factors','',colnames(deconvKorean))
colnames(deconvBelgian)=gsub('q05_spot_factors','',colnames(deconvBelgian))

#polishing for incorporating the deconvolutions

id_korean <- c()
barcode_korean <- c()
for (i in 1:length(deconvKorean$spot_id)) {
  elements <- strsplit(deconvKorean$spot_id[i], '_')[[1]]
  if (elements[1] == "Count") {
    id_korean <- c(id_korean, paste(elements[1:4], collapse = '_'))
    barcode_korean <- c(barcode_korean, elements[5])
  } else {
    id_korean <- c(id_korean, paste(elements[1:3], collapse = '_'))
    barcode_korean <- c(barcode_korean, elements[4])
  }
}

id_belgian <- c()
barcode_belgian <- c()
for (i in 1:length(deconvBelgian$spot_id)) {
  elements <- strsplit(deconvBelgian$spot_id[i], '_')[[1]]
  if (elements[1] == "Count") {
    id_belgian <- c(id_belgian, paste(elements[1:4], collapse = '_'))
    barcode_belgian <- c(barcode_belgian, elements[5])
  } else {
    id_belgian <- c(id_belgian, paste(elements[1:3], collapse = '_'))
    barcode_belgian <- c(barcode_belgian, elements[4])
  }
}

deconvKorean$ST.exp<-id_korean
deconvBelgian$ST.exp<-id_belgian

deconvKorean$ST.barcode<-barcode_korean
deconvBelgian$ST.barcode<-barcode_belgian


library(dplyr)
v573_1_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN048_A121573_Rep1')
v573_1_MD<-v573_1_MD[,2:39]
v573_1_MD<-v573_1_MD[v573_1_MD$ST.barcode%in%rownames(V573_1@meta.data),]
V573_1<-V573_1[,v573_1_MD$ST.barcode]
rownames(v573_1_MD)=rownames(V573_1@meta.data)
colnames(v573_1_MD)=gsub('$','_KoreanDeconv',colnames(v573_1_MD))
V573_1@meta.data<-cbind(V573_1@meta.data,v573_1_MD)


library(dplyr)
v573_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN048_A121573_Rep2')
v573_2_MD<-v573_2_MD[,2:39]
v573_2_MD<-v573_2_MD[v573_2_MD$ST.barcode%in%rownames(V573_2@meta.data),]
V573_2<-V573_2[,v573_2_MD$ST.barcode]
rownames(v573_2_MD)=rownames(V573_2@meta.data)
colnames(v573_2_MD)=gsub('$','_KoreanDeconv',colnames(v573_2_MD))
V573_2@meta.data<-cbind(V573_2@meta.data,v573_2_MD)


library(dplyr)
v371_1_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN048_A416371_Rep1')
v371_1_MD<-v371_1_MD[,2:39]
v371_1_MD<-v371_1_MD[v371_1_MD$ST.barcode%in%rownames(V371_1@meta.data),]
V371_1<-V371_1[,v371_1_MD$ST.barcode]
rownames(v371_1_MD)=rownames(V371_1@meta.data)
colnames(v371_1_MD)=gsub('$','_KoreanDeconv',colnames(v371_1_MD))
V371_1@meta.data<-cbind(V371_1@meta.data,v371_1_MD)


library(dplyr)
v371_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN048_A416371_Rep2')
v371_2_MD<-v371_2_MD[,2:39]
v371_2_MD<-v371_1_MD[v371_2_MD$ST.barcode%in%rownames(V371_2@meta.data),]
V371_2<-V371_1[,v371_2_MD$ST.barcode]
rownames(v371_2_MD)=rownames(V371_2@meta.data)
colnames(v371_2_MD)=gsub('$','_KoreanDeconv',colnames(v371_2_MD))
V371_2@meta.data<-cbind(V371_2@meta.data,v371_2_MD)


library(dplyr)
v763_1_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN123_A551763_Rep1')
v763_1_MD<-v763_1_MD[,2:39]
v763_1_MD<-v763_1_MD[v763_1_MD$ST.barcode%in%rownames(V763_1@meta.data),]
V763_1<-V763_1[,v763_1_MD$ST.barcode]
rownames(v763_1_MD)=rownames(V763_1@meta.data)
colnames(v763_1_MD)=gsub('$','_KoreanDeconv',colnames(v763_1_MD))
V763_1@meta.data<-cbind(V763_1@meta.data,v763_1_MD)


library(dplyr)
v763_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN124_A551763_Rep2')
v763_2_MD<-v763_2_MD[,2:39]
v763_2_MD<-v763_2_MD[v763_2_MD$ST.barcode%in%rownames(V763_2@meta.data),]
V763_2<-V763_2[,v763_2_MD$ST.barcode]
rownames(v763_2_MD)=rownames(V763_2@meta.data)
colnames(v763_2_MD)=gsub('$','_KoreanDeconv',colnames(v763_2_MD))
V763_2@meta.data<-cbind(V763_2@meta.data,v763_2_MD)


library(dplyr)
v688_1_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN123_A595688_Rep1')
v688_1_MD<-v688_1_MD[,2:39]
v688_1_MD<-v688_1_MD[v688_1_MD$ST.barcode%in%rownames(V688_1@meta.data),]
V688_1<-V688_1[,v688_1_MD$ST.barcode]
rownames(v688_1_MD)=rownames(V688_1@meta.data)
colnames(v688_1_MD)=gsub('$','_KoreanDeconv',colnames(v688_1_MD))
V688_1@meta.data<-cbind(V688_1@meta.data,v688_1_MD)


library(dplyr)
v688_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN124_A595688_Rep2')
v688_2_MD<-v688_2_MD[,2:39]
v688_2_MD<-v688_2_MD[v688_2_MD$ST.barcode%in%rownames(V688_2@meta.data),]
V688_2<-V688_2[,v688_2_MD$ST.barcode]
rownames(v688_2_MD)=rownames(V688_2@meta.data)
colnames(v688_2_MD)=gsub('$','_KoreanDeconv',colnames(v688_2_MD))
V688_2@meta.data<-cbind(V688_2@meta.data,v688_2_MD)


library(dplyr)
v015_1_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN123_A798015_Rep1')
v015_1_MD<-v015_1_MD[,2:39]
v015_1_MD<-v015_1_MD[v015_1_MD$ST.barcode%in%rownames(V015_1@meta.data),]
V015_1<-V015_1[,v015_1_MD$ST.barcode]
rownames(v015_1_MD)=rownames(V015_1@meta.data)
colnames(v015_1_MD)=gsub('$','_KoreanDeconv',colnames(v015_1_MD))
V015_1@meta.data<-cbind(V015_1@meta.data,v015_1_MD)


library(dplyr)
v015_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN124_A798015_Rep2')
v015_2_MD<-v015_2_MD[,2:39]
v015_2_MD<-v015_2_MD[v015_2_MD$ST.barcode%in%rownames(V015_2@meta.data),]
V015_2<-V015_2[,v015_2_MD$ST.barcode]
rownames(v015_2_MD)=rownames(V015_2@meta.data)
colnames(v015_2_MD)=gsub('$','_KoreanDeconv',colnames(v015_2_MD))
V015_2@meta.data<-cbind(V015_2@meta.data,v015_2_MD)

library(dplyr)
v797_1_MD<-deconvKorean%>%
  filter(ST.exp=='SN123_A938797_Rep1')
v797_1_MD<-v797_1_MD[,2:39]
v797_1_MD<-v797_1_MD[v797_1_MD$ST.barcode%in%rownames(V797_1@meta.data),]
V797_1<-V797_1[,v797_1_MD$ST.barcode]
rownames(v797_1_MD)=rownames(V797_1@meta.data)
colnames(v797_1_MD)=gsub('$','_KoreanDeconv',colnames(v797_1_MD))
V797_1@meta.data<-cbind(V797_1@meta.data,v797_1_MD)


library(dplyr)
v797_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN124_A938797_Rep2')
v797_2_MD<-v797_2_MD[,2:39]
v797_2_MD<-v797_2_MD[v797_2_MD$ST.barcode%in%rownames(V797_2@meta.data),]
V797_2<-V797_2[,v797_2_MD$ST.barcode]
rownames(v797_2_MD)=rownames(V797_2@meta.data)
colnames(v797_2_MD)=gsub('$','_KoreanDeconv',colnames(v797_2_MD))
V797_2@meta.data<-cbind(V797_2@meta.data,v797_2_MD)


library(dplyr)
v838_1_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN84_A120838_Rep1')
v838_1_MD<-v838_1_MD[,2:39]
v838_1_MD<-v838_1_MD[v838_1_MD$ST.barcode%in%rownames(V838_1@meta.data),]
V838_1<-V838_1[,v838_1_MD$ST.barcode]
rownames(v838_1_MD)=rownames(V838_1@meta.data)
colnames(v838_1_MD)=gsub('$','_KoreanDeconv',colnames(v838_1_MD))
V838_1@meta.data<-cbind(V838_1@meta.data,v838_1_MD)

library(dplyr)
v838_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN84_A120838_Rep2')
v838_2_MD<-v838_2_MD[,2:39]
v838_2_MD<-v838_2_MD[v838_2_MD$ST.barcode%in%rownames(V838_2@meta.data),]
V838_2<-V838_2[,v838_2_MD$ST.barcode]
rownames(v838_2_MD)=rownames(V838_2@meta.data)
colnames(v838_2_MD)=gsub('$','_KoreanDeconv',colnames(v838_2_MD))
V838_2@meta.data<-cbind(V838_2@meta.data,v838_2_MD)


{
  library(dplyr)
  v573_1_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN048_A121573_Rep1')
  v573_1_MD<-v573_1_MD[,2:43]
  v573_1_MD<-v573_1_MD[v573_1_MD$ST.barcode%in%rownames(V573_1@meta.data),]
  V573_1<-V573_1[,v573_1_MD$ST.barcode]
  rownames(v573_1_MD)=rownames(V573_1@meta.data)
  colnames(v573_1_MD)=gsub('$','_BelgianDeconv',colnames(v573_1_MD))
  V573_1@meta.data<-cbind(V573_1@meta.data,v573_1_MD)
  
  
  library(dplyr)
  v573_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN048_A121573_Rep2')
  v573_2_MD<-v573_2_MD[,2:43]
  v573_2_MD<-v573_2_MD[v573_2_MD$ST.barcode%in%rownames(V573_2@meta.data),]
  V573_2<-V573_2[,v573_2_MD$ST.barcode]
  rownames(v573_2_MD)=rownames(V573_2@meta.data)
  colnames(v573_2_MD)=gsub('$','_BelgianDeconv',colnames(v573_2_MD))
  V573_2@meta.data<-cbind(V573_2@meta.data,v573_2_MD)
  
  
  library(dplyr)
  v371_1_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN048_A416371_Rep1')
  v371_1_MD<-v371_1_MD[,2:43]
  v371_1_MD<-v371_1_MD[v371_1_MD$ST.barcode%in%rownames(V371_1@meta.data),]
  V371_1<-V371_1[,v371_1_MD$ST.barcode]
  rownames(v371_1_MD)=rownames(V371_1@meta.data)
  colnames(v371_1_MD)=gsub('$','_BelgianDeconv',colnames(v371_1_MD))
  V371_1@meta.data<-cbind(V371_1@meta.data,v371_1_MD)
  
  
  library(dplyr)
  v371_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN048_A416371_Rep2')
  v371_2_MD<-v371_2_MD[,2:43]
  v371_2_MD<-v371_1_MD[v371_2_MD$ST.barcode%in%rownames(V371_2@meta.data),]
  V371_2<-V371_1[,v371_2_MD$ST.barcode]
  rownames(v371_2_MD)=rownames(V371_2@meta.data)
  colnames(v371_2_MD)=gsub('$','_BelgianDeconv',colnames(v371_2_MD))
  V371_2@meta.data<-cbind(V371_2@meta.data,v371_2_MD)
  
  
  library(dplyr)
  v763_1_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN123_A551763_Rep1')
  v763_1_MD<-v763_1_MD[,2:43]
  v763_1_MD<-v763_1_MD[v763_1_MD$ST.barcode%in%rownames(V763_1@meta.data),]
  V763_1<-V763_1[,v763_1_MD$ST.barcode]
  rownames(v763_1_MD)=rownames(V763_1@meta.data)
  colnames(v763_1_MD)=gsub('$','_BelgianDeconv',colnames(v763_1_MD))
  V763_1@meta.data<-cbind(V763_1@meta.data,v763_1_MD)
  
  
  library(dplyr)
  v763_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN124_A551763_Rep2')
  v763_2_MD<-v763_2_MD[,2:43]
  v763_2_MD<-v763_2_MD[v763_2_MD$ST.barcode%in%rownames(V763_2@meta.data),]
  V763_2<-V763_2[,v763_2_MD$ST.barcode]
  rownames(v763_2_MD)=rownames(V763_2@meta.data)
  colnames(v763_2_MD)=gsub('$','_BelgianDeconv',colnames(v763_2_MD))
  V763_2@meta.data<-cbind(V763_2@meta.data,v763_2_MD)
  
  
  library(dplyr)
  v688_1_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN123_A595688_Rep1')
  v688_1_MD<-v688_1_MD[,2:43]
  v688_1_MD<-v688_1_MD[v688_1_MD$ST.barcode%in%rownames(V688_1@meta.data),]
  V688_1<-V688_1[,v688_1_MD$ST.barcode]
  rownames(v688_1_MD)=rownames(V688_1@meta.data)
  colnames(v688_1_MD)=gsub('$','_BelgianDeconv',colnames(v688_1_MD))
  V688_1@meta.data<-cbind(V688_1@meta.data,v688_1_MD)
  
  
  library(dplyr)
  v688_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN124_A595688_Rep2')
  v688_2_MD<-v688_2_MD[,2:43]
  v688_2_MD<-v688_2_MD[v688_2_MD$ST.barcode%in%rownames(V688_2@meta.data),]
  V688_2<-V688_2[,v688_2_MD$ST.barcode]
  rownames(v688_2_MD)=rownames(V688_2@meta.data)
  colnames(v688_2_MD)=gsub('$','_BelgianDeconv',colnames(v688_2_MD))
  V688_2@meta.data<-cbind(V688_2@meta.data,v688_2_MD)
  
  
  library(dplyr)
  v015_1_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN123_A798015_Rep1')
  v015_1_MD<-v015_1_MD[,2:43]
  v015_1_MD<-v015_1_MD[v015_1_MD$ST.barcode%in%rownames(V015_1@meta.data),]
  V015_1<-V015_1[,v015_1_MD$ST.barcode]
  rownames(v015_1_MD)=rownames(V015_1@meta.data)
  colnames(v015_1_MD)=gsub('$','_BelgianDeconv',colnames(v015_1_MD))
  V015_1@meta.data<-cbind(V015_1@meta.data,v015_1_MD)
  
  
  library(dplyr)
  v015_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN124_A798015_Rep2')
  v015_2_MD<-v015_2_MD[,2:43]
  v015_2_MD<-v015_2_MD[v015_2_MD$ST.barcode%in%rownames(V015_2@meta.data),]
  V015_2<-V015_2[,v015_2_MD$ST.barcode]
  rownames(v015_2_MD)=rownames(V015_2@meta.data)
  colnames(v015_2_MD)=gsub('$','_BelgianDeconv',colnames(v015_2_MD))
  V015_2@meta.data<-cbind(V015_2@meta.data,v015_2_MD)
  
  library(dplyr)
  v797_1_MD<-deconvBelgian%>%
    filter(ST.exp=='SN123_A938797_Rep1')
  v797_1_MD<-v797_1_MD[,2:43]
  v797_1_MD<-v797_1_MD[v797_1_MD$ST.barcode%in%rownames(V797_1@meta.data),]
  V797_1<-V797_1[,v797_1_MD$ST.barcode]
  rownames(v797_1_MD)=rownames(V797_1@meta.data)
  colnames(v797_1_MD)=gsub('$','_BelgianDeconv',colnames(v797_1_MD))
  V797_1@meta.data<-cbind(V797_1@meta.data,v797_1_MD)
  
  
  library(dplyr)
  v797_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN124_A938797_Rep2')
  v797_2_MD<-v797_2_MD[,2:43]
  v797_2_MD<-v797_2_MD[v797_2_MD$ST.barcode%in%rownames(V797_2@meta.data),]
  V797_2<-V797_2[,v797_2_MD$ST.barcode]
  rownames(v797_2_MD)=rownames(V797_2@meta.data)
  colnames(v797_2_MD)=gsub('$','_BelgianDeconv',colnames(v797_2_MD))
  V797_2@meta.data<-cbind(V797_2@meta.data,v797_2_MD)
  
  
  library(dplyr)
  v838_1_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN84_A120838_Rep1')
  v838_1_MD<-v838_1_MD[,2:43]
  v838_1_MD<-v838_1_MD[v838_1_MD$ST.barcode%in%rownames(V838_1@meta.data),]
  V838_1<-V838_1[,v838_1_MD$ST.barcode]
  rownames(v838_1_MD)=rownames(V838_1@meta.data)
  colnames(v838_1_MD)=gsub('$','_BelgianDeconv',colnames(v838_1_MD))
  V838_1@meta.data<-cbind(V838_1@meta.data,v838_1_MD)
  
  library(dplyr)
  v838_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN84_A120838_Rep2')
  v838_2_MD<-v838_2_MD[,2:43]
  v838_2_MD<-v838_2_MD[v838_2_MD$ST.barcode%in%rownames(V838_2@meta.data),]
  V838_2<-V838_2[,v838_2_MD$ST.barcode]
  rownames(v838_2_MD)=rownames(V838_2@meta.data)
  colnames(v838_2_MD)=gsub('$','_BelgianDeconv',colnames(v838_2_MD))
  V838_2@meta.data<-cbind(V838_2@meta.data,v838_2_MD)
}


#We use the Korean deconvolution as a major groundtruth 

process_C2L_Korean<-function(seurat_object,df=FALSE){
  if(df==FALSE){
    mt <- seurat_object@meta.data
    cell_type_columns <- mt[, stringr::str_detect(pattern = "^(?!ST\\.exp_)(?!ST\\.barcode_).+_KoreanDeconv$", string = colnames(mt))]
    row_totals <- rowSums(cell_type_columns)
    cell_type_percentages <- sweep(cell_type_columns, 1, row_totals, "/") * 100
    colnames(cell_type_percentages) <- paste0(colnames(cell_type_percentages), "_percentageKorean")
    mt <- cbind(mt, cell_type_percentages)
    seurat_object@meta.data <- mt
    return(seurat_object)
  }else{
    cell_type_columns <- seurat_object[, stringr::str_detect(pattern = "^(?!spot.id_)(?!ST\\.exp_)(?!ST\\.barcode_).+_KoreanDeconv$", string = colnames(seurat_object))]
    row_totals <- rowSums(cell_type_columns)
    cell_type_percentages <- sweep(cell_type_columns, 1, row_totals, "/") * 100
    colnames(cell_type_percentages) <- paste0(colnames(cell_type_percentages), "_percentageKorean")
    seurat_object <- cbind(seurat_object, cell_type_percentages)
    return(seurat_object)
  }
}

V015_1<-process_C2L_Korean(V015_1)
V015_2<-process_C2L_Korean(V015_2)

V371_1<-process_C2L_Korean(V371_1)
V371_2<-process_C2L_Korean(V371_2)

V573_1<-process_C2L_Korean(V573_1)
V573_2<-process_C2L_Korean(V573_2)

V688_1<-process_C2L_Korean(V688_1)
V688_2<-process_C2L_Korean(V688_2)

V763_1<-process_C2L_Korean(V763_1)
V763_2<-process_C2L_Korean(V763_2)

V797_1<-process_C2L_Korean(V797_1)
V797_2<-process_C2L_Korean(V797_2)

V838_1<-process_C2L_Korean(V838_1)
V838_2<-process_C2L_Korean(V838_2)

colnames(deconvKorean)=gsub('$','_KoreanDeconv',colnames(deconvKorean))
deconvKorean<-process_C2L_Korean(deconvKorean,df = T)


V015_1<-NormalizeData(V015_1)
V015_2<-NormalizeData(V015_2)
V371_1<-NormalizeData(V371_1)
V371_2<-NormalizeData(V371_2)
V573_1<-NormalizeData(V573_1)
V573_2<-NormalizeData(V573_2)
V688_1<-NormalizeData(V688_1)
V688_2<-NormalizeData(V688_2)
V763_1<-NormalizeData(V763_1)
V763_2<-NormalizeData(V763_2)
V797_1<-NormalizeData(V797_1)
V797_2<-NormalizeData(V797_2)
V838_1<-NormalizeData(V838_1)
V838_2<-NormalizeData(V838_2)




rownames(V015_1@assays$Spatial@features)->rownames(V015_1@assays$Spatial@layers$counts)
rownames(V015_1@assays$Spatial@features)->rownames(V015_1@assays$Spatial@layers$data)
rownames(V015_2@assays$Spatial@features)->rownames(V015_2@assays$Spatial@layers$counts)
rownames(V015_2@assays$Spatial@features)->rownames(V015_2@assays$Spatial@layers$data)
rownames(V371_1@assays$Spatial@features)->rownames(V371_1@assays$Spatial@layers$counts)
rownames(V371_1@assays$Spatial@features)->rownames(V371_1@assays$Spatial@layers$data)
rownames(V371_2@assays$Spatial@features)->rownames(V371_2@assays$Spatial@layers$counts)
rownames(V371_2@assays$Spatial@features)->rownames(V371_2@assays$Spatial@layers$data)
rownames(V573_1@assays$Spatial@features)->rownames(V573_1@assays$Spatial@layers$counts)
rownames(V573_1@assays$Spatial@features)->rownames(V573_1@assays$Spatial@layers$data)
rownames(V573_2@assays$Spatial@features)->rownames(V573_2@assays$Spatial@layers$counts)
rownames(V573_2@assays$Spatial@features)->rownames(V573_2@assays$Spatial@layers$data)
rownames(V688_1@assays$Spatial@features)->rownames(V688_1@assays$Spatial@layers$counts)
rownames(V688_1@assays$Spatial@features)->rownames(V688_1@assays$Spatial@layers$data)
rownames(V688_2@assays$Spatial@features)->rownames(V688_2@assays$Spatial@layers$counts)
rownames(V688_2@assays$Spatial@features)->rownames(V688_2@assays$Spatial@layers$data)
rownames(V763_1@assays$Spatial@features)->rownames(V763_1@assays$Spatial@layers$counts)
rownames(V763_1@assays$Spatial@features)->rownames(V763_1@assays$Spatial@layers$data)
rownames(V763_2@assays$Spatial@features)->rownames(V763_2@assays$Spatial@layers$counts)
rownames(V763_2@assays$Spatial@features)->rownames(V763_2@assays$Spatial@layers$data)
rownames(V797_1@assays$Spatial@features)->rownames(V797_1@assays$Spatial@layers$counts)
rownames(V797_1@assays$Spatial@features)->rownames(V797_1@assays$Spatial@layers$data)
rownames(V797_2@assays$Spatial@features)->rownames(V797_2@assays$Spatial@layers$counts)
rownames(V797_2@assays$Spatial@features)->rownames(V797_2@assays$Spatial@layers$data)
rownames(V838_1@assays$Spatial@features)->rownames(V838_1@assays$Spatial@layers$counts)
rownames(V838_1@assays$Spatial@features)->rownames(V838_1@assays$Spatial@layers$data)
rownames(V838_2@assays$Spatial@features)->rownames(V838_2@assays$Spatial@layers$counts)
rownames(V838_2@assays$Spatial@features)->rownames(V838_2@assays$Spatial@layers$data)



#We adopt the same gating strategy as we did for Wang Long et al. - 20% and 0.85 quantile 

deconvKorean1<-deconvKorean[deconvKorean$Myofibroblasts_KoreanDeconv+deconvKorean$Stromal.1_KoreanDeconv+deconvKorean$Stromal.2_KoreanDeconv+deconvKorean$Stromal.3_KoreanDeconv>quantile(deconvKorean$Myofibroblasts_KoreanDeconv+deconvKorean$Stromal.1_KoreanDeconv+deconvKorean$Stromal.2_KoreanDeconv+deconvKorean$Stromal.3_KoreanDeconv,0.85),]
deconvKorean1<-deconvKorean1[deconvKorean1$Myofibroblasts_KoreanDeconv_percentageKorean+deconvKorean1$Stromal.1_KoreanDeconv_percentageKorean+deconvKorean1$Stromal.2_KoreanDeconv_percentageKorean+deconvKorean1$Stromal.3_KoreanDeconv_percentageKorean>0.2,]


deconvKorean1<-deconvKorean1%>%
  nest_by(ST.exp_KoreanDeconv)


V573_1_fibro<-V573_1[,deconvKorean1$data[[1]]$ST.barcode_KoreanDeconv]
V573_2_fibro<-V573_2[,deconvKorean1$data[[2]]$ST.barcode_KoreanDeconv]


V371_1_fibro<-V371_1[,deconvKorean1$data[[3]]$ST.barcode_KoreanDeconv]
V371_2_fibro<-V371_2[,deconvKorean1$data[[4]]$ST.barcode_KoreanDeconv]


V763_1_fibro<-V763_1[,deconvKorean1$data[[5]]$ST.barcode_KoreanDeconv]
V763_2_fibro<-V763_2[,deconvKorean1$data[[8]]$ST.barcode_KoreanDeconv]


V688_1_fibro<-V688_1[,deconvKorean1$data[[6]]$ST.barcode_KoreanDeconv]
V688_2_fibro<-V688_2[,deconvKorean1$data[[9]]$ST.barcode_KoreanDeconv]


V015_1_fibro<-V015_1[,deconvKorean1$data[[7]]$ST.barcode_KoreanDeconv]
V015_2_fibro<-V015_2[,deconvKorean1$data[[10]]$ST.barcode_KoreanDeconv]


V797_1_fibro<-V797_1[,deconvKorean1$data[[14]]$ST.barcode_KoreanDeconv]
V797_2_fibro<-V797_2[,deconvKorean1$data[[11]]$ST.barcode_KoreanDeconv]


V838_1_fibro<-V838_1[,deconvKorean1$data[[12]]$ST.barcode_KoreanDeconv]
V838_2_fibro<-V838_2[,deconvKorean1$data[[13]]$ST.barcode_KoreanDeconv]


load("~/Report_OSR_Candiolo_2024/data/signfinal.Robj")
sign$EMR=NULL
sign$iCAF_Cords=NULL
sign$mCAF_Cords=NULL
sign$dCAF_Cords=NULL
sign$apCAF_Cords=NULL
sign$pericyte_Cords=NULL
sign$vCAFrCAF_Cords=NULL
sign$tCAF_Cords=NULL
sign$mrCAF=NULL

library(reticulate)
use_python('/home/carlo/.local/share/r-miniconda/envs/giotto_env/bin/python', required = TRUE)

library(circlize)
library(data.table)
process_PAGEanalysis1 <- function(seurat_object, signatures_all, selected_signatures,only_fibro=T) {
  if(only_fibro==T){signatures_all <- signatures_all[!names(signatures_all) %in% c("CCM", "EMRM")]}
  raw_exprs <- seurat_object@assays$Spatial@layers$counts
  rownames(raw_exprs) <- rownames(seurat_object@assays$Spatial@features@.Data)
  colnames(raw_exprs) <- rownames(seurat_object@meta.data)
  spatial_locs <- as.data.table(GetTissueCoordinates(seurat_object)[,1:2])
  colnames(spatial_locs) <- c("x", "y")
  
  myGiottoObj <- createGiottoObject(raw_exprs = raw_exprs, spatial_locs = spatial_locs)
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
    min_overlap_genes = 2,
    output_enrichment = c("original", "zscore")
  )
  
  # Run PAGE enrichment analysis with p-values
  myGiottoObj_pval <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 2,
    p_value = TRUE,
    reverse_log_scale = FALSE
  )
  
  # Extract p-value results
  pval_df <- as.data.frame(myGiottoObj_pval@spatial_enrichment$PAGE)
  zscore_df<-as.data.frame(myGiottoObj_initial@spatial_enrichment$PAGE)
  
  zscore_df <- as.data.frame(zscore_df)
  colnames(zscore_df) <- paste0(colnames(zscore_df), "PAGE")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, zscore_df)
  
  # Identify significantly enriched spots - transformed values corresponds to pval < 0.05
  significant_spots <- pval_df > 1.301
  
  # Store significant enrichment information in Seurat metadata
  significant_df <- as.data.frame(significant_spots)
  colnames(significant_df) <- paste0(colnames(significant_df), "PAGE_significant")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, significant_df)
  
  return(seurat_object)
  
}


library(Giotto)

#we are going to perform a PAGE enrichment test and a permutation test on a subset of the original dataset.
#if the dataset is too small, permutation test may fail

#we are going to consider results that succeed in the test (indicated by #!! )

#we acknowledge this is the crucial part of the analysis. Even though trying with different cutoffs to select more fibroblasts-enriched spots, the result does not change

#V573_1_fibro <- process_PAGEanalysis1(seurat_object = V573_1_fibro, signatures_all = sign)
#V573_2_fibro <- process_PAGEanalysis1(seurat_object = V573_2_fibro, signatures_all = sign)
#V371_1_fibro <- process_PAGEanalysis1(seurat_object = V371_1_fibro, signatures_all = sign)
#V371_2_fibro <- process_PAGEanalysis1(seurat_object = V371_2_fibro, signatures_all = sign) 
#V763_1_fibro <- process_PAGEanalysis1(seurat_object = V763_1_fibro, signatures_all = sign) 
#V763_2_fibro <- process_PAGEanalysis1(seurat_object = V763_2_fibro, signatures_all = sign)
V688_1_fibro <- process_PAGEanalysis1(seurat_object = V688_1_fibro, signatures_all = sign) #!!
#V688_2_fibro <- process_PAGEanalysis1(seurat_object = V688_2_fibro, signatures_all = sign) 
V015_1_fibro <- process_PAGEanalysis1(seurat_object = V015_1_fibro, signatures_all = sign) #!!
V015_2_fibro <- process_PAGEanalysis1(seurat_object = V015_2_fibro, signatures_all = sign) #!!
V797_1_fibro <- process_PAGEanalysis1(seurat_object = V797_1_fibro, signatures_all = sign) #!!
V797_2_fibro <- process_PAGEanalysis1(seurat_object = V797_2_fibro, signatures_all = sign) #!!
#V838_1_fibro <- process_PAGEanalysis1(seurat_object = V838_1_fibro, signatures_all = sign) 
V838_2_fibro <- process_PAGEanalysis1(seurat_object = V838_2_fibro, signatures_all = sign) #!!


GSM7058756_C1_fibroenriched@images$slice1@boundaries$centroids@cells=rownames(GSM7058756_C1_fibroenriched@meta.data)
GSM7058757_C2_fibroenriched@images$slice1@boundaries$centroids@cells=rownames(GSM7058757_C2_fibroenriched@meta.data)
GSM7058758_C3_fibroenriched@images$slice1@boundaries$centroids@cells=rownames(GSM7058758_C3_fibroenriched@meta.data)
GSM7058759_C4_fibroenriched@images$slice1@boundaries$centroids@cells=rownames(GSM7058759_C4_fibroenriched@meta.data)


#function change slightly for Wang Long et al. data 
process_PAGEanalysis1 <- function(seurat_object, signatures_all, selected_signatures,only_fibro=T) {
  if(only_fibro==T){signatures_all <- signatures_all[!names(signatures_all) %in% c("CCM", "EMRM")]}
  raw_exprs <- seurat_object@assays$RNA@layers$counts
  rownames(raw_exprs) <- rownames(seurat_object@assays$RNA@features@.Data)
  colnames(raw_exprs) <- rownames(seurat_object@meta.data)
  spatial_locs <- as.data.table(GetTissueCoordinates(seurat_object)[,1:2])
  colnames(spatial_locs) <- c("x", "y")
  
  myGiottoObj <- createGiottoObject(raw_exprs = raw_exprs, spatial_locs = spatial_locs)
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
    min_overlap_genes = 2,
    output_enrichment = c("original", "zscore")
  )
  
  # Run PAGE enrichment analysis with p-values
  myGiottoObj_pval <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 2,
    p_value = TRUE,
    reverse_log_scale = FALSE
  )
  
  # Extract p-value results
  pval_df <- as.data.frame(myGiottoObj_pval@spatial_enrichment$PAGE)
  zscore_df<-as.data.frame(myGiottoObj_initial@spatial_enrichment$PAGE)
  
  zscore_df <- as.data.frame(zscore_df)
  colnames(zscore_df) <- paste0(colnames(zscore_df), "PAGE")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, zscore_df)
  
  # Identify significantly enriched spots (p-value < 0.05)
  significant_spots <- pval_df > 1.301
  
  # Store significant enrichment information in Seurat metadata
  significant_df <- as.data.frame(significant_spots)
  colnames(significant_df) <- paste0(colnames(significant_df), "PAGE_significant")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, significant_df)
  
  return(seurat_object)
  
}
GSM7058756_C1_fibroenriched<-process_PAGEanalysis1(seurat_object = GSM7058756_C1_fibroenriched,signatures_all = sign)
GSM7058757_C2_fibroenriched<-process_PAGEanalysis1(seurat_object = GSM7058757_C2_fibroenriched,signatures_all = sign)
#GSM7058758_C3_fibroenriched<-process_PAGEanalysis1(seurat_object = GSM7058758_C3_fibroenriched,signatures_all = sign)
GSM7058759_C4_fibroenriched<-process_PAGEanalysis1(seurat_object = GSM7058759_C4_fibroenriched,signatures_all = sign)

#let's extract only PAGE values 

V688_1_fibro_md<- V688_1_fibro@meta.data[,121:130]
V015_1_fibro_md<- V015_1_fibro@meta.data[,121:130]
V015_2_fibro_md<- V015_2_fibro@meta.data[,121:130]
V797_1_fibro_md<-V797_1_fibro@meta.data[,121:130]
V797_2_fibro_md<-V797_2_fibro@meta.data[,121:130]
V838_2_fibro_md<-V838_2_fibro@meta.data[,121:130]

GSM7058756_C1_fibroenriched_md<-GSM7058756_C1_fibroenriched@meta.data[,77:86]
GSM7058757_C2_fibroenriched_md<-GSM7058757_C2_fibroenriched@meta.data[,77:86]
GSM7058759_C4_fibroenriched_md<-GSM7058759_C4_fibroenriched@meta.data[,77:86]


# define the Fisher's z-transformation and its inverse
fisher_z <- function(r) {
  return(0.5 * log((1 + r) / (1 - r)))
}

inverse_fisher_z <- function(z) {
  return((exp(2 * z) - 1) / (exp(2 * z) + 1))
}

# list of correlation matrices and corresponding number of spots
correlation_matrices <- list(cor(V688_1_fibro_md), 
                             cor(V015_1_fibro_md), 
                             cor(V015_2_fibro_md), 
                             cor(V797_1_fibro_md),
                             cor(V797_2_fibro_md),
                             cor(V838_2_fibro_md),
                             cor(GSM7058756_C1_fibroenriched_md),
                             cor(GSM7058757_C2_fibroenriched_md),
                             cor(GSM7058759_C4_fibroenriched_md)) 

num_spots <- c(dim(V688_1_fibro_md)[1], 
               dim(V015_1_fibro_md)[1], 
               dim(V015_2_fibro_md)[1], 
               dim(V797_1_fibro_md)[1],
               dim(V797_2_fibro_md)[1], 
               dim(V838_2_fibro_md)[1],
               dim(GSM7058756_C1_fibroenriched_md)[1],
               dim(GSM7058757_C2_fibroenriched_md)[1],
               dim(GSM7058759_C4_fibroenriched_md)[1]
) 

# Fisher's z-transformation to each matrix
z_matrices <- lapply(correlation_matrices, function(mat) {
  apply(mat, c(1, 2), fisher_z)
})

# Weight the transformed values by the number of spots
weighted_z_matrices <- mapply(function(z_matrix, n_spots) {
  z_matrix * n_spots
}, z_matrices, num_spots, SIMPLIFY = FALSE)

# Compute the sum of weighted z-values and total number of spots
sum_weighted_z_matrix <- Reduce("+", weighted_z_matrices)
total_spots <- sum(num_spots)

# Compute the weighted average of z-values
mean_weighted_z_matrix <- sum_weighted_z_matrix / total_spots

# Apply the inverse Fisher's z-transformation to the averaged matrix
merged_correlation_matrix <- apply(mean_weighted_z_matrix, c(1, 2), inverse_fisher_z)

# Ensure the diagonal remains 1 (since it's a correlation matrix)
diag(merged_correlation_matrix) <- 1

# Perform hierarchical clustering
row_model <- hclust(dist(merged_correlation_matrix))
col_model <- hclust(t(dist(merged_correlation_matrix)))



GSM7058757_C2$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',GSM7058757_C2$Stromal_C2L_percentage>GSM7058757_C2$Tumoral_C2L_percentage)))

p3<-SpatialFeaturePlot(GSM7058757_C2_fibroenriched,features = 'MF1PAGE',image.alpha = 0)+ggtitle('MF1 PAGE Enrichment - GSM7058757_C2')
p4<-SpatialFeaturePlot(GSM7058757_C2_fibroenriched,features = 'iCAFPAGE',image.alpha = 0)+ggtitle('iCAF PAGE Enrichment - GSM7058757_C2')
p5<-SpatialFeaturePlot(GSM7058757_C2_fibroenriched,features = 'mCAFPAGE',image.alpha = 0)+ggtitle('mCAF PAGE Enrichment - GSM7058757_C2')
p6<-SpatialFeaturePlot(GSM7058757_C2_fibroenriched,features = 'vCAFPAGE',image.alpha = 0)+ggtitle('vCAF PAGE Enrichment - GSM7058757_C2')
p7<-SpatialDimPlot(GSM7058757_C2,group.by = 'Tissue_Classification',image.alpha = 0)

save(p3,file = 'p3.Robj')
save(p4,file = 'p4.Robj')
save(p5,file = 'p5.Robj')
save(p6,file = 'p6.Robj')
save(p7,file = 'p7.Robj')


library(ComplexHeatmap)
# Create a Heatmap object with row and column dendrograms
p8 <- Heatmap(merged_correlation_matrix, 
              name = "Merged\nCorrelation\ncorrected for n.spots",
              col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
              cluster_rows = row_model,
              cluster_columns = col_model,
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_dend_side = "left",
              column_dend_side = "top",
              row_names_side = "left",
              column_names_side = "bottom",
              heatmap_legend_param = list(title = "Merged\nCorrelation\ncorrected\nfor n.spots"),
              row_dend_width = unit(4, "cm"),
              column_dend_height = unit(4, "cm"))

save(p8,file = 'p8.Robj')


  
V688_1_fibro_pval<- V688_1_fibro@meta.data[,132:141]
V015_1_fibro_pval<- V015_1_fibro@meta.data[,132:141]
V015_2_fibro_pval<- V015_2_fibro@meta.data[,132:141]
V797_1_fibro_pval<-V797_1_fibro@meta.data[,132:141]
V797_2_fibro_pval<-V797_2_fibro@meta.data[,132:141]
V838_2_fibro_pval<-V838_2_fibro@meta.data[,132:141]

GSM7058756_C1_fibroenriched_pval<-GSM7058756_C1_fibroenriched@meta.data[,88:97]
GSM7058757_C2_fibroenriched_pval<-GSM7058757_C2_fibroenriched@meta.data[,88:97]
GSM7058759_C4_fibroenriched_pval<-GSM7058759_C4_fibroenriched@meta.data[,88:97]

calculate_hypergeometric_pvals <- function(df) {
  categories <- colnames(df)
  comb_matrix <- combn(categories, 2, simplify = FALSE)
  
  pval_matrix <- matrix(0, nrow = length(categories), ncol = length(categories))
  colnames(pval_matrix) <- categories
  rownames(pval_matrix) <- categories
  
  N <- nrow(df)
  
  for (comb in comb_matrix) {
    cat1 <- comb[1]
    cat2 <- comb[2]
    
    K <- sum(df[[cat1]])
    n <- sum(df[[cat2]])
    k <- sum(df[[cat1]] & df[[cat2]])
    
    pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    
    pval_matrix[cat1, cat2] <- pval
    pval_matrix[cat2, cat1] <- pval
  }
  
  diag(pval_matrix) <- 1  # Diagonal elements should be 1 (or another suitable value, as self-comparisons are not meaningful)
  
  return(pval_matrix)
}


hyp_V688_1<-calculate_hypergeometric_pvals(V688_1_fibro_pval)
hyp_V015_1<-calculate_hypergeometric_pvals(V015_1_fibro_pval)
hyp_V015_2<-calculate_hypergeometric_pvals(V015_2_fibro_pval)
hyp_V797_1<-calculate_hypergeometric_pvals(V797_1_fibro_pval)
hyp_V797_2<-calculate_hypergeometric_pvals(V797_2_fibro_pval)
hyp_V838_2<-calculate_hypergeometric_pvals(V838_2_fibro_pval)
hyp_C1<-calculate_hypergeometric_pvals(GSM7058756_C1_fibroenriched_pval)
hyp_C2<-calculate_hypergeometric_pvals(GSM7058757_C2_fibroenriched_pval)
hyp_C4<-calculate_hypergeometric_pvals(GSM7058759_C4_fibroenriched_pval)


# Function to combine p-values using Fisher's method
combine_pvals <- function(pval_matrices) {
  categories <- rownames(pval_matrices[[1]])
  combined_pval_matrix <- matrix(0, nrow = length(categories), ncol = length(categories))
  colnames(combined_pval_matrix) <- categories
  rownames(combined_pval_matrix) <- categories
  
  for (i in 1:length(categories)) {
    for (j in 1:length(categories)) {
      pvals <- sapply(pval_matrices, function(x) x[i, j])
      if (all(pvals == 1)) {
        combined_pval <- 1
      } else {
        chisq_stat <- -2 * sum(log(pvals))
        combined_pval <- pchisq(chisq_stat, df = 2 * length(pvals), lower.tail = FALSE)
      }
      combined_pval_matrix[i, j] <- combined_pval
    }
  }
  
  return(combined_pval_matrix)
}


# Combine p-values from all datasets
combined_pval_matrix <- combine_pvals(list(hyp_V688_1,
                                           hyp_V015_1,
                                           hyp_V015_2,
                                           hyp_V797_1,
                                           hyp_V797_2,
                                           hyp_V838_2,
                                           hyp_C1,
                                           hyp_C2,
                                           hyp_C4))
combined_pval_df <- as.data.frame(combined_pval_matrix)
combined_pval_df


# Perform hierarchical clustering
row_model <- hclust(dist(combined_pval_df))
col_model <- hclust(t(dist(combined_pval_df)))

library(ComplexHeatmap)
# Create a Heatmap object with row and column dendrograms
p9 <- Heatmap(-log(combined_pval_df+5e-324), 
              name = "Combined\nenrichment\nscore\nvia Fisher's\nMethod",
              col = colorRamp2(c(0, 200, 600), c("white","#E94B3C",'#6C1413')),
              cluster_rows = row_model,
              cluster_columns = col_model,
              show_row_dend = TRUE,
              show_column_dend = T,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_dend_side = "left",
              column_dend_side = "top",
              row_names_side = "left",
              column_names_side = "bottom",
              heatmap_legend_param = list(title = "Combined\nenrichment\nscore\nvia Fisher's\nMethod"),
              row_dend_width = unit(3, "cm"),
              column_dend_height = unit(4, "cm"))

save(p9,file = 'p9.Robj')


#Now we have to do the entire analysis again with more signatures for supplementary data 


load("~/Report_OSR_Candiolo_2024/data/signfinal.Robj")
sign$EMR=NULL


process_PAGEanalysis1 <- function(seurat_object, signatures_all, selected_signatures,only_fibro=T) {
  if(only_fibro==T){signatures_all <- signatures_all[!names(signatures_all) %in% c("CCM", "EMRM")]}
  raw_exprs <- seurat_object@assays$Spatial@layers$counts
  rownames(raw_exprs) <- rownames(seurat_object@assays$Spatial@features@.Data)
  colnames(raw_exprs) <- rownames(seurat_object@meta.data)
  spatial_locs <- as.data.table(GetTissueCoordinates(seurat_object)[,1:2])
  colnames(spatial_locs) <- c("x", "y")
  
  myGiottoObj <- createGiottoObject(raw_exprs = raw_exprs, spatial_locs = spatial_locs)
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
    min_overlap_genes = 2,
    output_enrichment = c("original", "zscore")
  )
  
  # Run PAGE enrichment analysis with p-values
  myGiottoObj_pval <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 2,
    p_value = TRUE,
    reverse_log_scale = FALSE
  )
  
  # Extract p-value results
  pval_df <- as.data.frame(myGiottoObj_pval@spatial_enrichment$PAGE)
  zscore_df<-as.data.frame(myGiottoObj_initial@spatial_enrichment$PAGE)
  
  zscore_df <- as.data.frame(zscore_df)
  colnames(zscore_df) <- paste0(colnames(zscore_df), "PAGE")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, zscore_df)
  
  # Identify significantly enriched spots - transformed values corresponds to pval < 0.05
  significant_spots <- pval_df > 1.301
  
  # Store significant enrichment information in Seurat metadata
  significant_df <- as.data.frame(significant_spots)
  colnames(significant_df) <- paste0(colnames(significant_df), "PAGE_significant")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, significant_df)
  
  return(seurat_object)
  
}
#V573_1_fibro <- process_PAGEanalysis1(seurat_object = V573_1_fibro, signatures_all = sign)
#V573_2_fibro <- process_PAGEanalysis1(seurat_object = V573_2_fibro, signatures_all = sign)
#V371_1_fibro <- process_PAGEanalysis1(seurat_object = V371_1_fibro, signatures_all = sign)
#V371_2_fibro <- process_PAGEanalysis1(seurat_object = V371_2_fibro, signatures_all = sign) 
#V763_1_fibro <- process_PAGEanalysis1(seurat_object = V763_1_fibro, signatures_all = sign) 
#V763_2_fibro <- process_PAGEanalysis1(seurat_object = V763_2_fibro, signatures_all = sign)
V688_1_fibro <- process_PAGEanalysis1(seurat_object = V688_1_fibro, signatures_all = sign) #!!
#V688_2_fibro <- process_PAGEanalysis1(seurat_object = V688_2_fibro, signatures_all = sign) 
V015_1_fibro <- process_PAGEanalysis1(seurat_object = V015_1_fibro, signatures_all = sign) #!!
V015_2_fibro <- process_PAGEanalysis1(seurat_object = V015_2_fibro, signatures_all = sign) #!!
V797_1_fibro <- process_PAGEanalysis1(seurat_object = V797_1_fibro, signatures_all = sign) #!!
V797_2_fibro <- process_PAGEanalysis1(seurat_object = V797_2_fibro, signatures_all = sign) #!!
#V838_1_fibro <- process_PAGEanalysis1(seurat_object = V838_1_fibro, signatures_all = sign) 
V838_2_fibro <- process_PAGEanalysis1(seurat_object = V838_2_fibro, signatures_all = sign) #!!


GSM7058756_C1_fibroenriched@images$slice1@boundaries$centroids@cells=rownames(GSM7058756_C1_fibroenriched@meta.data)
GSM7058757_C2_fibroenriched@images$slice1@boundaries$centroids@cells=rownames(GSM7058757_C2_fibroenriched@meta.data)
GSM7058758_C3_fibroenriched@images$slice1@boundaries$centroids@cells=rownames(GSM7058758_C3_fibroenriched@meta.data)
GSM7058759_C4_fibroenriched@images$slice1@boundaries$centroids@cells=rownames(GSM7058759_C4_fibroenriched@meta.data)

process_PAGEanalysis1 <- function(seurat_object, signatures_all, selected_signatures,only_fibro=T) {
  if(only_fibro==T){signatures_all <- signatures_all[!names(signatures_all) %in% c("CCM", "EMRM")]}
  raw_exprs <- seurat_object@assays$RNA@layers$counts
  rownames(raw_exprs) <- rownames(seurat_object@assays$RNA@features@.Data)
  colnames(raw_exprs) <- rownames(seurat_object@meta.data)
  spatial_locs <- as.data.table(GetTissueCoordinates(seurat_object)[,1:2])
  colnames(spatial_locs) <- c("x", "y")
  
  myGiottoObj <- createGiottoObject(raw_exprs = raw_exprs, spatial_locs = spatial_locs)
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
    min_overlap_genes = 2,
    output_enrichment = c("original", "zscore")
  )
  
  # Run PAGE enrichment analysis with p-values
  myGiottoObj_pval <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 2,
    p_value = TRUE,
    reverse_log_scale = FALSE
  )
  
  # Extract p-value results
  pval_df <- as.data.frame(myGiottoObj_pval@spatial_enrichment$PAGE)
  zscore_df<-as.data.frame(myGiottoObj_initial@spatial_enrichment$PAGE)
  
  zscore_df <- as.data.frame(zscore_df)
  colnames(zscore_df) <- paste0(colnames(zscore_df), "PAGE")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, zscore_df)
  
  # Identify significantly enriched spots (p-value < 0.05)
  significant_spots <- pval_df > 1.301
  
  # Store significant enrichment information in Seurat metadata
  significant_df <- as.data.frame(significant_spots)
  colnames(significant_df) <- paste0(colnames(significant_df), "PAGE_significant")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, significant_df)
  
  return(seurat_object)
  
}
GSM7058756_C1_fibroenriched<-process_PAGEanalysis1(seurat_object = GSM7058756_C1_fibroenriched,signatures_all = sign)
GSM7058757_C2_fibroenriched<-process_PAGEanalysis1(seurat_object = GSM7058757_C2_fibroenriched,signatures_all = sign)
#GSM7058758_C3_fibroenriched<-process_PAGEanalysis1(seurat_object = GSM7058758_C3_fibroenriched,signatures_all = sign)
GSM7058759_C4_fibroenriched<-process_PAGEanalysis1(seurat_object = GSM7058759_C4_fibroenriched,signatures_all = sign)

V688_1_fibro_md_sup<- V688_1_fibro@meta.data[,143:160]
V015_1_fibro_md_sup<- V015_1_fibro@meta.data[,143:160]
V015_2_fibro_md_sup<- V015_2_fibro@meta.data[,143:160]
V797_1_fibro_md_sup<-V797_1_fibro@meta.data[,143:160]
V797_2_fibro_md_sup<-V797_2_fibro@meta.data[,143:160]
V838_2_fibro_md_sup<-V838_2_fibro@meta.data[,143:160]

GSM7058756_C1_fibroenriched_md_sup<-GSM7058756_C1_fibroenriched@meta.data[,99:116]
GSM7058757_C2_fibroenriched_md_sup<-GSM7058757_C2_fibroenriched@meta.data[,99:116]
GSM7058759_C4_fibroenriched_md_sup<-GSM7058759_C4_fibroenriched@meta.data[,99:116]

correlation_matrices <- list(cor(V688_1_fibro_md_sup), 
                             cor(V015_1_fibro_md_sup), 
                             cor(V015_2_fibro_md_sup), 
                             cor(V797_1_fibro_md_sup),
                             cor(V797_2_fibro_md_sup),
                             cor(V838_2_fibro_md_sup),
                             cor(GSM7058756_C1_fibroenriched_md_sup),
                             cor(GSM7058757_C2_fibroenriched_md_sup),
                             cor(GSM7058759_C4_fibroenriched_md_sup)) 

num_spots <- c(dim(V688_1_fibro_md_sup)[1], 
               dim(V015_1_fibro_md_sup)[1], 
               dim(V015_2_fibro_md_sup)[1], 
               dim(V797_1_fibro_md_sup)[1],
               dim(V797_2_fibro_md_sup)[1], 
               dim(V838_2_fibro_md_sup)[1],
               dim(GSM7058756_C1_fibroenriched_md_sup)[1],
               dim(GSM7058757_C2_fibroenriched_md_sup)[1],
               dim(GSM7058759_C4_fibroenriched_md_sup)[1]
) 

z_matrices <- lapply(correlation_matrices, function(mat) {
  apply(mat, c(1, 2), fisher_z)
})

# Weight the transformed values by the number of spots
weighted_z_matrices <- mapply(function(z_matrix, n_spots) {
  z_matrix * n_spots
}, z_matrices, num_spots, SIMPLIFY = FALSE)

# Compute the sum of weighted z-values and total number of spots
sum_weighted_z_matrix <- Reduce("+", weighted_z_matrices)
total_spots <- sum(num_spots)

# Compute the weighted average of z-values
mean_weighted_z_matrix <- sum_weighted_z_matrix / total_spots

# Apply the inverse Fisher's z-transformation to the averaged matrix
merged_correlation_matrix <- apply(mean_weighted_z_matrix, c(1, 2), inverse_fisher_z)

# Ensure the diagonal remains 1 (since it's a correlation matrix)
diag(merged_correlation_matrix) <- 1

# Perform hierarchical clustering
row_model <- hclust(dist(merged_correlation_matrix))
col_model <- hclust(t(dist(merged_correlation_matrix)))



library(ComplexHeatmap)
# Create a Heatmap object with row and column dendrograms
p10 <- Heatmap(merged_correlation_matrix, 
              name = "Merged\nCorrelation\ncorrected for n.spots",
              col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
              cluster_rows = row_model,
              cluster_columns = col_model,
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_dend_side = "left",
              column_dend_side = "top",
              row_names_side = "left",
              column_names_side = "bottom",
              heatmap_legend_param = list(title = "Merged\nCorrelation\ncorrected\nfor n.spots"),
              row_dend_width = unit(4, "cm"),
              column_dend_height = unit(4, "cm"))

save(p10,file = 'p10.Robj')


V688_1_fibro_pval_sup<- V688_1_fibro@meta.data[,162:179]
V015_1_fibro_pval_sup<- V015_1_fibro@meta.data[,162:179]
V015_2_fibro_pval_sup<- V015_2_fibro@meta.data[,162:179]
V797_1_fibro_pval_sup<-V797_1_fibro@meta.data[,162:179]
V797_2_fibro_pval_sup<-V797_2_fibro@meta.data[,162:179]
V838_2_fibro_pval_sup<-V838_2_fibro@meta.data[,162:179]

GSM7058756_C1_fibroenriched_pval_sup<-GSM7058756_C1_fibroenriched@meta.data[,118:135]
GSM7058757_C2_fibroenriched_pval_sup<-GSM7058757_C2_fibroenriched@meta.data[,118:135]
GSM7058759_C4_fibroenriched_pval_sup<-GSM7058759_C4_fibroenriched@meta.data[,118:135]

calculate_hypergeometric_pvals <- function(df) {
  categories <- colnames(df)
  comb_matrix <- combn(categories, 2, simplify = FALSE)
  
  pval_matrix <- matrix(0, nrow = length(categories), ncol = length(categories))
  colnames(pval_matrix) <- categories
  rownames(pval_matrix) <- categories
  
  N <- nrow(df)
  
  for (comb in comb_matrix) {
    cat1 <- comb[1]
    cat2 <- comb[2]
    
    K <- sum(df[[cat1]])
    n <- sum(df[[cat2]])
    k <- sum(df[[cat1]] & df[[cat2]])
    
    pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    
    pval_matrix[cat1, cat2] <- pval
    pval_matrix[cat2, cat1] <- pval
  }
  
  diag(pval_matrix) <- 1  # Diagonal elements should be 1 (or another suitable value, as self-comparisons are not meaningful)
  
  return(pval_matrix)
}


hyp_V688_1<-calculate_hypergeometric_pvals(V688_1_fibro_pval_sup)
hyp_V015_1<-calculate_hypergeometric_pvals(V015_1_fibro_pval_sup)
hyp_V015_2<-calculate_hypergeometric_pvals(V015_2_fibro_pval_sup)
hyp_V797_1<-calculate_hypergeometric_pvals(V797_1_fibro_pval_sup)
hyp_V797_2<-calculate_hypergeometric_pvals(V797_2_fibro_pval_sup)
hyp_V838_2<-calculate_hypergeometric_pvals(V838_2_fibro_pval_sup)
hyp_C1<-calculate_hypergeometric_pvals(GSM7058756_C1_fibroenriched_pval_sup)
hyp_C2<-calculate_hypergeometric_pvals(GSM7058757_C2_fibroenriched_pval_sup)
hyp_C4<-calculate_hypergeometric_pvals(GSM7058759_C4_fibroenriched_pval_sup)


combined_pval_matrix <- combine_pvals(list(hyp_V688_1,
                                           hyp_V015_1,
                                           hyp_V015_2,
                                           hyp_V797_1,
                                           hyp_V797_2,
                                           hyp_V838_2,
                                           hyp_C1,
                                           hyp_C2,
                                           hyp_C4))
combined_pval_df <- as.data.frame(combined_pval_matrix)
combined_pval_df

row_model <- hclust(dist(combined_pval_df),method = 'mcquitty')
col_model <- hclust(t(dist(combined_pval_df)),method = 'mcquitty')

library(ComplexHeatmap)
p11 <- Heatmap(-log(combined_pval_df+5e-324), 
              name = "Combined\nenrichment\nscore\nvia Fisher's\nMethod",
              col = colorRamp2(c(0, 400, 800), c("white","#E94B3C",'#6C1413')),
              cluster_rows = row_model,
              cluster_columns = col_model,
              show_row_dend = TRUE,
              show_column_dend = T,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_dend_side = "left",
              column_dend_side = "top",
              row_names_side = "left",
              column_names_side = "bottom",
              heatmap_legend_param = list(title = "Combined\nenrichment\nscore\nvia Fisher's\nMethod"),
              row_dend_width = unit(3, "cm"),
              column_dend_height = unit(4, "cm"))

save(p11,file = 'p11.Robj')


