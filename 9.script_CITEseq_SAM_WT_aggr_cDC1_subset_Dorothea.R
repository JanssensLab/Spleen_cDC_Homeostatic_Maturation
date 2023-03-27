# Script for TF analysis (DoRothEA) on SCT assay of cDC1 subset of WT_aggr object

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library('dorothea')
library('tidyr')
library('pheatmap')
library('tibble')

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/")

sampleName <- "SAM2and3_WT_subset" 
sampleFolder<-"SAM2and3_WT/"

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')
source('~/VIB/DATA/Roos/Daan 1/script_featurePlots.R')

##### Read object
seuratObj_filtered<-readRDS(file=paste0(sampleFolder,"results_subset/Robjects/seuratObj_filteredADT_",sampleName,"_clint_2021.rds"))

seuratObj<-UpdateSeuratObject(seuratObj_filtered)
rm(seuratObj_filtered)
gc()

## Dorothea TF list
dir.create(paste0(sampleFolder,"results_subset/Dorothea/"))
mm<-dorothea_mm
mm_Abcg1<-mm[which(mm$target == "Abcg1"),]

write.csv(mm_Abcg1,paste0(sampleFolder,"results_subset/Dorothea/Mm_Abcg1_Dorothea_",sampleName,".csv"))

mm_Nr1h2<-mm[which(mm$tf == "Nr1h2"),]
write.csv(mm_Nr1h2,paste0(sampleFolder,"results_subset/Dorothea/Mm_Nr1h2_tf_Dorothea_",sampleName,".csv"))

## Dorothea analysis
# Holland et al. (2020) showed that clustering the cells based on their TF activity profiles can also be very interesting. 
# Indeed, clustering cells using TF activity computed with VIPER and DoRothEA performs better than using the expression level of the same TFs. 
# In addition, it brings complementary information to the clusters based on transcriptomics profiles.
# 
# Here, we first run VIPER on DoRothEA’s regulons to obtain TFs activity, by using the wrapper function run_viper(). 
# This function can deal with different input types such as matrix, dataframe, ExpressionSet or even Seurat objects. 
# In case of a seurat object the function returns the same seurat object with an additonal assay called dorothea 
# containing the TF activities in the slot data.

## We read Dorothea Regulons for mouse:
dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_mouse %>%
  dplyr::filter(confidence %in% c("A","B","C"))

## Extra documentation
regulon_Srebf <- regulon %>% dplyr::filter(tf %in% "Srebf")
regulon_Srebf <- regulon %>% dplyr::filter(tf %in% "Srebf2")
regulon_Srebf1 <- regulon %>% dplyr::filter(tf %in% "Srebf1")
regulon_Srebf1_target <- regulon %>% dplyr::filter(target %in% "Srebf1")
regulon_Srebf2_target <- regulon %>% dplyr::filter(target %in% "Srebf2")
regulon_filtered<-rbind(regulon_Srebf1,regulon_Srebf,regulon_Srebf1_target,regulon_Srebf2_target)

write.csv(regulon_filtered,paste0(sampleFolder,"results_subset/Dorothea/Mm_Srebf_Dorothea_",sampleName,".csv"))

##Remove superfluous data temporarily (issue with memory)
seuratObj@assays$ADT<-NULL
seuratObj@assays$RNA<-NULL

## We compute Viper Scores 
seuratObj <- run_viper(seuratObj, regulon, assay_key = "SCT",
                  options = list(method = "scale", minsize = 4, 
                                 eset.filter = FALSE, cores = 4, 
                                 verbose = FALSE))

# We then apply Seurat to cluster the cells following the same protocol as usual but using TF activity scores.

## We compute the Nearest Neighbours to perform cluster
DefaultAssay(object = seuratObj) <- "dorothea"
seuratObj <- ScaleData(seuratObj)
seuratObj <- RunPCA(seuratObj, features = rownames(seuratObj), verbose = FALSE)
seuratObj <- FindNeighbors(seuratObj, dims = 1:20, verbose = FALSE)
seuratObj <- FindClusters(seuratObj, resolution = 0.8, verbose = FALSE)

seuratObj <- RunUMAP(seuratObj, dims = 1:20, umap.method = "uwot", metric = "cosine")

## Save clusters
seuratObj@meta.data$dorothea_clusters<-seuratObj@meta.data$dorothea_snn_res.0.8

## Add data back in after heavy analysis
seuratObj@assays$ADT<-seuratObj_filtered@assays$ADT
seuratObj@assays$RNA<-seuratObj_filtered@assays$RNA

## Save object
saveRDS(seuratObj, paste0(sampleFolder,"results_subset/Robjects/seuratObj_Dorothea_",sampleName,".rds"))

## Marker analysis 
## New dorothea clusters
seuratObj.markers <- FindAllMarkers(seuratObj, only.pos = TRUE, min.pct = 0.10, 
                               logfc.threshold = 0.25, verbose = FALSE)
table(seuratObj.markers$cluster)
saveRDS(seuratObj.markers, file=paste0(sampleFolder,"results_subset/Robjects/Dorothea_markersList_Dorotheaclus_",sampleName,".rds"))

### Create list with markers
totalNrDorotheaclusters_Dorotheaclus<-max(as.numeric(names(table(seuratObj.markers$cluster))))
totalNrDorotheaclusters_DorotheaclusPlusOne<-totalNrDorotheaclusters_Dorotheaclus+1
DorotheamarkersList_Dorotheaclus<-list()

for(i in 1:totalNrDorotheaclusters_DorotheaclusPlusOne){
  clusterNr<-i-1
  
  tmp<-seuratObj.markers[seuratObj.markers$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  DorotheamarkersList_Dorotheaclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(DorotheamarkersList_Dorotheaclus)<-paste0("Dorotheacluster",0:totalNrDorotheaclusters_Dorotheaclus)

### Write to Excel
library('openxlsx')
write.xlsx(DorotheamarkersList_Dorotheaclus, file =paste0(sampleFolder, "results_subset/Dorothea/DorotheamarkersList_Dorotheaclus_",sampleName,".xlsx"))

## Final SCT annotation
Idents(seuratObj)<-seuratObj@meta.data$final_annotation2021
seuratObj.markers.2 <- FindAllMarkers(seuratObj, only.pos = TRUE, min.pct = 0.10, 
                                    logfc.threshold = 0.25, verbose = FALSE, assay = "dorothea")
table(seuratObj.markers.2$cluster)
saveRDS(seuratObj.markers.2, file=paste0(sampleFolder,"results_subset/Robjects/Dorothea_markersList_SCTclus_",sampleName,".rds"))

### Create list with markers
totalNrDorotheaclusters_SCTclus<-names(table(seuratObj.markers.2$cluster))
DorotheamarkersList_SCTclus<-list()

for(i in totalNrDorotheaclusters_SCTclus){
  tmp<-seuratObj.markers.2[seuratObj.markers.2$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  DorotheamarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
library('openxlsx')
write.xlsx(DorotheamarkersList_SCTclus, file =paste0(sampleFolder, "results_subset/Dorothea/DorotheamarkersList_SCTclus_",sampleName,".xlsx"))

################################################################################

## Create plots
U1<-DimPlot(seuratObj, reduction = "umap", group.by = "dorothea_clusters", label = TRUE, label.size = 8, pt.size = 0.5)
U2<-DimPlot(seuratObj, reduction = "umap", group.by = "final_annotation2021", label = TRUE, label.size = 5, pt.size = 0.5)
U3<-DimPlot(seuratObj, reduction = "SCT_umap", group.by = "dorothea_clusters", label = TRUE, label.size = 8, pt.size = 0.5)
U4<-DimPlot(seuratObj, reduction = "SCT_umap", group.by = "final_annotation2021", label = TRUE, label.size = 5, pt.size = 0.5)

DefaultAssay(object = seuratObj) <- "dorothea"
F1<-FeaturePlot(object = seuratObj, features = c("Srebf1", "Srebf2","Nr1h2","Nr1h3"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
F2<-FeaturePlot(object = seuratObj, features = c("Srebf1", "Srebf2","Nr1h2","Nr1h3"), cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
F3<-FeaturePlot(object = seuratObj, features = "Nr1h2", cols = c("yellow", "red"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
F4<-FeaturePlot(object = seuratObj, features = c("Fos", "Jun","Junb"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
F5<-FeaturePlot(object = seuratObj, features = c("Fos", "Jun","Junb"), cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)

## Save plots
pdf(file=paste0(sampleFolder,"results_subset/Dorothea/Dim_and_feature_plot_Dorothea_",sampleName,".pdf"), height = 14, width = 20)
U1
U2
F1
U3
U4
F2
dev.off()

pdf(file=paste0(sampleFolder,"results_subset/Dorothea/Nr1h2_Dorothea_TF_activity_",sampleName,".pdf"), height = 7, width = 10)
F3
dev.off()

pdf(file=paste0(sampleFolder,"results_subset/Dorothea/Dim_and_feature_plot_extra_ERK_Dorothea_",sampleName,".pdf"), height = 14, width = 20)
U1
U2
F4
U3
U4
F5
dev.off()

##################################################################

## Read object
seuratObj <- readRDS(paste0(sampleFolder,"results_subset/Robjects/seuratObj_Dorothea_",sampleName,".rds"))

## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df <- GetAssayData(seuratObj, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()

## We create a data frame containing the cells and their clusters
CellsOrigClusters <- data.frame(cell = names(Idents(seuratObj)), 
                            cell_type = as.character(Idents(seuratObj)),
                            check.names = F)

CellsDorotheaClusters <- data.frame(cell = rownames(seuratObj@meta.data), 
                            cell_type = as.character(seuratObj@meta.data$dorothea_clusters),
                            check.names = F)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_Origclusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsOrigClusters)

viper_scores_Dorotheaclusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsDorotheaClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_Origclusters_scores <- viper_scores_Origclusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

summarized_viper_Dorotheaclusters_scores <- viper_scores_Dorotheaclusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

# For visualization purposes, we select the 20 most variable TFs across clusters according to our scores.
## We select the 20 most variable TFs. (20*6 populations = 120)
highly_variable_tfs_Origclusters <- summarized_viper_Origclusters_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(180, var) %>%
  distinct(tf)

highly_variable_tfs_Dorotheaclusters <- summarized_viper_Dorotheaclusters_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(180, var) %>%
  distinct(tf)

## We prepare the data for the plot
summarized_viper_Origclusters_scores_df <- summarized_viper_Origclusters_scores %>%
  semi_join(highly_variable_tfs_Origclusters, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 
summarized_viper_Origclusters_scores_df<-summarized_viper_Origclusters_scores_df[c(5,6,1,3,2,4),]
rownames(summarized_viper_Origclusters_scores_df)[1]<-"Pre-cDC1s"
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

## Order on expr Late mature cDC1s
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df))
order_tfs<-order_tfs[order(order_tfs$`Late mature cDC1s`, decreasing = T),]

my_breaks_Origclusters <- c(seq(min(summarized_viper_Origclusters_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_Origclusters_scores_df)/palette_length, 
                   max(summarized_viper_Origclusters_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap_Origclusters <- pheatmap(order_tfs,fontsize=14, #t(summarized_viper_Origclusters_scores_df)
                       fontsize_row = 10, cluster_cols = F, cluster_rows = F,
                       color=my_color, breaks = my_breaks_Origclusters, 
                       main = "DoRothEA (Orig clusters)", angle_col = 45,
                       treeheight_col = 0,  border_color = NA) 

summarized_viper_Dorotheaclusters_scores_df <- summarized_viper_Dorotheaclusters_scores %>%
  semi_join(highly_variable_tfs_Dorotheaclusters, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks_Dorothea_clusters <- c(seq(min(summarized_viper_Dorotheaclusters_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_Dorotheaclusters_scores_df)/palette_length, 
                   max(summarized_viper_Dorotheaclusters_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap_Dorotheaclusters <- pheatmap(t(summarized_viper_Dorotheaclusters_scores_df),fontsize=14, 
                                    fontsize_row = 10,
                                    color=my_color, breaks = my_breaks_Dorothea_clusters, 
                                    main = "DoRothEA (Dorothea clusters)", angle_col = 45,
                                    treeheight_col = 0,  border_color = NA) 

pdf(file=paste0(sampleFolder,"results_subset/Dorothea/Viper_heatmap_Orig_clusters_",sampleName,".pdf"), height = 7, width = 10)
viper_hmap_Origclusters
dev.off()

pdf(file=paste0(sampleFolder,"results_subset/Dorothea/Viper_heatmap_Dorothea_clusters_",sampleName,".pdf"), height = 7, width = 10)
viper_hmap_Dorotheaclusters
dev.off()