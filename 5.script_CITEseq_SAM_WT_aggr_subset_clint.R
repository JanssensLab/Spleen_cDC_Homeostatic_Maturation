# Script for further processing of cDC1 subset of WT_aggr object

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/")

sampleName <- "SAM2and3_WT_subset" #Change for this analysis!!!
sampleFolder<-"SAM2and3_WT/"

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')
source('~/VIB/DATA/Roos/Daan 1/script_featurePlots.R')

##### Create new clusters
##new clusters
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['SCT_tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['SCT_umap']]@cell.embeddings, stringsAsFactors = F)

Idents(seuratObj)<-seuratObj@meta.data$SCT_clusters
seuratObj@meta.data$sliced_clusters<-seuratObj@meta.data$SCT_clusters

#Use ADT clusters for removing contaminants from SCT clusters
#ADT Cluster 8: Sirpa + CD4 + ESAM (CD172a) -> Mig cDC2s!!
#ADT Cluster 9: TCRb + CD8b + Ly49D + CD127 -> NK cells and T cells
#ADT Cluster 10: random? Not clustered on RNA -> contaminants??
#ADT cluster 11: HSCs??
#ADT cluster 12: CD21 CD79b IgD -> B cells

# Split off cDC2s (ADT cluster!!!!!)
Idents(seuratObj)<-seuratObj@meta.data$ADT_clusters #Change annotation!
wantedCells<-WhichCells(seuratObj, idents = 8)
colorSomeCells(clusterMatrix, umapTable, wantedCells)

Idents(seuratObj)<-seuratObj@meta.data$sliced_clusters #Return to SCT clusters of new Ident!!
seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 18)
DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 8)
seuratObj@meta.data$sliced_clusters<-Idents(seuratObj) #Save new Ident in sliced clusters!!

# Split off NK cells and T cells (ADT cluster!!!!!)
Idents(seuratObj)<-seuratObj@meta.data$ADT_clusters #Change annotation!
wantedCells<-WhichCells(seuratObj, idents = 9)
colorSomeCells(clusterMatrix, umapTable, wantedCells)

Idents(seuratObj)<-seuratObj@meta.data$sliced_clusters #Return to SCT clusters of new Ident!!
seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 19)
DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 8)
seuratObj@meta.data$sliced_clusters<-Idents(seuratObj) #Save new Ident in sliced clusters!!

# Split off cluster 10 ADT (random?) (ADT cluster!!!!!)
Idents(seuratObj)<-seuratObj@meta.data$ADT_clusters #Change annotation!
wantedCells<-WhichCells(seuratObj, idents = 10)
colorSomeCells(clusterMatrix, umapTable, wantedCells)

Idents(seuratObj)<-seuratObj@meta.data$sliced_clusters #Return to SCT clusters of new Ident!!
seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 20)
DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 8)
seuratObj@meta.data$sliced_clusters<-Idents(seuratObj) #Save new Ident in sliced clusters!!

# Split off HSCs (ADT cluster!!!!!)
Idents(seuratObj)<-seuratObj@meta.data$ADT_clusters #Change annotation!
wantedCells<-WhichCells(seuratObj, idents = 11)
colorSomeCells(clusterMatrix, umapTable, wantedCells)

Idents(seuratObj)<-seuratObj@meta.data$sliced_clusters #Return to SCT clusters of new Ident!!
seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 21)
DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 8)
seuratObj@meta.data$sliced_clusters<-Idents(seuratObj) #Save new Ident in sliced clusters!!

# Split up B cells (ADT cluster!!!!!)
Idents(seuratObj)<-seuratObj@meta.data$ADT_clusters #Change annotation!
wantedCells<-WhichCells(seuratObj, idents = 12)
colorSomeCells(clusterMatrix, umapTable, wantedCells)

Idents(seuratObj)<-seuratObj@meta.data$sliced_clusters #Return to SCT clusters of new Ident!!
seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 22)
DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 8)
seuratObj@meta.data$sliced_clusters<-Idents(seuratObj) #Save new Ident in sliced clusters!!

##Save new clustering in seuratobject
seuratObj@meta.data$sliced_clusters <- seuratObj@active.ident #Sliced clustering
seuratObj@meta.data$annotated_clusters <- seuratObj@active.ident #Sliced clustering

seuratObj@meta.data$sliced_clusters<- factor(seuratObj@meta.data$sliced_clusters,sort(as.numeric(levels(seuratObj@meta.data$sliced_clusters)))) #reorder levels
seuratObj@meta.data$annotated_clusters<- seuratObj@meta.data$sliced_clusters
DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 8, group.by = "sliced_clusters")
DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 8)

##### Read object
seuratObj <- readRDS(file=paste0(sampleFolder,"results_subset/Robjects/seuratObj_sliced_",sampleName,"_clint.rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"results_subset/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"results_subset/Robjects/seuratObj_sliced_",sampleName,"_clint.rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

################################################################################
########## PLOTS
################################################################################
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['SCT_tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['SCT_umap']]@cell.embeddings, stringsAsFactors = F)

## New function
drawUMI_mitoPlot_new<-function(coordsTable, reductionType, clusterMatrix, columnName, titleInfo){
  
  columnNr<-which(colnames(clusterMatrix)==columnName)
  
  p <- ggplot()+
    geom_point(aes(x=sctTSNE_1,y=sctTSNE_2, colour=clusterMatrix[,columnNr]), data=coordsTable, size=2, shape=20) 
  
  if(reductionType=="umap"){
    p <- ggplot()+
      geom_point(aes(x=sctUMAP_1,y=sctUMAP_2, colour=clusterMatrix[,columnNr]), data=coordsTable, size=2, shape=20) 
  }
  
  p<-p +
    scale_colour_gradientn(colours = c("darkblue","cyan","green","yellow","orange","darkred")) +
    ggtitle(paste0(titleInfo," (",reductionType,")")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  
  return(p)
}

########## UMI plot ##########
p1<-drawUMI_mitoPlot_new(tsneTable, 'tsne', clusterMatrix, 'nCount_RNA',"UMI")
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'nCount_RNA',"UMI")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results_subset/QC/1_UMI_",sampleName,".png"), width = 20)
ggsave(p2, file=paste0(sampleFolder,"results_subset/QC/1_UMI_",sampleName,"_new.png"), width = 10)

########## mito.genes plot ##########
p1<-drawUMI_mitoPlot_new(tsneTable, 'tsne', clusterMatrix, 'subsets_Mito_percent',"mito")
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'subsets_Mito_percent',"mito")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results_subset/QC/2_percMito_",sampleName,".png"), width = 20)
ggsave( p2,  file=paste0(sampleFolder,"results_subset/QC/2_percMito_",sampleName,"_new.png"), width = 10)

########## PCA plot ##########
# pdf(file=paste0(sampleFolder,"results_subset/QC/13a_PCA_",sampleName,".pdf"), width=10)
DimPlot(object = seuratObj, reduction = "SCT_pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "SCT_pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "SCT_pca", dims = c(1,3))
# dev.off()

U1<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, label.size = 4, group.by = "SCT_clusters") +
  labs(title = "SCT clusters on SCT UMAP")
U2<-DimPlot(seuratObj, reduction = "ADT_umap", label = T, repel = T, label.size = 4, group.by = "ADT_clusters") +
  labs(title = "ADT clusters on ADT UMAP")

#ADT clustering on SCT UMAP and vice versa
U3<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, label.size = 4, group.by = "ADT_clusters") +
  labs(title = "ADT clusters on SCT UMAP")
U4<-DimPlot(seuratObj, reduction = "ADT_umap", label = T, repel = T, label.size = 4, group.by = "SCT_clusters") +
  labs(title = "SCT clusters on ADT UMAP")


T1<-DimPlot(seuratObj, reduction = "SCT_tsne", label = T, repel = T, label.size = 4, group.by = "SCT_clusters") +
  labs(title = "SCT clusters on SCT tSNE")
T2<-DimPlot(seuratObj, reduction = "ADT_tsne", label = T, repel = T, label.size = 4, group.by = "ADT_clusters") +
  labs(title = "ADT clusters on ADT tSNE")

#ADT clustering on SCT tSNE and vice versa
T3<-DimPlot(seuratObj, reduction = "SCT_tsne", label = T, repel = T, label.size = 4, group.by = "ADT_clusters") +
  labs(title = "ADT clusters on SCT tSNE")
T4<-DimPlot(seuratObj, reduction = "ADT_tsne", label = T, repel = T, label.size = 4, group.by = "SCT_clusters") +
  labs(title = "SCT clusters on ADT tSNE")

################################################################################
########## ANNOTATION
################################################################################
dir.create(paste0(sampleFolder,"results_subset/Annotation/"))

pdf(file=paste0(sampleFolder,"results_subset/Annotation/1_annotation_",sampleName,".pdf"), width = 15)
grid.arrange(U1, U3, ncol=2)
grid.arrange(U2, U4, ncol=2)
grid.arrange(T1, T3, ncol=2)
grid.arrange(T2, T4, ncol=2)
dev.off()


################################################################################
########## COLOR CELLS ACCORDING TO ADT on SCT + vice versa
################################################################################
DimPlot(seuratObj, reduction = "SCT_umap", label = F, group.by="SCT_clusters")
DimPlot(seuratObj, reduction = "SCT_umap", label = F, group.by="ADT_clusters")

#ADT clusters
pdf(file=paste0(sampleFolder,"results_subset/Annotation/1_annotation_Color_ADT_clusters_on_SCT_UMAP_",sampleName,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$ADT_clusters))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$ADT_clusters==i),])))
  C1<-C1+ggtitle(paste0("ADT_cluster_",i))
  print(C1)
}
dev.off()

#SCT clusters
pdf(file=paste0(sampleFolder,"results_subset/Annotation/1_annotation_Color_SCT_clusters_on_SCT_UMAP_",sampleName,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$SCT_clusters))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$SCT_clusters==i),])))
  C1<-C1+ggtitle(paste0("SCT_cluster_",i))
  print(C1)
}
dev.off()

#########################################################################################################################################

########################################
##### all clusters vs all clusters
########################################

dir.create(paste0(sampleFolder,"results_subset/Marker_lists"))

Idents(seuratObj)<-seuratObj@meta.data$ADT_clusters

### Find ADTmarkers for every ADT cluster compared to all remaining cells, report only the positive ones
ADTMarkers_ADTclus <- FindAllMarkers(seuratObj, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_ADTclus$cluster)
saveRDS(ADTMarkers_ADTclus, file=paste0(sampleFolder,"results_subset/Robjects/ADTmarkersList_ADTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerADTcluster']]<-paste0(table(ADTMarkers_ADTclus$cluster)," ADT markers for ADT cluster ",rownames(table(ADTMarkers_ADTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrADTclusters_ADTclus<-max(as.numeric(names(table(ADTMarkers_ADTclus$cluster))))
totalNrADTclusters_ADTclusPlusOne<-totalNrADTclusters_ADTclus+1
ADTmarkersList_ADTclus<-list()

for(i in 1:totalNrADTclusters_ADTclusPlusOne){
  clusterNr<-i-1
  
  tmp<-ADTMarkers_ADTclus[ADTMarkers_ADTclus$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  ADTmarkersList_ADTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(ADTmarkersList_ADTclus)<-paste0("ADTcluster",0:totalNrADTclusters_ADTclus)

### Write to Excel
library('openxlsx')
write.xlsx(ADTmarkersList_ADTclus, file =paste0(sampleFolder, "results_subset/Marker_lists/ADTmarkersList_ADTclus_",sampleName,".xlsx"))

########################################

### Find SCTmarkers for every ADT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_ADTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_ADTclus$cluster)
saveRDS(SCTMarkers_ADTclus, file=paste0(sampleFolder,"results_subset/Robjects/SCTmarkersList_ADTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerADTcluster']]<-paste0(table(SCTMarkers_ADTclus$cluster)," SCT markers for ADT cluster ",rownames(table(SCTMarkers_ADTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrSCTclusters_ADTclus<-max(as.numeric(names(table(SCTMarkers_ADTclus$cluster))))
totalNrSCTclusters_ADTclusPlusOne<-totalNrSCTclusters_ADTclus+1
SCTmarkersList_ADTclus<-list()

for(i in 1:totalNrSCTclusters_ADTclusPlusOne){
  clusterNr<-i-1
  
  tmp<-SCTMarkers_ADTclus[SCTMarkers_ADTclus$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  SCTmarkersList_ADTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(SCTmarkersList_ADTclus)<-paste0("ADTcluster",0:totalNrSCTclusters_ADTclus)

### Write to Excel
library('openxlsx')
write.xlsx(SCTmarkersList_ADTclus, file =paste0(sampleFolder, "results_subset/Marker_lists/SCTmarkersList_ADTclus_",sampleName,".xlsx"))

########################################

### Find ADTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
Idents(seuratObj)<-seuratObj@meta.data$SCT_clusters

ADTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_SCTclus$cluster)
saveRDS(ADTMarkers_SCTclus, file=paste0(sampleFolder,"results_subset/Robjects/ADTmarkersList_SCTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerSCTcluster']]<-paste0(table(ADTMarkers_SCTclus$cluster)," ADT markers for SCT cluster ",rownames(table(ADTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrADTclusters_SCTclus<-max(as.numeric(names(table(ADTMarkers_SCTclus$cluster))))
totalNrADTclusters_SCTclusPlusOne<-totalNrADTclusters_SCTclus+1
ADTmarkersList_SCTclus<-list()

for(i in 1:totalNrADTclusters_SCTclusPlusOne){
  clusterNr<-i-1
  
  tmp<-ADTMarkers_SCTclus[ADTMarkers_SCTclus$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  ADTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(ADTmarkersList_SCTclus)<-paste0("SCTcluster",0:totalNrADTclusters_SCTclus)

### Write to Excel
library('openxlsx')
write.xlsx(ADTmarkersList_SCTclus, file =paste0(sampleFolder, "results_subset/Marker_lists/ADTmarkersList_SCTclus_",sampleName,".xlsx"))

########################################

### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results_subset/Robjects/SCTmarkersList_SCTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTcluster']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset/Robjects/diagnostics_",sampleName,"_clint.rds"))

### Create list with markers
totalNrSCTclusters_SCTclus<-max(as.numeric(names(table(SCTMarkers_SCTclus$cluster))))
totalNrSCTclusters_SCTclusPlusOne<-totalNrSCTclusters_SCTclus+1
SCTmarkersList_SCTclus<-list()

for(i in 1:totalNrSCTclusters_SCTclusPlusOne){
  clusterNr<-i-1
  
  tmp<-SCTMarkers_SCTclus[SCTMarkers_SCTclus$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  SCTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(SCTmarkersList_SCTclus)<-paste0("SCTcluster",0:totalNrSCTclusters_SCTclus)

### Write to Excel
library('openxlsx')
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results_subset/Marker_lists/SCTmarkersList_SCTclus_",sampleName,".xlsx"))

################################################################################
########################################
##### Markers annotated clusters
########################################
seuratObj@meta.data$annotated_clusters <- seuratObj@active.ident 
levels(seuratObj@meta.data$annotated_clusters) <- c("Res cDC1s 1","Res cDC1s 5","Res cDC1s 2","Res cDC1s 3","Res cDC1s 4",'Prolif Res cDC1s 1',"Ccl3+Ccl4+ Res cDC1s",
                                                    "Mig cDC1s",'Prolif Res cDC1s 2',"Cxcl9+Cxcl10+ Res cDC1s",'Prolif Res cDC1s 4',"S100a+Ccr2+CD300LG+ Res cDC1s",
                                                    "Lower quality Res cDC1s 1",'Prolif Res cDC1s 3',"Lower UMI Res cDC1s","Lower quality Res cDC1s 2",
                                                    "IFN1+ Res cDC1s","Lower quality Mig cDC1s",'Doublet cDC1s/cDC2s',"Doublets cDC1s/NKcells", "Unknown contaminating cells",
                                                    "CD41+ Precursor cells","Doublets cDC1s/Bcells")
U_annot<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)
ggsave(U_annot, file=paste0(sampleFolder,"results_subset/Annotation/2_UMAP_annotated_",sampleName,".png"), height = 10, width = 20, dpi = "retina")

Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters

### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results_subset/Robjects/SCTmarkersList_SCTclus_",sampleName,"_annotated.rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTclusterannotated']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

### Create list with markers
totalNrSCTclusters_SCTclus<-names(table(SCTMarkers_SCTclus$cluster))
SCTmarkersList_SCTclus<-list()

for(i in totalNrSCTclusters_SCTclus){
  tmp<-SCTMarkers_SCTclus[SCTMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC

  SCTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
library('openxlsx')
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results_subset/Marker_lists/SCTmarkersList_SCTclus_",sampleName,"_annotated.xlsx"))

######################################################

### Make heatmap for annotated clusters
SCTMarkers_SCTclus<- readRDS(file=paste0(sampleFolder,"results_subset/Robjects/SCTmarkersList_SCTclus_",sampleName,"_annotated.rds"))

## Perform on a subset -> better view of smaller clusters!!
seuratObj.small <- subset(seuratObj, downsample = 300)

## Heatmap SCT markers on SCT clusters
top10 <- SCTMarkers_SCTclus %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
D1<-DoHeatmap(seuratObj.small, features = top10$gene, group.by = "annotated_clusters") + NoLegend()
ggsave(D1, file=paste0(sampleFolder, "results_subset/Heatmaps/Heatmap_SCTmarkersList_Annotatedclus_",sampleName,".png"), 
       height = 20, width = 12, dpi = "retina")

pdf(file=paste0(sampleFolder, "results_subset/Heatmaps/Heatmap_SCTmarkersList_Annotatedclus_",sampleName,".pdf"), 
    height = 25, width = 25)
DoHeatmap(seuratObj.small, features = top10$gene, group.by = "annotated_clusters") + NoLegend()
dev.off()


########################################
#### Subset seuratobject
########################################

##Check precursors!!
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "CD41+ Precursor cells"))
Precursors<-WhichCells(seuratObj, idents = "CD41+ Precursor cells")
seuratObj@assays$SCT@data[c("Itgax","Flt3","Siglech","Cd19","Cd3d","H2-K1","H2-Eb1"),Precursors]
seuratObj@assays$ADT@data[c("CD11c","CD135","Siglec","CD19","CD3e","IA-IE", "H2Kb"),Precursors] 

# Remove the bad clusters!!
seuratObjNew<-subset(seuratObj, idents = c("Res cDC1s 1","Res cDC1s 5","Res cDC1s 2","Res cDC1s 3","Res cDC1s 4",
                                           'Prolif Res cDC1s 1',"Ccl3+Ccl4+ Res cDC1s","Mig cDC1s",'Prolif Res cDC1s 2',
                                           "Cxcl9+Cxcl10+ Res cDC1s",'Prolif Res cDC1s 4',"S100a+Ccr2+CD300LG+ Res cDC1s",
                                           'Prolif Res cDC1s 3',"IFN1+ Res cDC1s"))

seuratObjNew@meta.data$annotated_clusters<-factor(seuratObjNew@meta.data$annotated_clusters, 
                                           levels=c("S100a+Ccr2+CD300LG+ Res cDC1s","Res cDC1s 1","Res cDC1s 2",
                                                    "Res cDC1s 3","Res cDC1s 4","Ccl3+Ccl4+ Res cDC1s","Res cDC1s 5",
                                                    'Prolif Res cDC1s 1','Prolif Res cDC1s 2','Prolif Res cDC1s 3','Prolif Res cDC1s 4',
                                                    "Cxcl9+Cxcl10+ Res cDC1s","IFN1+ Res cDC1s","Mig cDC1s"))

Idents(seuratObjNew)<-seuratObjNew@meta.data$annotated_clusters

U2 <- DimPlot(seuratObjNew, reduction = "SCT_umap", label = T, repel = T, label.size = 4)
ggsave(U2, file=paste0(sampleFolder,"results_subset/Annotation/2_UMAP_annotated_cleaned_",sampleName,".png"), height = 10, width = 20, dpi = "retina")

##Replace seuratObj with clean version for further steps
seuratObj<-seuratObjNew
rm(seuratObjNew)
gc()

##########################################################
######### Automatic clustering annotation 
##########################################################

DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, group.by = "SCT_snn_res.0.3", label.size = 4)
DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, group.by = "SCT_snn_res.0.8", label.size = 4)
DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)

seuratObj@meta.data$final_annotation<-as.character(seuratObj@meta.data$annotated_clusters)

#Change annotation
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="Prolif Res cDC1s 1"),"final_annotation"] <- "Proliferating cDC1s"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="Prolif Res cDC1s 2"),"final_annotation"] <- "Proliferating cDC1s"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="Prolif Res cDC1s 3"),"final_annotation"] <- "Proliferating cDC1s"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="Prolif Res cDC1s 4"),"final_annotation"] <- "Proliferating cDC1s"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="S100a+Ccr2+CD300LG+ Res cDC1s"),"final_annotation"] <- "pre-cDC1s"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="Mig cDC1s"),"final_annotation"] <- "Late mature cDC1s"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="Res cDC1s 1"),"final_annotation"] <- "Early immature cDC1s"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="Res cDC1s 2"),"final_annotation"] <- "Early immature cDC1s"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="Res cDC1s 3"),"final_annotation"] <- "Late immature cDC1s"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="Res cDC1s 4"),"final_annotation"] <- "Late immature cDC1s"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="IFN1+ Res cDC1s"),"final_annotation"] <- "Early mature cDC1s"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="Cxcl9+Cxcl10+ Res cDC1s"),"final_annotation"] <- "Early mature cDC1s"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="Res cDC1s 5"),"final_annotation"] <- "Unknown immature cDC1s 2"
seuratObj@meta.data[which(seuratObj@meta.data$final_annotation=="Ccl3+Ccl4+ Res cDC1s"),"final_annotation"] <- "Unknown immature cDC1s"

DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, group.by = "final_annotation", label.size = 4)
Idents(seuratObj)<-seuratObj@meta.data$final_annotation

## Further cleaning performed based on UMAP coordinates (CellSelector and cutoffs for various clusters)
## CellSelector example Prolif cDC1s
U1 <- DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 4)
seuratObj <- CellSelector(U1, object=seuratObj, ident="Proliferating cDC1s")

## UMAP example Unknown immature cDC1s
clusterMatrix<-seuratObj@meta.data
umapTable<-as.data.frame(seuratObj[['SCT_umap']]@cell.embeddings, stringsAsFactors = F)
colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, idents = "Unknown immature cDC1s"))
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., sctUMAP_1 < -3)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = "Unknown immature cDC1s"))
colorSomeCells(clusterMatrix, umapTable, wantedCells)
seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = "Late immature cDC1s")
DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 3)

#Fix issue naming and factorization
seuratObj@meta.data$final_annotation<-as.character(Idents(seuratObj))
seuratObj@meta.data$final_annotation<-as.factor(seuratObj@meta.data$final_annotation)
Idents(seuratObj)<-seuratObj@meta.data$final_annotation
DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 3)

##########################################################################################

## Get cell cycle assignments (updated 2021)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Convert with BioMart
# Basic function to convert human to mouse gene names 
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

## Cell cycle scoring

# We assign scores in the CellCycleScoring function, which stores S and G2/M scores in object meta data, 
# along with the predicted classification of each cell in either G2M, S or G1 phase. 
# CellCycleScoring can also set the identity of the Seurat object to the cell-cycle phase 
# by passing set.ident = TRUE (the original identities are stored as old.ident). 
# Please note that Seurat does not use the discrete classifications (G2M/G1/S) in downstream cell cycle regression. 
# Instead, it uses the quantitative scores for G2M and S phase. 
# However, we provide our predicted classifications in case they are of interest.

seuratObj <- CellCycleScoring(seuratObj, s.features = m.s.genes, g2m.features = m.g2m.genes, assay= "SCT", set.ident = F)

# view cell cycle scores and phase assignments 2021!!!
head(seuratObj[[]])

U_CC<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, group.by = "Phase", label.size = 4)
ggsave(U_CC, file=paste0(sampleFolder,"results_subset/Cell_cycle/2_UMAP_Cell_cycle_phases_",sampleName,"_2021.png"), height = 10, width = 20, dpi = "retina")

pdf(file = paste0(sampleFolder,"results_subset/Cell_cycle/2_UMAP_Cell_cycle_phases_",sampleName,"_2021.pdf"),height = 7, width = 10)
U_CC
dev.off()

###########################################################################################################################

## New transition marker genes January 2021
seuratObj@meta.data$final_annotation2021<-seuratObj@meta.data$final_annotation

levels(seuratObj@meta.data$final_annotation2021) <- c("Early immature cDC1s","Early mature cDC1s","Late immature cDC1s",
                                                      "Late mature cDC1s","pre-cDC1s","Proliferating cDC1s",     
                                                      "Late immature cDC1s","Late immature cDC1s")

seuratObj@meta.data$final_annotation2021<-factor(seuratObj@meta.data$final_annotation2021,
                                             as.character(c("pre-cDC1s","Early immature cDC1s","Proliferating cDC1s",
                                                            "Late immature cDC1s","Early mature cDC1s","Late mature cDC1s"))) #reorder levels

###########################################################################################################################

#Unannotated UMAP for supplementary paper
pdf(file=paste0("Paper/2.Supplementary_paper_unannotated_",sampleName,"_2021.pdf"), height = 10, width = 15)
DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, group.by = "sliced_clusters", label.size = 6)
dev.off()

# Detailed annotated UMAP for supplementary paper
pdf(file=paste0("Paper/2.Supplementary_paper_detail_annotated_",sampleName,"_2021.pdf"), height = 10, width = 17)
DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)
dev.off()

###########################################################################################################################

library("RColorBrewer")

## Reorder again
seuratObj@meta.data$final_annotation2021<-factor(seuratObj@meta.data$final_annotation2021,
                                                 as.character(c("pre-cDC1s","Proliferating cDC1s","Early immature cDC1s",
                                                                "Late immature cDC1s","Early mature cDC1s","Late mature cDC1s"))) #reorder levels

#Save as pdf
pdf(file=paste0("Paper/1b.Figure1_paper_annotated_",sampleName,"_2021_optionb.pdf"), height = 10, width = 15)
DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, group.by = "final_annotation2021", label.size = 5,
        cols =c(brewer.pal(n = 6, name = "Set3")[1],"Yellow",brewer.pal(n = 6, name = "Set3")[c(3:6)])) 
dev.off()

Idents(seuratObj)<-seuratObj@meta.data$final_annotation2021

## Dotplots
Colors_dotplot<-c("#071AE5","#F50635")

wantedGenes<-c("S100a10","Ccr2","Vim","Fcer1g","Anxa2","Birc5","Top2a","Stmn1","Mki67","Cdca8",
               "Cd8a","Manf","Creld2","Pdia6","Hspa5","Apoe","Apol7c","Ccl4","Egr3","Nr4a2","Cadm1",
               "Cxcl10","Cxcl9","Nfkbia","Cd40","Isg15","Cd63","Ccr7","Fscn1","Il12b","Ccl5")
wantedGenes_ADT<-c("CD300LG","CD226","CD62L","CD43","CD11b-mh","CD8a","CD39","CD38","IA-IE","CD207-mh",
                   "CD103","XCR1","CD80","CD274","CD86","CD63","CD83", "ESAM", "CD197","CD278","CD107a")

wantedGenes<-rev(wantedGenes)
wantedGenes_ADT<-rev(wantedGenes_ADT)

pdf(file=paste0("Paper/1c.Figure1_paper_Dotplot_population_markers_",sampleName,"_2021.pdf"), width = 10, height = 6)
DotPlot(seuratObj, features = wantedGenes, group.by = "final_annotation2021", cols = Colors_dotplot, col.min = -2, col.max = 2) + RotatedAxis()
dev.off()

pdf(file=paste0("Paper/1c.Figure1_paper_Dotplot_population_ADTmarkers_",sampleName,"_2021.pdf"), width = 10, height = 6)
DotPlot(seuratObj, assay = "ADT", features = wantedGenes_ADT, group.by = "final_annotation2021", cols = Colors_dotplot, col.min = -2, col.max = 2) + RotatedAxis()
dev.off()

#############################
### Extra detail clusters 2021
############################
preDC_to_EarlyImm<-FindMarkers(seuratObj, ident.1 = "Early immature cDC1s", ident.2 = "pre-cDC1s", logfc.threshold = 0.20, only.pos = FALSE)
preDC_to_EarlyImm_ADT<-FindMarkers(seuratObj, ident.1 = "Early immature cDC1s", ident.2 = "pre-cDC1s", logfc.threshold = 0.20, only.pos = FALSE, assay = "ADT")
EarlyImm_to_LateImm<-FindMarkers(seuratObj, ident.1 = "Late immature cDC1s", ident.2 = "Early immature cDC1s", logfc.threshold = 0.20, only.pos = FALSE)
EarlyImm_to_LateImm_ADT<-FindMarkers(seuratObj, ident.1 = "Late immature cDC1s", ident.2 = "Early immature cDC1s", logfc.threshold = 0.20, only.pos = FALSE, assay = "ADT")
LateImm_to_EarlyMat<-FindMarkers(seuratObj, ident.1 = "Early mature cDC1s", ident.2 = "Late immature cDC1s", logfc.threshold = 0.20, only.pos = FALSE)
LateImm_to_EarlyMat_ADT<-FindMarkers(seuratObj, ident.1 = "Early mature cDC1s", ident.2 = "Late immature cDC1s", logfc.threshold = 0.20, only.pos = FALSE, assay = "ADT")
EarlyMat_to_LateMat<-FindMarkers(seuratObj, ident.1 = "Late mature cDC1s", ident.2 = "Early mature cDC1s",  logfc.threshold = 0.20, only.pos = FALSE)
EarlyMat_to_LateMat_ADT<-FindMarkers(seuratObj, ident.1 = "Late mature cDC1s", ident.2 = "Early mature cDC1s", logfc.threshold = 0.20, only.pos = FALSE, assay = "ADT")


##### Create list
listDEgenesExtra<-tibble::lst(preDC_to_EarlyImm, preDC_to_EarlyImm_ADT, EarlyImm_to_LateImm,
                              EarlyImm_to_LateImm_ADT, LateImm_to_EarlyMat, LateImm_to_EarlyMat_ADT,
                              EarlyMat_to_LateMat, EarlyMat_to_LateMat_ADT)

##Add geneSymbol in column (for the export)
listDEgenesExtra<-lapply(listDEgenesExtra,function(x){x<-cbind(x,'gene'=rownames(x))})
##Filter on adj.P-value
listDEgenesExtra<-lapply(listDEgenesExtra, function(x){dplyr::filter(x, p_val_adj<0.01)})
##Add score
listDEgenesExtra<-lapply(listDEgenesExtra, function(x){rbind(x[x$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/(pct.2+0.001)*avg_logFC),
                                                             x[x$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/(pct.1+0.001)*avg_logFC))})
##Sort on logFC
listDEgenesExtra<-lapply(listDEgenesExtra,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(listDEgenesExtra,file=paste0(sampleFolder,"results_subset/Robjects/detailClusters_sliced_",sampleName,"_final_2021.rds"))

##write to Excel
library('openxlsx')
write.xlsx(listDEgenesExtra, paste0(sampleFolder,"results_subset/Marker_lists/detailClusters_sliced_",sampleName,"_final_2021.xlsx"))

##########################################################################################

### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
# Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters3
SCTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results_subset/Robjects/SCTmarkersList_SCTclus_",sampleName,"_final_2021.rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTcluster_final2021']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

### Create list with markers
totalNrSCTclusters_SCTclus<-names(table(SCTMarkers_SCTclus$cluster))
SCTmarkersList_SCTclus<-list()

for(i in totalNrSCTclusters_SCTclus){
  tmp<-SCTMarkers_SCTclus[SCTMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  SCTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
library('openxlsx')
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results_subset/Marker_lists/SCTmarkersList_SCTclus_",sampleName,"_final_2021.xlsx"))


### Find ADTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
ADTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_SCTclus$cluster)
saveRDS(ADTMarkers_SCTclus, file=paste0(sampleFolder,"results_subset/Robjects/ADTmarkersList_SCTclus_",sampleName,"_final_2021.rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerSCTcluster_final2021']]<-paste0(table(ADTMarkers_SCTclus$cluster)," ADT markers for SCT cluster ",rownames(table(ADTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

### Create list with markers
totalNrADTclusters_SCTclus<-names(table(ADTMarkers_SCTclus$cluster))
ADTmarkersList_SCTclus<-list()

for(i in totalNrADTclusters_SCTclus){
  tmp<-ADTMarkers_SCTclus[ADTMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  ADTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
library('openxlsx')
write.xlsx(ADTmarkersList_SCTclus, file =paste0(sampleFolder, "results_subset/Marker_lists/ADTmarkersList_SCTclus_",sampleName,"_final_2021.xlsx"))

##########################################################################################

# Markers unannotated UMAP
library(future)
plan("multiprocess", workers = 6)

Idents(seuratObj)<-seuratObj@meta.data$sliced_clusters

### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results_subset/Robjects/SCTmarkersList_SCTclus_",sampleName,"_unannotated_2021.rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTcluster_unannotated2021']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

### Create list with markers
totalNrSCTclusters_SCTclus<-names(table(SCTMarkers_SCTclus$cluster))
SCTmarkersList_SCTclus<-list()

for(i in totalNrSCTclusters_SCTclus){
  tmp<-SCTMarkers_SCTclus[SCTMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  SCTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
library('openxlsx')
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results_subset/Marker_lists/SCTmarkersList_SCTclus_",sampleName,"_unannotated_2021.xlsx"))

### Find ADTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
ADTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_SCTclus$cluster)
saveRDS(ADTMarkers_SCTclus, file=paste0(sampleFolder,"results_subset/Robjects/ADTmarkersList_SCTclus_",sampleName,"_unannotated_2021.rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerSCTcluster_unannotated2021']]<-paste0(table(ADTMarkers_SCTclus$cluster)," ADT markers for SCT cluster ",rownames(table(ADTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

### Create list with markers
totalNrADTclusters_SCTclus<-names(table(ADTMarkers_SCTclus$cluster))
ADTmarkersList_SCTclus<-list()

for(i in totalNrADTclusters_SCTclus){
  tmp<-ADTMarkers_SCTclus[ADTMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  ADTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
library('openxlsx')
write.xlsx(ADTmarkersList_SCTclus, file =paste0(sampleFolder, "results_subset/Marker_lists/ADTmarkersList_SCTclus_",sampleName,"_unannotated_2021.xlsx"))

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"results_subset/Robjects/seuratObj_final_",sampleName,"_clint_2021.rds"))

##### Read object
seuratObj<-readRDS(file=paste0(sampleFolder,"results_subset/Robjects/seuratObj_final_",sampleName,"_clint_2021.rds"))

##########################################################################################

## Remove ADTs not expressed

## Load in ADT panel
library(openxlsx)
library(dplyr)
library(tidyverse)

Panel_Niels<-read.xlsx("CITE-seq_antibody_panel.xlsx")
Other_panel<-read.xlsx("PanelSAM_Clint.xlsx")
ADT_names<-rownames(seuratObj@assays$ADT)

intersect(Other_panel$AB_name,Panel_Niels$Target)
ADT_remove<-setdiff(Other_panel$AB_name,Panel_Niels$Target) #Missing 2! 19 instead of 21
setdiff(Panel_Niels$Target,Other_panel$AB_name)
setdiff(Panel_Niels$Target,intersect(Other_panel$AB_name,Panel_Niels$Target))

Test1<-sort(c(Panel_Niels$Target,rep("X",21)))
Test2<-sort(Other_panel$AB_name)

Testdf<-cbind(Test1,Test2)
Extra_ADT<-c("CD278.1","CD309-A0553") #2 missing ADT

rownames(Other_panel)<-make.unique(Other_panel$AB_name)
Other_panel_filtered<-Other_panel[!(rownames(Other_panel) %in% ADT_remove),] #Initial filter to 171

rownames(Other_panel_filtered)<-(Other_panel_filtered$whitelist_name)
Other_panel_filtered2<-Other_panel_filtered[!(rownames(Other_panel_filtered) %in% Extra_ADT),] #Subsequent filter to 169

intersect(Other_panel$AB_name,ADT_names)
intersect(Other_panel$whitelist_name,ADT_names)
intersect(Other_panel_filtered2$whitelist_name,ADT_names) #This is the list to filter on!!!! Matches seuratobj names!

seuratObj_filtered<-seuratObj
seuratObj_filtered@assays$ADT@counts<-seuratObj_filtered@assays$ADT@counts[Other_panel_filtered2$whitelist_name,]
seuratObj_filtered@assays$ADT@data<-seuratObj_filtered@assays$ADT@data[Other_panel_filtered2$whitelist_name,]
seuratObj_filtered@assays$ADT@scale.data<-seuratObj_filtered@assays$ADT@scale.data[intersect(seuratObj_filtered@assays[["ADT"]]@var.features,Other_panel_filtered2$whitelist_name),]
seuratObj_filtered@assays$ADT@var.features<-intersect(seuratObj_filtered@assays[["ADT"]]@var.features,Other_panel_filtered2$whitelist_name)
seuratObj_filtered@assays$ADT@meta.features<-seuratObj_filtered@assays$ADT@meta.features[Other_panel_filtered2$whitelist_name,]

##### Save object
saveRDS(seuratObj_filtered, file=paste0(sampleFolder,"results_subset/Robjects/seuratObj_filteredADT_",sampleName,"_clint_2021.rds"))

##### Read object
seuratObj_filtered<-readRDS(file=paste0(sampleFolder,"results_subset/Robjects/seuratObj_filteredADT_",sampleName,"_clint_2021.rds"))

##########################################################################################

## Create diet object for online tool (15/03/23)
seuratObj_filtered_diet<-DietSeurat(seuratObj_filtered, counts = T, data = T, scale.data = F,
                                    assays = c("SCT","ADT"), dimreducs = "SCT_umap", graphs = NULL)

## Metadata columns
DimPlot(seuratObj_filtered_diet, reduction = "SCT_umap", label = T, repel = T, group.by = "SCT_clusters", label.size = 5,
        cols =gg_color_hue(14))
        # cols =c(brewer.pal(n = 6, name = "Set3")[1],"Yellow",brewer.pal(n = 6, name = "Set3")[c(3:6)]))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

All<-c(levels(as.factor(seuratObj_filtered_diet$orig.ident)),levels(as.factor(seuratObj_filtered_diet$Phase)),
  levels(as.factor(seuratObj_filtered_diet$SCT_clusters)), #Remove ADT (also numbers like SCT!)
  levels(as.factor(seuratObj_filtered_diet$annotated_clusters)),levels(as.factor(seuratObj_filtered_diet$final_annotation2021)))
Filtered_wrong_order<-c(levels(as.factor(seuratObj_filtered_diet$orig.ident)),levels(as.factor(seuratObj_filtered_diet$Phase)),
  levels(as.factor(as.character(seuratObj_filtered_diet$SCT_clusters))),
  levels(as.factor(as.character(seuratObj_filtered_diet$annotated_clusters))),levels(as.factor(seuratObj_filtered_diet$final_annotation2021)))
Metadata_info<-c(intersect(All,Filtered_wrong_order),levels(as.factor(as.character(seuratObj_filtered_diet$ADT_clusters))))
Color_info<-c(gg_color_hue(2),gg_color_hue(3),gg_color_hue(14),gg_color_hue(14),
              c(brewer.pal(n = 6, name = "Set3")[1],"Yellow",brewer.pal(n = 6, name = "Set3")[c(3:6)]),
              gg_color_hue(8))
Metadata_column<-c(rep("orig.ident",2),rep("Phase",3),rep("SCT_clusters",14),rep("annotated_clusters",14),
                   rep("final_annotation2021",6),rep("ADT_clusters",8))
Info_Kevin<-as.data.frame(cbind(Metadata_info,Color_info,Metadata_column))

write.xlsx(Info_Kevin, file =paste0(sampleFolder, "results_subset/Annotation/Info_Kevin_cDC1_object",sampleName,".xlsx"))

## Change to name metadata columns tool (performed after upload to tool Kevin, but before upload to GEO)
# "Annotation", "Annotation (detailed)", "Sample", "Phase", "Clustering cDNA", "Clustering ADT"
seuratObj_filtered_diet$Annotation<-seuratObj_filtered_diet$final_annotation2021
seuratObj_filtered_diet$Annotation_detailed<-seuratObj_filtered_diet$annotated_clusters
seuratObj_filtered_diet$Sample<-seuratObj_filtered_diet$orig.ident
seuratObj_filtered_diet$Clustering_cDNA<-seuratObj_filtered_diet$SCT_clusters
seuratObj_filtered_diet$Clustering_ADT<-seuratObj_filtered_diet$ADT_clusters

##### Save object
saveRDS(seuratObj_filtered_diet, file=paste0(sampleFolder,"results_subset/Robjects/seuratObj_paper_diet",sampleName,"_2023.rds"))

##### Read object
seuratObj_filtered_diet<-readRDS(file=paste0(sampleFolder,"results_subset/Robjects/seuratObj_paper_diet",sampleName,"_2023.rds"))

##########################################################################################

###################
# Feature plots
###################
library(viridis)

## Pdf format
features<-c("Ccr2", "Fcer1g", "Klf4" , "Cd24a", "Pdia3", "Pdia6", "Hspa5", "Calr","Apol7c","Apoe","Fdps",
            "Nr4a2","Nr4a3","Egr3","Egr1","Dnase1l3","Clec9a","Ucp2","Cxcl10", "Cxcl9", "Iigp1", "Ifi47", "Gbp2","Gbp5",
            "Cd63", "Fscn1", "Il4i1", "Socs2")


pdf(file=paste0(sampleFolder,"results_subset/Feature_plots/Paper/Feature_plots/FigS4A_Feature_plot_paper_28_markers_",sampleName,"_blue_grey.pdf"), height = 7, width = 10)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features =feature, cols = c("grey", "blue"), 
                  reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) #+
    # scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

## Extra plots ADT Victor 
features<-c("CD62L", "XCR1" , "CD103", "CD207-mh", "ESAM", "CD197")


pdf(file=paste0(sampleFolder,"results_subset/Feature_plots/Paper/Feature_plots/Fig4C_Feature_plot_paper_6_ADTmarkers_",sampleName,".pdf"), height = 7, width = 10)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features = features, cols = c("grey", "blue"), 
                  reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) +
    scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

## Extra markers paper discussion 
features<-c("Abcg1", "Apol10b", "Birc5" , "Ccr7")

pdf(file=paste0(sampleFolder,"results_subset/Feature_plots/Paper/Feature_plots/FigS4A_Feature_plot_paper_4_extra_markers_",sampleName,"_blue_grey.pdf"), height = 7, width = 10)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features =feature, cols = c("grey", "blue"), 
                  reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) #+
  # scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

## Extra markers paper 2
features<-c("Dnajc3","Creld2","Hspa1a","Hspa1b","Ccl3","Ccl4","Ifit2")

pdf(file=paste0(sampleFolder,"results_subset/Feature_plots/Paper/Feature_plots/FigS4A_Feature_plot_paper_7_extra_markers_",sampleName,"_blue_grey.pdf"), height = 7, width = 10)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features =feature, cols = c("grey", "blue"), 
                  reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) #+
  # scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

## Check LXR genes
library("RColorBrewer")
Colset<-brewer.pal(n = 6, name = "Set3")

levels(Idents(seuratObj_filtered))[1]<-"Pre-cDC1s"

F1<-FeaturePlot(object = seuratObj_filtered, features = c("Nr1h2", "Nr1h3"), cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
V1<-VlnPlot(object = seuratObj_filtered, features = c("Nr1h2", "Nr1h3"), cols = Colset)

pdf(file=paste0(sampleFolder,"results_subset/Feature_plots/Feature_and_Violin_plot_LXR_markers_",sampleName,"_blue_grey.pdf"), height = 14, width = 20)
F1/V1
dev.off()

## Check CCR7 and Dnase1l3 genes 
library("RColorBrewer")
Colset<-brewer.pal(n = 6, name = "Set3")

levels(Idents(seuratObj_filtered))[1]<-"Pre-cDC1s"

F1<-FeaturePlot(object = seuratObj_filtered, features = c("Ccr7", "CD197"), cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
V1<-VlnPlot(object = seuratObj_filtered, features = c("Ccr7", "CD197"), cols = Colset)

F2<-FeaturePlot(object = seuratObj_filtered, features = c("Dnase1l3"), cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
V2<-VlnPlot(object = seuratObj_filtered, features = c("Dnase1l3"), cols = Colset)

pdf(file=paste0(sampleFolder,"results_subset/Feature_plots/Feature_and_Violin_plot_rebuttal_",sampleName,"_blue_grey.pdf"), height = 14, width = 20)
F1/V1
dev.off()

pdf(file=paste0(sampleFolder,"results_subset/Feature_plots/Feature_and_Violin_plot_rebuttal2_",sampleName,"_blue_grey.pdf"), height = 14, width = 10)
F2/V2
dev.off()

## Check SREBF genes 
F1<-FeaturePlot(object = seuratObj_filtered, features = c("Srebf1", "Srebf2"), cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
V1<-VlnPlot(object = seuratObj_filtered, features = c("Srebf1", "Srebf2"))

pdf(file=paste0(sampleFolder,"results_subset/Feature_plots/Feature_and_Violin_plot_SREBF_markers_",sampleName,"_blue_grey.pdf"), height = 14, width = 20)
F1/V1
dev.off()

##########################################################################################

## Extra markers Srebf2 for paper 
features<-c("Abca1","Acly","Acot7",'Fdft1',"Hmgcr",'Insig1','Ldlr','Lrp1','Mttp',
            'Npc1l1','Pcsk9','Pcyt2','Pon1','Sqle','Stard4','Cebpa',"Hmgcs1","Fdps")

pdf(file=paste0(sampleFolder,"results_subset/Feature_plots/Paper/Feature_plots/FigS4A_Feature_plot_paper_Srebp_markers_",sampleName,"_blue_grey.pdf"), height = 7, width = 10)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features =feature, cols = c("grey", "blue"), 
                  reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) #+
  # scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

# PER 36
plots1 <- FeaturePlot(seuratObj, features = features[1:16],order=T, pt.size=0.3,min.cutoff = "q2", max.cutoff = "q98",combine = FALSE)
plots1 <- lapply(X = plots1, FUN = function(x) x + theme(text = element_text(size = 12),plot.title = element_text(size = 20), axis.text = element_text(size = 5)))

pdf(file = paste0("SAM2and3_WT/results_subset/Dorothea/SREBF2_Featureplots_q2-q98_per16.pdf"), width = 30, height = 30)
CombinePlots(plots = plots1, ncol=4)
dev.off()

## Extra markers Nr1h2 for paper 
features<-c("Abca1","Abcg1","Bhlhe40","Srebf1","Fasn","Pparg","Cdc42ep4",
            "Dnajb12","Eif1","Mnt","Ncor2","Sqstm1","Tbc1d1")

pdf(file=paste0(sampleFolder,"results_subset/Feature_plots/Paper/Feature_plots/FigS4A_Feature_plot_paper_Nr1h2_targets_",sampleName,"_blue_grey.pdf"), height = 7, width = 10)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features =feature, cols = c("grey", "blue"), 
                  reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) #+
  # scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

## Actual markers for paper 
features<-c("Hmgcr",'Ldlr',"Fdps")

pdf(file=paste0(sampleFolder,"results_subset/Feature_plots/Paper/Feature_plots/FigS4A_Feature_plot_paper_Cholesterol_markers_",sampleName,"_blue_grey.pdf"), height = 7, width = 10)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features =feature, cols = c("grey", "blue"), 
                  reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) #+
  # scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

##########################################################################################

## Detailed marker heatmaps
SCTMarkers_SCTclus<-readRDS(file=paste0(sampleFolder,"results_subset/Robjects/SCTmarkersList_SCTclus_",sampleName,"_final_2021.rds"))
ADTMarkers_SCTclus<-readRDS(file=paste0(sampleFolder,"results_subset/Robjects/ADTmarkersList_SCTclus_",sampleName,"_final_2021.rds"))

### Calculate scores for top10 lists!!!!
totalNrSCTclusters_SCTclus<-names(table(SCTMarkers_SCTclus$cluster))
SCTmarkersList_SCTclus<-list()

for(i in totalNrSCTclusters_SCTclus){

  tmp<-SCTMarkers_SCTclus[SCTMarkers_SCTclus$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  SCTmarkersList_SCTclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

## Heatmap SCT markers on SCT clusters
library("RColorBrewer")
Colset<-brewer.pal(n = 6, name = "Set3")

# top10 <- SCTMarkers_SCTclus %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top5_1 <- SCTmarkersList_SCTclus$`pre-cDC1s` %>% top_n(n = 5, wt = score)
top5_2 <- SCTmarkersList_SCTclus$`Proliferating cDC1s` %>% top_n(n = 5, wt = score)
top10_3 <- SCTmarkersList_SCTclus$`Early immature cDC1s` %>% top_n(n = 10, wt = score)
top21_4 <- SCTmarkersList_SCTclus$`Late immature cDC1s` %>% top_n(n = 21, wt = score)
top10_5 <- SCTmarkersList_SCTclus$`Early mature cDC1s` %>% top_n(n = 10, wt = score)
top10_6 <- SCTmarkersList_SCTclus$`Late mature cDC1s` %>% top_n(n = 10, wt = score)
top10<-rbind(top5_1,top5_2,top10_3,top21_4,top10_5,top10_6)
top10<-top10[-grep("Gm",top10$gene),] #Filter Gm genes

D1<-DoHeatmap(seuratObj, assay = "SCT", features = top10$gene, size = 3,
              draw.lines = T, group.colors = Colset) + 
  scale_color_manual(values = Colset) +
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill", na.value = "white") 

pdf(file=paste0(sampleFolder, "results_subset/Heatmaps/Heatmap_SCTmarkersList_SCTclus_",sampleName,"_2021.pdf"), 
    height = 15, width = 15)
D1
dev.off()

##########################################################################################
dir.create(paste0(sampleFolder,"results_subset/Rebuttal/"))

## Rebuttal average expression heatmap (use SCT markers from heatmap above)
Idents(seuratObj)

cluster.averages <- AverageExpression(seuratObj, return.seurat = TRUE)
cluster.averages

## Create heatmap
library("RColorBrewer")
Colset<-brewer.pal(n = 6, name = "Set3")

H1 <- DoHeatmap(cluster.averages, assay = "SCT", 
                features = top10$gene, size = 3,
                draw.lines = F, group.colors = Colset) + 
  scale_color_manual(values = Colset) +
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill")

pdf(file=paste0(sampleFolder,"results_subset/Rebuttal/Fig4A_Rebuttal_heatmap_maturation_average_population_expression_",sampleName,".pdf"), width = 15, height = 15)
H1
dev.off()

##############

## Rebuttal violin plots
V1<-VlnPlot(object = seuratObj, features = c("Ccr2", "Birc5","Hspa5","Nr4a2","Cxcl10","Cd63"),
            cols = Colset, log = F, ncol = 2)

pdf(file=paste0(sampleFolder,"results_subset/Rebuttal/Fig4B_Rebuttal_Violin_plots_",sampleName,".pdf"), height = 20, width = 14)
V1
dev.off()
