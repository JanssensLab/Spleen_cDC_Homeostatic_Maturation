# Script for further processing of WT_aggr object

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

sampleName <- "SAM2and3_WT"
sampleFolder<-paste0(sampleName,"/")

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')
source('~/VIB/DATA/Roos/Daan 1/script_featurePlots.R')

## Changes to seuratObject
##### Create new clusters
##new clusters
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['SCT_tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['SCT_umap']]@cell.embeddings, stringsAsFactors = F)

Idents(seuratObj)<-seuratObj@meta.data$SCT_clusters
seuratObj@meta.data$sliced_clusters<-seuratObj@meta.data$SCT_clusters

################################################################################

################################################################################
########## PLOTS
################################################################################
clusterMatrix<-seuratObj@meta.data
# logTable<-as.matrix(seuratObj[['RNA']]@data)
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

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results/QC/11a_UMI.png"), width = 20)

########## mito.genes plot ##########
p1<-drawUMI_mitoPlot_new(tsneTable, 'tsne', clusterMatrix, 'subsets_Mito_percent',"mito")
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'subsets_Mito_percent',"mito")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results/QC/11b_percMito.png"), width = 20)

## Check annotations
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
dir.create(paste0(sampleFolder,"results/Annotation/"))

pdf(file=paste0(sampleFolder,"results/Annotation/1_annotation_",sampleName,".pdf"), width = 15)
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
pdf(file=paste0(sampleFolder,"results/Annotation/1_annotation_Color_ADT_clusters_on_SCT_UMAP_",sampleName,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$ADT_clusters))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$ADT_clusters==i),])))
  C1<-C1+ggtitle(paste0("ADT_cluster_",i))
  print(C1)
}
dev.off()

#SCT clusters
pdf(file=paste0(sampleFolder,"results/Annotation/1_annotation_Color_SCT_clusters_on_SCT_UMAP_",sampleName,".pdf"), width = 15)
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

dir.create(paste0(sampleFolder,"results/Marker_lists"))

Idents(seuratObj)<-seuratObj@meta.data$ADT_clusters

### Find ADTmarkers for every ADT cluster compared to all remaining cells, report only the positive ones
ADTMarkers_ADTclus <- FindAllMarkers(seuratObj, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_ADTclus$cluster)
saveRDS(ADTMarkers_ADTclus, file=paste0(sampleFolder,"results/Robjects/ADTmarkersList_ADTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerADTcluster']]<-paste0(table(ADTMarkers_ADTclus$cluster)," ADT markers for ADT cluster ",rownames(table(ADTMarkers_ADTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
write.xlsx(ADTmarkersList_ADTclus, file =paste0(sampleFolder, "results/Marker_lists/ADTmarkersList_ADTclus_",sampleName,".xlsx"))

########################################

### Find SCTmarkers for every ADT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_ADTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_ADTclus$cluster)
saveRDS(SCTMarkers_ADTclus, file=paste0(sampleFolder,"results/Robjects/SCTmarkersList_ADTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerADTcluster']]<-paste0(table(SCTMarkers_ADTclus$cluster)," SCT markers for ADT cluster ",rownames(table(SCTMarkers_ADTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
write.xlsx(SCTmarkersList_ADTclus, file =paste0(sampleFolder, "results/Marker_lists/SCTmarkersList_ADTclus_",sampleName,".xlsx"))

########################################

### Find ADTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
Idents(seuratObj)<-seuratObj@meta.data$SCT_clusters

ADTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_SCTclus$cluster)
saveRDS(ADTMarkers_SCTclus, file=paste0(sampleFolder,"results/Robjects/ADTmarkersList_SCTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerSCTcluster']]<-paste0(table(ADTMarkers_SCTclus$cluster)," ADT markers for SCT cluster ",rownames(table(ADTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
write.xlsx(ADTmarkersList_SCTclus, file =paste0(sampleFolder, "results/Marker_lists/ADTmarkersList_SCTclus_",sampleName,".xlsx"))

########################################

### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results/Robjects/SCTmarkersList_SCTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTcluster']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results/Marker_lists/SCTmarkersList_SCTclus_",sampleName,".xlsx"))


########################################################################################################################

DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 4)

### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results/Robjects/SCTmarkersList_SCTclus_",sampleName,"_annotated.rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTclusterannotated']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

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
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results/Marker_lists/SCTmarkersList_SCTclus_",sampleName,"_annotated.xlsx"))


##############################################################################################################################

##### Transpose annotation full object and clean up
##### Read Full aggr object
seuratObj_Full <- readRDS(file="SAM_aggr/results/Robjects/seuratObj_sliced_SAM_aggr_clint.rds")

# Get overlap of cells (no loss of cells)
Cells_complete<-rownames(seuratObj_Full@meta.data[which(seuratObj_Full@meta.data$orig.ident=="SAM2"|seuratObj_Full@meta.data$orig.ident=="SAM3"),])
Cells_subset<-colnames(seuratObj@assays$RNA)

#Subset SeuratObjFull
seuratObj_Full<-seuratObj_Full[,Cells_complete]

#Adjust cell suffix to match WT seuratobj
Cells_complete[grep("-2",Cells_complete)]<-gsub("-2","-1",Cells_complete[grep("-2",Cells_complete)]) #13049
Cells_complete[grep("-3",Cells_complete)]<-gsub("-3","-2",Cells_complete[grep("-3",Cells_complete)]) #17190

Cells_WT<-intersect(Cells_complete,Cells_subset)

#Can't rename features seuratObj -> revert back cell names for full object
Cells_WT_Full<-Cells_WT
Cells_WT_Full[grep("-2",Cells_WT_Full)]<-gsub("-2","-3",Cells_WT_Full[grep("-2",Cells_WT_Full)]) #17190->17189
Cells_WT_Full[grep("-1",Cells_WT_Full)]<-gsub("-1","-2",Cells_WT_Full[grep("-1",Cells_WT_Full)]) #13049->13035

# Change to character otherwise issue with factor
seuratObj@meta.data$annotated_clusters<-seuratObj@active.ident
seuratObj@meta.data[,"annotated_clusters"]<-as.character(seuratObj@meta.data[,"annotated_clusters"])
seuratObj@meta.data[Cells_WT,"annotated_clusters"]<-as.character(seuratObj_Full@meta.data[Cells_WT_Full,"annotated_clusters"])
seuratObj@meta.data[,"annotated_clusters"]<-as.factor(seuratObj@meta.data[,"annotated_clusters"])

# Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters
U_aggr<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, label.size = 4, group.by = "annotated_clusters")
T_aggr<-DimPlot(seuratObj, reduction = "SCT_tsne", label = T, repel = T, label.size = 4, group.by = "annotated_clusters")

pdf(file=paste0(sampleFolder,"results/Annotation/1_annotation_full_aggregate_",sampleName,".pdf"), width = 15)
U_aggr
T_aggr
dev.off()

## Annotation notes: use full aggr as template!! More cells, so more detailed.
seuratObj@meta.data$annotated_clusters_full<-seuratObj@meta.data$annotated_clusters

## Further cleaning performed based on UMAP coordinates (CellSelector and cutoffs for various clusters)
clusterMatrix<-seuratObj@meta.data
umapTable<-as.data.frame(seuratObj[['SCT_umap']]@cell.embeddings, stringsAsFactors = F)

#UMAP coordinates
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., sctUMAP_2 > 0)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = "NK cells"))
colorSomeCells(clusterMatrix, umapTable, wantedCells)
seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = "Doublets cDC1s/RBCs")

#Cellselector
U1 <- DimPlot(seuratObj, reduction = "SCT_umap", label = T, label.size = 4)
seuratObj <- CellSelector(U1, object=seuratObj, ident="Mig cDC1s")

seuratObj@meta.data$annotated_clusters_clean <- seuratObj@active.ident #Sliced clustering

Order_cellpop<-data.frame(table(seuratObj@meta.data$annotated_clusters_clean))
Order_cellpop<-Order_cellpop[order(Order_cellpop$Freq, decreasing = T),]

seuratObj@meta.data$annotated_clusters_clean<-factor(seuratObj@meta.data$annotated_clusters_clean,as.character(Order_cellpop$Var1)) #reorder levels

U_annot<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, label.size = 4, group.by = "annotated_clusters_clean")
ggsave(U_annot, file=paste0(sampleFolder,"results/Annotation/2_UMAP_annotated_",sampleName,".png"), height = 15, width = 20, dpi = "retina")

## Annotated figure paper:
# Remove the bad clusters!! 
seuratObjNew<-subset(seuratObj, idents = c("Res cDC1s","Prolif Res cDC1s","Res cDC2s","T cells","Cd209a+ESAM- Res cDC2s","Mig cDC1s",             
                                           "NK cells","Cxcl9+ Res cDC1s","Ccr2+ Res cDC1s","Mig cDC2s","Prolif Res cDC2s",
                                           "B cells","RBCs","Doublets cDC1s/RBCs","pDCs","MCs","NFs","Plasma cells",
                                           "MFs","Basophils"))

# Rename clusters
levels(Idents(seuratObjNew))<-c("Resident cDC1s","Proliferating Resident cDC1s","ESAM+ Resident cDC2s","T cells",
                                "ESAM- Resident cDC2s","Migratory cDC1s","NK cells","Cxcl9+ Resident cDC1s",
                                "Ccr2+ Resident cDC1s","Migratory cDC2s","Proliferating Resident cDC2s",
                                "B cells","RBCs","Doublets cDC1s/RBCs","pDCs","Monocytes","Neutrophils","Plasma cells",
                                "Macrophages","Basophils")

# Save as pdf
## Bad clusters (low Qual + definite doublets) removed and clusters ordered according to number of cells (decreasing)
pdf(file=paste0(sampleFolder,"results/Annotation/4_Figure1a_paper_annotated_",sampleName,".pdf"), height = 13, width = 15)
DimPlot(seuratObjNew, reduction = "SCT_umap", label = T, repel = T, label.size = 5)
dev.off()

pdf(file=paste0("Paper/1a.Figure1_paper_annotated_",sampleName,".pdf"), height = 13, width = 15)
DimPlot(seuratObjNew, reduction = "SCT_umap", label = T, repel = T, label.size = 5)
dev.off()
