# Script for further processing of cDC2 subset of WT_aggr object

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

sampleName <- "SAM2and3_WT_subset2_v2" 
sampleFolder<-"SAM2and3_WT/"

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')
source('~/VIB/DATA/Roos/Daan 1/script_featurePlots.R')

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

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results_subset2_v2/QC/1_UMI_",sampleName,".png"), width = 20)

########## mito.genes plot ##########
p1<-drawUMI_mitoPlot_new(tsneTable, 'tsne', clusterMatrix, 'subsets_Mito_percent',"mito")
p2<-drawUMI_mitoPlot_new(umapTable, 'umap', clusterMatrix, 'subsets_Mito_percent',"mito")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results_subset2_v2/QC/2_percMito_",sampleName,".png"), width = 20)


########## PCA plot ##########
# pdf(file=paste0(sampleFolder,"results_subset2_v2/QC/13a_PCA_",sampleName,".pdf"), width=10)
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
dir.create(paste0(sampleFolder,"results_subset2_v2/Annotation/"))

pdf(file=paste0(sampleFolder,"results_subset2_v2/Annotation/1_annotation_",sampleName,".pdf"), width = 15)
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
pdf(file=paste0(sampleFolder,"results_subset2_v2/Annotation/1_annotation_Color_ADT_clusters_on_SCT_UMAP_",sampleName,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$ADT_clusters))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$ADT_clusters==i),])))
  C1<-C1+ggtitle(paste0("ADT_cluster_",i))
  print(C1)
}
dev.off()

#SCT clusters
pdf(file=paste0(sampleFolder,"results_subset2_v2/Annotation/1_annotation_Color_SCT_clusters_on_SCT_UMAP_",sampleName,".pdf"), width = 15)
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

dir.create(paste0(sampleFolder,"results_subset2_v2/Marker_lists"))

Idents(seuratObj)<-seuratObj@meta.data$ADT_clusters

### Find ADTmarkers for every ADT cluster compared to all remaining cells, report only the positive ones
ADTMarkers_ADTclus <- FindAllMarkers(seuratObj, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_ADTclus$cluster)
saveRDS(ADTMarkers_ADTclus, file=paste0(sampleFolder,"results_subset2_v2/Robjects/ADTmarkersList_ADTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerADTcluster']]<-paste0(table(ADTMarkers_ADTclus$cluster)," ADT markers for ADT cluster ",rownames(table(ADTMarkers_ADTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset2_v2/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
write.xlsx(ADTmarkersList_ADTclus, file =paste0(sampleFolder, "results_subset2_v2/Marker_lists/ADTmarkersList_ADTclus_",sampleName,".xlsx"))

########################################

### Find SCTmarkers for every ADT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_ADTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_ADTclus$cluster)
saveRDS(SCTMarkers_ADTclus, file=paste0(sampleFolder,"results_subset2_v2/Robjects/SCTmarkersList_ADTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerADTcluster']]<-paste0(table(SCTMarkers_ADTclus$cluster)," SCT markers for ADT cluster ",rownames(table(SCTMarkers_ADTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset2_v2/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
write.xlsx(SCTmarkersList_ADTclus, file =paste0(sampleFolder, "results_subset2_v2/Marker_lists/SCTmarkersList_ADTclus_",sampleName,".xlsx"))

########################################

### Find ADTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
Idents(seuratObj)<-seuratObj@meta.data$SCT_clusters

ADTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "ADT", only.pos = TRUE)
table(ADTMarkers_SCTclus$cluster)
saveRDS(ADTMarkers_SCTclus, file=paste0(sampleFolder,"results_subset2_v2/Robjects/ADTmarkersList_SCTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['ADTmarkersPerSCTcluster']]<-paste0(table(ADTMarkers_SCTclus$cluster)," ADT markers for SCT cluster ",rownames(table(ADTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset2_v2/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
write.xlsx(ADTmarkersList_SCTclus, file =paste0(sampleFolder, "results_subset2_v2/Marker_lists/ADTmarkersList_SCTclus_",sampleName,".xlsx"))

########################################

### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
SCTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results_subset2_v2/Robjects/SCTmarkersList_SCTclus_",sampleName,".rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTcluster']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset2_v2/Robjects/diagnostics_",sampleName,"_clint.rds"))

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
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results_subset2_v2/Marker_lists/SCTmarkersList_SCTclus_",sampleName,".xlsx"))

########################################################################################################################

########################################
##### Markers annotated clusters
########################################

seuratObj@meta.data$annotated_clusters <- seuratObj@active.ident 
levels(seuratObj@meta.data$annotated_clusters) <- c("ESAM+ Res cDC2s","ESAM+ Res cDC2s","ESAM- Res cDC2s","Mig cDC2s",
                                                    "Ccl3+Ccl4+ESAM+ Res cDC2s","Sox4+ESAM- Res cDC2s",'Prolif Res cDC2s 1','Prolif Res cDC2s 2',
                                                    "Cd209a+CD301a+ESAM- Res cDC2s","Ifitm1+ESAM+ Res cDC2s","Hsp+ Res cDC2s",
                                                    "Unknown Res cDC2s","Mig cDC2s")
U_annot<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, group.by = "annotated_clusters", label.size = 4)
ggsave(U_annot, file=paste0(sampleFolder,"results_subset2_v2/Annotation/2_UMAP_annotated_",sampleName,".png"), height = 10, width = 20, dpi = "retina")


seuratObj@meta.data$minimalist_annotation <- seuratObj@active.ident 
levels(seuratObj@meta.data$minimalist_annotation) <- c("ESAM+ Res cDC2s","ESAM+ Res cDC2s","ESAM- Res cDC2s","Mig cDC2s",
                                                      "ESAM+ Res cDC2s","ESAM- Res cDC2s",'Prolif Res cDC2s','Prolif Res cDC2s',
                                                      "ESAM- Res cDC2s","ESAM+ Res cDC2s","Hsp+ Res cDC2s",
                                                      "Unknown Res cDC2s","Mig cDC2s")
U_annot<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, group.by = "minimalist_annotation", label.size = 4)
ggsave(U_annot, file=paste0(sampleFolder,"results_subset2_v2/Annotation/2_UMAP_minim_annotated_",sampleName,".png"), height = 10, width = 20, dpi = "retina")


### Find SCTmarkers for every SCT cluster compared to all remaining cells, report only the positive ones
Idents(seuratObj) <- seuratObj@meta.data$annotated_clusters

SCTMarkers_SCTclus <- FindAllMarkers(seuratObj, assay = "SCT", only.pos = TRUE)
table(SCTMarkers_SCTclus$cluster)
saveRDS(SCTMarkers_SCTclus, file=paste0(sampleFolder,"results_subset2_v2/Robjects/SCTmarkersList_SCTclus_",sampleName,"_annotated.rds"))

### Add to diagnostics
diagnostics[['SCTmarkersPerSCTclusterannotated']]<-paste0(table(SCTMarkers_SCTclus$cluster)," SCT markers for SCT cluster ",rownames(table(SCTMarkers_SCTclus$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset2_v2/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

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
write.xlsx(SCTmarkersList_SCTclus, file =paste0(sampleFolder, "results_subset2_v2/Marker_lists/SCTmarkersList_SCTclus_",sampleName,"_annotated.xlsx"))

########################################################################################################################

########################################################
####### Create final object: subset + reannotate #######
########################################################

## Final annotation paper 2021
Idents(seuratObj)<-seuratObj@meta.data$minimalist_annotation 

## Remove unknown cluster: low quality markers + Mito markers!!!!!!!!!!!
seuratObj<-subset(seuratObj, idents = c("ESAM+ Res cDC2s","ESAM- Res cDC2s","Mig cDC2s","Prolif Res cDC2s","Hsp+ Res cDC2s"))

## Refactor levels
levels(Idents(seuratObj)) <- c("ESAM+ Resident cDC2s","ESAM- Resident cDC2s","Migratory cDC2s","Proliferating cDC2s","ESAM- Resident cDC2s")

seuratObj@meta.data$final_annotation2021<-Idents(seuratObj)

## Split off clear pre-cDC2s (based on Vim, Ccr2, S100a10 expression!!)
colorSomeCells(clusterMatrix, umapTable, 
               WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$annotated_clusters=="Cd209a+CD301a+ESAM- Res cDC2s"),])))

seuratObj@meta.data$final_annotation2021<-as.character(seuratObj@meta.data$final_annotation2021)
seuratObj@meta.data[WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$annotated_clusters=="Cd209a+CD301a+ESAM- Res cDC2s"),])),"final_annotation2021"] <- "pre-cDC2s"

## Reorder levels
seuratObj@meta.data$final_annotation2021<-factor(seuratObj@meta.data$final_annotation2021,
                                                 as.character(c("pre-cDC2s","Proliferating cDC2s","ESAM- Resident cDC2s","ESAM+ Resident cDC2s",
                                                                "Migratory cDC2s"))) #reorder levels

Idents(seuratObj)<-seuratObj@meta.data$final_annotation2021

## Change names annotation 
levels(Idents(seuratObj)) <- c("pre-cDC2s","Proliferating cDC2s","ESAM- immature cDC2s","ESAM+ immature cDC2s","Mature cDC2s")
seuratObj@meta.data$final_annotation2021<-Idents(seuratObj)

## Create annotated UMAP
U_annot<-DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, label.size = 5)

pdf(file = paste0(sampleFolder,"results_subset2_v2/Annotation/FigureS7_paper_annotated_",sampleName,"_2021.pdf"), height = 7, width = 10)
U_annot
dev.off()

## Create pdf with important markers
library(viridis)
features<-c("Pdia3", "Pdia6", "Hspa5" , "Calr", "Apol7c", "Apoe", "Dnase1l3","Cxcl10", "Cxcl9", "Iigp1", "Ifi47",
            "Isg15","CD62L")

pdf(file=paste0(sampleFolder,"results_subset2_v2/Feature_plots/Paper/Feature_plots/FigureS7_Feature_plot_paper_all_markers_",sampleName,"_blue_grey.pdf"), height = 7, width = 10)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features =feature, cols = c("grey", "blue"), 
                  reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) #+
    # scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

## Create pdf with extra marker 
features<-c("Dnajc3")

pdf(file=paste0(sampleFolder,"results_subset2_v2/Feature_plots/Paper/Feature_plots/FigureS7_Feature_plot_paper_extra_marker_",sampleName,"_blue_grey.pdf"), height = 7, width = 10)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features =feature, cols = c("grey", "blue"), 
                  reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) #+
  # scale_color_viridis(option = "C")
  print(F1)
}
dev.off()
########################################

##### Read object
seuratObj <- readRDS(file=paste0(sampleFolder,"results_subset2_v2/Robjects/seuratObj_",sampleName,"_final2021.rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"results_subset2_v2/Robjects/diagnostics_",sampleName,"_final2021.rds"))

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"results_subset2_v2/Robjects/seuratObj_",sampleName,"_final2021.rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"results_subset2_v2/Robjects/diagnostics_",sampleName,"_final2021.rds"))


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

## Filter seperate parts of seuratObj!
seuratObj_filtered<-seuratObj
seuratObj_filtered@assays$ADT@counts<-seuratObj_filtered@assays$ADT@counts[Other_panel_filtered2$whitelist_name,]
seuratObj_filtered@assays$ADT@data<-seuratObj_filtered@assays$ADT@data[Other_panel_filtered2$whitelist_name,]
seuratObj_filtered@assays$ADT@scale.data<-seuratObj_filtered@assays$ADT@scale.data[intersect(seuratObj_filtered@assays[["ADT"]]@var.features,Other_panel_filtered2$whitelist_name),]
seuratObj_filtered@assays$ADT@var.features<-intersect(seuratObj_filtered@assays[["ADT"]]@var.features,Other_panel_filtered2$whitelist_name)
seuratObj_filtered@assays$ADT@meta.features<-seuratObj_filtered@assays$ADT@meta.features[Other_panel_filtered2$whitelist_name,]

##### Save object
saveRDS(seuratObj_filtered, file=paste0(sampleFolder,"results_subset2_v2/Robjects/seuratObj_filteredADT_",sampleName,"_final2021.rds"))

##### Read object
seuratObj_filtered<-readRDS(file=paste0(sampleFolder,"results_subset2_v2/Robjects/seuratObj_filteredADT_",sampleName,"_final2021.rds"))

##########################################################################################

## Create diet object for online tool 
seuratObj_filtered_diet<-DietSeurat(seuratObj_filtered, counts = T, data = T, scale.data = F,
                                    assays = c("SCT","ADT"), dimreducs = "SCT_umap", graphs = NULL)

## Metadata columns
# "orig.ident" "SCT_clusters"  "ADT_clusters" "annotated_clusters" "final_annotation2021" 
DimPlot(seuratObj_filtered_diet, reduction = "SCT_umap", label = T, repel = T, group.by = "final_annotation2021", label.size = 5,
        cols = gg_color_hue(5))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

All<-c(levels(as.factor(seuratObj_filtered_diet$orig.ident)),
       levels(as.factor(seuratObj_filtered_diet$SCT_clusters)), #Remove ADT (also numbers like SCT!)
       levels(as.factor(seuratObj_filtered_diet$annotated_clusters)),levels(as.factor(seuratObj_filtered_diet$final_annotation2021)))
Filtered_wrong_order<-c(levels(as.factor(seuratObj_filtered_diet$orig.ident)),
                        levels(as.factor(as.character(seuratObj_filtered_diet$SCT_clusters))),
                        levels(as.factor(as.character(seuratObj_filtered_diet$annotated_clusters))),levels(as.factor(seuratObj_filtered_diet$final_annotation2021)))
Metadata_info<-c(intersect(All,Filtered_wrong_order),levels(as.factor(as.character(seuratObj_filtered_diet$ADT_clusters))))
Color_info<-c(gg_color_hue(2),gg_color_hue(12),gg_color_hue(10),
              gg_color_hue(5), gg_color_hue(8))
Metadata_column<-c(rep("orig.ident",2),rep("SCT_clusters",12),rep("annotated_clusters",10),
                   rep("final_annotation2021",5),rep("ADT_clusters",8))
Info_Kevin<-as.data.frame(cbind(Metadata_info,Color_info,Metadata_column))

write.xlsx(Info_Kevin, file =paste0(sampleFolder, "results_subset2_v2/Annotation/Info_Kevin_cDC2_object",sampleName,".xlsx"))

## Change to name metadata columns tool (performed after upload to tool Kevin, but before upload to GEO)
# "Annotation", "Annotation (detailed)", "Sample", "Phase", "Clustering cDNA", "Clustering ADT"
seuratObj_filtered_diet$Annotation<-seuratObj_filtered_diet$final_annotation2021
seuratObj_filtered_diet$Annotation_detailed<-seuratObj_filtered_diet$annotated_clusters
seuratObj_filtered_diet$Sample<-seuratObj_filtered_diet$orig.ident
seuratObj_filtered_diet$Clustering_cDNA<-seuratObj_filtered_diet$SCT_clusters
seuratObj_filtered_diet$Clustering_ADT<-seuratObj_filtered_diet$ADT_clusters

##### Save object
saveRDS(seuratObj_filtered_diet, file=paste0(sampleFolder,"results_subset2_v2/Robjects/seuratObj_paper_diet",sampleName,"_2023.rds"))

##### Read object
seuratObj_filtered_diet<-readRDS(file=paste0(sampleFolder,"results_subset2_v2/Robjects/seuratObj_paper_diet",sampleName,"_2023.rds"))

##########################################################################################

## Check Dnase1l3 rebuttal 
F2<-FeaturePlot(object = seuratObj_filtered, features = c("Dnase1l3"), cols = c("grey", "blue"), 
                reduction = "SCT_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
V2<-VlnPlot(object = seuratObj_filtered, features = c("Dnase1l3"))

pdf(file=paste0(sampleFolder,"results_subset2_v2/Feature_plots/Feature_and_Violin_plot_rebuttal_",sampleName,"_blue_grey.pdf"), height = 14, width = 10)
F2/V2
dev.off()

###########################