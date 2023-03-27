# Script for further trajectory inference on ADT assay of cDC1 subset of WT_aggr object

library(dyno)
library(tidyverse)
library(Seurat)

setwd("~/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/")

sampleName <- "SAM2and3_WT_subset"
sampleFolder<-"SAM2and3_WT/"

########################################
##### Some variables
########################################
listLabels<-list('SAM2','SAM3')

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')
source('~/VIB/DATA/Roos/Daan 1/script_featurePlots.R')

##### Read object
seuratObj<-readRDS(file=paste0(sampleFolder,"results_subset/Robjects/seuratObj_final_",sampleName,"_clint_2021.rds")) #New annotation 2021, but no other changes!!
diagnostics <- readRDS(file=paste0(sampleFolder,"results_subset/Robjects/diagnostics_sliced_",sampleName,"_clint.rds"))

DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, label.size = 4)

#########################################################################################

#Create cluster matrix
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['SCT_tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['SCT_umap']]@cell.embeddings, stringsAsFactors = F)

##Find first cell (based on expression markers: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3858649/)
C1<-WhichCells(object = seuratObj, expression = adt_CD83 == 0)
C2<-WhichCells(object = seuratObj, expression = adt_CD80 == 0)
C3<-WhichCells(object = seuratObj, expression = adt_CD86 == 0)
C4<-WhichCells(object = seuratObj, expression = `adt_IA-IE` < 0.3)
C5<-WhichCells(object = seuratObj, expression = Sell > 0.5)
C6<-WhichCells(object = seuratObj, expression = adt_CD62L > 0.5)
C7<-WhichCells(object = seuratObj, idents = "pre-cDC1s")

Cells<-intersect(intersect(intersect(intersect(intersect(intersect(C1,C2),C3),C4),C5),C6),C7)
colorSomeCells(clusterMatrix, umapTable, Cells)

First_cell<-Cells

#Convert seuratobj to Dyno object
object_counts <- Matrix::t(as(as.matrix(seuratObj@assays$ADT@counts), 'sparseMatrix'))
object_expression <- Matrix::t(as(as.matrix(seuratObj@assays$ADT@data), 'sparseMatrix'))

object_dyn<- wrap_expression(
  counts = object_counts, 
  expression = object_expression
)

#Add info of first cell 
#PAGA_Tree
object_dyn <- add_prior_information(
  object_dyn,
  start_id = First_cell
)

# #Add dimred
Cell_names<-rownames(object_counts)

#Add annotation info
Annotation<-as.matrix(seuratObj@meta.data$final_annotation2021) #2021!!
rownames(Annotation)<-Cell_names

object_dyn <- add_grouping(
  object_dyn,
  Annotation
)

## Reorder levels
object_dyn$grouping<-factor(object_dyn$grouping,
                            as.character(c("pre-cDC1s","Proliferating cDC1s","Early immature cDC1s",
                                           "Late immature cDC1s","Early mature cDC1s","Late mature cDC1s")))
object_dyn$group_ids<-levels(object_dyn$grouping)
  
#Start Dyno pipeline
guidelines <- guidelines_shiny(object_dyn)

# # Reproduces the guidelines as created in the shiny app
# answers <- dynguidelines::answer_questions(
#   multiple_disconnected = FALSE, 
#   expect_topology = TRUE, 
#   expected_topology = "linear", 
#   n_cells = 18479, 
#   n_features = 190, 
#   time = "16h", 
#   memory = "60GB", 
#   prior_information = "start_id", 
#   docker = FALSE
# )
# guidelines <- dynguidelines::guidelines(answers = answers) 

# Select methods proposed by guidelines (run in loop??)
methods_selected <- guidelines$methods_selected

# Reproducibility
set.seed(42)

# Try out with first method proposed
model <- infer_trajectory(object_dyn, methods_selected[4], verbose=T) #give_priors = c("start_id", "groups_id", "dimred"),

Dimred<-as.matrix(seuratObj[['SCT_umap']]@cell.embeddings)

object_dyn <- add_dimred(
  object_dyn,
  Dimred
)

#Check results in model
model$milestone_network

head(model$progressions, 10)

##Choose method!!!
TI_method<-"scorpius_ADT"

###############################################################

## Root the trajectory (add arrow)
model<-add_root(model, root_milestone_id = "milestone_begin")

###############################################################

## Save model
saveRDS(model, file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/Minim_Non-regressed/Scorpius_ADT/TI_model_",TI_method,"_2021.rds"))

model<-readRDS(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/Minim_Non-regressed/Scorpius_ADT/TI_model_scorpius_ADT_2021.rds"))

################################################################
library(RColorBrewer)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## Rearrange colors
Colset<-brewer.pal(n = 12, name = "Set3")[c(3,5,4,6,1,2)]

#Plot grouping
P1<-plot_dimred(
  model,
  dimred = Dimred,
  color_cells = "grouping",
  expression_source = object_dyn$expression, 
  grouping = object_dyn$grouping) + 
  ggtitle("Cell grouping") + #scale_fill_manual(values =gg_color_hue(length(levels(seuratObj@meta.data$annotated_clusters))))
  scale_fill_manual(values = Colset)

ggsave(P1, file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/",TI_method,"_annotation.png"), width = 15, height = 15)

pdf(paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/",TI_method,"_annotation_2021.pdf"), width = 15, height = 15)
plot_dimred(
  model,
  dimred = Dimred,
  color_cells = "grouping",
  expression_source = object_dyn$expression, 
  grouping = object_dyn$grouping) + 
  ggtitle("Cell grouping") + #scale_fill_manual(values =gg_color_hue(length(levels(seuratObj@meta.data$annotated_clusters))))
  scale_fill_manual(values = Colset)
dev.off()

#Plot expression
P2<-plot_dimred(
  model, 
  dimred = Dimred,
  expression_source = object_dyn$expression,
  feature_oi = "Ccr7"
) + ggtitle("Feature expression")

ggsave(P2, file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/",TI_method,"_Ccr7_expression_new2.png"), width = 15, height = 15)

# #Coloring by milestone
P3<-plot_dimred(model, dimred = Dimred, label_milestones = T) + ggtitle("Cell ordering")
ggsave(P3, file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/",TI_method,"_milestones_new2.png"), width = 15, height = 15)

#Pseudotime
P4<-plot_dimred(model, "pseudotime", dimred = Dimred, pseudotime = calculate_pseudotime(model)) + ggtitle("Pseudotime")

P4[["layers"]][[2]][["aes_params"]][["colour"]]<-"#ff0000"
P4[["layers"]][[3]][["aes_params"]][["colour"]]<-"#ff0000"
P4[["layers"]][[4]][["aes_params"]][["colour"]]<-"#ff0000"

pdf(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/",TI_method,"_pseudotime_2021.pdf"), width = 15, height = 15)
P4
dev.off()
