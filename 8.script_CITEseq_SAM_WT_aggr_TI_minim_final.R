# Script for further trajectory inference on SCT assay of cDC1 subset of WT_aggr object

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
object_counts <- Matrix::t(as(as.matrix(seuratObj@assays$SCT@counts), 'sparseMatrix'))
object_expression <- Matrix::t(as(as.matrix(seuratObj@assays$SCT@data), 'sparseMatrix'))

object_dyn<- wrap_expression(
  counts = object_counts, 
  expression = object_expression
)

#Add info of first cell (first attempt) 
#PAGA_Tree
object_dyn <- add_prior_information(
  object_dyn,
  start_id = First_cell
)

# #Add dimred
Cell_names<-rownames(object_counts)


#Add annotation info
Annotation<-as.matrix(seuratObj@meta.data$final_annotation2021) 
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

####################################################

#Start Dyno pipeline
guidelines <- guidelines_shiny(object_dyn)

# # Reproduces the guidelines as created in the shiny app
# answers <- dynguidelines::answer_questions(
#   multiple_disconnected = TRUE, 
#   expect_topology = NULL, 
#   expected_topology = NULL, 
#   n_cells = 18479, 
#   n_features = 18752, 
#   memory = "60GB", 
#   prior_information = c("start_id", "groups_id", "dimred"), 
#   docker = FALSE
# )
# guidelines <- dynguidelines::guidelines(answers = answers)

# Select methods proposed by guidelines (run in loop??)
methods_selected <- guidelines$methods_selected

# Reproducibility
set.seed(42)

# Try out with first method proposed
model <- infer_trajectory(object_dyn, first(methods_selected), verbose=T) #give_priors = c("start_id", "groups_id", "dimred"),
model_paga_tree<-model

Dimred<-as.matrix(seuratObj[['SCT_umap']]@cell.embeddings)

object_dyn <- add_dimred(
  object_dyn,
  Dimred
)

#Check results in model
model$milestone_network

head(model$progressions, 10)

##Choose method!!!
TI_method<-"paga_tree"

###############################################################

## Save model
saveRDS(model, file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/Minim_Non-regressed/TI_model_",TI_method,".rds"))

model<-readRDS(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/Minim_Non-regressed/PAGA_Tree/TI_model_paga_tree.rds"))

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

pdf(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/",TI_method,"_SCT_annotation_2021.pdf"), width = 15, height = 15)
P1
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

pdf(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/",TI_method,"_SCT_pseudotime_2021.pdf"), width = 15, height = 15)
P4
dev.off()

####################################################################

### Extra: Transform model! Simplify the trajectory
model_new<-model 

##Progressions
model2<- model$progressions %>% filter(to == 1) %>% filter(from == 3)
model2<- rbind(model2,model$progressions %>% filter(to == 2) %>% filter(from == 1))
model2<- rbind(model2,model$progressions %>% filter(to == 7) %>% filter(from == 2))
model2<- rbind(model2,model$progressions %>% filter(to == 8) %>% filter(from == 7))
model2<- rbind(model2,model$progressions %>% filter(to == 13) %>% filter(from == 8))

model_new$progressions<-model2

##Subset cell ids
cells_subset<-model2$cell_id

model_new$cell_ids<-model2$cell_id

cells_subset_extra<-paste0(cells_subset,".1") ##extra

##milestone_percentages (MP)
model2_MP<-model$milestone_percentages
model2_MP<-as.matrix(model2_MP)
rownames(model2_MP)<-model2_MP[,1]
model2_MP_filtered<-model2_MP[grep(cells_subset[1],rownames(model2_MP)),]
for (i in 2:length(cells_subset)) {
  model2_MP_filtered<-rbind(model2_MP_filtered, model2_MP[grep(cells_subset[i],rownames(model2_MP)),])
}
# model2_MP<-model2_MP[grep(cells_subset[i],rownames(model2_MP)),]
rownames(model2_MP_filtered)<-seq(1:nrow(model2_MP_filtered))
model2_MP_df<-as.data.frame(model2_MP_filtered)

model_new$milestone_percentages<-model2_MP_df

## Add Dimred
model_new$dimred<-Dimred[cells_subset,]

## Add Dimred_milestones
model_new$dimred_milestones<-model_new$dimred_milestones[c(1,2,3,7,8,13),]

## Add milestone_ids
model_new$milestone_ids<-c(1,2,3,7,8,13)

## Add milestone_network
model_new$milestone_network<-model$milestone_network[c(5,3,2,10,11),] 

##Add segment_progressions
model2_drsp<- model$dimred_segment_progressions %>% filter(to == 1) %>% filter(from == 3)
model2_drsp<- rbind(model2_drsp,model$dimred_segment_progressions %>% filter(to == 2) %>% filter(from == 1))
model2_drsp<- rbind(model2_drsp,model$dimred_segment_progressions %>% filter(to == 7) %>% filter(from == 2))
model2_drsp<- rbind(model2_drsp,model$dimred_segment_progressions %>% filter(to == 8) %>% filter(from == 7))
model2_drsp<- rbind(model2_drsp,model$dimred_segment_progressions %>% filter(to == 13) %>% filter(from == 8))

model_new$dimred_segment_progressions<-model2_drsp

## Add segment_points
model_new$dimred_segment_points<-model_new$dimred_segment_points[c(7:17,30:35,78:109),] #Filter
model_new$dimred_segment_points<-model_new$dimred_segment_points[c(12:17,8:11,1:7,18:49),] #Reorder

## Make linear: simplify trajectory!!
model_new<-simplify_trajectory(model_new, allow_self_loops = FALSE)

####################################################

## New
saveRDS(model_new,file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/Minim_Non-regressed/PAGA_Tree/TI_model_paga_tree_filtered_2021.rds"))
model_new<-readRDS(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Trajectory_inference/Minim_Non-regressed/PAGA_Tree/TI_model_paga_tree_filtered_2021.rds"))
