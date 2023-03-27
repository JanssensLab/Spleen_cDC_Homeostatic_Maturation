## Script for creating Bulk RNA-Seq heatmaps of tolerogenic mature vs immature splenic cDC1s 

library(dplyr)
library(sqldf)
library(ggplot2)

########################################
##### Functions
########################################

###Normalize per gene
normalizePerGene<-function(expMatrix){
  resultMatrix<-t(apply(expMatrix, 1, function(x)(x-min(x))/(max(x)-min(x))))
  
  return(resultMatrix)
}


################################################################################
######### LOAD DATA AND CREATE SUBSETS
################################################################################

getwd()
setwd("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Victor/")

###### eBayes ##########
listDEgenes<-readRDS(file="results/Robjects/listDEgenes_clint_oct2019.rds")

## Between mig and res ##
## WT
WT_mig_vs_res <- listDEgenes$WT_mig_vs_res$gene
length(WT_mig_vs_res) ##1557

# WT_mig_vs_res_top_30 <- head(listDEgenes$WT_mig_vs_res[order(listDEgenes$WT_mig_vs_res$adj.P.Val),], n=30) #Top 30 according to P-value
# WT_mig_vs_res_top_30 <- WT_mig_vs_res_top_30[order(WT_mig_vs_res_top_30$logFC, decreasing = T),"gene"] #Ordered by decreasing logFC
# 
# WT_mig_vs_res_ordered <- listDEgenes$WT_mig_vs_res[order(listDEgenes$WT_mig_vs_res$adj.P.Val),]  #Order all according to P-value
# WT_mig_vs_res_top_15 <- head(WT_mig_vs_res_ordered[which(WT_mig_vs_res_ordered$logFC > 0),], n=15) #Top 15 with logFC>0
# WT_mig_vs_res_bottom_15 <- head(WT_mig_vs_res_ordered[which(WT_mig_vs_res_ordered$logFC < 0),], n=15) #Top 15 with logFC<0
# WT_mig_vs_res_top_bottom_15 <- rbind(WT_mig_vs_res_top_15,WT_mig_vs_res_bottom_15) #Combine two dataframes
# WT_mig_vs_res_top_bottom_15 <- WT_mig_vs_res_top_bottom_15[order(WT_mig_vs_res_top_bottom_15$logFC, decreasing = T),"gene"] #Extract genes
# 
WT_mig_vs_res_logFC_top100 <- head(listDEgenes$WT_mig_vs_res, n=180) #Ordered by decreasing logFC
# WT_mig_vs_res_logFC_bottom100 <- tail(listDEgenes$WT_mig_vs_res, n=20) #Ordered by decreasing logFC
# 
# WT_mig_vs_res_top_bottom_100 <- rbind(WT_mig_vs_res_logFC_top100,WT_mig_vs_res_logFC_bottom100) #Combine two dataframes
# WT_mig_vs_res_top_bottom_100 <- WT_mig_vs_res_top_bottom_100[order(WT_mig_vs_res_top_bottom_100$logFC, decreasing = T),"gene"] #Extract genes

WT_mig_vs_res_logFC_top100_new<-head(WT_mig_vs_res_logFC_top100[order(WT_mig_vs_res_logFC_top100$logFC, decreasing = T),"gene"],100) #14/10/20: Only UP genes

##### expTable #####
expTable<-readRDS(file="results/Robjects/expTable_afterCombat_clint.rds")
dim(expTable)
# 11843    36

## Calculate mean ExpTable
colsMeanWT_Mig<-c(grep("fl_mig",colnames(expTable)))
colsMeanIre1KO_Mig<-grep("IRE1KO_mig",colnames(expTable))
colsMeanXpb1KO_Mig<-grep("XBP1KO_mig",colnames(expTable))

colsMeanWT_Res<-c(grep("fl_res",colnames(expTable)))
colsMeanIre1KO_Res<-grep("IRE1KO_res",colnames(expTable))
colsMeanXpb1KO_Res<-grep("XBP1KO_res",colnames(expTable))

WT_mean_Mig<-apply(expTable[,colsMeanWT_Mig],1,mean)
Ire1KO_mean_Mig<-apply(expTable[,colsMeanIre1KO_Mig],1,mean)
Xbp1KO_mean_Mig<-apply(expTable[,colsMeanXpb1KO_Mig],1,mean)

WT_mean_Res<-apply(expTable[,colsMeanWT_Res],1,mean)
Ire1KO_mean_Res<-apply(expTable[,colsMeanIre1KO_Res],1,mean)
Xbp1KO_mean_Res<-apply(expTable[,colsMeanXpb1KO_Res],1,mean)

expTable_mean<-cbind(WT_mean_Mig,WT_mean_Res,Ire1KO_mean_Mig,Ire1KO_mean_Res,Xbp1KO_mean_Mig,Xbp1KO_mean_Res)
colnames(expTable_mean)<-c('WT_Mig','WT_Res','Ire1KO_Mig','Ire1KO_Res','Xbp1KO_Mig','Xbp1KO_Res')

################################################################################
######### GET SAMPLES AND GENES
################################################################################

#### Get samples ####

## WT_mig_vs_res 
colsWT_res<-c(grep("XBP1fl_res",colnames(expTable)),grep("XBP1flIRE1fl_res",colnames(expTable)))
colsWT_mig<-c(grep("XBP1fl_mig",colnames(expTable)),grep("XBP1flIRE1fl_mig",colnames(expTable)))
wantedSamples<-c(colsWT_res,colsWT_mig)

#### Get genes ####

## WT_mig_vs_res 
# wantedGenes<-WT_mig_vs_res #eBayes all
# wantedGenes<-WT_mig_vs_res_top_30 #eBayes genes top 30
# wantedGenes<-WT_mig_vs_res_top_bottom_15 #eBayes genes top 15 (pos and neg logFC)
# wantedGenes<-efflux_DE #Victor efflux overlap Oct 2020
# wantedGenes<-effluxGenes #Victor efflux overlap Oct 2020
wantedGenes<-WT_mig_vs_res_logFC_top100_new #eBayes genes top UP 100 genes Oct 2020

################################################################################
######### Normalize and order
################################################################################

#### Normalize  ####
## All samples
expProfiles<-expTable[wantedGenes,wantedSamples] 

expProfilesNorm<-normalizePerGene(expProfiles)
gem<-apply(expProfiles,1,mean)
expProfilesNorm2<-expProfiles-gem

## All WT
colnames(expProfilesNorm2) <- c(rep("CCR7- cDC1s",8), 
                                rep("CCR7+ cDC1s",8))

################################################################################
######### Create heatmap
################################################################################

#### Prepare heatmap ####
library('pheatmap')
library('grid')
myColorPalette<-c("#08306b", "#08326e", "#083573", "#083876", "#083b7b", "#083d7f", "#084083", "#084387", "#08468c", "#084990", 
                  "#084b94", "#0a4f97", "#0b5199", "#0e559d", "#0f579f", "#125ba2", "#135da4", "#1661a7", "#1864aa", "#1967ad", 
                  "#1c6ab0", "#1f6db1", "#2372b4", "#2675b6", "#2a79b8", "#2e7ebb", "#3181bc", "#3685bf", "#3888c1", "#3d8dc3", 
                  "#4090c5", "#4794c7", "#4c98c9", "#539dcb", "#5ba1ce", "#60a5d0", "#67aad2", "#6cadd4", "#74b2d6", "#79b5d8", 
                  "#81badb", "#8bbfde", "#99c7e2", "#a8d0e6", "#b2d5e9", "#c1dded", "#cbe3f0", "#daebf4", "#e4f0f7", "#f3f8fb", 
                  "#fffefd", "#fff8f5", "#feefe9", "#fee9e1", "#fee0d4", "#fedacc", "#fdd1c0", "#fdcbb8", "#fdc2ac", "#fcb9a0", 
                  "#fcb398", "#fcab8f", "#fca68a", "#fc9e82", "#fc997c", "#fc9174", "#fc8c6e", "#fb8466", "#fb7d5d", "#fb7758", 
                  "#fb7050", "#f96a4d", "#f66348", "#f45e45", "#f15640", "#ee4e3b", "#ec4937", "#ea4133", "#e83c2f", "#e5342a", 
                  "#e32f27", "#dd2c25", "#d92924", "#d32622", "#cd2220", "#c9201f", "#c31d1d", "#bf1a1c", "#b9171a", "#b51419", 
                  "#ae1117", "#a91016", "#a00e15", "#980c14", "#920a13", "#890812", "#840711", "#7b0510", "#75030f", "#6d010e")

paletteLength<-100

#### Create heatmap ####
myBreaks <- c(seq(min(expProfilesNorm2), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(expProfilesNorm2)/paletteLength, max(expProfilesNorm2), length.out=floor(paletteLength/2)))

## Test Oct 2020 annotated heatmaps
p<-pheatmap(as.matrix(expProfilesNorm2),color=myColorPalette,scale="none", cluster_rows=F, cluster_cols=F,
            border_color = "gray25", cellwidth = 11, cellheight = 8, gaps_col = c(8), fontsize=8, #gaps_row = c(342)
            show_rownames = T, show_colnames = T, breaks = myBreaks)
# ggsave(p,file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Victor/results/heatmaps/migvsres/2020/DE_heatmap.png"),dpi=300, height = 5, limitsize = FALSE)


# use this function to make row or column names bold
# parameters:
#   mat: the matrix passed to pheatmap
#   rc_fun: either rownames or colnames
#   rc_names: vector of names that should appear in boldface
library(tidyverse)
make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}

# At the gene level, well-known markers of dendritic cell maturation were found (Cd70, H2-M2, Cd40, …), 
# as well as typical genes associated with cell migration and cytoskeletal rearrangements (Ccr7, Fscn1, Arc, Synpo2, Marcksl1,…) 
# or intracellular vesicle transport (Cd63, Mreg…) (Suppl. 1C). 
# Several top DE genes appeared associated with the TGF pathway, such as Scube3, a TGFRII ligand {Wu:2011cg}, 
# or the integrin Itgb8 crucial for activation of TGF in Tregs {Travis:2007fk}, which -according to literature {Ardouin:2016gz}-, 
# are specifically induced in homeostatic conditions. Together with genes such as Ccl22 (a chemokine needed for communication 
# of dendritic cells with CCR4+ Tregs {Rapp:2019gy}) or Aldh1a2, catalyzing retinoic acid synthesis, this couples homeostatic 
# maturation with maintaining and inducing peripheral Treg formation. In addition, several genes were found known to be linked with 
# immunosuppression such as Socs2, IL4i1, Arg2 or Cd274, but remarkably not all of them are unique to the homeostatic maturation program 
# {Ardouin:2016gz}(Suppl. 1C). 


siggenes_DCmat<-c("Cd70", "H2-M2") #, "Cd40"
siggenes_migration<-c("Ccr7", "Fscn1", "Arc", "Synpo2") #, "Marcksl1"
siggenes_vesicle<-c("Cd63", "Mreg","Pla2gf4","Kif26b")
siggenes_Tgfb<-c("Scube3","Itgb8","Cdkn2b","Foxh1") #,"Aldh1a2"
siggenes_imm_suppr<-c("Oprd1", "Jag1","Socs2","Il4i1","Arg2","Adora2a","Ccl22") #,"Cd274"

# siggenes_imm_suppr %in% rownames(expProfilesNorm2)

### Row metadata
categories_top100<-c("DC Maturation","DC Migration","Organelle Trafficking","TGFb Signalling","Immunosuppression")

rowMeta_df_top100 <- data.frame(Function = rep("Other", 100), 
                                stringsAsFactors = F,
                                row.names = rownames(as.matrix(expProfilesNorm2)))
for (gene_v in siggenes_DCmat) rowMeta_df_top100[rownames(rowMeta_df_top100) == gene_v, "Function"] <- "DC Maturation"
for (gene_v in siggenes_migration) rowMeta_df_top100[rownames(rowMeta_df_top100) == gene_v, "Function"] <- "DC Migration"
for (gene_v in siggenes_vesicle) rowMeta_df_top100[rownames(rowMeta_df_top100) == gene_v, "Function"] <- "Organelle Trafficking"
for (gene_v in siggenes_Tgfb) rowMeta_df_top100[rownames(rowMeta_df_top100) == gene_v, "Function"] <- "TGFb Signalling"
for (gene_v in siggenes_imm_suppr) rowMeta_df_top100[rownames(rowMeta_df_top100) == gene_v, "Function"] <- "Immunosuppression"

### Make color list
library(RColorBrewer)
colors_top100 <- c(brewer.pal(5, "Set1"),"#A6A6A6")
colors_top100 <- colors_top100[c(1:length(categories_top100), 6)]
names(colors_top100) <- c(categories_top100, "Other")
annColors_top100 <- list("Function" = colors_top100)

p<-pheatmap(as.matrix(expProfilesNorm2),color=myColorPalette,scale="none", cluster_rows=F, cluster_cols=F,
         border_color = "gray25", cellwidth = 11, cellheight = 8, gaps_col = c(8), fontsize=8, #gaps_row = c(342)
         show_rownames = T, show_colnames = T, breaks = myBreaks,
         labels_row = make_bold_names(as.matrix(expProfilesNorm2), rownames, c(siggenes_DCmat, siggenes_migration, siggenes_vesicle,
                                                                               siggenes_Tgfb,siggenes_imm_suppr)),
         annotation_row = rowMeta_df_top100, annotation_colors = annColors_top100, main = "Top 100 DE genes Migratory vs Resident cDC1s")

ggsave(p,file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Victor/results/heatmaps/migvsres/2020/Annotated_heatmap_Mig_vs_Res_top100.png"),dpi=300, height = 15, limitsize = FALSE)

cairo_pdf(file=paste0("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Victor/results/heatmaps/migvsres/2020/Annotated_heatmap_Mig_vs_Res_top100.pdf"), height = 15)
pheatmap(as.matrix(expProfilesNorm2),color=myColorPalette,scale="none", cluster_rows=F, cluster_cols=F,
         border_color = "gray25", cellwidth = 11, cellheight = 8, gaps_col = c(8), fontsize=8, #gaps_row = c(342)
         show_rownames = T, show_colnames = T, breaks = myBreaks,
         labels_row = make_bold_names(as.matrix(expProfilesNorm2), rownames, c(siggenes_DCmat, siggenes_migration, siggenes_vesicle,
                                                                               siggenes_Tgfb,siggenes_imm_suppr)),
         annotation_row = rowMeta_df_top100, annotation_colors = annColors_top100, main = "Top 100 DE genes Migratory vs Resident cDC1s")
dev.off()
