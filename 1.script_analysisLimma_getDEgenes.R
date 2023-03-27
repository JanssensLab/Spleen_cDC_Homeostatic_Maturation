## Script for processing Bulk RNA-Seq data
## Script includes other conditions besides just immature and tolerogenic mature splenic cDC1s
## Those conditions/comparisons were however beyond the scope of the manuscript and data was not shown.

library("edgeR")
library("limma")
library("ggplot2")
library("dplyr")

########################################
##### Functions
########################################

###Get DE genes
getDEgenes<-function(expMatrix, pValCutOff, logFCcutOff){
  topgenes<-expMatrix[expMatrix$adj.P.Val<pValCutOff,]
  genes_up<-topgenes[topgenes$logFC>logFCcutOff,]
  genes_down<-topgenes[topgenes$logFC< -logFCcutOff,]
  ##Sort genes on logFC
  genes_up<-genes_up[order(genes_up$logFC, decreasing=TRUE),]
  genes_down<-genes_down[order(genes_down$logFC, decreasing=TRUE),]
  genes_de_sorted<-rbind(genes_up, genes_down)
  
  return(genes_de_sorted)
}

###Normalize per gene
normalizePerGene<-function(expMatrix){
  resultMatrix<-t(apply(expMatrix, 1, function(x)(x-min(x))/(max(x)-min(x))))
  
  return(resultMatrix)
}

###Get ggplot colors
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

################################################################################
######### LOAD DATA
################################################################################

getwd()
setwd("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Victor/")
combatFolder<-"no_combat/"

##### Load raw counts
countData<-read.table(file="outputHTSeqCount/all_counts.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE)
dim(countData)
# 24421    36

### Load meta data
colData<-read.table(file="documentation/metadata.txt",sep="\t", header=TRUE, stringsAsFactors=TRUE)
rownames(colData)<-colData$fileName
colData<-colData[,-1]
dim(colData)
# 36  2

### Rename
cbind(colnames(countData), rownames(colData))
colnames(countData)<-rownames(colData)
dim(countData)
# 24421    36

################################################################################
########## CREATE OBJECT
################################################################################

y <- DGEList(counts = countData)

################################################################################
########## FILTER DATA
################################################################################

##### Filter low count genes
## always work with count-per-million (CPM) instead of raw counts
## Usually a gene is required to have a count of 5-10 in a library to be considered expressed in that library
## Imagine that the lowest lib size is around 6 million reads => threshold is set on CPM>1
## But what if lib size is 20 million? Here CPM of 1 means 20 counts. Do we still use 5 counts and then the CPM cut off
## will be 0.25 or is the threshold of 5 counts a cut off for lib sizes of around 5-6 million? Then we need to put the
## cut off on 20 for lib sizes around 20 million and so use a CPM of 1.
## Threshold needs to be true in at least x samples. x is always the lowest number of replicates.
## for example: 3 samples with each 2 replicates => x set on 2
## => This ensures that a gene will be retained if it is only expressed in both replicates of a certain group

## Do filtering
yNoFilter<-y
myCpm<-cpm(y)

keep = rowSums(cpm(y)>1) >= 4
y = y[keep,]
dim(y)
##11843
dim(yNoFilter)
##24421

##### Reset lib sizes
y$samples$lib.size = colSums(y$counts)
y$samples

################################################################################
########## NORMALISATION
################################################################################

##### Scale normalisation
yNoNormalisation<-y
y <- calcNormFactors(y)

##### Create colors for the samples
colorsUnique<-ggplotColours(8)
replicates<-c(4,4,4,4,4,4,6,6)

theColors<-c()
for(i in 1:length(colorsUnique)){
  theColors<-c(theColors, rep(colorsUnique[i],replicates[i]))
}

cbind(colData, theColors)

##### MDS-plot
plotMDS(y, cex=0.9, col=theColors)

plotMDS(y, pch=c(rep(1,8),rep(16,8),rep(6,8),rep(17,12)), col=as.factor(colData$condition))

#capture par settings, then add space to the right
opar <- par(no.readonly = TRUE)
par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 6))
#draw plot and legend
plotMDS(y, pch=c(rep(1,8),rep(16,8),rep(6,8),rep(17,12)), col=as.factor(colData$condition))
legend(par("usr")[2], par("usr")[4], c("Xbp1fl","XBP1flIRE1fl","XBP1ko","XBP1koIRE1ko"), pch = c(1,16,6,17), 
       col = c("black","red","green","blue"), bty = "n")
#set par back to original
par(opar)

################################################################################
########## LOG2 TRANSFORMATION
################################################################################

#### Create design matrix
TS <- paste(colData$cell, colData$condition, sep=".")
TS <- factor(TS, levels=unique(TS))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
design

##### Do voom
png(file="results/plots/meanVariancePlot.png",width = 1515, height = 1138)
v <- voom(y, design, plot = TRUE) # Transforms count data to logCPM + estimates Mean-Variance relationship (weights)
dev.off()

expTable<-v$E

##### Normalised counts
countData_norm<-cpm(y, normalized.lib.sizes=TRUE, log=FALSE)

##################################################
########## CORRECT WITH COMBAT
##################################################
library('sva')

combatFolder<-"with_combat/"

## 2 Samples performed on a different day
combatM<-colnames(expTable)
combatM<-cbind(combatM, c(rep("1",length(colnames(expTable)))))
combatM<-cbind(combatM, as.character(colData$condition))
test<-as.data.frame(combatM, stringsAsFactors=F)
rownames(test)<-test[,1]
test<-test[,-1]
test[grep('_(7|8)',rownames(test)),1]<-"2"
combatM<-test
colnames(combatM)<-c("batch","subgroups")
combatM

modcombat<-model.matrix(~subgroups, data=combatM)
expTable_combat<-ComBat(dat=expTable, batch=as.numeric(combatM$batch), mod=modcombat, par.prior=TRUE)

expTable<-expTable_combat

### MDS plot
plotMDS(expTable, cex=0.9, col=theColors)

### MDS plot version 2
png(file=paste0("results/plots/",combatFolder,"4b_MDSplot.png"), width = 1140, height = 850)
#capture par settings, then add space to the right
opar <- par(no.readonly = TRUE)
par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 6))
#draw plot and legend
plotMDS(expTable, pch=c(rep(1,8),rep(16,8),rep(6,8),rep(17,12)), col=as.factor(colData$condition))
legend(par("usr")[2], par("usr")[4], c("Xbp1fl","XBP1flIRE1fl","XBP1ko","XBP1koIRE1ko"), pch = c(1,16,6,17), 
       col = c("black","red","green","blue"), bty = "n")
#set par back to original
par(opar)
dev.off()


### Write results
write.table(expTable, "results/expTableCombat_clint.txt", sep="\t")
saveRDS(expTable, file="results/Robjects/expTable_afterCombat_clint.rds")

################################################################################
########## CHECK FILTERING AND NORMALISATION
################################################################################

#################### BARPLOT ####################

##### Barplot lib sizes raw counts
png(file="results/plots/1a_barplot_beforeNorm.png", width = 1515, height = 1138)
par(mar = c(9,3,3,1)) #more margin: bottom, left, top, right
bp<-barplot(yNoFilter$samples$lib.size*1e-6,axisnames=FALSE,main="Barplot lib sizes of raw counts",ylab="Library size (millions)")
axis(1, labels=rownames(yNoFilter$samples), at = bp, las=2, cex.axis=0.8)
dev.off()

##### Barplot lib sizes normalised counts
# countData_norm<-cpm(y, normalized.lib.sizes=TRUE, log=FALSE)

png(file="results/plots/1b_barplot_afterNorm.png", width = 1515, height = 1138)
par(mar = c(9,4,3,1)) #more margin: bottom, left, top, right
bp<-barplot(colSums(countData_norm)*1e-6,axisnames=FALSE,main="Barplot lib sizes of normalised counts",ylab="Library size (millions)")
axis(1, labels=colnames(countData_norm), at = bp, las=2, cex.axis=0.7)
dev.off()

#################### BOXPLOT ####################
col <- rainbow(nrow(colData))

y2<-y
y2$samples$norm.factors<-1
y2$counts[,1]<-ceiling(y2$counts[,1]*0.05)
y2$counts[,2]<-y2$counts[,2]*5

par(mfrow=c(1,2), mar = c(12,4,3,1)) #more margin: bottom, left, top, right
png(file="results/plots/2a_boxplot_beforeNorm.png", width = 1515, height = 1138)
lcpm<-cpm(y2, log=T)
boxplot(lcpm, las=2, main="Unnormalised data", ylab="log-cpm", col=col)
dev.off()

png(file="results/plots/2b_boxplot_afterNorm.png", width = 1515, height = 1138)
y2<-calcNormFactors(y2)
lcpm<-cpm(y2, log=T)
boxplot(lcpm, las=2, main="Normalised data", ylab="log-cpm", col=col)
dev.off()

#################### DENSITY PLOT ####################
col <- topo.colors(nrow(colData))

png(file="results/plots/3a_densityPlot.png", width = 1515, height = 1138)
par(mfrow=c(1,2))
### Plot log2-CPM values of each sample before filtering
theCpmNoFilter<-cpm(yNoFilter, log=TRUE)

plot(density(theCpmNoFilter[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="raw data", xlab="log cpm")
abline(v=0, lty=3)
ci = 2
for (i in 2:nrow(colData)){
  print(paste0(i,". Add line for sample ",colnames(theCpmNoFilter)[i]))
  
  den <- density(theCpmNoFilter[,i])
  lines(den$x, den$y, col=col[ci], lwd=2)
  ci = ci+1
}
# legend("topright", rownames(y$samples), text.col=col, cex=0.6, bty="n")

### Plot log2-CPM values of each sample after filtering (and normalisation)
theCpm<-cpm(y, log=TRUE)

plot(density(theCpm[,1]), col=col[1], lwd=2, ylim=c(0,0.25), las=2, main="filtered data", xlab="log cpm")
abline(v=0, lty=3)
ci = 2
for (i in 2:nrow(colData)){
  print(paste0(i,". Add line for sample ",colnames(theCpm)[i]))
  
  den <- density(theCpm[,i])
  lines(den$x, den$y, col=col[ci], lwd=2)
  ci = ci+1
}
dev.off()

par(mfrow=c(1,1))

#################### HISTOGRAM OF EXPTABLE ####################

png(file="results/plots/3b_histogramFiltering.png", width = 1515, height = 1138)
par(mfrow=c(1,2))
### Histogram
hist(expTable)

### Density plot
plot(density(expTable[,1]), col=col[1], lwd=2, ylim=c(0,0.25), las=2, main="Density of expTable", xlab="log2")
abline(v=0, lty=3)
ci = 2
for (i in 2:nrow(colData)){
  print(paste0("Add line for sample ",colnames(expTable)[i]))
  
  den <- density(expTable[,i])
  lines(den$x, den$y, col=col[ci], lwd=2)
  ci = ci+1
}
# legend("topright", rownames(y$samples), text.col=col, cex=0.6, bty="n")
dev.off()


par(mfrow=c(1,1))

################################################################################
########## PCA
################################################################################
library("rgl")

### Calculate variance
variance<-apply(expTable, 1, var)
varianceSorted<-sort(variance, decreasing=TRUE, index.return=TRUE)
### Get top 15%
numberOfGenes<-0.15*length(variance)
indexTopVariance<-varianceSorted$ix[1:numberOfGenes]
matrixPCAtmp<-expTable[indexTopVariance,]

### Prepare PCA-plot
pca<-prcomp(scale(t(matrixPCAtmp)))
matrixPCA<-cbind(pca$x[,1],pca$x[,2],pca$x[,3])
PCAcolors<-theColors

PoV <- pca$sdev^2/sum(pca$sdev^2)
summary(pca)

### Draw PCA-plot, all labels
pcaPlot<-plot3d(matrixPCA, main="",col=PCAcolors,pch=21, type="s", radius=2, legend=TRUE, xlab=paste0("pc1 (",round(PoV[1]*100,2),"%)"), 
                ylab=paste0("pc2 (",round(PoV[2]*100,2),"%)"), zlab=paste0("pc3 (",round(PoV[3]*100,2),"%)"))
text3d(x=matrixPCA[,1], y=(matrixPCA[,2]-2), z=(matrixPCA[,3]), rownames(matrixPCA) ,col=PCAcolors, cex=1)
# grid3d(c("x+", "y-", "z"))

rgl.viewpoint(0, 0)
rgl.snapshot(paste0("results/plots/",combatFolder,"4_pca_view1.png"))
rgl.viewpoint(35, 0)
rgl.snapshot(paste0("results/plots/",combatFolder,"4_pca_view2.png"))

### Save 3D
dirToSave<-paste0(getwd(),"/results/plots/",combatFolder)
writeWebGL(dir = dirToSave, filename = file.path(dirToSave, "4_pca.html"),
          template = system.file(file.path("WebGL", "template.html"), package = "rgl"),
          snapshot = TRUE, font = "Arial")


### Create nice plot
rgl.snapshot(paste0("results/plots/",combatFolder,"4_pca_withGrid.png"))

################################################################################
########## CORRELATION HEATMAP SAMPLES
################################################################################
library("RColorBrewer")
library("gplots")

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

##heatmap 1: based on distance
distsRL <- dist(t(expTable),method="euclidean")
hc <- hclust(distsRL,method="ward.D")

pdf(paste0("results/plots/",combatFolder,"4_corrSamples_distance_clint.pdf"))
heatmap.2(as.matrix(distsRL),
          Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col=rev(hmcol),margin=c(13, 13), cexRow=0.6, cexCol=0.6)
dev.off()

##heatmap 2: based on correlation
cm=cor(expTable)

distsRL <- as.dist(1-cor(expTable))
hc <- hclust(distsRL,method="ward.D")

pdf(paste0("results/plots/",combatFolder,"4_corrSamples_correlation.pdf"))
heatmap.2(cm, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col=hmcol, margin=c(13, 13), cexRow=0.6, cexCol=0.6)
dev.off()



################################################################################
########## GET DE GENES
################################################################################

head(design)

#### Fit linear model on data
fit <- lmFit(expTable, design)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=migcDC.XBP1ko-(migcDC.Xbp1fl+migcDC.XBP1flIRE1fl)/2,
                             group2=migcDC.XBP1koIRE1ko-(migcDC.Xbp1fl+migcDC.XBP1flIRE1fl)/2,
                             group3=rescDC.XBP1ko-(rescDC.Xbp1fl+rescDC.XBP1flIRE1fl)/2, 
                             group4=rescDC.XBP1koIRE1ko-(rescDC.Xbp1fl+rescDC.XBP1flIRE1fl)/2,
                             group5=migcDC.XBP1koIRE1ko-migcDC.XBP1ko, 
                             group6=migcDC.XBP1ko-rescDC.XBP1ko, 
                             group7=migcDC.XBP1koIRE1ko-rescDC.XBP1koIRE1ko, 
                             group8=(migcDC.Xbp1fl+migcDC.XBP1flIRE1fl)/2-(rescDC.Xbp1fl+rescDC.XBP1flIRE1fl)/2, #Only this comparison is used in the manuscript!!
                             group9=migcDC.XBP1ko-migcDC.XBP1koIRE1ko, 
                             group10=(migcDC.Xbp1fl+migcDC.XBP1flIRE1fl)/2-migcDC.XBP1ko, 
                             group11=(migcDC.Xbp1fl+migcDC.XBP1flIRE1fl)/2-migcDC.XBP1koIRE1ko,
                             group12=rescDC.XBP1koIRE1ko-rescDC.XBP1ko,
                             Diff=(migcDC.XBP1koIRE1ko-rescDC.XBP1koIRE1ko)-(((migcDC.Xbp1fl+migcDC.XBP1flIRE1fl)/2)-((rescDC.Xbp1fl+rescDC.XBP1flIRE1fl)/2)),
                             levels=design)
cont.matrix
fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Quick list of DE genes
summa.fit <- decideTests(fit.eb, p.value = 0.05, lfc = 1)
summary(summa.fit)

summa.fit <- decideTests(fit.eb, p.value = 0.01, lfc = 1)
summary(summa.fit)

volcanoplot(fit.eb,coef=1,highlight=100, names = rownames(fit.eb))
volcanoplot(fit.eb,coef=2,highlight=100, names = rownames(fit.eb))

########################################
##### Extract DE genes
########################################

##### Put DE genes in list
coef<-c(1:13)
listDEgenes<-list()
listDEgenes_moreStrict<-list()
listallgenes<-list()

for(i in 1:length(coef)){
  allGenes<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=coef[i])
  DEgenes<-getDEgenes(allGenes,0.05,1)
  DEgenes_strict<-DEgenes[DEgenes$AveExpr >0,]
  listallgenes[[i]]<-allGenes
  listDEgenes[[i]]<-DEgenes
  listDEgenes_moreStrict[[i]]<-DEgenes_strict
}

names(listDEgenes)<-c("migDC_Xbp1KO_vs_WT","migDC_Xbp1Ire1KO_vs_WT","resDC_Xbp1KO_vs_WT","resDC_Xbp1Ire1KO_vs_WT",
                      "migDC_Xbp1Ire1KO_vs_Xbp1KO", "Xbp1KO_mig_vs_res","Xbp1Ire1KO_mig_vs_res","WT_mig_vs_res",
                      "migDC_Xbp1KO_vs_Xbp1Ire1KO", "migDC_WT_vs_Xbp1KO","migDC_WT_vs_Xbp1Ire1KO",
                      "resDC_Xbp1Ire1KO_vs_Xbp1KO", "migcDC_XBP1koIRE1ko_vs_rescDC_XBP1koIRE1ko_VS_migcDC_WT_vs_rescDC_WT")

names(listDEgenes_moreStrict)<-names(listDEgenes)
names(listallgenes)<-names(listDEgenes)

##### Get numbers of DE genes
lapply(listDEgenes,dim)

#$WT_mig_vs_res
#[1] 1557    6

########################################
##### Write results
########################################

### Add geneSymbol in column (for the export)
listDEgenes<-lapply(listDEgenes,function(x){dplyr::mutate(x,'gene'=rownames(x))})
listallgenes<-lapply(listallgenes,function(x){dplyr::mutate(x,'gene'=rownames(x))})

### Write results
library('openxlsx')
write.xlsx(listDEgenes[c(1:8,12:13)], file = "results/summary_DEgenes_clint_oct2019.xlsx")
write.xlsx(listallgenes[c(1:8,12:13)], file = "results/summary_allgenes_clint_june2020.xlsx")

# ### Save results
saveRDS(listDEgenes, file="results/Robjects/listDEgenes_clint_oct2019.rds")
saveRDS(listallgenes, file="results/Robjects/listallgenes_clint_june2020.rds")

### Read results
listDEgenes<-readRDS(file="results/Robjects/listDEgenes_clint_oct2019.rds")
listallgenes<-readRDS(file="results/Robjects/listallgenes_clint_june2020.rds")
expTable<-readRDS(file="results/Robjects/expTable_afterCombat_clint.rds")

## Update 2023 for upload to GEO
countData<-read.table(file="outputHTSeqCount/all_counts.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE)
expTable<-readRDS(file="results/Robjects/expTable_afterCombat_clint.rds")

# Check names
all(colnames(countData) == colnames(expTable))

# Filter samples to those used in paper
cbind(colnames(countData)[grep("XBP1fl",colnames(countData))], colnames(expTable)[grep("XBP1fl",colnames(expTable))])
countData_paper<-countData[,grep("XBP1fl",colnames(countData))] #Without non-WT samples
expTable_paper<-expTable[,grep("XBP1fl",colnames(expTable))] #Without non-WT samples
cbind(colnames(countData_paper), colnames(expTable_paper))

write.table(countData_paper, "results/GEO/Bulk_RNA_seq_raw_counts_cDC1_tolerogenic_maturation.txt", sep="\t")
write.table(expTable_paper, "results/GEO/Bulk_RNA_seq_normalized_counts_cDC1_tolerogenic_maturation.txt", sep="\t")
