# script created by KSB, 06.06.18
# Perform DE analysing using RUVs vs a traditional linear modelling approach


# Load dependencies and set input paths --------------------------

# Load dependencies:
library(edgeR)
library(plyr)
library(NineteenEightyR)
library(RColorBrewer)
library(biomaRt)
library(ggpubr)
library(ggplot2)
library(ggsignif)
library(pheatmap)
library(viridis)
library(gplots)
library(circlize)
library(ComplexHeatmap)
library(Homo.sapiens)
require(VennDiagram)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(eulerr)
library(reshape2)

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/"

# Set output directory and create it if it does not exist:
outputdir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/RUVvsLinearModel/"

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir)
}

# set up colour palette. The "wes" palette will be used for island and other statistical information, whereas NineteenEightyR will be used for batch information
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
# set up colour palette for batch
batch.col=electronic_night(n=3)
dev.off()


# Begin Analysis --------------------------------------------------------------------------------------------------

# Load the RUVs-corrected set1 object
load(paste0(inputdir, "RUVs_Setup/set1_RUVsCorrectedObject.Rda"))
# load the filtered, normalised y DGE list object
load(paste0(inputdir, "dataPreprocessing/indoRNA.read_counts.TMM.filtered.Rda"))
# lcpm
load(paste0(inputdir, "dataPreprocessing/indoRNA.logCPM.TMM.filtered.Rda"))

# define sample names
samplenames <- as.character(y$samples$samples)
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
samplenames <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)

# define variables
batch=y$samples$batch
Island=y$samples$Island

# get which samples are replicates
load(paste0(inputdir, "dataPreprocessing/covariates.Rda"))
allreps=covariates[,1][which(covariates$replicate)]
allreps=unique(allreps)
allreplicated=as.factor(samplenames %in% allreps)

# get which samples are replicates - we need this information to label replicates in the PCA plots
allreplicated=as.factor(samplenames %in% allreps)

## RUVs ----------------------------------------

# load z, the normalised count matrix
load(paste0(inputdir, "DE_Island/RUVs/z_UQNormalised.Rda"))

# create a new variable for blocking using sample IDs
z$samples$ind <- samplenames

# setup RUV design
design.RUV <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set1))
colnames(design.RUV)=gsub("Island", "", colnames(design.RUV))
colnames(design.RUV)[3]="Mappi"
contr.matrix.RUVs <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design.RUV))
v.RUV <- voom(z, design.RUV, plot=FALSE)
dupcor <- duplicateCorrelation(v.RUV, design.RUV, block=z$samples$ind)
# run voom a second time
vDup <- voom(z, design.RUV, plot=TRUE, block=z$samples$ind, correlation=dupcor$consensus)
dupcor <- duplicateCorrelation(vDup, design.RUV, block=z$samples$ind) # get warning message: Too much damping - convergence tolerance not achievable

# fit linear models
voomDupVfit <- lmFit(vDup, design.RUV, block=z$samples$ind, correlation=dupcor$consensus)
voomDupVfit <- contrasts.fit(voomDupVfit, contrasts=contr.matrix.RUVs)
efit.RUV <- eBayes(voomDupVfit, robust=T)

# for now, we'll continue sticking with a k of 5. Let's see how many DE genes we get setting a lfc of 1 and 0.01 pvalue
dt <- decideTests(efit.RUV, p.value=0.01, lfc=1)
summary(dt)

#       SMBvsMTW SMBvsMPI MTWvsMPI
#Down          8       58       87
#NotSig    12939    12761    12729
#Up           28      156      159


# now let's do the same on a regular limma LM pipeline ----------------------------------------------------------------------------------------------------

# set up design
y$samples$Age[which(is.na(y$samples$Age) == T)]=45
design.lm <- model.matrix(~0 + Island + y$samples$Age + batch + y$samples$RIN + y$sample$CD8T + y$sample$CD4T + y$sample$NK + y$sample$Bcell + y$sample$Mono + y$sample$Gran)
colnames(design.lm)=gsub("Island", "", colnames(design.lm))
#rename columns to exclude spaces and unrecognised characters
colnames(design.lm)[c(3,4,7:13)]=c("Mappi","Age","RIN", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

# set up blocking identifiers
y$samples$ind <- samplenames

# set up contrast matrix using nested design
contr.matrix.limma <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design.lm))
v.lm <- voom(y, design.lm, plot=FALSE)
dupcor <- duplicateCorrelation(v.lm, design.lm, block=y$samples$ind)
vDup.lm <- voom(y, design.lm, plot=TRUE, block=y$samples$ind, correlation=dupcor$consensus)
dupcor <- duplicateCorrelation(vDup.lm, design.lm, block=y$samples$ind) # get warning message: Too much damping - convergence tolerance not achievable

# fit linear models
vfit.lm <- lmFit(vDup.lm, design.lm, block=y$samples$ind, correlation=dupcor$consensus)
vfit.lm <- contrasts.fit(vfit.lm, contrasts=contr.matrix.limma)
efit.lm <- eBayes(vfit.lm, robust=T)

dt.lm <- decideTests(efit.lm, p.value=0.01, lfc=1)
summary(dt.lm)

#       SMBvsMTW SMBvsMPI MTWvsMPI
#Down          6       87       96
#NotSig    12940    12662    12748
#Up           29      226      131


# LM vs RUV similiarities --------------------------------

# Get which genes overlap between RUVs and linear model
pdf(paste0(outputdir,"VennDiagram_LMvsRUVsComparison.pdf"), height=4, width=12)
# sapply(1:3, function(x) vennDiagram(cbind(dt[,x],dt.lm[,x]), circle.col=c("red","blue"), names=c("RUVs", "LinearModel"), main=colnames(dt.lm)[x]))
par(mfrow=  c(1,3),mar = c(0,0,0,0)) 
for (x in c(1:3)){
    vennDiagram(cbind(dt[,x],dt.lm[,x]), circle.col=c("red","blue"), names=c("RUVs", "LinearModel"))
    title(colnames(dt.lm)[x], line=-6)
}
dev.off()

for (i in 1:3){
	commonGenes=dt[,i] & dt.lm[,i]
	commonGenes=which(commonGenes == TRUE)
	commonGenes.results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=names(commonGenes), filters="ensembl_gene_id")
	write.table(commonGenes.results,file=paste0(outputdir,"allCommonGenes_",colnames(dt)[i],".txt"))
}

# test top genes -----------------------------------

# transform ensembl IDs to entrez IDs to be compatible with human c2 dataset (below)
ensembl_genes=rownames(y)
entrez=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "description"),values= ensembl_genes,mart=ensembl.mart.90)
y$entrezID=entrez[match(rownames(y), entrez[,1]), 2]

# gene set testing with Camera
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata")) 
idx <- ids2indices(Hs.c2,id=y$entrezID) 

# analysis of top gene pathways through Camera
pdf(paste0(outputdir,"Top10CameraGeneSets_VennDiagram.pdf"), height=12, width=14)
for (i in 1:3){
    camera.lm.matrix=camera(v.lm,idx,design.lm,contrast=contr.matrix.limma[,i])
    camera.RUV.matrix=camera(v.RUV,idx,design.RUV,contrast=contr.matrix.RUVs[,i])
    LM=rownames(camera.lm.matrix[1:100,])
    RUV=rownames(camera.RUV.matrix[1:100,])
    intersect=length(intersect(LM, RUV))
    fit2 <- euler(c(LM=length(LM)-intersect,RUV=length(RUV)-intersect,"LM&RUV"=intersect))
    assign(colnames(efit.lm)[i],plot(fit2, fills = c("dodgerblue4","darkgoldenrod1"),edges = FALSE,fontsize = 8,quantities = list(fontsize = 10, col="white"), alpha=0.8, main=colnames(efit.lm)[i], cex=1, counts=T))
    write.table(camera.lm.matrix, file=paste0(outputdir,"CAMERA_LM_",colnames(efit.lm)[i],".txt"))
    write.table(camera.lm.matrix, file=paste0(outputdir,"CAMERA_RUV_",colnames(efit.lm)[i],".txt"))
}
grid.arrange(SMBvsMTW, SMBvsMPI, MTWvsMPI, ncol=3, widths=c(1.1,1,1.1))
# you can also plot this with cowplot: plot_grid(SMBvsMTW, SMBvsMPI, MTWvsMPI,rel_heights = c(1/2, 1/4, 1/4), align="h", ncol=3). Source: https://stackoverflow.com/questions/36198451/specify-widths-and-heights-of-plots-with-grid-arrange
dev.off()

# look at the top genes using top table
for (i in 1:3){
    topTable.lm <- topTable(efit.RUV, coef=i, p.value=0.01, n=Inf, sort.by="P")
    topTable.RUV <- topTable(efit.lm, coef=i, p.value=0.01, n=Inf, sort.by="P")
    LM=topTable.lm$SYMBOL[1:100]
    RUV=topTable.RUV$SYMBOL[1:100]
    intersect=length(intersect(LM, RUV))
    fit2 <- euler(c(LM=length(LM)-intersect,RUV=length(RUV)-intersect,"LM&RUV"=intersect))
    assign(colnames(efit.lm)[i],plot(fit2, fills = c("dodgerblue4","darkgoldenrod1"),edges = FALSE,fontsize = 8,quantities = list(fontsize = 10, col="white"), alpha=0.8, main=colnames(efit.lm)[i], cex=1, counts=T))
}
grid.arrange(SMBvsMTW, SMBvsMPI, MTWvsMPI, ncol=3, widths=c(1.1,1,1))

# plot pca and associations of both -------------------------------------------

# rename column names of lcpm to sample names (so that they are shorter and easier to read)
colnames(lcpm)=samplenames

# first set up batch-corrected lm data
design <- model.matrix(~0 + Island)
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[(3)]=c("Mappi")
# here, we're using log2-cpm values of the TMM-normalised data. This seems to be suffiecient (i.e., voom-corrected output is not necessary). See link here: https://support.bioconductor.org/p/76837/.
batch.corrected.lcpm <- removeBatchEffect(lcpm, batch=batch, covariates = cbind(y$samples$Age, y$samples$RIN, y$sample$CD8T, y$sample$CD4T, y$sample$NK, y$sample$Bcell, y$sample$Monoy$sample$Gran),design=design)

# This is our PCA plotting function. We might need it later to explore differnet dimensions of the PCA
plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(dataToPca), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    for (i in 1:9){
        pca_axis1=i
        pca_axis2=i+1
        plot(pca$x[,pca_axis1], pca$x[,pca_axis2], col=speciesCol, pch=namesPch, cex=2, xlab=paste0("PC", pca_axis1, " (", round(pca.var[pca_axis1]*100, digits=2), "% of variance)"), ylab=paste0("PC", pca_axis2, " (", round(pca.var[pca_axis2]*100, digits=2), "% of variance)", sep=""), main=name)
        # points(pca$x[,pca_axis1][which(allreplicated==T)], pca$x[,pca_axis2][which(allreplicated==T)], col="black", pch=8, cex=2)
        # text(pca$x[,pca_axis1][which(allreplicated==T)], pca$x[,pca_axis2][which(allreplicated==T)], labels=samplenames[which(allreplicated==T)], pos=3)
        legend(legend=unique(sampleNames), pch=16, x="bottomright", col=unique(speciesCol), cex=0.6, title=name, border=F, bty="n")
        legend(legend=unique(as.numeric(y$samples$batch)), "topright", pch=unique(as.numeric(y$samples$batch)) + 14, title="Batch", cex=0.6, border=F, bty="n")
    }

    return(pca)
}

# PCA association function
pc.assoc <- function(pca.data){
    all.pcs <- data.frame()
    for (i in 1:ncol(pca.data$x)){
        all.assoc <- vector()
        for (j in 1:ncol(all.covars.df)){
            test.assoc <- anova(lm(pca.data$x[,i] ~ all.covars.df[,j]))[1,5]
            all.assoc <- c(all.assoc, test.assoc)
        }
        single.pc <- c(i, all.assoc)
        all.pcs <- rbind(all.pcs, single.pc)
    }
    names(all.pcs) <- c("PC", colnames(all.covars.df))

    print ("Here are the relationships between PCs and some possible covariates")
    print (all.pcs)
    return (all.pcs)
}

# Prepare covariate matrix
# assign covariate names
# subtract variables we don't need
subtract=c("group", "norm.factors", "samples")
# get index of unwanted variables
subtract=which(colnames(y$samples) %in% subtract)
covariate.names = colnames(y$samples)[-subtract]
for (name in covariate.names){
 assign(name, y$samples[[paste0(name)]])
}

# get rid of covariates we aren't interested in
covariate.names=covariate.names[grep("lib.size|ID|microscopy.pos|PCR.pos|fract.pfpx.reads|replicate|ind",covariate.names, invert=T)]
# Prepare covariate matrix
all.covars.df <- y$samples[,covariate.names]

# Let's just look at the forst pca for now and make a function to do this (without any labels to make the visualisation of batch clear)
one.dimension=function(data) {
    pca <- prcomp(t(data), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    plot(pca$x[,1], pca$x[,2], col=batch.col[as.numeric(batch)], pch=as.numeric(y$samples$batch) + 14, cex=2, xlab=paste0("PC", 1, " (", round(pca.var[1]*100, digits=2), "% of variance)"), ylab=paste0("PC", 2, " (", round(pca.var[2]*100, digits=2), "% of variance)", sep=""))
    legend(legend=unique(batch), pch=unique(as.numeric(y$samples$batch)) + 14, x="bottomright", col=unique(batch.col[as.numeric(batch)]), cex=0.6, title="batch", border=F, bty="n")
}

# Compare the first PCA of the LM-corrected method to the RUVs corrected method
pdf(paste0(outputdir,"PCA_RUVvsLM_FirstDimension.pdf"), height=8, width=15)
par(mfrow=c(1,2))
one.dimension(batch.corrected.lcpm)
title("Limma-LM")
# here, we're looking at the RUVs-corrected, log2-transformed output
one.dimension(cpm(normCounts(set1), log=T))
title("RUVs")
dev.off()

# We also want to look at a heatmap of the significant covariates compared to each dimension of the PCA

# first for limma-corrected data
name="batch"
pcaresults.limma <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=batch.col[as.numeric(batch)],namesPch=as.numeric(y$samples$batch) + 14,sampleNames=batch)
all.pcs.limma <- pc.assoc(pcaresults.limma)
all.pcs.limma$Variance <- pcaresults.limma$sdev^2/sum(pcaresults.limma$sdev^2)
# we can save the table to look at the association in all of the dimensions
write.table(all.pcs.limma, file=paste0(outputdir,"PC_associations_limmaCorrectedData.txt"))

# now for the RUVs-corrected data
pcaresults.RUVs <- plot.pca(dataToPca=cpm(normCounts(set1), log=T), speciesCol=batch.col[as.numeric(batch)],namesPch=as.numeric(y$samples$batch) + 14,sampleNames=batch)
all.pcs.RUVs <- pc.assoc(pcaresults.RUVs)
all.pcs.RUVs$Variance <- pcaresults.RUVs$sdev^2/sum(pcaresults.RUVs$sdev^2)
# save the associations for the RUVs data
write.table(all.pcs.RUVs, file=paste0(outputdir,"PC_associations_RUVsCorrectedData.txt"))

# set up heatmap information
# first limma
limma.pc=all.pcs.limma[1:5,covariate.names]
limma.hm=melt(limma.pc)
limma.hm$Dimension=rep(1:5,ncol(limma.pc))
colnames(limma.hm)[c(1,2)]=c("Covariate","Association")

# RUVs
RUV.pc=all.pcs.RUVs[1:5,covariate.names]
RUV.hm=melt(RUV.pc)
RUV.hm$Dimension=rep(1:5,ncol(RUV.pc))
colnames(RUV.hm)[c(1,2)]=c("Covariate","Association")

# set up colour palette
hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')
# Save output of both plots
limma=ggplot(data = limma.hm, aes(x = Covariate, y = Dimension)) + geom_tile(aes(fill = Association))  + scale_fill_gradientn(colours = hm.palette(100)) + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Limma") 
RUV=ggplot(data = RUV.hm, aes(x = Covariate, y = Dimension)) + geom_tile(aes(fill = Association))  + scale_fill_gradientn(colours = hm.palette(100)) + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("RUVs") 
# use ggarange to plot both heatmpas
pdf(paste0(outputdir,"Heatmap_SigCovar_RUVvsLM.pdf"), height=7, width=15)
ggarrange(limma,RUV)
dev.off()

# Compare p-value distributions of both -----------------------------------------------

pdf(paste0(outputdir,"VolcanoPlot_RUVvsLM.pdf"), height=10, width=15)
layout(matrix(c(1,2,3,4,5,6), nrow=3,byrow = F))
for (i in 1:ncol(efit.lm)){
    plot(efit.lm$coef[,i], -log10(as.matrix(efit.lm$p.value)[,i]), pch=20, main=paste("LM",colnames(efit)[i],sep="\n"), xlab="log2FoldChange", ylab="-log10(pvalue)")
    points(efit.lm$coef[,i][which(names(efit.lm$coef[,i]) %in% hkControls)], -log10(as.matrix(efit.lm$p.value)[,i][which(names(efit.lm$coef[,i]) %in% hkControls)]) , pch=20, col=4, xlab="log2FoldChange", ylab="-log10(pvalue)")
    legend("topleft", "genes", "hk genes",fill=4)
    abline(v=c(-1,1))
}
for (i in 1:ncol(efit.RUV)){
    plot(efit.RUV$coef[,i], -log10(as.matrix(efit.RUV$p.value)[,i]), pch=20, main=paste("RUV",colnames(efit)[i],sep="\n"), xlab="log2FoldChange", ylab="-log10(pvalue)")
    points(efit.RUV$coef[,i][which(names(efit.RUV$coef[,i]) %in% hkControls)], -log10(as.matrix(efit.RUV$p.value)[,i][which(names(efit.RUV$coef[,i]) %in% hkControls)]) , pch=20, col=2, xlab="log2FoldChange", ylab="-log10(pvalue)")
    legend("topleft", "genes", "hk genes",fill=2)
    abline(v=c(-1,1))
}
dev.off()

# See which genes are not in common and why ----------------------------------------------------------------------

# set up not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

# see which genes aren't common and why this might be
pdf(paste0(outputdir,"genesNotInCommon_RUVsVsLM.pdf"), height=10, width=15)
for (i in 1:ncol(efit.lm)){
    RUVs=dt[,i][which(abs(dt[,i])==T)]
    lm=dt.lm[,i][which(abs(dt.lm[,i])==T)]
    notInLM=names(RUVs[which(names(RUVs) %!in% names(lm))])
    notInRUV=names(lm[which(names(lm) %!in% names(RUVs))])
    topTable.lm <- topTable(efit.lm, coef=i, n=Inf)
    topTable.RUV <- topTable(efit, coef=i, n=Inf)

    plot(topTable.lm[notInLM,]$adj.P.Val, abs(topTable.lm[notInLM,]$logFC), pch=16, main=paste("Top DE Genes Not in Common", colnames(efit)[i], sep="\n"), xlab="adjusted p.value", ylab="logFC")
    points(topTable.RUV[notInRUV,]$adj.P.Val, abs(topTable.RUV[notInRUV,]$logFC), col="red", pch=16)
    text(topTable.lm[notInLM,]$adj.P.Val, abs(topTable.lm[notInLM,]$logFC), labels=topTable.lm[notInLM,]$SYMBOL, pos=1)
    text(topTable.RUV[notInRUV,]$adj.P.Val, abs(topTable.RUV[notInRUV,]$logFC), labels=topTable.RUV[notInRUV,]$SYMBOL, col="red", pos=2)
    legend("topright", c("genes not in lm", "genes not in RUVs"), pch=16, col=c("black", "red"))
    abline(v=0.01, lty=2)
    abline(h=1, lty=2)
}
dev.off()









