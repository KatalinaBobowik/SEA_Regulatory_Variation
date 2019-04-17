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
library(RUVSeq)


# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/"
housekeepingdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/BatchEffects/"
# Set output directory and create it if it does not exist:
outputdir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/RUVs_Setup/"

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir)
}

# set up colour palette. The "wes" palette will be used for island and other statistical information, whereas NineteenEightyR will be used for batch information
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
# set up colour palette for batch
batch.col=electronic_night(n=3)
lightcolors=c("thistle", "darkblue")
dev.off()

# Begin Analysis --------------------------------------------------------------------------------------------------

# load the filtered, normalised y DGE list object
load(paste0(inputdir, "dataPreprocessing/indoRNA.read_counts.TMM.filtered.Rda"))

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

# RUVSeq works with a genes-by-samples numeric matrix or a SeqExpressionSet object containing read counts. Let's set up a SeqExpressionSet with our counts matrix
set <- newSeqExpressionSet(as.matrix(y$counts), phenoData = data.frame(Island, row.names=colnames(y)))

# Set up PCA function. This is the hacked code of the plotPCA package in RUVseq. We'll need this to make various PCA plots
hackedPCA=function(object,pc1,pc2,pch,color){
    logTransform <- apply(log(object+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
    s <- svd(logTransform)
    percent <- s$d^2/sum(s$d^2)*100
    labs <- sapply(seq_along(percent), function(i) {
            paste("PC ", i, " (", round(percent[i], 2), "%)", sep="")
    })
    plot(s$u[,pc1], s$u[,pc2], xlab=labs[pc1], ylab=labs[pc2], pch=pch, col=color)
    points(s$u[,pc1][which(allreplicated==T)], s$u[,pc2][which(allreplicated==T)], col="black", pch=8, cex=1)
    text(s$u[,pc1][which(allreplicated==T)], s$u[,pc2][which(allreplicated==T)], labels=samplenames[which(allreplicated==T)], pos=1, cex=0.8)
}

# Normalisation can be performed using "median","upper", or "full", however when passing to edgeR's normalisation method (below), the only options are "TMM","RLE", and "upperquartile". In order to keep consistency, we'll go ahead and choose the "upper" method, since it's in both.
# check normalisation before and after performing upper qiuartile normalisation
pdf(paste0(outputdir,"Normalisation_beforeandAfterUQnorm.pdf"), height=15, width=15)
par(mfrow=c(2,2))
plotRLE(set, outline=FALSE, ylim=c(-2,2), col=batch.col[batch], xaxt='n', main="No Normalisation")
legend(legend=unique(as.numeric(y$samples$batch)), "topright", pch=15, title="Batch", cex=0.6, border=F, bty="n", col=unique(batch.col))
# get PCA of batch and check how closely replictes sit to one another
hackedPCA(object=counts(set), pc1=1, pc2=2, pch=as.numeric(allreplicated)+15, col=batch.col[batch])
title("No Normalisation")
legend(legend=unique(as.numeric(batch)), "topright", pch=15, title="Batch", cex=0.6, border=F, bty="n", col=unique(batch.col))

# normalise with UQ normalisation
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, col=batch.col[batch], xaxt='n', main="Upper Quartile Normalisation", ylim=c(-2, 2))
legend(legend=unique(as.numeric(y$samples$batch)), "topright", pch=15, title="Batch", cex=0.6, border=F, bty="n", col=unique(batch.col[as.numeric(batch)]))
hackedPCA(object=normCounts(set), pc1=1, pc2=2, pch=as.numeric(allreplicated)+15, col=batch.col[batch])
title("Upper Quartile Normalisation")
legend(legend=unique(as.numeric(batch)), "topright", pch=15, title="Batch", cex=0.6, border=F, bty="n", col=unique(batch.col))
dev.off()


# First, construct a matrix specifying the replicates. 
# create matrix filled with -1s 
replicates=matrix(-1, nrow=length(allreps), ncol=3)
rownames(replicates)=unique(samplenames[samplenames %in% allreps])
for (i in 1:nrow(replicates)){
    replicates[i,1:length(grep(rownames(replicates)[i], samplenames))] = grep(rownames(replicates)[i], samplenames)
}

# set up genes used for RUVs
genes <- rownames(y)

# explore which value of k is the best by looking at variation in RLE and PCA plots
pdf(paste0(outputdir,"RUVsNormalisation_choosingK.pdf"), height=15, width=12)
par(mar=c(4.1,4.1,3.1,1.1), mfrow=c(6,2))
for (i in c(1:6)){
    set1 <- RUVs(set, genes, k=i, replicates)
    plotRLE(set1, main=paste0("k = ",i), col=as.numeric(batch), outline=FALSE, xaxt="n")
    # plot PCA to see how closely replicates sit
    hackedPCA(normCounts(set1),pc1=1,pc2=2,pch=as.numeric(allreplicated)+15,color=lightcolors[allreplicated]) 
    title(paste0("k = ",i))
}
dev.off()


# Set up list of housekeeping genes as controls (from Eisenberg and Levanon, 2003) for volcano plots
housekeeping=read.table(paste0(housekeepingdir,"Housekeeping_ControlGenes.txt"), as.is=T, header=F)
# if this is broken, use host = "uswest.ensembl.org"
ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'www.ensembl.org', ensemblRedirect = FALSE)
biomart.results.table <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=housekeeping, filters="hgnc_symbol")
hkGenes=as.vector(biomart.results.table[,1])
hkControls=hkGenes[which(hkGenes %in% rownames(y$counts))]

# look at the p-value distributions and choose the best k
deGenes=vector()
k=vector()
allGenes=list()
for (j in 1:3){
    counter=0
    pdf(paste0(outputdir,"pvalueDist_choosingK_RUVs_",j,".pdf"))
    for (i in c(1:6)){
        counter=counter+1
        set1 <- RUVs(set, genes, k=i, replicates)
        design <- cbind(model.matrix(~0 + Island), pData(set1)[2:(i+1)])
        colnames(design)=gsub("Island", "", colnames(design))
        colnames(design)[3]="Mappi"
        z <- DGEList(counts=counts(set1), group=Island)
        # here, we use uppserquarile normalisation since "upper" is used for between lane normalisation of our 'set' object
        z <- calcNormFactors(z, method="upperquartile")
        contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
        v <- voom(z, design, plot=FALSE, normalize.method = "quantile")
        vfit <- lmFit(v, design)
        vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
        efit <- eBayes(vfit)
        topTable <- topTable(efit, coef=j, p.value=0.01, lfc=1, n=Inf)
        deGenes[counter]=nrow(topTable)
        k[counter]=i
        # get p-value distribution from raw pvalues
        hist(efit$p.value[,j], main=paste(colnames(efit)[j],i,sep="\n"), ylim=c(0,max(table(round(efit$p.value[,j], 1)))+1000), xlab="p-value")
        # get volcano plots
        plot(efit$coef[,j], -log10(as.matrix(efit$p.value)[,j]), pch=20, main=paste(colnames(efit)[j],i,sep="\n"), xlab="log2FoldChange", ylab="-log10(pvalue)")
        points(efit$coef[,j][which(names(efit$coef[,j]) %in% hkControls)], -log10(as.matrix(efit$p.value)[,j][which(names(efit$coef[,j]) %in% hkControls)]) , pch=20, col=4, xlab="log2FoldChange", ylab="-log10(pvalue)")
        legend("topleft", c("genes", "hk genes"),fill=c("black",4))
        abline(v=c(-1,1), lty=2)
    }
    dev.off()
    # pdf(paste("numberDeGenes_choosingK_RUVs_",colnames(efit)[j],".pdf"))
    # plot(k, deGenes, main=colnames(efit)[j])
    allGenes[[colnames(efit)[j]]]=deGenes
    # dev.off()
}

# plot the number of DE genes for all population comparisons and all varying degrees of k
pdf(paste0(outputdir,"numberDeGenes_choosingK_RUVs_.pdf"), height=5, width=15)
par(mfrow=c(1,3))
for (i in c(1:3)){
    plot(k, allGenes[[i]], main=names(allGenes)[i], ylab="n DE Genes")
}
dev.off()

# applying RUVs ----------------------------------------------------------------------------

# we'll set k to 5 for now
set1 <- RUVs(set, genes, k=5, replicates)

# write table of values of k, the unwnated factors of variation
write.table(pData(set1), file=paste0(outputdir,"RUVs_unwantedFactorsOfVariation_Table.txt"), quote=F, sep="\t")
# save the set1 object
save(set1, file = paste0(outputdir, "set1_RUVsCorrectedObject.Rda"))

# now plot all of the covariates and see how the batch correction worked

# assign covariate names
# subtract variables we don't need
subtract=c("group", "norm.factors", "samples")
# get index of unwanted variables
subtract=which(colnames(y$samples) %in% subtract)
covariate.names = colnames(y$samples)[-subtract]
for (name in covariate.names){
 assign(name, y$samples[[paste0(name)]])
}

# Age, RIN, and library size need to be broken up into chunks for easier visualisation of trends (for instance in Age, we want to see high age vs low age rather than the effect of every single age variable)
for (name in c("Age","RIN","lib.size")){
  assign(name, cut(as.numeric(as.character(y$samples[[paste0(name)]])), breaks=5))
}

# assign names to covariate names so you can grab individual elements by name
names(covariate.names)=covariate.names

# assign factor variables
factorVariables=c(colnames(Filter(is.factor,y$samples))[which(colnames(Filter(is.factor,y$samples)) %in% covariate.names)], "Age", "lib.size", "RIN")
numericVariables=colnames(Filter(is.numeric,y$samples))[which(colnames(Filter(is.numeric,y$samples)) %in% covariate.names)] %>% subset(., !(. %in% factorVariables))

# get rid of covariates we aren't interested in
covariate.names=covariate.names[grep("lib.size|ID|microscopy.pos|PCR.pos|fract.pfpx.reads|replicate",covariate.names, invert=T)]
# Prepare covariate matrix
all.covars.df <- y$samples[,covariate.names]

# PCA plotting function
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

# Plot PCA
for (name in factorVariables){
  if (nlevels(get(name)) < 26){
    pdf(paste0(outputdir,"pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=cpm(normCounts(set1), log=T), speciesCol=as.numeric(get(name)),namesPch=as.numeric(y$samples$batch) + 14,sampleNames=get(name))
    dev.off()
  } else {
    pdf(paste0(outputdir,"pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=cpm(normCounts(set1), log=T), speciesCol=as.numeric(get(name)),namesPch=20,sampleNames=get(name))
    dev.off()
  }
}

# plot batch
pdf(paste0(outputdir,"pcaresults_batch.pdf"))
name="batch"
pcaresults <- plot.pca(dataToPca=cpm(normCounts(set1), log=T), speciesCol=batch.col[as.numeric(batch)],namesPch=as.numeric(y$samples$batch) + 14,sampleNames=batch)
dev.off()
  
# plot blood
for (name in numericVariables){
    initial <- cut(get(name), breaks = seq(min(get(name), na.rm=T), max(get(name), na.rm=T), len = 80),include.lowest = TRUE)
    bloodCol <- colorRampPalette(c("blue", "red"))(79)[initial]
    pdf(paste0(outputdir,"pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=cpm(normCounts(set1), log=T), speciesCol=bloodCol,namesPch=as.numeric(y$samples$batch) + 14,sampleNames=get(name))
    dev.off()
}


# Get PCA associations
all.pcs <- pc.assoc(pcaresults)
all.pcs$Variance <- pcaresults$sdev^2/sum(pcaresults$sdev^2)

# plot pca covariates association matrix to illustrate any potential confounding and evidence for batches
pdf(paste0(outputdir,"significantCovariates_AnovaHeatmap.pdf"))
pheatmap(all.pcs[1:5,covariate.names], cluster_col=F, col= colorRampPalette(brewer.pal(11, "RdYlBu"))(100), cluster_rows=F, main="Significant Covariates \n Anova")
dev.off()

# Write out the covariates:
write.table(all.pcs, file=paste0(outputdir,"pca_covariates_blood.txt"), col.names=T, row.names=F, quote=F, sep="\t")