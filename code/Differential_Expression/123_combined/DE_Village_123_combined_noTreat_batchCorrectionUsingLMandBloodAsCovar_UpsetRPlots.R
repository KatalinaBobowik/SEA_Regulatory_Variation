# script created by KSB, 08.08.18
# Perform DE analysing relationship between villages

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
library(EnsDb.Hsapiens.v86)
library(wesanderson)

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/dataPreprocessing/"
housekeepingdir="/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/BatchEffects/"

# Set output directory and create it if it does not exist:
outputdir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Village/"

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir)
}

# Load colour schemes:
# Load colour schemes:
mappi <- wes_palette("Zissou1", 20, type = "continuous")[20]
mentawai <- wes_palette("Zissou1", 20, type = "continuous")[1]
sumba <- wes_palette("Zissou1", 20, type = "continuous")[11]

smb_mtw <- wes_palette("Darjeeling1", 9, type = "continuous")[3]
smb_mpi <- wes_palette("Darjeeling1", 9, type = "continuous")[7]
mtw_mpi <- "darkorchid4"

# set up colour palette for batch
batch.col=electronic_night(n=3)
village.col=c("#EBCC2A","chocolate","chocolate","#3B9AB2","#F21A00","chocolate","chocolate","chocolate","#78B7C5","orange","chocolate")

dev.off()


# Load log CPM matrix and y object:
# lcpm
load(paste0(inputdir, "indoRNA.logCPM.TMM.filtered.Rda"))
# y DGE list object
load(paste0(inputdir, "indoRNA.read_counts.TMM.filtered.Rda"))

# Removing heteroscedascity with voom and fitting linear models -----------------------------------------------------------------------

# First, remove samples that have less than ten individuals per village
table(y$samples$Sampling.Site)
# Sampling.Site
#   Anakalung    Bilarenge    Hupu Mada      Madobag        Mappi  Padira Tana 
#          20            1            5           17           21            3 
# Patiala Bawa        Rindi     Taileleu        Wunga   Wura Homba 
#           1            5           32           17            1 

# remove Bilarenge, Hupu Mada, Padira Tana, Patiala Bawa, Rindi, and Wura Homba
y=y[,-grep("Bilarenge|Hupu Mada|Padira Tana|Patiala Bawa|Rindi|Wura Homba", y$samples$Sampling.Site)]
# drop unused levels
y$samples=droplevels(y$samples)

# set up designs
design <- model.matrix(~0 + y$samples$Sampling.Site + y$samples$Age + y$samples$batch + y$samples$RIN + y$samples$CD8T + y$samples$CD4T + y$samples$NK + y$samples$Bcell + y$samples$Mono + y$samples$Gran)
# Get rid of 'Sampling Site' and 'y$samples' from column names
colnames(design)=gsub("Sampling.Site", "", colnames(design))
colnames(design)=gsub("[\\y$]", "", colnames(design))
colnames(design)=gsub("samples", "", colnames(design))

# set up contrast matrix
#contr.matrix <- makeContrasts(ANKvsMDB=Anakalung-Madobag, ANKvsMPI=Anakalung-Mappi, ANKvsTLL=Anakalung-Taileleu, ANKvsWNG=Anakalung-Wunga, MDBvsMPI=Madobag-Mappi, MDBvsTLL=Madobag-Taileleu, MDBvsWNG=Madobag-Wunga, MPIvsTLL=Mappi-Taileleu, MPIvsWNG=Mappi-Wunga, TLLvsWNG=Taileleu-Wunga, levels=colnames(design))
contr.matrix <- makeContrasts(ANKvsMDB=Anakalung-Madobag, ANKvsMPI=Anakalung-Mappi, ANKvsTLL=Anakalung-Taileleu, ANKvsWNG=Anakalung-Wunga, MDBvsMPI=Madobag-Mappi, MDBvsTLL=Madobag-Taileleu, WNGvsMDB=Wunga-Madobag, TLLvsMPI=Taileleu-Mappi, WNGvsMPI=Wunga-Mappi, WNGvsTLL=Wunga-Taileleu, levels=colnames(design)) # Contrasts are ordered in the same order as the island ones, in case we want to look at directional effects

# now go ahead with voom normalisation
# Using duplicate correlation and blocking -----------------------------------------------------

# First, we need to perform voom normalisation
v <- voom(y, design, plot=F)

# create a new variable for blocking using sample IDs
# define sample names
samplenames <- as.character(y$samples$samples)
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
samplenames <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)

y$samples$ind <- samplenames

# Estimate the correlation between the replicates.
# Information is borrowed by constraining the within-block corre-lations to be equal between genes and by using empirical Bayes methods to moderate the standarddeviations between genes 
dupcor <- duplicateCorrelation(v, design, block=y$samples$ind)

# run voom a second time with the blocking variable and estimated correlation
# The  vector y$samples$ind indicates the  two  blocks  corresponding  to  biological  replicates
pdf(paste0(outputdir,"Limma_voomDuplicateCorrelation_TMMNormalisation.pdf"), height=8, width=12)
par(mfrow=c(1,2))
vDup <- voom(y, design, plot=TRUE, block=y$samples$ind, correlation=dupcor$consensus)
dupcor <- duplicateCorrelation(vDup, design, block=y$samples$ind) # get warning message: Too much damping - convergence tolerance not achievable

# With duplicate correction and blocking:
# the inter-subject correlation is input into the linear model fit
voomDupVfit <- lmFit(vDup, design, block=y$samples$ind, correlation=dupcor$consensus)
voomDupVfit <- contrasts.fit(voomDupVfit, contrasts=contr.matrix)
voomDupEfit <- eBayes(voomDupVfit, robust=T)

plotSA(voomDupEfit, main="Mean-variance trend elimination with duplicate correction")
dev.off()

# save voom and efit object
save(vDup, file = paste0(outputdir, "vDup_Village.Rda"))
save(voomDupEfit, file = paste0(outputdir, "voomDupEfit_Village.Rda"))

# get top genes using toptable
allDEresults <- list()

for(i in 1:ncol(voomDupEfit)){
    allDEresults[[i]] <- topTable(voomDupEfit, coef=i, n=Inf, sort.by="p")
}

# get number of DE genes at differnet thresholds
# noLFC
summary(decideTests(voomDupEfit, method="separate", adjust.method = "BH", p.value = 0.01))
#       ANKvsMDB ANKvsMPI ANKvsTLL ANKvsWNG MDBvsMPI MDBvsTLL WNGvsMDB TLLvsMPI
#Down         59     1167      475        1      391      111      567     2112
#NotSig    12844    10821    12066    12973    12027    12504    12005     9005
#Up           72      987      434        1      557      360      403     1858
#       WNGvsMPI WNGvsTLL
#Down       2126      690
#NotSig     8875    11756
#Up         1974      529

# LFC of 0.5
summary(decideTests(voomDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
#       ANKvsMDB ANKvsMPI ANKvsTLL ANKvsWNG MDBvsMPI MDBvsTLL WNGvsMDB TLLvsMPI
#Down         19      323       78        0      160        8      149      719
#NotSig    12903    12108    12687    12974    12482    12922    12657    11815
#Up           53      544      210        1      333       45      169      441
#       WNGvsMPI WNGvsTLL
#Down        725      139
#NotSig    11344    12607
#Up          906      229

# LFC of 1
summary(decideTests(voomDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1))
#       ANKvsMDB ANKvsMPI ANKvsTLL ANKvsWNG MDBvsMPI MDBvsTLL WNGvsMDB TLLvsMPI
#Down          5       50        9        0       52        2       17      121
#NotSig    12948    12717    12912    12974    12798    12969    12908    12723
#Up           22      208       54        1      125        4       50      131
#       WNGvsMPI WNGvsTLL
#Down        129        8
#NotSig    12543    12914
#Up          303       53

for (i in 1:ncol(efit)){
    write.table(allDEresults[[i]], file=paste0(outputdir,"topTable.voomNoNorm.tmm.filtered.dup_corrected.village.", colnames(contr.matrix)[i], ".txt"))
}

# Without duplicate correlation -------------------------------------------------

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust=T)

pdf(file=paste0(edaoutput, "voomNoNorm.tmm.filtered.indoRNA.no_dup_correction.mean-variance-trend.village.pdf"))
    plotSA(efit, main="Mean-variance trend elimination without duplicate correction")
dev.off()

# get top genes using toptable
allDEresultsNoDup <- list()
for(i in 1:ncol(efit)){
    allDEresultsNoDup[[i]] <- topTable(efit, coef=i, n=Inf, sort.by="p")
}

# no LFC threshold
summary(decideTests(efit, method="separate", adjust.method = "BH", p.value = 0.01))
#       ANKvsMDB ANKvsMPI ANKvsTLL ANKvsWNG MDBvsMPI MDBvsTLL WNGvsMDB TLLvsMPI
#Down         95     1568      557        1      503      118      675     2186
#NotSig    12785    10184    11936    12973    11857    12514    11792     8939
#Up           95     1223      482        1      615      343      508     1850
#       WNGvsMPI WNGvsTLL
#Down       2459      867
#NotSig     8388    11403
#Up         2128      705

# LFC of 0.05
summary(decideTests(efit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5))
#Down         23      435       95        1      189        9      177      728
#NotSig    12887    11911    12644    12973    12434    12923    12587    11790
#Up           65      629      236        1      352       43      211      457
#       WNGvsMPI WNGvsTLL
#Down        816      176
#NotSig    11194    12493
#Up          965      306

# LFC of 1
summary(decideTests(efit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1))
#       ANKvsMDB ANKvsMPI ANKvsTLL ANKvsWNG MDBvsMPI MDBvsTLL WNGvsMDB TLLvsMPI
#Down          3       61        9        1       61        2       21      121
#NotSig    12948    12674    12903    12973    12782    12969    12884    12727
#Up           24      240       63        1      132        4       70      127
#       WNGvsMPI WNGvsTLL
#Down        146       13
#NotSig    12509    12894
#Up          320       68

for (i in 1:ncol(efit)){
    write.table(allDEresultsNoDup[[i]], file=paste0(outputdir,"topTable.voomNoNorm.tmm.filtered.not_dup_corrected.village.", colnames(contr.matrix)[i], ".txt"))
}

for (i in 1:ncol(efit)){
    withWithout <- join(allDEresults[[i]], allDEresultsNoDup[[i]], by="ENSEMBL")
    print(paste0("Spearman correlation between methods in ", colnames(contr.matrix)[i], ":"))
    print(cor(withWithout[,6], withWithout[,12], method="spearman"))
}

#[1] "Spearman correlation between methods in ANKvsMDB:"
#[1] 0.950506
#[1] "Spearman correlation between methods in ANKvsMPI:"
#[1] 0.9495114
#[1] "Spearman correlation between methods in ANKvsTLL:"
#[1] 0.942682
#[1] "Spearman correlation between methods in ANKvsWNG:"
#[1] 0.9459589
#[1] "Spearman correlation between methods in MDBvsMPI:"
#[1] 0.958493
#[1] "Spearman correlation between methods in MDBvsTLL:"
#[1] 0.9442371
#[1] "Spearman correlation between methods in WNGvsMDB:"
#[1] 0.9507572
#[1] "Spearman correlation between methods in TLLvsMPI:"
#[1] 0.9445926
#[1] "Spearman correlation between methods in WNGvsMPI:"
#[1] 0.9471946
#[1] "Spearman correlation between methods in WNGvsTLL:"
#[1] 0.9451869


# QC after fitting linear models --------------------------------------------------------------------------------------

# check to see p-value distribution is normal
pdf("PvalueDist_NotAdjusted.pdf", height=15, width=10)
par(mfrow=c(5,2))
for (i in 1:ncol(efit)){
    hist(efit$p.value[,i], main=colnames(efit)[i], ylim=c(0,max(table(round(efit$p.value[,i], 1)))+1000), xlab="p-value")
}
dev.off()

# check p-value distribution for adjusted p-values
pdf("PvalueDist_Adjusted.pdf", height=15, width=10)
par(mfrow=c(5,2))
for (i in 1:ncol(efit)){
    topTable <- topTable(efit, coef=i, n=Inf)
    hist(topTable$adj.P.Val, main=colnames(efit)[i], xlab="p-value")
}
dev.off()

# Verify that control housekeeping genes are not significantly DE. Set up list of housekeeping genes as controls (from Eisenberg and Levanon, 2003)
housekeeping=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/BatchEffects/Housekeeping_ControlGenes.txt", as.is=T, header=F)
# if this is broken, use host = "uswest.ensembl.org"
ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'www.ensembl.org', ensemblRedirect = FALSE)
biomart.results.table <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=housekeeping, filters="hgnc_symbol")
hkGenes=as.vector(biomart.results.table[,1])
hkControls=hkGenes[which(hkGenes %in% rownames(y$counts))]

# Volcano plot with points of housekeeping genes
pdf("VolcanoPlots.pdf", height=15, width=10)
par(mfrow=c(5,2))
for (i in 1:ncol(efit)){
    plot(efit$coef[,i], -log10(as.matrix(efit$p.value)[,i]), pch=20, main=colnames(efit)[i], xlab="log2FoldChange", ylab="-log10(pvalue)")
    points(efit$coef[,i][which(names(efit$coef[,i]) %in% hkControls)], -log10(as.matrix(efit$p.value)[,i][which(names(efit$coef[,i]) %in% hkControls)]) , pch=20, col=4, xlab="log2FoldChange", ylab="-log10(pvalue)")
    legend("topleft", "genes", "hk genes",fill=4)
    abline(v=c(-1,1))
}
dev.off()

# PCA visualisation after correction and association with covariates ------------------------------------------------------------

# let's also visualise how our PCAs look after limma correction by using removeBatcheffect. Help on design of removeBatcheffects was given by the lovely John Blischak.
design <- model.matrix(~0 + y$samples$Sampling.Site)
colnames(design)=gsub("Sampling.Site", "", colnames(design))
colnames(design)=gsub("[\\y$]", "", colnames(design))
colnames(design)=gsub("samples", "", colnames(design))

# Remove unused columns from lcpm matrix
lcpm=lcpm[,which(colnames(lcpm) %in% colnames(y))]
batch.corrected.lcpm <- removeBatchEffect(lcpm, batch=y$samples$batch, covariates = cbind(y$samples$Age, y$samples$RIN, y$samples$CD8T, y$samples$CD4T, y$samples$NK, y$samples$Bcell, y$samples$Mono, y$samples$Gran), design=design)

# rename column names of lcpm to sample names (so that they are shorter and easier to read)
colnames(lcpm)=samplenames

# get which samples are replicates
load(paste0(inputdir, "covariates.Rda"))
allreps=covariates[,1][which(covariates$replicate)]
allreps=unique(allreps)
allreplicated=as.factor(samplenames %in% allreps)

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

# Prepare covariate matrix
all.covars.df <- y$samples[,covariate.names]

# Plot PCA
for (name in factorVariables){
  if (nlevels(get(name)) < 26){
    pdf(paste0(outputdir,"pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=as.numeric(get(name)),namesPch=as.numeric(y$samples$batch) + 14,sampleNames=get(name))
    dev.off()
  } else {
    pdf(paste0(outputdir,"pcaresults_",name,"_RemoveBatchEffect.pdf"))
    pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=as.numeric(get(name)),namesPch=20,sampleNames=get(name))
    dev.off()
  }
}

# plot batch
pdf(paste0(outputdir,"pcaresults_batch.pdf"))
name="batch"
pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=batch.col[as.numeric(batch)],namesPch=as.numeric(y$samples$batch) + 14,sampleNames=batch)
dev.off()
  
# plot numeric variables
for (name in numericVariables){
    initial = .bincode(get(name), breaks=seq(min(get(name), na.rm=T), max(get(name), na.rm=T), len = 80),include.lowest = TRUE)
    bloodCol <- colorRampPalette(c("blue", "red"))(79)[initial]
    pdf(paste0(outputdir,"pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=bloodCol,namesPch=as.numeric(y$samples$batch) + 14,sampleNames=get(name))
    legend(legend=c("High","Low"), pch=16, x="bottomright", col=c(bloodCol[which.max(get(name))], bloodCol[which.min(get(name))]), cex=0.6, title=name, border=F, bty="n")
    legend(legend=unique(as.numeric(y$samples$batch)), "topright", pch=unique(as.numeric(y$samples$batch)) + 14, title="Batch", cex=0.6, border=F, bty="n")
    dev.off()
}

# Get PCA associations
all.pcs <- pc.assoc(pcaresults)
all.pcs$Variance <- pcaresults$sdev^2/sum(pcaresults$sdev^2)
write.table(all.pcs, file=paste0(outputdir,"pca_covariates_blood_RNASeqDeconCell.txt"), col.names=T, row.names=F, quote=F, sep="\t")

# plot pca covariates association matrix to illustrate any remaining confounding/batch
pdf(paste0(outputdir,"significantCovariates_AnovaHeatmap.pdf"))
pheatmap(log(all.pcs[1:5,covariate.names]), cluster_col=F, col= colorRampPalette(brewer.pal(11, "RdYlBu"))(100), cluster_rows=F, main="Significant Covariates \n Anova")
dev.off()

# Summary and visualisation of gene trends ---------------------------------------------------------------------------

dt <- decideTests(efit, p.value=0.01, lfc=1)

# plot log2 fold change between islands
pdf(paste0(outputdir,"log2FC_IslandComparisons_pval01_dupCor.pdf"))
# note 'p.value' is the cutoff value for adjusted p-values
reordered.efit=voomDupEfit[,c("ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsMPI", "WNGvsMPI", "MDBvsMPI", "TLLvsMPI", "ANKvsWNG", "MDBvsTLL")]
col=c(rep(smb_mtw,4), rep(smb_mpi, 2), rep(mtw_mpi, 2), sumba, mentawai)
topTable <- topTable(reordered.efit, coef=1, n=Inf, p.value=0.01)
plot(density(topTable$logFC), col=col[1], xlim=c(-3,3), ylim=c(0,1.5), main="LogFC Density", xlab="LogFC", ylab="Density")
abline(v=c(-1,-0.5,0.5,1), lty=3)
counter=0
for (i in 2:ncol(reordered.efit)){
    counter=counter+1
    topTable <- topTable(reordered.efit, coef=i, n=Inf, p.value=0.01)
    lines(density(topTable$logFC), col=col[i], xlim=c(-3,3), ylim=c(0,1.5))
}
legend(x="topright", bty="n", col=col, legend=colnames(voomDupEfit), lty=1, lwd=2)
dev.off()

# graphical representation of DE results through MD plot
pdf("MD_Plots_pval01_lfc1.pdf", height=15, width=10)
par(mfrow=c(5,2))
for (i in 1:ncol(voomDupEfit)){
    o <- which(names(voomDupEfit$Amean) %in% names(which(abs(dt[,i])==1)))
    x <- voomDupEfit$Amean
    z <- voomDupEfit$coefficients[,i]
    t=which(names(voomDupEfit$coefficients[,i]) %in% names(which(abs(dt[,i])==1)))
    G <- voomDupEfit$genes[which(abs(dt[,i])==1),]$SYMBOL
    plotMD(voomDupEfit, column=i, status=dt[,i], main=colnames(voomDupEfit)[i], hl.col=c("blue","red"), values=c(-1,1))
    abline(h=c(1,-1), lty=2)
    legend(legend=paste(names(summary(dt)[,i]), summary(dt)[,i], sep="="), x="bottomright", border=F, bty="n")
    text(x[o], z[t], labels=G)
}
dev.off()


# We can also look at the top ten DE genes with a heatmap of logCPM values for the top 100 genes. Each gene (or row) is scaled so that mean expression is zero and the standard deviation is one (we're using 'E' from the voom object which is a numeric matrix of normalized expression values on the log2 scale). Samples with relatively high expression of a given gene are marked in red and samples with relatively low expression are marked in blue. Lighter shades and white represent genes with intermediate expression levels. Samples and genes are reordered by the method of hierarchical clustering
# first, make a heatmap of all top genes in one pdf

# reset ensemble row names to gene symbols
rownames(vDup$E)=v$genes$SYMBOL

# set up annotation
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

df1=data.frame(village = as.character(y$samples$Sampling.Site))
df2=data.frame(batch = as.numeric(batch))
ha1 = HeatmapAnnotation(df = df1, col = list(village = c("Anakalung" =  1, "Madobag" = 2, "Mappi" = 3, "Taileleu"=4, "Wunga"=5)))

for (i in 1:ncol(efit)){
    topTable <- topTable(efit, coef=i, p.value=0.01, lfc=1, n=Inf, sort.by="p")
    index <- which(v$genes$ENSEMBL %in% topTable$ENSEMBL[1:10])
    pdf(paste0("Heatmap_",colnames(efit)[i],".pdf"), height=15, width=15)
    draw(Heatmap(t(scale(t(v$E[index,]))), col=col_fun, column_title = colnames(efit)[i], top_annotation = ha1, show_row_names = T, show_heatmap_legend = F, show_column_names = F, name = "Z-Score"),show_annotation_legend = FALSE,newpage=F)
    dev.off()
}

# We can also make individual pdfs of the top genes
village1=c("Anakalung","Madobag","Mappi","Taileleu","Wunga")
village2=c("Anakalung","Madobag","Mappi","Taileleu","Wunga")

counter=0
for (i1 in village1){
    village2=village2[-1]
    for (i2 in village2){
        counter=counter+1
        topTable <- topTable(efit, coef=counter, p.value=0.01, lfc=1, n=Inf, sort.by="p")
        index <- which(v$genes$ENSEMBL %in% topTable$ENSEMBL[1:10])
        df=data.frame(village = as.character(y$samples$Sampling.Site[grep(paste(i1,i2,sep="|"), y$samples$Sampling.Site)]))
        ha =  HeatmapAnnotation(df = df, col = list(village = c("Anakalung" =  1, "Madobag" = 2, "Mappi" = 3, "Taileleu"=4, "Wunga" = 5)))
        pdf(paste("HeatmapTopGenes",i1,i2,".pdf",sep="_"), height=10, width=15)
        draw(Heatmap(t(scale(t(v$E[index,])))[,grep(paste(i1,i2,sep="|"), y$samples$Sampling.Site)], col=col_fun, column_title = colnames(efit)[counter], top_annotation = ha, show_row_names = T, show_heatmap_legend = T, show_column_names = T, name = "Z-Score"),show_annotation_legend = TRUE,newpage=F)
        dev.off()
    }
}

# Look at which genes are in common using UpsetR

byVillages <- decideTests(voomDupEfit, method="separate", adjust.method = "BH", p.value = 0.01)
byVillages05 <- decideTests(voomDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=0.5)
byVillages1 <- decideTests(voomDupEfit, method="separate", adjust.method = "BH", p.value = 0.01, lfc=1)

pdf(paste0(outputdir, "UpsetR_SamplingSiteComparison_by_village_allfcs.pdf"), width=12)
    upset(as.data.frame(abs(byVillages)), sets = c("ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsMPI", "WNGvsMPI", "MDBvsMPI", "TLLvsMPI", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,4), rep(smb_mpi, 2), rep(mtw_mpi, 2), sumba, mentawai), nintersects=50,  order.by = "freq", keep.order=T)
    upset(as.data.frame(abs(byVillages05)), sets = c("ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsMPI", "WNGvsMPI", "MDBvsMPI", "TLLvsMPI", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,4), rep(smb_mpi, 2), rep(mtw_mpi, 2), sumba, mentawai), nintersects=50,  order.by = "freq", keep.order=T)
    upset(as.data.frame(abs(byVillages1)), sets = c("ANKvsMDB", "ANKvsTLL", "WNGvsMDB", "WNGvsTLL", "ANKvsMPI", "WNGvsMPI", "MDBvsMPI", "TLLvsMPI", "ANKvsWNG", "MDBvsTLL"), sets.bar.color = c(rep(smb_mtw,4), rep(smb_mpi, 2), rep(mtw_mpi, 2), sumba, mentawai), nintersects=50,  order.by = "freq", keep.order=T)
dev.off()

# first see which logFC threshold is best
logFC.df=matrix(nrow=3,ncol=ncol(efit))
colnames(logFC.df)=colnames(efit)
counter=0
for (number in c(0,0.5,1)){
    counter=counter+1
    dt <- decideTests(efit, p.value=0.01, lfc=number)
    values=summary(abs(dt))[3,]
    logFC.df[counter,]=values
}
# add in column specifying logFC
logFC.df=cbind(logFC = c(0,0.5,1), logFC.df)
write.table(logFC.df, file="logFC_thresholds.txt", sep="\t", quote=F)




