# by KSB, 09.22.17

# create an RNASeq pipeline for the Yamagishi et al 2014 reads from the paper: Interactive transcriptome analysis of malaria patients and infecting Plasmodium falciparum. 


# load packages
library(Rsubread)
library(RColorBrewer)
library(edgeR)
library(Homo.sapiens)
library(limma)
library(ensembldb)
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(readr)
library(openxlsx)
library(pheatmap)
library(devtools)
library(ggbiplot)
library(biomaRt)
library(biomartr)
library(gplots)
library(sva)
library(magrittr)
library(dendextend)
library(qvalue)
library(rowr)
library(reshape2)
library(RUVSeq)
library(doParallel)
library(car)
library(ggpubr)
library(GO.db)
library(goseq)
library(ggplot2)
library(ggsignif)
library(wesanderson)
library(treemap)
library(NineteenEightyR)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(vioplot)
library(ReactomePA)
library(EGSEA)

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Sick/"
housekeepingdir="/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/BatchEffects/"

# set colour palette
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
fractExpr.pal=electronic_night(n=3)
dev.off()

# Set output directory
outputdir="/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/Yamagishi/Rarefaction_Subset/"

# read in count files for sick samples
files.sick=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Yamagishi", pattern="Filter", full.names=T)
files.controls=list.files(path="/Users/katalinabobowik/Documents/UniMelb_PhD/Projects/Yamagishi/HealthyControls", pattern="Filter", full.names=T)

# set up DGE matrix combining healthy and sick samples
y <- readDGE(c(files.sick, files.controls), columns=c(1,3)) 
# Organising gene annotation using bioconductor's Homo.sapiens package, version 1.3.1 (Based on genome:  hg19)
geneid <- rownames(y)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENSEMBL")
# Check for and remove duplicated gene IDs, then add genes dataframe to DGEList object
genes <- genes[!duplicated(genes$ENSEMBL),]
y$genes <- genes
dim(y)
# [1] 27413   147

# assign healthy and control samples
y$samples$diseaseStatus[grep("Controls", colnames(y))]="control"
y$samples$diseaseStatus[grep("Controls", colnames(y), invert = T)]="malaria"
# make disease status into factor
y$samples$diseaseStatus=as.factor(y$samples$diseaseStatus)

# Trim file names into shorter sample names and apply to column names
samplenames <- sapply(1:length(colnames(y)), function(x) tail(strsplit(colnames(y)[x],"_")[[1]],1))
colnames(y) <- samplenames

# Get initial statistics before pre-processing
# Visualise library size
pdf(paste0(outputdir,"librarysizeYamagishi_preFiltering.pdf"), height=10, width=15)
barplot(y$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=as.numeric(factor(y$samples$diseaseStatus)), names=colnames(y), las=3, ylim=c(0,max(y$samples$lib.size*1e-6)+10))
legend(x="topright", col=unique(as.numeric(factor(y$samples$diseaseStatus))), legend=c("malaria", "controls"), pch=15, cex=0.8)
dev.off()

# Total number of genes
pdf(paste0(outputdir,"totalGenesYamagishi_preFiltering.pdf"), height=10, width=15)
barplot(apply(y$counts, 2, function(c)sum(c!=0)), ylab="n Genes", cex.names=0.75, col=as.numeric(factor(y$samples$diseaseStatus)), names=colnames(y), las=3, ylim=c(0,max(apply(y$counts, 2, function(c)sum(c!=0)))+2000))
legend(x="topright", col=unique(as.numeric(factor(y$samples$diseaseStatus))), legend=c("malaria", "controls"), pch=15, cex=0.8)
dev.off()

# Load in Yamagishi et al 2014 supplementary table 11 (with patient information)
sup11.sick=read.xlsx("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/Yamagishi/Supplemental_Table_11.xlsx", sheet=1)
sup11.controls=read.xlsx("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/Yamagishi/Supplemental_Table_11.xlsx", sheet=3)
sup11.all=rbind(sup11.sick, sup11.controls)

# load in SRA run table info for sick samples and controls
sra.sick=read.delim("/Users/katalinabobowik/Documents/Singapore_StemCells/Projects/Sumba/Papers/SupplementaryMaterials/SraRunTable.txt", as.is=T, strip.white=T)
sra.controls=read.delim("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/Yamagishi/SraRunTable_Controls.txt", as.is=T, strip.white=T)
sra.sick=sra.sick[,c("Library_Name_s","Run_s")]
sra.controls=sra.controls[,c("Library_Name","Run")]
colnames(sra.sick)=colnames(sra.controls)
sra.all=rbind(sra.sick,sra.controls)

# Create empty matrix and fill in with patient names, Age, Sex, and SRA Run ID
sample_summary=matrix("NA",nrow=nrow(sup11.all),ncol=7)
sample_summary[,1:5]=c(sup11.all$Patient_Name,sup11.all$Age, sup11.all$Sex, sup11.all[,"%Pf_tags"], sup11.all$From)
colnames(sample_summary)=c("Patient_Name", "Age", "Sex", "PF_Tags","From", "sample_ID", "SRA_ID")
# convert to data frame
sample_summary=as.data.frame(sample_summary)

# There is a discrepancy between the SRA table and Supplementary 11 table- one says "malaria7#09" and one says "malaria7#009". We need to change "malaria7#09" to "malaria7#009"
sample_summary$Patient_Name=gsub("malaria7#09", "malaria7#009", sample_summary$Patient_Name)

# Match patient names with the patient names in the SRA table, then grab the corresponding SRA run IDs
sra.all[match(sample_summary$Patient_Name, sra.all$Library_Name),"Run"]
# we seem to have NA values in the dataframe. Which ones are these?
sample_summary$Patient_Name[which(is.na(match(sample_summary$Patient_Name, sra.all$Library_Name)))]
# [1] "malaria11#1" "malaria11#2" "malaria11#3" "malaria11#4" "malaria11#5"
# [6] "malaria11#6" "malaria11#7" "malaria11#8" "malaria11#9"

# when checking the sra run table sheet, these samples have a 0 in front of them. Let's take this out to keep sample names the same.
sra.all$Library_Name=gsub("malaria11#0","malaria11#", sra.all$Library_Name)
sample_summary[,"SRA_ID"]=sra.all[match(sample_summary$Patient_Name, sra.all$Library_Name),"Run"]

# Assign gender, age, location, and PF load
y$samples$gender <- as.factor(sample_summary[match(colnames(y), sample_summary$SRA_ID),"Sex"])
y$samples$age <- as.numeric(sample_summary[match(colnames(y), sample_summary$SRA_ID),"Age"])
y$samples$location <- as.factor(sample_summary[match(colnames(y), sample_summary$SRA_ID),"From"])
y$samples$PFload <- as.numeric(as.character(sample_summary[match(colnames(y), sample_summary$SRA_ID),"PF_Tags"]))

# let's also add our own PF load made by getting the total number of unmapped reads that mapped to the combined PFPX genome
alignedPFPXreads=read.table(paste0(inputdir,"Malaria_summary_table_Yamagishi.txt"), header=T)
colnames(alignedPFPXreads)[1]="SampleID"
alignedPFPXreads$SampleID=gsub("_Controls","",alignedPFPXreads$SampleID) %>% gsub("_Sick","",.)
identical(alignedPFPXreads$SampleID, colnames(y))
# [1] TRUE
y$samples$alignedPFPX=alignedPFPXreads$fract.reads.pfpx.yam

#assign covariates to variables
gender <- addNA(y$samples$gender)
age <- addNA(cut(as.numeric(as.character(y$samples$age)), c(0,15,30,70), labels=c("0-15", "15-30", "30-70")))
lib.size <- cut(as.numeric(as.character(y$samples$lib.size)), c(0,1000000,10000000,20000000,40000000, 60000000), labels=c("0-1","1-10", "10-20", "20-40", "40-60"))
PFload=addNA(cut(y$samples$PFload, c(-1,20,40,60,80), labels=c("0-20","20-40", "40-60", "60-80")))
diseaseStatus=y$samples$diseaseStatus
location=addNA(y$samples$location)
alignedPFPX=y$samples$alignedPFPX

# assign covariate names
covariate.names=c("gender", "age", "lib.size", "PFload","diseaseStatus","location")

png(paste0(outputdir,"rarefactionCurves.png"))
for (name in covariate.names) {
  plot(1:length(y$counts[,1]), cumsum(sort(y$counts[,1], decreasing=T)/sum(y$counts[,1])), log="x", type="n", xlab="Number of genes", ylab="Fraction of reads pool", ylim=c(0,1), main=name) ## initialize the plot area
  counter=0
  for (sample in colnames(y)){
    counter=counter+1
    lines(1:length(y$counts[,sample]), cumsum(sort(y$counts[,sample], decreasing=T)/sum(y$counts[,sample])), lwd=2, col=as.numeric(get(name))[counter])
  }
  levels=levels(get(name))
  levels[which(is.na(levels))] = "NA"
  legend(x="bottomright", bty="n", col=1:length(levels(get(name))), legend=levels, lty=1, lwd=2)
}
dev.off()

# get rid of globing genes
globin_genes = c("HBA1", "HBA2", "HBB", "MB", "CYGB", "NGB")
remove_globin_index = grep(paste("^",globin_genes,"$",sep="",collapse="|"), y$genes$SYMBOL)
y=y[-c(remove_globin_index),, keep.lib.sizes=FALSE]

# wait, what? Why do the rarefactions curve look like this?? It seems like some samples are expressing the majority of their genes at less than ten genes!
# let's filter the data (to get rid of outlier samples) and then look at the rarefaction curves again.

# Data pre-processing ------------------------------------------------------------------------

# Filter out samples with library size <5 million and samples with no gender information and assign covariate names
y=y[,which(y$samples$lib.size >= 9000000)]
dim(y)
# [1] 27413  124

minFractExpr <- data.frame(Min=integer(), FirstQ=integer(), Median=integer(),Mean=integer(),ThirdQ=integer(),Max=integer())
counter=0
for (sample in colnames(y)) {
  counter=counter + 1
  cumulativeSum=cumsum(sort(y$counts[,sample], decreasing=T)/sum(y$counts[,sample]))
  summaryCumulativeSum=summary(cumulativeSum)
  minFractExpr[counter,] <- summaryCumulativeSum
}
rownames(minFractExpr)=colnames(y)

# order by minimum number of genes expresses
minFractExpr=minFractExpr[order(minFractExpr$Min),]
susSamples=rownames(tail(minFractExpr, n=10))
normalSamples=rownames(head(minFractExpr, n=10))

# add these in as a group to the dataframe
y$samples$FractExpr=rep("Intermediary", ncol(y))
y$samples[susSamples,]$FractExpr="High"
y$samples[normalSamples,]$FractExpr="Low"
y$samples$FractExpr=as.factor(y$samples$FractExpr)

# reassign covariate names
covariate.names=c(covariate.names, "FractExpr")

# reassign covariates to variables
gender <- addNA(y$samples$gender)
age <- addNA(cut(as.numeric(as.character(y$samples$age)), c(0,15,30,70), labels=c("0-15", "15-30", "30-70")))
lib.size <- cut(as.numeric(as.character(y$samples$lib.size)), c(0,1000000,10000000,20000000,40000000, 60000000), labels=c("0-1","1-10", "10-20", "20-40", "40-60"))
PFload=addNA(cut(y$samples$PFload, c(-1,20,40,60,80), labels=c("0-20","20-40", "40-60", "60-80")))
diseaseStatus=y$samples$diseaseStatus
location=addNA(y$samples$location)
alignedPFPX=y$samples$alignedPFPX
FractExpr=y$samples$FractExpr

png(paste0(outputdir,"rarefactionCurves_postfiltering.png"))
for (name in covariate.names) {
  plot(1:length(y$counts[,1]), cumsum(sort(y$counts[,1], decreasing=T)/sum(y$counts[,1])), log="x", type="n", xlab="Number of genes", ylab="Fraction of reads pool", ylim=c(0,1), main=name) ## initialize the plot area
  counter=0
  for (sample in colnames(y)){
    counter=counter+1
    lines(1:length(y$counts[,sample]), cumsum(sort(y$counts[,sample], decreasing=T)/sum(y$counts[,sample])), lwd=2, col=as.numeric(get(name))[counter])
  }
  levels=levels(get(name))
  levels[which(is.na(levels))] = "NA"
  legend(x="bottomright", bty="n", col=1:length(levels(get(name))), legend=levels, lty=1, lwd=2)
}
dev.off()

# Visualise library size after filtering
pdf(paste0(outputdir,"librarysizeYamagishi_postFiltering.pdf"), height=10, width=15)
barplot(y$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=fractExpr.pal[as.numeric(factor(y$samples$FractExpr))], names=colnames(y), las=3, ylim=c(0,max(y$samples$lib.size*1e-6)+10))
legend(x="topright", col=fractExpr.pal[unique(as.numeric(factor(y$samples$FractExpr)))], legend=unique(factor(y$samples$FractExpr)), pch=15, cex=0.8)
dev.off()

# n genes
pdf(paste0(outputdir,"totalGenesYamagishi_postFiltering.pdf"), height=10, width=15)
barplot(apply(y$counts, 2, function(c)sum(c!=0)), ylab="n Genes", cex.names=0.75, col=fractExpr.pal[as.numeric(factor(y$samples$FractExpr))], names=colnames(y), las=3, ylim=c(0,max(apply(y$counts, 2, function(c)sum(c!=0)))+2000))
legend(x="topright", col=fractExpr.pal[unique(as.numeric(factor(y$samples$FractExpr)))], legend=unique(factor(y$samples$FractExpr)), pch=15, cex=0.8)
dev.off()

# How does this look like for the PF load?
pdf(paste0(outputdir,"totalPFPXload_fractExpr.pdf"), height=10, width=15)
barplot(y$samples$alignedPFPX, col=fractExpr.pal[as.numeric(factor(y$samples$FractExpr))])
legend(x="topright", col=fractExpr.pal[unique(as.numeric(factor(y$samples$FractExpr)))], legend=unique(factor(y$samples$FractExpr)), pch=15, cex=0.8)
dev.off()

# Transformation from the raw scale --------------------------------------------------------------------

# Transform raw counts onto a scale that accounts for library size differences. Here, we transform to CPM and log-CPM values (prior count for logCPM = 0.25). 
cpm <- cpm(y) 
lcpm <- cpm(y, log=TRUE)
susSamples.lcpm=lcpm[,susSamples]
normalSamples.lcpm=lcpm[,normalSamples]

# plot the most highly expressed gene for each of our ten suspicious samples
pdf(paste0(outputdir,"topGenes_FractExpr.pdf"), height=10, width=15)
for (i in 1:length(susSamples)){
  topExprGene=tail(susSamples.lcpm[order(susSamples.lcpm[,i]),susSamples[i]], n=10)
  names(topExprGene)=y[names(topExprGene),]$genes$SYMBOL
  barplot(rev(topExprGene),las=3,col="coral", )
}
dev.off()

pdf(paste0(outputdir,"topGenes_FractExpr_normal.pdf"), height=10, width=15)
for (i in 1:length(normalSamples)){
  topExprGene=tail(normalSamples.lcpm[order(normalSamples.lcpm[,i]),normalSamples[i]], n=10)
  names(topExprGene)=y[names(topExprGene),]$genes$SYMBOL
  barplot(rev(topExprGene),las=3,col="purple", )
}
dev.off()

# make a barplot of the most highly expressed gene 
top=vector()
for (i in 1:ncol(y)){
  topExprGene=tail(lcpm[order(lcpm[,i]),i], n=1)
  top=c(top,topExprGene)
}

pdf(paste0(outputdir,"topGenes_FractExpr_allSamples.pdf"), height=10, width=15)
barplot(top,las=3,col=fractExpr.pal[as.numeric(y$samples$FractExpr)])
legend(x="topright", col=fractExpr.pal[unique(as.numeric(factor(y$samples$FractExpr)))], legend=unique(factor(y$samples$FractExpr)), pch=15, cex=0.8)
dev.off()

# Remove genes that are lowly expressed -----------------------------------------

# get histogram of number of genes expressed at log2 cpm > 0.5 and 1 (before filtering)
pdf(paste0(outputdir,"historgram_nGenesExpressed.pdf"), height=10, width=15)
hist(rowSums(cpm>0.5), main= "n Genes expressed at cpm > 0.5 \n pre-filtering", xlab="samples", col=4)
hist(rowSums(cpm>1), main= "n Genes expressed at cpm > 1 \n pre-filtering", xlab="samples", col=4)
dev.off()

# Remove genes that are lowly expressed- a gene is only retained if it is expressed at log-CPM > 1 in at least half of the libraries
keep.exprs <- rowSums(cpm>1) >= (nrow(y$samples)*0.5)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y)
# [1] 11096   124

# Compare library sizes before and after removing lowly-expressed genes
nsamples <- ncol(y)
col <- as.numeric(factor(y$samples$diseaseStatus))
pdf(paste0(outputdir,"libraryDensity_afterFiltering.pdf"), height=10, width=15)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm)$y)), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", levels(gender), ncol=1, cex=0.6, text.col=unique(col), bty="n")

lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm)$y)+0.2), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(y), ncol=2, cex=0.6, text.col=col, bty="n")
dev.off()

# Normalise gene expression distributions (i.e., no bias introduced during sample preparation/sequencing)
y <- calcNormFactors(y, method = "TMM")

# Duplicate data, set normalisation back to 1, and plot difference between normalised and non-normalised data
y2 <- y
y2$samples$norm.factors <- 1
lcpm <- cpm(y2, log=TRUE)
pdf(paste0(outputdir,"NormalisedGeneExpressionDistribution_YamagishiHealthyvsControls.pdf"), height=15, width=15)
boxplot(lcpm, las=2, col=col, main="", cex.axis=0.75)
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
y2 <- calcNormFactors(y2)
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="", cex.axis=0.75)
title(main="B. Example: Normalised data",ylab="Log-cpm")
dev.off()

# get density plot after normalisation
pdf(paste0(outputdir,"libraryDensity_afterFilteringandNormalisation.pdf"), height=10, width=15)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm)$y)+0.2), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(y), ncol=2, cex=0.6, text.col=col, bty="n")
dev.off()

# Data exploration --------------------------------------------------------------------

# Reassign covariates
gender <- addNA(y$samples$gender)
age <- addNA(cut(as.numeric(as.character(y$samples$age)), c(0,15,30,70), labels=c("0-15", "15-30", "30-70")))
lib.size <- cut(as.numeric(as.character(y$samples$lib.size)), c(0,1000000,10000000,20000000,40000000, 60000000), labels=c("0-1","1-10", "10-20", "20-40", "40-60"))
PFload=addNA(cut(y$samples$PFload, c(-1,20,40,60,80), labels=c("0-20","20-40", "40-60", "60-80")))
diseaseStatus=y$samples$diseaseStatus
location=addNA(y$samples$location)
alignedPFPX=y$samples$alignedPFPX
FractExpr=y$samples$FractExpr

# plot MDS
pdf(paste0(outputdir,"MDSPLots.pdf"), height=10, width=15)
for (name in covariate.names) {
    plotMDS(lcpm, labels=get(name), col=as.numeric(get(name)))
    title(main=name)
}
dev.off()


# PCA plotting function
plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(dataToPca), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    for (i in 1:9){
        pca_axis1=i
        pca_axis2=i+1
        plot(pca$x[,pca_axis1], pca$x[,pca_axis2], col=speciesCol, pch=namesPch, cex=2, xlab=paste0("PC", pca_axis1, " (", round(pca.var[pca_axis1]*100, digits=2), "% of variance)"), ylab=paste0("PC", pca_axis2, " (", round(pca.var[pca_axis2]*100, digits=2), "% of variance)", sep=""), main=name)
       	legend(legend=unique(sampleNames), pch=16, x="bottomright", col=unique(speciesCol), cex=0.6, title=name, border=F, bty="n")
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
all.covars.df <- y$samples[,c(covariate.names, "alignedPFPX")] 

# Plot PCA
for (name in covariate.names){
    pdf(paste0(outputdir,"pcaresults_",name,".pdf"))
    pcaresults <- plot.pca(dataToPca=lcpm, speciesCol=as.numeric(get(name)),namesPch=15,sampleNames=get(name))
    dev.off()
}

# plot numeric variables
name="alignedPFPX"
initial = .bincode(alignedPFPX, breaks=seq(min(alignedPFPX, na.rm=T), max(alignedPFPX, na.rm=T), len = 80),include.lowest = TRUE)
bloodCol <- colorRampPalette(c("blue", "red"))(79)[initial]
pdf(paste0(outputdir,"pcaresults_alignedPlasmodiumLoad.pdf"))
pcaresults <- plot.pca(dataToPca=lcpm, speciesCol=bloodCol,namesPch=15,sampleNames=alignedPFPX)
legend(legend=c("High","Low"), pch=16, x="bottomright", col=c(bloodCol[which.max(alignedPFPX)], bloodCol[which.min(alignedPFPX)]), cex=0.6, title=name, border=F, bty="n")
dev.off()

# Get PCA associations
all.pcs <- pc.assoc(pcaresults)
all.pcs$Variance <- pcaresults$sdev^2/sum(pcaresults$sdev^2)

# plot pca covariates association matrix to illustrate any potential confounding and evidence for batches
pdf(paste0(outputdir,"significantCovariates_AnovaHeatmap.pdf"))
pheatmap(all.pcs[1:5,c(covariate.names,"alignedPFPX")], cluster_col=F, col= colorRampPalette(brewer.pal(11, "RdYlBu"))(100), cluster_rows=F, main="Significant Covariates \n Anova")
dev.off()

# Write out the covariates:
write.table(all.pcs, file=paste0(outputdir,"pca_covariates_significanceLevels.txt"), col.names=T, row.names=F, quote=F, sep="\t")

# look for sample outliers from PCA
pca.outliers.final=matrix(nrow=0, ncol=3)

for (i in 1:ncol(pcaresults$x)){
    pca.dim=c()
    outlier.sample=c()
    outlier.zscore=c()
    zscore=scale(pcaresults$x[,i])
    outliers=which(abs(zscore) >= 3)
    if (length(outliers) > 0){
        pca.dim=c(pca.dim, i)
        outlier.sample=c(outlier.sample, names(pcaresults$x[,i][outliers]))
        outlier.zscore=c(outlier.zscore, zscore[outliers])
        pca.outliers=matrix(c(rep(pca.dim,length(outlier.sample)), outlier.sample, outlier.zscore), nrow=length(outlier.sample), ncol=3)
        pca.outliers.final=rbind(pca.outliers.final, pca.outliers)
    }
}
colnames(pca.outliers.final)=c("Pca.dim", "Samples", "Z.score")
write.table(pca.outliers.final, file=paste0(outputdir,"sample_outliersInPCA.txt"), quote=F, row.names=F)

# Dissimilarity matrix with euclidean distances
pdf(paste0(outputdir,"SampleDistances_Euclidean.pdf"), height=8, width=15)
par(mar=c(6.1,4.1,4.1,2.1))
eucl.distance <- dist(t(lcpm), method = "euclidean")
eucl.cluster <- hclust(eucl.distance, method = "complete")
dend.eucl=as.dendrogram(eucl.cluster)
labels_colors(dend.eucl)=as.numeric(diseaseStatus)[order.dendrogram(dend.eucl)]
plot(dend.eucl, main="log2-CPM \n Euclidean Distances")
labels_colors(dend.eucl)=as.numeric(FractExpr)[order.dendrogram(dend.eucl)]
plot(dend.eucl, main="log2-CPM \n Euclidean Distances")
dev.off()

# Heatmap of lcpm distances
# set up annotation
df1=data.frame(DiseaseStatus = as.character(diseaseStatus))
df2=data.frame(fractExpr = as.character(FractExpr))
ha1 = HeatmapAnnotation(df = df1, col = list(DiseaseStatus = c("malaria" =  1, "control" = 2)))
ha2 = rowAnnotation(df = df2, col= list(fractExpr=c("High" = fractExpr.pal[1], "Intermediary" = fractExpr.pal[2], "Low" = fractExpr.pal[3])))

# plot
pdf(paste0(outputdir,"lcpmCorrelationHeatmaps.pdf"), height=10, width=15)
Heatmap(cor(lcpm,method="pearson"), col=magma(100), column_title = "Pearson Correlation \n log2-CPM", name="Corr Coeff", top_annotation = ha1, show_row_names = FALSE, column_names_gp=gpar(fontsize = 8)) + ha2
dev.off()

# DE analysis ------------------

# First, let's look at what the top genes are when comparing groups where the fraction of expressed genes is different
# Set up design matrix
design <- model.matrix(~0 + y$samples$FractExpr + y$samples$diseaseStatus + y$samples$alignedPFPX)
#design <- model.matrix(~0 + y$samples$FractExpr)
colnames(design)=gsub("diseaseStatus", "", colnames(design))
colnames(design)=gsub("[\\y$]", "", colnames(design))
colnames(design)=gsub("samples", "", colnames(design))

# set up contrast matrix
contr.matrix <- makeContrasts(HighExprVsLow=FractExprHigh - FractExprLow, HighExprVsInter=FractExprHigh - FractExprIntermediar, IntermediarExprVsLow=FractExprIntermediar - FractExprLow, levels=colnames(design))

# Remove heteroscedascity from count data
pdf(paste0(outputdir,"Limma_voom_TMM_cyclicLoess.pdf"))
v <- voom(y, design, plot=TRUE)
# fit linear models
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)
dev.off()

dt <- decideTests(efit,p.value=0.01,lfc=1)
summary(dt)

#       HighExprVsLow HighExprVsInter IntermediarExprVsLow
#Down              24               0                   23
#NotSig         10709           10808                11014
#Up               363             288                   59

# look at different DE thresholds
# test different logFC thresholds
logFC.df=matrix(nrow=3,ncol=1)
counter=0
for (number in c(0,0.5,1)){
    counter=counter+1
    dt <- decideTests(efit, p.value=0.01, lfc=number)
    logFC.df[counter,]=sum(abs(dt))
}
logFC.df=cbind(logFC = c(0,0.5,1), logFC.df)
write.table(logFC.df, file=paste0(outputdir,"logFC_thresholds.txt"))

# get top genes
toptable=topTable(efit, coef=1, p.value=0.01, n=Inf, sort.by="p")
write.table(toptable,file=paste0(outputdir,"topTable_fractExprGenes.txt"))

# Visual QC of duplicate correlation voom output after fitting linear models --------------------------------------------------------------------------------------

# check to see p-value distribution is normal
pdf(paste0(outputdir,"PvalueDist_NotAdjusted_dupCor.pdf"), height=15, width=10)
par(mfrow=c(3,1))
for (i in 1:ncol(efit)){
    hist(efit$p.value[,i], main=colnames(efit)[i], ylim=c(0,max(table(round(efit$p.value[,i], 1)))+1000), xlab="p-value")
}
dev.off()

# check p-value distribution for adjusted p-values
pdf(paste0(outputdir,"PvalueDist_Adjusted_dupCor.pdf"), height=15, width=10)
par(mfrow=c(3,1))
for (i in 1:ncol(efit)){
    topTable <- topTable(efit, coef=i, n=Inf)
    hist(topTable$adj.P.Val, main=colnames(efit)[i], ylim=c(0,1500), xlab="p-value")
}
dev.off()

# Verify that control housekeeping genes are not significantly DE. Set up list of housekeeping genes as controls (from Eisenberg and Levanon, 2003)
housekeeping=read.table(paste0(housekeepingdir,"Housekeeping_ControlGenes.txt"), as.is=T, header=F)
# if this is broken, use host = "uswest.ensembl.org"
ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'www.ensembl.org', ensemblRedirect = FALSE)
biomart.results.table <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=housekeeping, filters="hgnc_symbol")
hkGenes=as.vector(biomart.results.table[,1])
hkControls=hkGenes[which(hkGenes %in% rownames(y$counts))]

# Volcano plot with points of housekeeping genes
pdf(paste0(outputdir,"VolcanoPlots_fractExpr.pdf"), height=15, width=10)
par(mfrow=c(3,1))
for (i in 1:ncol(efit)){
    plot(efit$coef[,i], -log10(as.matrix(efit$p.value)[,i]), pch=20, main=colnames(efit)[i], xlab="log2FoldChange", ylab="-log10(pvalue)")
    points(efit$coef[,i][which(names(efit$coef[,i]) %in% hkControls)], -log10(as.matrix(efit$p.value)[,i][which(names(efit$coef[,i]) %in% hkControls)]) , pch=20, col=4, xlab="log2FoldChange", ylab="-log10(pvalue)")
    legend("topleft", "genes", "hk genes",fill=4)
    abline(v=c(-1,1))
}
dev.off()

# PCA visualisation after correction and association with covariates ------------------------------------------------------------

# let's also visualise how our PCAs look after limma correction by using removeBatcheffect. Help on design of removeBatcheffects was given by the lovely John Blischak.
design <- model.matrix(~0 + y$samples$FractExpr)
colnames(design)=gsub("[\\y$]", "", colnames(design))
colnames(design)=gsub("samples", "", colnames(design))
batch.corrected.lcpm <- removeBatchEffect(lcpm, covariates = cbind(y$samples$diseaseStatus, y$samples$alignedPFPX))

# Plot PCA
for (name in covariate.names){
    pdf(paste0(outputdir,"pcaresults_",name,"_removeBatchEffect.pdf"))
    pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm, speciesCol=as.numeric(get(name)),namesPch=15,sampleNames=get(name))
    dev.off()
}

# plot numeric variables
name="alignedPFPX"
initial = .bincode(alignedPFPX, breaks=seq(min(alignedPFPX, na.rm=T), max(alignedPFPX, na.rm=T), len = 80),include.lowest = TRUE)
bloodCol <- colorRampPalette(c("blue", "red"))(79)[initial]
pdf(paste0(outputdir,"pcaresults_alignedPlasmodiumLoad.pdf"))
pcaresults <- plot.pca(dataToPca=batch.corrected.lcpm,speciesCol=bloodCol,namesPch=15,sampleNames=alignedPFPX)
legend(legend=c("High","Low"), pch=16, x="bottomright", col=c(bloodCol[which.max(alignedPFPX)], bloodCol[which.min(alignedPFPX)]), cex=0.6, title=name, border=F, bty="n")
dev.off()

# Get PCA associations
all.pcs <- pc.assoc(pcaresults)
all.pcs$Variance <- pcaresults$sdev^2/sum(pcaresults$sdev^2)

# plot pca covariates association matrix to illustrate any potential confounding and evidence for batches
pdf(paste0(outputdir,"significantCovariates_AnovaHeatmap_removeBatchEffect.pdf"))
pheatmap(all.pcs[1:5,c(covariate.names,"alignedPFPX")], cluster_col=F, col= colorRampPalette(brewer.pal(11, "RdYlBu"))(100), cluster_rows=F, main="Significant Covariates \n Anova")
dev.off()

# Visual exploration of top DE genes ------------------------------------

# plot log2 fold change between islands
pdf(paste0(outputdir,"log2FC_FractExprComparisons_pval01.pdf"))
# note 'p.value' is the cutoff value for adjusted p-values
topTable <- topTable(efit, coef=1, n=Inf, p.value=0.01)
plot(density(topTable$logFC), col=9, xlim=c(-3,8), ylim=c(0,1), main="LogFC Density", xlab="LogFC", ylab="Density")
abline(v=c(-1,1), lty=3)
counter=0
for (i in 2:ncol(efit)){
    counter=counter+1
    topTable <- topTable(efit, coef=i, n=Inf, p.value=0.01)
    lines(density(topTable$logFC), col=9+counter, xlim=c(-2,8), ylim=c(0,1))
}
legend(x="topright", bty="n", col=9:11, legend=colnames(efit), lty=1, lwd=2)
dev.off()

# graphical representation of DE results through MD plot
pdf(paste0(outputdir,"MD_Plots_pval01_lfc1_fractExpr_highVsLow.pdf"), height=10, width=15)
o <- which(names(efit$Amean) %in% names(which(abs(dt[,1])==1)))
x <- efit$Amean
z <- efit$coefficients[,1]
t=which(names(efit$coefficients[,1]) %in% names(which(abs(dt[,1])==1)))
G <- efit$genes[which(abs(dt[,1])==1),]$SYMBOL
plotMD(efit, column=1, status=dt[,1], main=colnames(efit)[1], hl.col=c("blue","red"), values=c(-1,1))
abline(h=c(1,-1), lty=2)
legend(legend=paste(names(summary(dt)[,1]), summary(dt)[,1], sep="="), x="bottomright", border=F, bty="n")
text(x[o], z[t], labels=G)
dev.off()

# We can also look at the top ten DE genes with a heatmap of logCPM values for the top 100 genes. Each gene (or row) is scaled so that mean expression is zero and the standard deviation is one (we're using 'E' from the voom object which is a numeric matrix of normalized expression values on the log2 scale). Samples with relatively high expression of a given gene are marked in red and samples with relatively low expression are marked in blue. Lighter shades and white represent genes with intermediate expression levels. Samples and genes are reordered by the method of hierarchical clustering
# first, make a heatmap of all top genes in one pdf

# set up annotation
FractExprHigh=topTable(efit, coef=1, p.value=0.01, lfc=1, n=Inf, sort.by="p")
highVsLow <- FractExprHigh$ENSEMBL[1:100]
i <- which(v$genes$ENSEMBL %in% highVsLow)
mycol <- colorpanel(1000,"blue","white","red")

# plot all samples
pdf(paste0(outputdir,"heatmap_pval01_FDR1_FractExprs-HighvsLow.pdf"), height=15, width=15)
heatmap.2(v$E[i,], scale="row", labRow=v$genes$SYMBOL[i], labCol=y$samples$diseaseStatus,col=mycol, trace="none", density.info="none", margin=c(8,6), lhei=c(2,10), dendrogram="column")
dev.off()

# now just for the top ten samples
pdf(paste0(outputdir,"heatmap_pval01_FDR1_FractExprs-HighvsLow_100.pdf"), height=15, width=15)
heatmap.2(v$E[i,c(susSamples,normalSamples)], scale="row", labRow=v$genes$SYMBOL[i], labCol=c(rep("high",10), rep("low",10)),col=mycol, trace="none", density.info="none", margin=c(8,6), lhei=c(2,10), dendrogram="column")
dev.off()

# show the number of DE genes between all islands
pdf(paste0(outputdir,"vennDiagram_allSigDEGenes_pval01_FDR1_FractExpr.pdf"), height=15, width=15)
vennDiagram(dt[,1:3], circle.col=c(9,10,11))
dev.off()

# Enrichment analysis for Gene Ontology ----------------------------------------------------------------------------------------------------

# transform ensembl IDs to entrez IDs to be compatible with human c2 dataset (below)
ensembl_genes=rownames(y)
entrez=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "description"),values= ensembl_genes,mart=ensembl.mart.90)
y$entrezID=entrez[match(rownames(y), entrez[,1]), 2]

# gene set testing with Camera
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata")) 
idx <- ids2indices(Hs.c2,id=y$entrezID) 
camera.matrix=camera(v,idx,design,contrast=contr.matrix[,1])
write.table(camera.matrix, file=paste0(outputdir,"cameraMatrix_healtyvsSick.txt"))

# EGSEA -----------------------------------------------------------------------

# transform ensembl IDs to entrez IDs to be compatible with human c2 dataset (below)
v$genes$entrezID=entrez[match(rownames(y), entrez[,1]), 2]

#subset voom object to contain only genes with an Entrez Gene ID
v <- v[which(!is.na(v$genes$entrezID)),]
v <- v[-c(which(duplicated(v$genes$entrezID))),]
rownames(v) <- v$genes$entrezID
rownames(v$genes) <- v$genes$entrezID

#build gene-set collections to compute enrichment for
gs.annots = buildIdx(entrezIDs=v$genes$entrezID, species="human",msigdb.gsets=c("c2", "c5"), go.part = TRUE)

#generate a symbolsMap
symbolsMap = v$genes[, c(4, 2)]
colnames(symbolsMap) = c("FeatureID", "Symbols")
symbolsMap[, "Symbols"] = as.character(symbolsMap[, "Symbols"])

#establish the list of independent gene-set enrichment methods for EGSEA to use
baseMethods = egsea.base()[-2]

# Ensemble testing with EGSEA
gsa = egsea(voom.results=v, contrasts=contr.matrix, gs.annots=gs.annots, symbolsMap=symbolsMap, baseGSEAs=baseMethods, sort.by="med.rank", num.threads = 8, report = FALSE)

# get a summary of the top gene sets
plotSummaryHeatmap(gsa, gs.label="c2", hm.vals = "avg.logfc.dir",file.name=paste0(outputdir,"summary_heatmaps_c2"), sort.by="med.rank")
plotSummaryHeatmap(gsa, gs.label="kegg", hm.vals = "avg.logfc.dir",file.name=paste0(outputdir,"summary_heatmaps_kegg"), sort.by="med.rank")
plotBars(gsa, gs.label = "c2", contrast = 1, file.name=paste0(outputdir,"comparison-c2-bars"))

# get top gene sets
topc2=topSets(gsa, contrast = "comparison",sort.by="med.rank", gs.label="c2", number=20)
topkegg=topSets(gsa, contrast = "comparison",sort.by="med.rank", gs.label="kegg")

# let's make a heatmap of all of the top enriched pathways in the c2 dataset
for (i in 1:20){
    plotHeatmap(gsa, gene.set=topc2[i], gs.label="c2", contrast = 1, file.name = paste0(outputdir,"c2_ComparisonHeatmap_",topc2[i]))
    plotPathway(gsa, gene.set = topkegg[i],contrast = "comparison", gs.label = "kegg",file.name = paste0(outputdir,"kegg_ComparisonPathway_",topkegg[i]))
}

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")



# DE between healthy and sick samples ------------------------------------

outputidr="/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/Yamagishi/"
# First, let's look at what the top genes are when comparing groups where the fraction of expressed genes is different
# Set up design matrix
design <- model.matrix(~0 + y$samples$diseaseStatus + y$samples$FractExpr + y$samples$alignedPFPX)
colnames(design)=gsub("diseaseStatus", "", colnames(design))
colnames(design)=gsub("[\\y$]", "", colnames(design))
colnames(design)=gsub("samples", "", colnames(design))

contr.matrix <- makeContrasts(HealthyvsSick=control - malaria, levels=colnames(design))

# Remove heteroscedascity from count data
pdf(paste0(outputdir,"Limma_voom_TMM_cyclicLoess.pdf"))
v <- voom(y, design, plot=TRUE)
# fit linear models
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)
dev.off()

dt <- decideTests(efit,p.value=0.01,lfc=1)
summary(dt)

#       HealthyvsSick
#Down            1397
#NotSig          8712
#Up               987

# look at different DE thresholds
# test different logFC thresholds
logFC.df=matrix(nrow=3,ncol=1)
counter=0
for (number in c(0,0.5,1)){
    counter=counter+1
    dt <- decideTests(efit, p.value=0.01, lfc=number)
    logFC.df[counter,]=sum(abs(dt))
}
logFC.df=cbind(logFC = c(0,0.5,1), logFC.df)
write.table(logFC.df, file=paste0(outputdir,"logFC_thresholds.txt"))

# get top genes
toptable=topTable(efit, coef=1, p.value=0.01, n=Inf, sort.by="p")
write.table(toptable,file=paste0(outputdir,"topTable_fractExprGenes.txt"))

# Visual exploration of top DE genes ------------------------------------

# plot log2 fold change between islands
pdf(paste0(outputdir,"log2FC_FractExprComparisons_pval01.pdf"))
# note 'p.value' is the cutoff value for adjusted p-values
topTable <- topTable(efit, coef=1, n=Inf, p.value=0.01)
plot(density(topTable$logFC), col=9, xlim=c(-3,3), ylim=c(0,0.6), main="LogFC Density", xlab="LogFC", ylab="Density")
abline(v=c(-1,1), lty=3)
dev.off()

# We can also look at the top ten DE genes with a heatmap of logCPM values for the top 100 genes. Each gene (or row) is scaled so that mean expression is zero and the standard deviation is one (we're using 'E' from the voom object which is a numeric matrix of normalized expression values on the log2 scale). Samples with relatively high expression of a given gene are marked in red and samples with relatively low expression are marked in blue. Lighter shades and white represent genes with intermediate expression levels. Samples and genes are reordered by the method of hierarchical clustering
# first, make a heatmap of all top genes in one pdf

# set up annotation
malariaVsControl=topTable(efit, coef=1, p.value=0.01, lfc=1, n=Inf, sort.by="p")
malariaVscontrol.DEGenes <- malariaVsControl$ENSEMBL[1:10]
i <- which(v$genes$ENSEMBL %in% malariaVscontrol.DEGenes)
mycol <- colorpanel(1000,"blue","white","red")

# plot all samples
png(paste0(outputdir,"heatmap_pval01_FDR1_FractExprs-HighvsLow.png"))
heatmap.2(v$E[i,], scale="row", labRow=v$genes$SYMBOL[i], labCol=y$samples$diseaseStatus,col=mycol, trace="none", density.info="none", margin=c(8,6), lhei=c(2,10), dendrogram="column")
dev.off()


------------

# Enrichment analysis for Gene Ontology ----------------------------------------------------------------------------------------------------

ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'www.ensembl.org', ensemblRedirect = FALSE)

# transform ensembl IDs to entrez IDs to be compatible with human c2 dataset (below)
ensembl_genes=rownames(y)
entrez=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "description"),values= ensembl_genes,mart=ensembl.mart.90)
y$entrezID=entrez[match(rownames(y), entrez[,1]), 2]

# gene set testing with Camera
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata")) 
idx <- ids2indices(Hs.c2,id=y$entrezID) 
camera.matrix=camera(v,idx,design,contrast=contr.matrix[,1])
write.table(camera.matrix, file=paste0(outputdir,"cameraMatrix_healtyvsSick.txt"))







