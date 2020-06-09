# Test out ABIS for deconvoluting blood using RNASeq data from 123 Indonesian samples
# Code developed by Katalina Bobowik, 20.03.2019
# using the shiny app from the paper from Monaco et al, 2019: https://www.cell.com/cell-reports/fulltext/S2211-1247(19)30059-2


# load packages
library(magrittr)
library(ggplot2)
library(ggpubr)
library(reshape2)

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/"
outputdir= "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/bloodDecon/"
deconestimatedir= "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/"

# Read in data

ABIS = read.table(paste0(inputdir, "ABIS/predictedCellCounts_ABIS_scaled.txt"))
# deconCell = read.table(paste0(inputdir, "predictedCellCounts_DeconCell_Paper.txt"))
# here, I'm using the unscaled values instead of the scaled ones, since Limma does not like all blood porprtions to add up to %100 in the design matrix
deconCell = read.table(paste0(inputdir, "predictedCellCounts_DeconCell.txt"))
methylation = read.table(paste0(deconestimatedir,"indonesian_cell_counts_rough_estimate_new.txt"), header=T)
# Problem 1: replicates are uneven between methylation data and RNASeq data, so there is an uneven number of rows.
# Solution: get average of duplicates

# ABIS ----------------------------------------

samplenames <- rownames(ABIS)
samplenames <- sub("\\.","-", samplenames) %>% sub("\\.","-", .)  
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
samplenames <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)
ABIS$samplenames= samplenames
replicate = samplenames[which(duplicated(samplenames))]
replicate = unique(replicate)
newDeconEst = matrix(nrow = length(replicate), ncol = 6)
rownames(newDeconEst) = replicate
colnames(newDeconEst) = colnames(ABIS)[1:6]
for(id in replicate){
	replicatedRowAvg = as.numeric(colMeans(ABIS[which(ABIS$samplenames %in% id),c(1:6)]))
	newDeconEst[id,] = replicatedRowAvg
}
newDeconEst = as.data.frame(newDeconEst)
newDeconEst$samplenames = rownames(newDeconEst)

# remove duplicate rows and replace with avergaes
ABIS = ABIS[-which(ABIS$samplenames %in% replicate),]
ABIS = rbind(newDeconEst, ABIS)
# remove female (MPI-296)
ABIS = ABIS[-grep("MPI-296",ABIS$samplenames),]

# rename columns
colnames(ABIS)[1:6] = c("Gran","Bcell","CD4T","CD8T","NK","Mono")
ABIS$type = "ABIS"

# deconCEll ----------------------------------------

samplenames <- rownames(deconCell)
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
samplenames <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)
deconCell$samplenames= samplenames
replicate = samplenames[which(duplicated(samplenames))]
replicate = unique(replicate)
newDeconEst = matrix(nrow = length(replicate), ncol = 6)
rownames(newDeconEst) = replicate
colnames(newDeconEst) = colnames(deconCell)[1:6]
for(id in replicate){
	replicatedRowAvg = as.numeric(colMeans(deconCell[which(deconCell$samplenames %in% id),c(1:6)]))
	newDeconEst[id,] = replicatedRowAvg
}
newDeconEst = as.data.frame(newDeconEst)
newDeconEst$samplenames = rownames(newDeconEst)

# remove duplicate rows and replace with avergaes
deconCell = deconCell[-which(deconCell$samplenames %in% replicate),]
deconCell = rbind(newDeconEst, deconCell)
# remove female (MPI-296)
deconCell = deconCell[-grep("MPI-296",deconCell$samplenames),]

# rename columns
colnames(deconCell)[1:6] = c("Gran","Bcell","CD4T","CD8T","NK","Mono")
deconCell$type = "deconCell"

# Methylation ----------------------------------------

# remove female (MPI-296)
methylation = methylation[-grep("MPI-296",methylation$ID),]
methylation$type = "Houseman"
methylation[,1:6] = sapply(1:6, function(x) methylation[,x]*100)
# Make violin plots -----------------------------------

methylation$samplenames = methylation$ID
methylation$ID = NULL
df <- melt(rbind(methylation, ABIS, deconCell), id.vars = c("samplenames", "type"))

# set up pvalue matrix
cell.types=colnames(methylation[1:6])
topGenes.pvalue=matrix(nrow=6, ncol=3)
rownames(topGenes.pvalue)=cell.types
colnames(topGenes.pvalue)=c("abisVdeconCell","abisVHouseman","deconCellVHouseman")
for (i in cell.types){
    # get significant genes over a logFC of 1 for all Island comparisons
    cell = df[grep(i, df$variable),]
    topGenes.pvalue[i,"abisVdeconCell"] = t.test(cell[grep("ABIS", cell$type),"value"], cell[grep("deconCell", cell$type),"value"], paired = TRUE, alternative = "two.sided")$p.value
    topGenes.pvalue[i,"abisVHouseman"] = t.test(cell[grep("ABIS", cell$type),"value"], cell[grep("Houseman", cell$type),"value"], paired = TRUE, alternative = "two.sided")$p.value
    topGenes.pvalue[i,"deconCellVHouseman"] = t.test(cell[grep("deconCell", cell$type),"value"], cell[grep("Houseman", cell$type),"value"], paired = TRUE, alternative = "two.sided")$p.value
 }   
 
# make pvalues into scientific notation with max 3 digits
topGenes.pvalue=formatC(topGenes.pvalue, format="e", digits=2, drop0trailing=T)
# convert e notation to base 10 notation
topGenes.pvalue=sub("e", "x10^", topGenes.pvalue)

annotation_df <- data.frame(start=c("ABIS","ABIS", "deconCell"), end=c("deconCell","Houseman","Houseman"), label=topGenes.pvalue["CD8T",])

counter=0
for(cell in cell.types){
    counter=counter+1
    gene.df <- data.frame(df[which(df$variable==cell),])
    #annotation_df <- data.frame(start=c("Sumba","Sumba", "Mentawai"), end=c("Mentawai","West Papua","West Papua"), y=c(max(gene.df[,1]+4),max(gene.df[,1]+5),max(gene.df[,1]+6)), label=paste("limma p-value =",topGenes.pvalue[ensembl,],sep=" "))
    annotation_df <- data.frame(start=c("ABIS","ABIS", "deconCell"), end=c("deconCell","Houseman","Houseman"), y=c(max(gene.df[,"value"]+10),max(gene.df[,"value"]+12),max(gene.df[,"value"]+18)),label=topGenes.pvalue[cell,])
    assign(cell.types[counter], ggviolin(gene.df, x = "type", y = "value", fill="type", add=c("boxplot"),main=cell,add.params = c(list(fill = "white"), list(width=0.05))) + geom_signif(data=annotation_df,aes(xmin=start, xmax=end, annotations=label, y_position=y),textsize = 4, vjust = -0.2,manual=TRUE) + ylim(NA, max(gene.df[,"value"])+20) + scale_fill_brewer(palette="Dark2"))
}

pdf(paste0(outputdir,"deconEstIsland_Unscaled.pdf"), height = 10, width = 15)
ggarrange(CD8T,CD4T,NK,Bcell,Mono,Gran, ncol = 3, nrow = 2)
dev.off()


# 	, main=favGenes[counter], palette=1:3, add.params = c(list(fill = "white"), list(width=0.05))) + geom_signif(data=annotation_df,aes(xmin=start, xmax=end, annotations=label, y_position=y),textsize = 3, vjust = -0.2,manual=TRUE) + ylim(NA, max(gene.df[,1])+7))
# }
# p = ggplot(df, aes(x=variable, y=value, fill=type)) +
#   geom_boxplot()
# # Use brewer color palettes
# p + scale_fill_brewer(palette="Dark2") + theme_classic() + facet_grid(rows = vars(drv))

#write.table(cor(predicted.cellcounts, deconCell, method="spearman"), file=paste0(outputdir,"SpearmanCorrelation_scaledvsNotScaled.txt"))

