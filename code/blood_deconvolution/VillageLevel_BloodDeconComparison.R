# Test out ABIS for deconvoluting blood using RNASeq data from 123 Indonesian samples
# Code developed by Katalina Bobowik, 20.03.2019
# using the shiny app from the paper from Monaco et al, 2019: https://www.cell.com/cell-reports/fulltext/S2211-1247(19)30059-2


# load packages
library(magrittr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(wesanderson)

wes=wes_palette("Zissou1",n=5)
village.col=c("#EBCC2A","#3B9AB2","#F21A00","#78B7C5",wes[4])
village.col2=c(wes[1],wes[2],wes[3],wes[4],wes[5])

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/batchRemoval/"
outputdir= "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/bloodDecon/"
deconestimatedir= "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/indoRNA_SequencingFiles/"

# Read in data
ABIS = read.table(paste0(inputdir, "ABIS/predictedCellCounts_ABIS_scaled.txt"))
# here, I'm using the unscaled values instead of the scaled ones, since Limma does not like all blood porprtions to add up to %100 in the design matrix
# deconCell = read.table(paste0(inputdir, "predictedCellCounts_DeconCell_Paper.txt"))
deconCell = read.table(paste0(inputdir, "predictedCellCounts_DeconCell.txt"))
methylation = read.table(paste0(deconestimatedir,"indonesian_cell_counts_rough_estimate_new.txt"), header=T)

# Problem 1: replicates are uneven between methylation data and RNASeq data, so there is an uneven number of rows.
# Solution: get average of duplicates

# ABIS ----------------------------------------

samplenames <- rownames(ABIS)
samplenames <- sub("\\.","-", samplenames) %>% sub("\\.","-", .)  
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
village <- sapply(strsplit(samplenames, "[-.]"), `[`, 2)
village[grep("_", village)] = "MPI"
other=c("BKB", "HPM","PDT", "PTB", "RIN", "WHB")
village[which(village %in% other)] = "OTHR"
samplenames <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)

ABIS$samplenames= samplenames
ABIS$village = village
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
village.est <- sapply(strsplit(rownames(newDeconEst), "[-.]"), `[`, 2)
village.est[2] = "MPI"
newDeconEst$village = village.est

# remove duplicate rows and replace with avergaes
ABIS = ABIS[-which(ABIS$samplenames %in% replicate),]
ABIS = rbind(newDeconEst, ABIS)
# remove female (MPI-296)
ABIS = ABIS[-grep("MPI-296",ABIS$samplenames),]

# rename columns
colnames(ABIS)[1:6] = c("Gran","Bcell","CD4T","CD8T","NK","Mono")
ABIS$type = "ABIS"
ABIS = ABIS[-which(ABIS$village=="OTHR"),]

# deconCEll ----------------------------------------

samplenames <- rownames(deconCell)
samplenames <- sub("([A-Z]{3})([0-9]{3})", "\\1-\\2", samplenames)
village <- sapply(strsplit(samplenames, "[-.]"), `[`, 2)
village[grep("_", village)] = "MPI"
other=c("BKB", "HPM","PDT", "PTB", "RIN", "WHB")
village[which(village %in% other)] = "OTHR"
samplenames <- sapply(strsplit(samplenames, "[_.]"), `[`, 1)

deconCell$samplenames= samplenames
deconCell$village = village

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
village.est <- sapply(strsplit(rownames(newDeconEst), "[-.]"), `[`, 2)
village.est[2] = "MPI"
newDeconEst$village = village.est

# remove duplicate rows and replace with avergaes
deconCell = deconCell[-which(deconCell$samplenames %in% replicate),]
deconCell = rbind(newDeconEst, deconCell)
# remove female (MPI-296)
deconCell = deconCell[-grep("MPI-296",deconCell$samplenames),]

# rename columns
#colnames(deconCell)[1:6] = c("Gran","Bcell","CD4T","CD8T","NK","Mono")
colnames(deconCell)[1:6]=c("Granulocytes","B cells (CD19+)","CD4+ T cells","CD8+ T cells","NK cells (CD3- CD56+)","Monocytes (CD14+)")
deconCell$type = "deconCell"
deconCell = deconCell[-which(deconCell$village=="OTHR"),]

# Methylation ----------------------------------------

# remove female (MPI-296)
methylation = methylation[-grep("MPI-296",methylation$ID),]
mpi.id = which(lengths(strsplit(as.character(methylation$ID), "[-.]")) == 2)

village <- sapply(strsplit(as.character(methylation$ID), "[-.]"), `[`, 2)
village[mpi.id] = "MPI"
other=c("BKB", "HPM","PDT", "PTB", "RIN", "WHB")
village[which(village %in% other)] = "OTHR"
methylation$village = village
methylation$type = "Houseman"
methylation = methylation[-which(methylation$village=="OTHR"),]
methylation$samplenames = methylation$ID
methylation$ID = NULL

# ANOVA -----------------------------------

df <- as.data.frame(deconCell)

# Pst hoc ANOVA
cellTypeSummary = matrix(nrow=6, ncol=2)
colnames(cellTypeSummary) = c("Cell.type","ANOVA.Pvalue")
counter=0
for(celltype in colnames(df)[1:6]){
	counter=counter+1
	cell.summary=summary(aov(df[,celltype] ~ df$village))
	cellTypeSummary[counter,1]=celltype
	cellTypeSummary[counter,2]=cell.summary[[1]][["Pr(>F)"]][1]
	TUKEY <- TukeyHSD(x=aov(df[,celltype] ~ df$village), 'df$village', conf.level=0.95)
	write.table(TUKEY[[1]],file=paste0(outputdir,celltype,"TUKEYtest_unscaled.txt"))
	# write out ANOVA results
	write.table(data.frame(summary(aov(df[,celltype] ~ df$village))[[1]]), file=paste0(outputdir,celltype,"ANOVAtest_unscaled.txt"), quote=F, sep="\t", col.names=T)
}
write.table(cellTypeSummary, file=paste0(outputdir,"cellTypeSummaryANOVA_unscaled.txt"))

# Violin plots -----------------------------------

df2 = as.data.frame(deconCell[c(1:6,8)])
df2 = melt(df2)
df2$village = gsub("MPI","KOR", df2$village)
df2$village = factor(df2$village, levels = c("TLL", "MDB", "ANK", "WNG","KOR"))
CellType.adjustedPval = data.frame(cellTypeSummary[,1], p.adjust(cellTypeSummary[,2],  method="BH"))
CellType.adjustedPval = formatC(CellType.adjustedPval$p.adjust.cellTypeSummary...2...method....BH..,format = "e", digits = 2)
p.adjust(cellTypeSummary[,2],  method="BH")

# make labeller for facet titles
counter=0
cell_names <- list()
for (name in levels(df2$variable)) {
	counter=counter+1
	cell_names[[name]]=paste0(name," (",CellType.adjustedPval[counter],")")
}

cell_labeller <- function(variable,value){
  return(cell_names[value])
}

# plot
pdf(paste0(outputdir,"deconEstVillage_Unscaled.pdf"), height = 10, width = 14)
#ggviolin(df2, x = "village", y = "value", fill="village", palette = village.col2, add=c("boxplot"),add.params = c(list(fill = "white"), list(width=0.05))) + facet_grid(rows = vars(variable), scales = "free") + ylab("Cell type proportion")
ggviolin(df2, x = "village", y = "value", fill="village", palette = village.col2, add=c("boxplot"),add.params = c(list(fill = "white"), list(width=0.05))) + facet_wrap(~variable, scales = "free",labeller=cell_labeller) + ylab("Cell type proportion")
dev.off()
