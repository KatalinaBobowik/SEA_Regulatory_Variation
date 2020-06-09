# Script for creating barplot of number of reads at all different processing stages
# Code developed by Katalina Bobowik, 10.12.2018

# load library and colour palette
library(viridis)
library(Rcmdr)
library(DescTools)

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Paper_Figures")

# read in summary file and tidy up 
a=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/Thesis_Writing/Tables/Table2_TotalnReads_AllSamples.txt", header=T, sep="\t")
a=a[1:123,1:5]
rownames(a)=make.unique(as.character(a[,1]))
a$Sample = NULL

# divide all rows by the first column, which is the total number of reads
b=t(apply(a[,1:4], 1, function(x){x/x[1]}))
# now subtract each column percentage from the previous column in order to get the cumulative number of reads
for(column in c(1:3)){
	b[,column]=b[,column]-b[,column+1]
}

# plot both n reads and percentages onto one pdf
pdf("TotalandPercentageOfReads_allFilteringStages.pdf", height=15,width=18)
par(mar=c(6.1,4.1,10.1,2.1), xpd=T, mfrow=c(2,1))
# first plot total number of reads
barplot(as.matrix(t(a)*1e-6), col=viridis(4), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, border="white", las=3, ylim=c(0,200), main="Number of Reads \nAll filtering stages", ylab="Number of Reads (millions)")
legend(120,240, col=viridis(4), legend=colnames(a), pch=15, cex=1.5)
# now plot reads as a percentage
barplot(Rev(t(b),margin=1),col=Rev(viridis(4)),cex.axis=1.5, cex.main=1.5, cex.lab=1.5,border="white", las=3,main="Percentage of Reads \nAll filtering stages", ylab="Percentage of Reads")
#legend(130,1.3, col=viridis(4), legend=colnames(a), pch=15, cex=0.9)
dev.off()

# plot just n reads
pdf("TotalReads_allFilteringStages.pdf", height=10,width=18)
par(mar=c(8.1,6.1,10.1,2.1), xpd=T)
barplot(as.matrix(t(a)*1e-6), col=viridis(4), cex.main=2, cex.lab=2, cex.axis=2, border="white", las=3, ylim=c(0,200), main="Number of Reads \nAll filtering stages", ylab="Number of Reads (millions)")
legend(120,240, col=viridis(4), legend=colnames(a), pch=15, cex=1.5)
dev.off()

