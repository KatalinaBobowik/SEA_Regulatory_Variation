# script created by KSB, 08.08.18
# Perform DE analysing relationship between islands

# load dependencies: libraries, human count data, and data preprocessing

library(edgeR)
library(plyr)
library(NineteenEightyR)
library(RColorBrewer)
library(biomaRt)
library(ggplot2)
library(ggsignif)
library(EGSEA)
library(goseq)
library(ReactomePA)

# Set paths:
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/"
revigodir= "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/REVIGO/"

# Set output directory and create it if it does not exist:
outputdir <- "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Village/GeneSetAnalysis/"

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir)
}

# Load colour schemes:
wes=c("#3B9AB2", "#EBCC2A", "#F21A00", "#00A08A", "#ABDDDE", "#000000", "#FD6467","#5B1A18")
palette(c(wes, brewer.pal(8,"Dark2")))
# set up colour palette for batch
batch.col=electronic_night(n=3)
dev.off()

# load files -------------------------------------------------------------------------------------------------

# load y, the normalised and filtered counts object
load(paste0(inputdir, "dataPreprocessing/indoRNA.read_counts.TMM.filtered.Rda"))
# load in the efit object
load(paste0(inputdir, "DE_Village/voomDupEfit_Village.Rda"))
# load in the voom object
load(paste0(inputdir, "DE_Village/vDup_Village.Rda"))

# set up design matrix ----------------------------------------------------------------------------------------

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

# Enrichment analysis for Gene Ontology ----------------------------------------------------------------------------------------------------

ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'www.ensembl.org', ensemblRedirect = FALSE)

# transform ensembl IDs to entrez IDs to be compatible with human c2 dataset (below)
ensembl_genes=rownames(y)
entrez=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "description"),values= ensembl_genes,mart=ensembl.mart.90)
y$entrezID=entrez[match(rownames(y), entrez[,1]), 2]

# gene set testing with Camera
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata")) 
idx <- ids2indices(Hs.c2,id=y$entrezID) 
for (i in 1:ncol(contr.matrix)){
    camera.matrix=camera(vDup,idx,design,contrast=contr.matrix[,i])
    write.table(camera.matrix, file=paste0(outputdir,"cameraMatrix_",colnames(contr.matrix)[i],".txt"))
}

# gene set testing with goSeq
for(pop in 1:ncol(voomDupEfit)){
    for(pval in c(0.05, 0.01)){
        topTable <- topTable(voomDupEfit, coef=pop, n=Inf, p.value=pval, lfc=0.5)
        gene.vector=as.integer(rownames(y) %in% rownames(topTable))
        names(gene.vector)=rownames(y)

        # set the probability weighting fcuntion, i.e., implement a weight for each gene dependent on its length
        pwf=nullp(gene.vector,"hg19","ensGene")
        # use  the  default  method  to  calculate  the  over  and  under  expressed  GO categories among DE genes
        GO.wall=goseq(pwf,"hg19","ensGene",use_genes_without_cat=TRUE)

        # now let's interpret the results. First we need to apply a multiple hypothesis testing correction set at 5% (BH method)
        enriched.GO=GO.wall[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05,]
        write.table(enriched.GO, file=paste0(outputdir,"enrichedGOterms_lfc0.5_",colnames(voomDupEfit)[pop],pval,".txt"), quote=F, row.names=F)
        if(nrow(enriched.GO)>0){
            zz=file(paste0(outputdir,"topTen_enrichedGOterms_",pop,pval,".txt"), open="wt")
            sink(zz)
            sink(zz, type = "message")
            # get GO terms for top ten enriched GO terms - write output with the sink() function
            for(go in 1:length(enriched.GO$category)){
                print(GOTERM[[enriched.GO$category[go]]])
            }
        }
        sink(type = "message")
        sink()

        # KEGG pathway analysis
        en2eg=as.list(org.Hs.egENSEMBL2EG)
        # Get the mapping from Entrez 2 KEGG
        eg2kegg=as.list(org.Hs.egPATH)
        # Define a function which gets all unique KEGG IDs
        # associated with a set of Entrez IDs
        grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
        # Apply this function to every entry in the mapping from
        # ENSEMBL 2 Entrez to combine the two maps
        kegg=lapply(en2eg,grepKEGG,eg2kegg)

        # produce PWF as before
        pwf=nullp(gene.vector,"hg19","ensGene")
        KEGG=goseq(pwf,gene2cat=kegg, use_genes_without_cat=TRUE)
        enriched.GO.kegg=KEGG[p.adjust(KEGG$over_represented_pvalue, method="BH")<.05,]
        write.table(enriched.GO.kegg, file=paste0(outputdir,"enrichedGOkegg_lfc0.5_",pop,pval,".txt"))
    }
}

# Enrichment Visualisation ----------------------------------------------------------------------------------------------------

# Having played around with many visualisation options, it seems as though using Revigo is the best and most user-friendly. Using the output from the GoSeq enriched Go results (FDR <0.05), we can plug our GO IDs into Revigo (http://revigo.irb.hr/) and then output this as an R script.
# First save  output (remember that SMBvsMTW has no significantly enriched GO pathways so we'll only do this for SMBvsMPI and MTWvsMPI)  

# Revigo for MTW vs MPI (pval of 0.01 and LFC of 0.05 FDR)
#source(paste0(revigodir,"REVIGO_MTWvsMPI_LFC05_Pval01.r"))
# SMB vs MPI
#source(paste0(revigodir,"REVIGO_SMBvsMPI_LFC05_Pval01.r"))


# EGSEA -----------------------------------------------------------------------

# transform ensembl IDs to entrez IDs to be compatible with human c2 dataset (below)
ensembl_genes=rownames(y)
entrez=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "description"),values= ensembl_genes,mart=ensembl.mart.90)
vDup$genes$entrezID=entrez[match(rownames(y), entrez[,1]), 2]

#subset voom object to contain only genes with an Entrez Gene ID
vDup <- vDup[which(!is.na(vDup$genes$entrezID)),]
vDup <- vDup[-c(which(duplicated(vDup$genes$entrezID))),]
rownames(vDup) <- vDup$genes$entrezID
rownames(vDup$genes) <- vDup$genes$entrezID

#build gene-set collections to compute enrichment for
gs.annots = buildIdx(entrezIDs=vDup$genes$entrezID, species="human",msigdb.gsets=c("c2"), go.part = TRUE)

#generate a symbolsMap
symbolsMap = vDup$genes[, c(4, 2)]
colnames(symbolsMap) = c("FeatureID", "Symbols")
symbolsMap[, "Symbols"] = as.character(symbolsMap[, "Symbols"])

#establish the list of independent gene-set enrichment methods for EGSEA to use
baseMethods = egsea.base()[-2]

# Ensemble testing with EGSEA
gsa = egsea(voom.results=vDup, contrasts=contr.matrix, gs.annots=gs.annots, symbolsMap=symbolsMap, baseGSEAs=baseMethods, sort.by="med.rank", num.threads = 8, report = FALSE)

# get a summary of the top gene sets
plotSummaryHeatmap(gsa, gs.label="c2", hm.vals = "avg.logfc.dir",file.name=paste0(outputdir,"summary_heatmaps_c2"), sort.by="med.rank")
plotSummaryHeatmap(gsa, gs.label="kegg", hm.vals = "avg.logfc.dir",file.name=paste0(outputdir,"summary_heatmaps_kegg"), sort.by="med.rank")

# get top gene sets
topc2=topSets(gsa, contrast = "comparison",sort.by="med.rank", gs.label="c2", number=20)
topkegg=topSets(gsa, contrast = "comparison",sort.by="med.rank", gs.label="kegg")

# let's make a heatmap of all of the top enriched pathways in the c2 dataset
for (i in 1:20){
    plotHeatmap(gsa, gene.set=topc2[i], gs.label="c2", contrast = "comparison", file.name = paste0(outputdir,"c2_ComparisonHeatmap_",topc2[i]))
    plotPathway(gsa, gene.set = topkegg[i],contrast = "comparison", gs.label = "kegg",file.name = paste0(outputdir,"kegg_ComparisonPathway_",topkegg[i]))
}

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

