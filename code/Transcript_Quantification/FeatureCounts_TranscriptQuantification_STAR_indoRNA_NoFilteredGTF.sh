# created by KSB, 09.13.2017

# quantify the abundance of transcripts from the STAR alignment of Indonesian reads using FeatureCounts
# The Subread package was downlaoded from Sourceforge for Linux, version 1.5.3 (subread-1.5.3-Linux-x86_64.tar.gz)

#download GTF file from Ensembl and unzip
wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz
gunzip /vlsci/SG0008/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf.gz
 
# load modules
module load samtools/1.9

# filter the bam file so that only uniquely mapped reads remain
inputdir=/data/scratch/projects/punim0586/kat/101BPIndonesianSamplesMapped
outdir=/data/scratch/projects/punim0586/kat/UniquelyMappedIndo101BP
for file in ${inputdir}/*.bam; do
    filename=`basename $file`
    # "-q 255" = unique reads; -b output bam
    echo samtools view -b -q 255 $file '>' ${outdir}/uniquelyMapped_${filename}
done > /data/cephfs/punim0586/kbobowik/Sumba/FCNoFiltering.txt

### FeatureCounts pipeline -------------------------

# execute FeatureCounts 
# here are the following flags I used:
# -T: Number of the threads.  The value should be between 1 and 32.  1 by default.
# -s: Indicate if strand-specific read counting should be performed. Acceptable  values:  0  (unstranded),  1  (stranded)  and  2  (re-versely stranded).  0 by default.  For paired-end reads, strand of the first read is taken as the strand of the whole fragment. FLAG field is used to tell if a read is first or second read in a pair.
# -p: If specified, fragments (or templates) will be counted instead of reads.  This option is only applicable for paired-end reads.
# -a: Give the name of an annotation file
# - t: Specify the feature type.  Only rows which have the matched feature type in the provided GTF annotation file will be included for read counting.  ‘exon’ by default.
# -o: Give the name of the output file
# -g: Specify the attribute type used to group features (eg.  exons) into meta-features (eg.  genes)

# for uniquely mapped reads ------------------------------

# load required modules
module load subread
gtf=/data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf
inputdir=/data/scratch/projects/punim0586/kat/UniquelyMappedIndo101BP
outdir=/data/scratch/projects/punim0586/kat/UniquelyMappedIndo101BP_FC

for file in ${inputdir}/uniquelyMapped*.bam; do
  sample=`basename $file _Aligned.sortedByCoord.out.bam`
  ID=${sample##*_}
  featureCounts -T 6 -s 2 -p -a $gtf -t exon -g gene_id -o ${outdir}/${sample}.txt $file
  cat ${outdir}/${sample}.txt | awk '{ print $1 "\t" $6 "\t" $7 }' | tail -n +2 | sed -e "1s~${file}~Count~" > ${outdir}/Filter_GeneLengthCount_${sample}.txt
done

# for multi-mapping mapped reads ------------------------------

# load required modules
gtf=/data/cephfs/punim0586/shared/genomes/hg38/GTF_annotation/Homo_sapiens.GRCh38.90.gtf
inputdir=/data/scratch/projects/punim0586/kat/101BPIndonesianSamplesMapped
outdir=/data/scratch/projects/punim0586/kat/MultiMappingIndo1o1BP_FC

for file in ${inputdir}/*.bam; do
  sample=`basename $file _Aligned.sortedByCoord.out.bam`
  ID=${sample##*_}
  echo featureCounts -T 6 -s 2 -p -a $gtf -t exon -g gene_id -o ${outdir}/${sample}.txt $file
  cat ${outdir}/${sample}.txt | awk '{ print $1 "\t" $6 "\t" $7 }' | tail -n +2 | sed -e "1s~${file}~Count~" > ${outdir}/Filter_GeneLengthCount_${sample}.txt
done > /data/cephfs/punim0586/kbobowik/Sumba/FCNoFilteredGTF.txt

