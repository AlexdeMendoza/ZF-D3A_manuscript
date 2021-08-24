library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)


# RNA-seq libraries were processed using:
# fastp for trimming the reads
#fastp -i "$readR1" -I "$readR2" -o "$baseName"_trimmed_R1.fastq.gz -O "$baseName"_trimmed_R2.fastq.gz

# Align the reads using hisat2
#hisat2 --time \
#--rna-strandness RF \
#--threads "$cores" \
#-x "$indexPath" \
#-1 "$baseName"_trimmed_R1.fastq.gz -2 "$baseName"_trimmed_R2.fastq.gz \
#-S "$baseName".sam

# convert to bam, sort and index alignment
#samtools view -bSu "$baseName".sam | samtools sort -T "$baseName"_sorted - > "$baseName".bam

# generate bigwigs for the IGV browser using DeepTools2
# while read bam ; do bamCoverage -b $bam -o ${bam%%.bam}.bw -p 20 -bs 1 --normalizeUsing CPM -bl ERCC92.bed; done < list_of_merged_bam

# find the path to the bam files
bamFiles <- list.files(pattern = ".bam$", full.names = TRUE)

# find the path to the UCSC gene annotation of hg19 and the ERCC spike-in annotation
genePath <- "genes.gtf"
erccPath <- "ERCC92.gtf"


txdb <- makeTxDbFromGFF(file = genePath)
genes <- exonsBy(x = txdb, by = "gene")



ercc_txdb <- makeTxDbFromGFF(file = erccPath)
ercc <- exonsBy(x = ercc_txdb, by = "gene")

all_genes <- c(genes, ercc)

countTranscriptOverlaps <- function(bamFiles, grl=all_genes,
                                    mode="IntersectionNotEmpty"){
  
  # Specify the bam files and load into list
  bams <- BamFileList(bamFiles)
  
  # Load the transcripts data and format for summarizeOverlaps
  
  # Count overlaps
  olaps <- summarizeOverlaps(features = grl,
                             mode=mode,
                             singleEnd=FALSE,
                             reads = bams,
                             ignore.strand = TRUE,
                             BPPARAM=SerialParam())
  
  counts <- assays(olaps)$counts
  
  return(counts)
}

counts <- countTranscriptOverlaps(bamFiles = bamFiles, grl = all_genes,
                                  outFile = outFile)


colnames(counts) <- colnames(counts) %>% str_replace(. , pattern=".merge.bam", replace = "")

write.table(counts, file = "ZF_Counts_table.tsv", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
