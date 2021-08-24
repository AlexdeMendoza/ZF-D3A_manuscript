###############################################################
library(stringr)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(tidyr)
library(dplyr)
library(cowplot)
library(ChIPpeakAnno)
library(GenomicFeatures)

# ChIP-seq fastq were adapter and quality trimmed with fastp
#fastp -i "$R1" -I "$R2" -o "$prefix"_trim.fastq.gz -O "$prefix2"_trim.fastq.gz

# Mapped using bowtie
#(bowtie2 -q --threads "$cores" \
#  -x "$index" -X 2000 \
#  -1 "$prefix"_trim.fastq.gz -2 "$prefix2"_trim.fastq.gz | \
#  samtools view -bSu - | samtools sort -T "$prefix" - > "$prefix".bam) 2> "$prefix".log

# duplicate remove with sambamba
#while read line; do
# out=$(echo $line |perl -pe "s/\.bam//")
# sambamba markdup -t 20 -r $line ${out}.dedup.bam
#done < bams

## Remove blacklist regions
#blacklist="hg19_blacklist_regions_ENCFF001TDO.bed"

# Remove MAPQ < 30 and filter blacklist regions
#parallel samtools view -L "$blacklist" -U {.}.filtered.bam -b {} '>' {.}.blacklist.bam ::: *dedup.bam

# Peaks called with MACS2 as in 
#macs2 callpeak -t ChIP_rep1.dedup.filtered.bam -c Input.dedup.filtered.bam -f BAMPE -q 0.05 --down-sample -g hs -n ChIP_rep1

# Duplicates merged with IDR


#Import ChIP-seq peaks of ZF
read_bed_to_GRobject <- function(bedfile){
  dat <- fread(bedfile)
  gr <- GRanges(seqnames = Rle(dat$V1),
                ranges = IRanges(start = dat$V2, end = dat$V3))
  dat <- dat[,!c(1:3)]
  mcols(gr) <- dat
  
  return(gr)
}

zf_mut_gr <- read_bed_to_GRobject("ZF598_D3amut_FEER_P2A_GFP_HA_ChIP.nofilt.idr_peaks")
zf_wt_gr <- read_bed_to_GRobject("ZF598_D3awt_P2A_GFP_HA_ChIP.nofilt.idr_peaks")

# Import features from UCSC annotation
genePath <- "genes.gtf"
txdb <- makeTxDbFromGFF(file = genePath)
#promoters extracted by transcript, not gene, to take into account alternative promoters
promoters <- promoters(transcripts, upstream = 2000, downstream = 200)

# ZF-peak in promoter:
genes_with_ZFmut_promoter <- promoters[overlapsAny( promoters, zf_mut_gr)] %>% names()

# res_mut object from RNA-seq comparisons of DEseq2
# genes diff expressed in mutant that have a ZF peak in the promoter
res_mut %>% data.frame() %>% dplyr::filter(padj < 0.05) %>% .$gene_id %in% genes_with_ZFmut_promoter %>% table()
#FALSE  TRUE 
#4739   575 

# Plots ZF overlaps with promoters compared to noDox-mut vs Dox-mut differential expression

res_mut$peak_in_promoter <- ifelse(res_mut$gene_id %in% genes_with_ZFmut_promoter, "ZF_peak_in_promoter", "No_peak")
res_mut$direction <- ifelse(res_mut$log2FoldChange > 0, "upregulated", "downregulated")

ggplot(data.frame(res_mut), aes(y = -log10(padj), x = log2FoldChange, colour = peak_in_promoter)) + geom_point()

plot_proportion_df <- data.frame(res_mut) %>% filter(padj < 0.05) %>% dplyr::select(peak_in_promoter, direction) %>% table() %>% data.frame()

gg_peaks_promoters <- ggplot(data.frame(res_mut) %>% filter(padj < 0.05), aes(x = log2FoldChange,fill = direction) ) + geom_histogram() + theme_light() +
  facet_grid(peak_in_promoter~.) 
ggsave(gg_peaks_promoters, filename = "DEG_with_ZF_peak_promoter_mutant_hist.pdf")

gg_peaks_promoters <- ggplot(plot_proportion_df, aes(x = peak_in_promoter, y = Freq, fill = direction, label = Freq)) +
  geom_bar(position="stack", stat="identity") +
  geom_text(size = 2, color = "white", position = position_stack(vjust = 0.5)) +
  ggtitle("Differentially expressed genes (FDR < 0.05)\nClassified by having a ZF on its promoter") +
  theme_light()

ggsave(gg_peaks_promoters, filename = "DEG_with_ZF_peak_promoter_mutantTrue.pdf", width = 4, height = 3)



