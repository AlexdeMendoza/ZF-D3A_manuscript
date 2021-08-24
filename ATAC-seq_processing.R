library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(stringr)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(tidyr)
library(dplyr)
library(cowplot)
library(ChIPpeakAnno)

# Trim the nextera adapters from the fastqs with bbduk
#sh bbduk2.sh \
#in="$R1" \
#in2="$R2" \
#out="$prefix".tmp_R1.fq \
#out2="$prefix".tmp_R2.fq \
#rliteral=GCGATCGAGGACGGCAGATGTGTATAAGAGACAG,CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG \
#ktrim=r \
#mink=3 \
#threads="$cores" \
#overwrite=true

## bowtie2 alignment
#(bowtie2 -q --threads "$cores" -X2000 \
#  -x "$index" \
#  -1 "$prefix".tmp_R1.fq -2 "$prefix".tmp_R2.fq -S "$prefix".sam) 2> "$prefix".log

## Convert to BAM file
#sambamba view -S -f bam "$prefix".sam > "$prefix".bam
#sambamba sort -o "$prefix"_sorted.bam "$prefix".bam

# Create the bam file index
#sambamba index "$prefix"_sorted.bam

#Filter the chrM data
#samtools idxstats "$prefix"_sorted.bam | cut -f 1 | grep -v chrM | xargs samtools view -b "$prefix"_sorted.bam > "$prefix".bam

# Get the insert size metrics to generate plots
#samtools view -f66 "$prefix".bam | cut -f 9 | sed 's/^-//' > "$prefix".InsertSizeMetrics.txt

# Remove PCR duplicates
#parallel -j4 sambamba markdup --overflow-list-size 600000 --hash-table-size 2000000 -t 8 -r -p {} {.}_dedup.bam ::: *.merge.bam

# Remove blacklisted regions
#blacklisted_regions="hg19_blacklist_regions_ENCFF001TDO.bed"
#parallel -j8 samtools view -L $blacklisted_regions -U {.}.filtered.bam -b {} ::: *_dedup.bam

# shift insertion sites with DeepTools2 alignmentSieve with --ATACshift option

# Call peaks
#macs2 callpeak --nomodel -t "$prefix".bam -f BAM -n "$prefix" --keep-dup all --gsize hs

# Merge peaks from replicates with IDR

# Merge all peaks with bedtools merge
#cat *narrowPeaks | cut -f1,2,3 | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > ATACall_peaks_merged.bed

# Import all peaks to R
read_bed_to_GRobject <- function(bedfile){
  dat <- fread(bedfile)
  gr <- GRanges(seqnames = Rle(dat$V1),
                ranges = IRanges(start = dat$V2, end = dat$V3))
  #dat <- dat[,!c(1:3)]
  #mcols(gr) <- dat
  
  return(gr)
}

atac_gr <- read_bed_to_GRobject("ATACall_peaks_merge.bed")

# then take counts of ATAC peaks

bamFiles <- list.files(path = "/ZF598/ZF_ATACseq",
                       pattern = ".shifted.sorted.bam$", full.names = TRUE) %>% BamFileList()


# Count overlaps

olaps <- summarizeOverlaps(features = atac_gr,
                           reads = bamFiles,
                           ignore.strand = TRUE,
                           BPPARAM=SerialParam())

counts <- assays(olaps)$counts
rownames(counts) <- olaps@rowRanges %>% add_loci() %>% .$loci

write.table(counts, file = "ATACall_peaks.merge.readCounts.tsv", quote = F, sep = "\t")

################################
# Import ATAC-seq peaks and normalise using DEseq2 and compute differential ATAC peaks across replicates

atac_counts <- fread("ATACall_peaks.merge.readCounts.tsv.gz") %>% data.frame()
rownames(atac_counts) <- atac_counts$V1
atac_counts$V1 <- NULL
colnames(atac_counts) <- colnames(atac_counts) %>% str_remove(., pattern = ".merge_dedup.filtered.bam.shifted.sorted.bam")


coldata_atac <- data.frame(sample=colnames(atac_counts)) %>%
  mutate(condition = ifelse(grepl("D3awt_noDox", sample), "noDox_wt",
                            ifelse(grepl("D3amut", sample), "Dox_mutant", 
                                   ifelse(grepl("Dox_withdrawal", sample), "DoxWD_wt","Dox_wt")))) 

coldata_atac$condition <- factor(coldata_atac$condition, levels = c("noDox_wt","Dox_wt","DoxWD_wt","Dox_mutant"))

dds_atac <- DESeqDataSetFromMatrix(countData = atac_counts,
                                   colData = coldata_atac,
                                   design = ~ condition)

dds_atac <- estimateSizeFactors(dds_atac)

atac_norm_summary <- data.frame(loci = rownames(counts(dds_atac)), noDox_atac = rowMeans(counts(dds_atac, normalize = TRUE)[,dds_atac$condition == "noDox_wt"]),
                                Dox_atac = rowMeans(counts(dds_atac, normalize = TRUE)[,dds_atac$condition == "Dox_wt"]),
                                DoxWD_atac = rowMeans(counts(dds_atac, normalize = TRUE)[,dds_atac$condition == "DoxWD_wt"]),
                                Dox_mutant = rowMeans(counts(dds_atac, normalize = TRUE)[,dds_atac$condition == "Dox_mutant"]))

# Diff test on noDox vs Dox 
dds_atac_Dox <- DESeqDataSetFromMatrix(countData = atac_counts[,c("RL1919_ZF598_D3awt_noDox_ATAC_rep1","RL1920_ZF598_D3awt_ATAC_rep1","RL1923_ZF598_D3awt_noDox_ATAC_rep2","RL1924_ZF598_D3awt_Dox_ATAC_rep2")],
                                       colData = coldata_atac[grepl(coldata_atac$condition, pattern = "Dox_wt"),],
                                       design = ~ condition)
dds_atac_mutant <- DESeqDataSetFromMatrix(countData = atac_counts[,c("RL1919_ZF598_D3awt_noDox_ATAC_rep1","RL1923_ZF598_D3awt_noDox_ATAC_rep2","RL1922_ZF598_D3amut_Dox_ATAC_rep1","RL1926_ZF598_D3amut_DoxATAC_rep2")],
                                          colData = coldata_atac[grepl(coldata_atac$condition, pattern = "noDox_wt|Dox_mutant"),],
                                          design = ~ condition)
dds_atac_WD <- DESeqDataSetFromMatrix(countData = atac_counts[,c("RL1919_ZF598_D3awt_noDox_ATAC_rep1","RL1923_ZF598_D3awt_noDox_ATAC_rep2","RL1921_ZF598_D3awt_Dox_withdrawal_ATAC_rep1","RL1925_ZF598_D3awt_Dox_withdrawal_ATAC_rep2")],
                                      colData = coldata_atac[grepl(coldata_atac$condition, pattern = "noDox_wt|DoxWD_wt"),],
                                      design = ~ condition)



filter_and_perform_DEseq_on_peaks <- function(dds_obj, ctrl = "noDox_wt", trt = "Dox_wt"){
  dds_obj <- dds_obj[rowSums(counts(dds_obj)) >= 10,]
  dds_obj <- dds_obj[rowSums(counts(dds_obj)==0) <= ncol(counts(dds_obj))*0.5,]
  
  
  dds_obj <- estimateSizeFactors(dds_obj)
  
  rowData(dds_obj)$control_expr <- rowMeans(counts(dds_obj, normalize = TRUE)[,dds_obj$condition ==  ctrl ])
  rowData(dds_obj)$condition_expr <- rowMeans(counts(dds_obj, normalize = TRUE)[,dds_obj$condition == trt ])
  
  dds_obj <- DESeq(dds_obj)
  
  return(dds_obj)
}


dds_atac_Dox <- filter_and_perform_DEseq_on_peaks(dds_atac_Dox)
dds_atac_mutant <- filter_and_perform_DEseq_on_peaks(dds_atac_mutant, trt = "Dox_mutant")
dds_atac_WD <- filter_and_perform_DEseq_on_peaks(dds_atac_WD, trt = "DoxWD_wt")

obtain_results_print_stats <- function(dds_obj){
  res <- results(dds_obj)
  
  message(paste("DE ATAC peaks",sum(res$padj < fdrlevel.de, na.rm = TRUE)))
  
  df <- data.frame(res) %>% mutate(change_direction = ifelse(log2FoldChange < 0, "downregulated", "upregulated")) %>%
    filter(padj < fdrlevel.de) %>% .$change_direction %>% table()
  message(paste("Downregulated", df[1]))
  message(paste("Upregulated",df[2]))
  res$control_expr <- rowData(dds_obj)$control_expr
  res$condition_expr <- rowData(dds_obj)$condition_expr
  
  
  res <- res %>% na.omit()
  
  return(res)
}

res_atac_Dox <- obtain_results_print_stats(dds_atac_Dox)
#27270 DE ATAC-peaks
#downregulated   upregulated
#12770          14500
res_atac_mutant <- obtain_results_print_stats(dds_atac_mutant)
# 1 single peak!!!
res_atac_WD <- obtain_results_print_stats(dds_atac_WD)


ggMAplot_ATAC <- function(deseq_res, lab = "(WT noDox/Dox)"){
  deseq_res$gene_id <- rownames(deseq_res)
  ggplot(data.frame(deseq_res), aes(x = baseMean, y = log2FoldChange)) + 
    geom_point(size = 0.5, alpha = 0.5, mapping = aes(color = padj < 0.01)) +
    ylab(paste0("log Fold Change", lab)) +
    xlab("mean normalised ATAC signal") +
    geom_hline(yintercept=0, col="grey") +
    geom_hline(yintercept=1, col="grey", linetype="dashed") +
    geom_hline(yintercept=-1, col="grey", linetype="dashed") +
    annotate(geom="text", x=1e+03, y=-4, label=paste("Downregulated peaks =",nrow(deseq_res[deseq_res$log2FoldChange < 0 & deseq_res$padj < 0.01,])),color="black") +
    annotate(geom="text", x=1e+03, y=4, label=paste("Upregulated peaks =",nrow(deseq_res[deseq_res$log2FoldChange > 0 & deseq_res$padj < 0.01,])),color="black") +
    theme_minimal() + scale_x_log10()
}

ggMA_ATAC_noDox_Dox <- ggMAplot_ATAC(deseq_res = res_atac_Dox, lab = "(WT Dox vs NoDox)")
ggMA_ATAC_noDox_Mutant <- ggMAplot_ATAC(deseq_res = res_atac_mutant, lab = "(WT Mutant Dox vs NoDox)")
ggMA_ATAC_noDox_WD <- ggMAplot_ATAC(deseq_res = res_atac_WD, lab = "(WT Dox withdrawal vs noDox)")

ggsave(plot_grid(ggMA_ATAC_noDox_Dox, ggMA_ATAC_noDox_Mutant, ggMA_ATAC_noDox_WD, ncol = 1), 
       filename = "MAplots-ATAC.png", height = 10, width = 6)

