library(bsseq)
library(dmrseq)
library(dplyr)
library(stringr)
library(BiocParallel)
library(data.table)
BiocParallel::register(BiocParallel::SerialParam())

#### WGBS reads were trimmed using BBduk
#bbduk.sh \
#in="$readR1" in2="$readR2" \
#out="$baseR1"_trimmed.fastq.gz out2="$baseR2"_trimmed.fastq.gz \
#literal="$known_adapter" threads="$cores" overwrite=t -Xmx20g \
#overwrite=true &&

#### Merge overlapping read pairs
#bbmerge.sh in1="$baseR1"_trimmed.fastq.gz in2="$baseR2"_trimmed.fastq.gz qtrim=r \
#out="$baseR1"_merged.fastq.gz outu1="$baseR1"_unmerged.fastq.gz outu2="$baseR2"_unmerged.fastq.gz \
#-Xmx20g

#### Paired-end alignment with BS-Seeker2
#cat "$baseR1"_map_manifest | parallel --colsep="\t" -j"$files" bs_seeker2-align.py \
#--aligner=bowtie2 --bt2--end-to-end --bt2-p "$bt_cores" -e 300 -X 2000 \
#--temp_dir=/scratchfs/amendoza/tmp \
#-1 {1} -2 {2} -o {3} \
#-d "$indexBS" \
#-g "$genome"

#### Sort the output bam files
#while read i; do sambamba sort -t "$cores" $i; done < "$baseR1"_out &&
  
#### Clip the overlaps in the sorted BAM files using bamUtil
#parallel -j"$files" bam clipOverlap --in {} --out {.}_clipped.bam ::: 0{1..6}."$baseR1"_unmerged_split.sorted.bam &&

#### Merge the bam files
#sambamba merge -t "$cores" "$baseR1"_pairs.bam 0{1..6}."$baseR1"_unmerged_split.sorted_clipped.bam &&
  
#### Map the merged pairs
#paste "$baseR1"_merged_fq "$baseR1"_merged_read_bams | parallel -j"$files" --colsep="\t" bs_seeker2-align.py \
#--aligner=bowtie2 --bt2--end-to-end --bt2-p "$bt_cores" -e 400 \
#--temp_dir=/scratchfs/amendoza/tmp \
#-i {1} -o {2} \
#-d "$indexBS" \
#-g "$genome" &&
  
#### Sort the output bam files
#while read i; do sambamba sort -t "$cores" $i; done < "$baseR1"_merged_read_bams &&
  
#### Merge the single end map bam files
#sambamba merge -t "$cores" "$baseR1"_merged_reads.bam 0{1..6}."$baseR1"_merged_split.sorted.bam &&
  
#### Merge all bam files
#sambamba merge -t "$cores" "$baseR1".bam "$baseR1"_merged_reads.bam "$baseR1"_pairs.bam &&
  
### Remove the PCR duplicates for the paired reads
#sambamba markdup -r -t "$cores" -p --tmpdir=/scratchfs/tmp "$i" "$i"_dedup.bam;

### Call methylation for each chromosome (parallel) with CGmapTools
#parallel -j "$cores" cgmaptools convert bam2cgmap \
#--bam {} --genome "$genome" -o {.} ::: "$prefix"_chr*_all_dedup_temp.bam

### Each methylome was read into bsseq using this function

# Path to a CGmap file from BSseeker2
path <- args[1]

make_bsseq_obj <- function(CGmap_path, context="CG"){
  
  gc()
  
  id <- basename(CGmap_path) %>% str_replace(pattern = ".CGmap.gz", replacement = "")
  
  load_CGmap <- function(CGmap_path, forceZipped = FALSE) {
    header <- c("chr", "base", "position", "triContext",
                "diContext", "mC", "C_reads", "CT_reads")
    if(forceZipped == TRUE){
      dat <- data.table::fread(paste('gzip -dc ', path),
                               header = TRUE, col.names = header)
      return(dat)
    }
    ext <- tools::file_ext(CGmap_path)
    if(ext == "gz"){
      dat <- data.table::fread(paste('gzip -dc ', CGmap_path),
                               header = TRUE, col.names = header)
    } else {
      dat <- data.table::fread(CGmap_path, sep = "\t",
                               header = TRUE, col.names = header)
    }
    return(dat)
  }
  
  message(str_c("Loading ", CGmap_path, "..."))
  dat <- load_CGmap(CGmap_path = CGmap_path)
  
  # subset context
  dat <- dat[dat$diContext == context, ]
  
  # Add strand
  dat$base <- ifelse(test = dat$base == "C", yes = "+", no = "-")
  
  #Load data into BSseq object
  message(str_c("Making BSseq object for ", id, "..."))
  
  bs_obj <- BSseq(chr = dat$chr, pos=dat$position,
                  sampleNames = id,
                  M = matrix(dat$C_reads),
                  Cov = matrix(dat$CT_reads))
  
  # Add strand
  strand(bs_obj) <- dat$base
  
  # Collapse strand for CG
  if (context=="CG"){
    bs_obj <- strandCollapse(bs_obj)
  }
  
  message("Saving BSseq object...")
  saveRDS(object = bs_obj, file = str_c(id, "_", context, "_bsseq_obj.Rds"))
  message("Done!")
}

make_bsseq_obj(CGmap_path=path, context="CG")


##########################################
## Detecting DMRs

### Perform DMR analysis for noDox vs Dox
obj_files <- c("ZF598_D3awt_P2A_GFP_noDox_MethylC_rep1.CG_bsseq_obj.Rds",
               "ZF598_D3awt_P2A_GFP_noDox_MethylC_rep2.CG_bsseq_obj.Rds",
               "ZF598_D3awt_P2A_GFP_Dox1000_MethylC_rep1.CG_bsseq_obj.Rds",
               "ZF598_D3awt_P2A_GFP_Dox1000_MethylC_rep2.CG_bsseq_obj.Rds")

add_loci <- function(gr){
  
  loci <- str_c(seqnames(gr), start(gr), sep=":") %>%
    str_c(end(gr), sep = "-")
  
  gr$loci <- loci
  
  return(gr)
  
}

#--------Call DMRs

# Load the data
obj_list <- lapply(X = obj_files, readRDS)
obj_list <- bsseq::combineList(x = obj_list)


pData(obj_list)$condition <- factor(c(rep("NoDox", times=2), rep("Dox", times=2)))
pData(obj_list)$Replicate <- c(1:2,1:2)

# Remove CpG with no coverage
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(obj_list, type="Cov")==0) == 0)
obj_list <- obj_list[loci.idx, ]
obj_list <- obj_list[seqnames(granges(obj_list)) != "chrY", ]

# Save object for easy retreval
saveRDS(object = obj_list,
        file = "NoDox_Dox.rds")

# calling DMRs
regions <- dmrseq(obj_list, testCovariate = "condition", bpSpan = 500, maxGap = 500, maxPerms = 20, chrsPerChunk = 1)

T <- getCoverage(obj_list, regions = regions, type="Cov", what = "perRegionTotal")
C <- getCoverage(obj_list, regions = regions, type="M", what = "perRegionTotal")
regions$delta <- rowMeans(C[,c(1,2)]) - rowMeans(C[,c(3,4)])

saveRDS(object = regions,
        file = "NoDox_Dox_dmrseq_DMRs.rds")
write.table(data.frame(regions), file = "NoDox_Dox_dmrseq_DMRs.tsv", quote = F, sep = "\t", row.names = F)


##########################################
## Detecting mCG levels on DMRs

# Making the "all samples" file: 
obj_files <- c("ZF598_D3awt_P2A_GFP_noDox_MethylC_rep1.CG_bsseq_obj.Rds",
               "ZF598_D3awt_P2A_GFP_noDox_MethylC_rep2.CG_bsseq_obj.Rds",
               "ZF598_D3awt_P2A_GFP_Dox1000_MethylC_rep1.CG_bsseq_obj.Rds",
               "ZF598_D3awt_P2A_GFP_Dox1000_MethylC_rep2.CG_bsseq_obj.Rds",
               "ZF598_D3awt_P2A_GFP_Dox1000_withdraw7d_MethylC_rep1.CG_bsseq_obj.Rds",
               "ZF598_D3awt_P2A_GFP_Dox1000_withdraw7d_MethylC_rep2.CG_bsseq_obj.Rds",
               "ZF598_D3amut_P2A_GFP_noDox_MethylC_rep1.CG_bsseq_obj.Rds",
               "ZF598_D3amut_P2A_GFP_Dox1000_MethylC_rep1.CG_bsseq_obj.Rds")


# Read a Bs_seq boject from .Rds file

obj_list <- lapply(X = obj_files, readRDS)
bs_all <- bsseq::combineList(x = obj_list)


pData(bs_all)$condition <- factor(c(rep("noDox", times=2), rep("Dox", times=2),
                                    rep("DoxWD", times=2), "Mut", "MutDox"))
pData(bs_all)$Replicate <- c(1:2,1:2,1:2,1,1)

saveRDS(object = bs_all, file = "all_methylomes_bsseq.rds")

# filter DMRs with FDR < 0.1 and delta mCG > 20%
regions_filt <- regions[regions$qval < 0.1 & abs(regions$delta) > 0.2]

get_weighted_mC <- function(gr, bs = bs_all){
  T <- getCoverage(bs, regions = gr, type="Cov", what = "perRegionTotal")
  C <- getCoverage(bs, regions = gr, type="M", what = "perRegionTotal")
  C <- C/T
  rownames(C) <- add_loci(gr) %>% .$loci
  colnames(C) <- paste0(pData(bs)$condition,"_rep",pData(bs)$Replicate)
  return(C)
}

dmrs_01 <- get_weighted_mC(gr = regions_filt)

#### Counts CpGs per DMR and weighted average per sample

library("BSgenome.Hsapiens.UCSC.hg19")
library(GenomicFeatures)

seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, regions_filt)
CpG_count <- dinucleotideFrequency(seqs)[ ,"CG"]
regions_filt$CpG <- CpG_count

regions_filt$NoDox_mCG <- rowMeans(dmrs_01[,c(1,2)])
regions_filt$Dox_mCG <- rowMeans(dmrs_01[,c(3,4)])

saveRDS(regions_filt, file = "dmrseq_NoDox_Dox_filteredDMRs.rds")
write.table(dmrs_01, file = "noDox_Dox_dmrs_mCGvalues.tsv", quote = F, sep = "\t")


##########################################
######### Call UMRs with MethylSeekR
library("BSgenome.Hsapiens.UCSC.hg19")
library(MethylSeekR)
library(bsseq)
library(dplyr)
bs_noDox1 <- readRDS("ZF598_D3awt_P2A_GFP_noDox_MethylC_rep1.CG_bsseq_obj.Rds")
bs_noDox2 <- readRDS("ZF598_D3awt_P2A_GFP_noDox_MethylC_rep2.CG_bsseq_obj.Rds")
bs_noDox_both <- bsseq::combine(bs_noDox1, bs_noDox2)
pData(bs_noDox_both)$group <- c("A","A")
bs_noDox_both <- collapseBSseq(bs_noDox_both, group = c("A", "A"))
sLengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
gr_methyl <- granges(bs_noDox_both) 
gr_methyl <- gr_methyl[seqnames(gr_methyl) != "chrL"]
seqlevels(gr, pruning.mode="coarse") <- names(sLengths[1:23])
seqlengths(gr) <- sLengths[1:23]
PMDsegments.gr <- segmentPMDs(m=gr, chr.sel="chr22", seqLengths=sLengths, num.cores=1)
savePMDSegments(PMDs=PMDsegments.gr,
                GRangesFilename="PMDs_noDox.gr.rds", TableFilename="PMDs_noDox.tab")
UMRLMRsegments.gr <- segmentUMRsLMRs(m=gr, meth.cutoff=0.5,
                                     nCpG.cutoff=4, PMDs=PMDsegments.gr,
                                     num.cores=1, myGenomeSeq=Hsapiens,
                                     seqLengths=sLengths)
saveUMRLMRSegments(segs=UMRLMRsegments.gr,
                   GRangesFilename="UMRsLMRs_noDox.gr.rds", TableFilename="UMRsLMRs_noDox.tab")


#### Calculate mCG on UMRs
bs_all <- readRDS("all_methylomes_bsseq.rds")
UMRs <- readRDS("UMRsLMRs_noDox.gr.rds")

get_weighted_mC_cov_CpG <- function(gr, bs = bs_all){
  T <- getCoverage(bs, regions = gr, type="Cov", what = "perRegionTotal")
  C <- getCoverage(bs, regions = gr, type="M", what = "perRegionTotal")
  C <- C/T
  CpGs <- getCoverage(bs, regions = gr, type = "Cov", what = "perBase")
  CpGs <- lapply(CpGs, FUN = length) %>% unlist()
  CpGs <- CpGs/nrow(bs@colData)
  mean_coverage <- getCoverage(bs, regions = gr, type = "Cov", what = "perRegionAverage")
  mean_coverage <- rowMeans(mean_coverage)
  rownames(C) <- add_loci(gr) %>% .$loci
  colnames(C) <- paste0(pData(bs)$condition,"_rep",pData(bs)$Replicate)
  C <- cbind(C,CpGs,mean_coverage)
  return(C)
}

add_loci <- function(gr){
  
  loci <- str_c(seqnames(gr), start(gr), sep=":") %>%
    str_c(end(gr), sep = "-")
  
  gr$loci <- loci
  
  return(gr)
  
}

umrs <- get_weighted_mC_cov_CpG(gr = UMRs)
write.table(umrs, file = "UMRsLMRs_all_mCG.tsv", quote = F, sep = "\t")

##########################################
# Calculate mCG levels on ZF-D3A-peaks
read_bed_to_GRobject <- function(bedfile){
  dat <- fread(bedfile)
  gr <- GRanges(seqnames = Rle(dat$V1),
                ranges = IRanges(start = dat$V2, end = dat$V3))
  #dat <- dat[,!c(1:3)]
  #mcols(gr) <- dat
  
  return(gr)
}

ZF_gr <- read_bed_to_GRobject("ZF598_D3awt_P2A_GFP_HA_ChIP.nofilt.idr_peaks.bed")

zf_peaks <- get_weighted_mC_cov_CpG(gr = ZF_gr)

write.table(zf_peaks, file = "ZF_WT_Peaks_all_mCG.tsv", quote = F, sep = "\t")

##########################################
# Calculate mCG levels on ATAC-peaks

atac_gr <- read_bed_to_GRobject("ATACall_peaks_merge.bed")

atac_gr <- get_weighted_mC_cov_CpG(gr = atac_gr)

write.table(atac_gr, file = "ATACall_peaks_merge.mCG.tsv", quote = F, sep = "\t")


##########################################
# Plots mCG level differences in scatterplots (Figure 3a)

dmrs_mCG <- fread("noDox_Dox_dmrs_mCGvalues.tsv") %>% data.frame() %>% 
  dplyr::rename(loci = V1) %>% 
  mutate(noDox = (noDox_rep1+noDox_rep2)/2, Dox = (Dox_rep1+Dox_rep2)/2,DoxWD = (DoxWD_rep1+DoxWD_rep2)/2,) %>% 
  dplyr::select(loci,noDox,Dox,DoxWD,Mut_rep1,MutDox_rep1)


Lab.palette <- colorRampPalette(c("navyblue", "yellow", "red"), space = "Lab")
pdf(file = "Heatmaps_DMRs.pdf", height = 4, width = 10)
par(mfrow=c(1,3))
smoothScatter(dmrs_mCG$noDox, dmrs_mCG$Dox, colramp = Lab.palette, xlab = "mCG/CG ZF_D3A_wt noDox", ylab = "mCG/CG ZF_D3A_wt Dox", main = "DMRs (dmrseq FDR < 0.1)" )
abline(coef = c(0,1), col="grey", lwd=3, lty=2)
smoothScatter(dmrs_mCG$noDox, dmrs_mCG$DoxWD, colramp = Lab.palette, xlab = "mCG/CG ZF_D3A_wt noDox", ylab = "mCG/CG ZF_D3A_wt Dox Withdrawal")
abline(coef = c(0,1), col="grey", lwd=3, lty=2)
smoothScatter(dmrs_mCG$noDox, dmrs_mCG$MutDox_rep1, colramp = Lab.palette, xlab = "mCG/CG ZF_D3A_wt noDox", ylab = "mCG/CG ZF_D3A_mutant Dox")
abline(coef = c(0,1), col="grey", lwd=3, lty=2)
dev.off()

