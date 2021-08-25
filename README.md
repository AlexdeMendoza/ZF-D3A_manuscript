# ZF-D3A


Files and code used for de Mendoza, A., Nguyen, T., Ford, E. et al. "Limited repressive capacity of promoter DNA methylation revealed through epigenome manipulation". 


**Code**
* RNAseq_processing.R
* RNAseq_analysis.R
* WGBS_processing_and_DMRs.R
* PromoterDMR_vs_Txn.R
* ChIP-seq_processing.R
* ATAC-seq_processing.R

**Files**
* `ZF_RNAseq_SpikeDesign.tsv`: ERCC-spike in design for the RNA-seq samples.
* `DMRs.tsv.gz`: DMRs obtained from dmrseq comparing noDox vs Dox and filtered by FDR < 0.1.
* `genes.gtf.gz`: Version of UCSC hg19 annotation file used in the study.
* `PMDs_noDox.bed.gz`: Partially Methylated Domains for noDox obtained with MethylSeekR.
* `UMRsLMRs_noDox.bed.gz`: Unmethylated and Lowly Methylated Regions for noDox obtained with MethylSeekR.
* `ZF_Counts_table.tsv.gz`: mRNA-seq counts table.
* `ZF598_D3awt_P2A_GFP_HA_ChIP.nofilt.idr_peaks.gz`: ZF-D3A-wt ChIP-seq peaks merged with IDR.
* `ZF598_D3amut_FEER_P2A_GFP_HA_ChIP.nofilt.idr_peaks.gz`: ZF-D3A-mut ChIP-seq peaks merged with IDR.
* `hg19_ZF598_scan.bed.gz`: Motif scan for the ZF-D3A motif obtained using HOMER2.
* `ATACall_peaks_merge.bed.gz`: ATAC-seq peaks for all samples merged.


