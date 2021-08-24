library(stringr)
library(ggplot2)
library(reshape2)
library(edgeR)
library(DESeq2)
library(limma)
library(tidyr)
library(dplyr)
library(cowplot)
library(ChIPpeakAnno)


# read output file from RNA_processing.R 
counts <- read.table("ZF_Counts_table.tsv")

ERCC_counts <- counts[grepl(rownames(counts), pattern = "ERCC-"),]
gene_counts <- counts[!grepl(rownames(counts), pattern = "ERCC-"),]


# metadata trung:
mdat <- read.table("ZF_RNAseq_SpikeDesign.tsv", sep = "\t", header = T) %>%
  mutate(sample = str_replace(ZF_samples, pattern = "_rep\\d+", replacement = ""))
mdat$sample <- factor(mdat$sample, levels = c("WT_noDox","WT_Dox","WT_Dox_wd", "Mut_noDox","Mut_Dox"))

pDat <- data.frame(Cells = mdat$Cells, ERCC_spike = mdat$Type_of_Spike.in, sample = mdat$sample)

y <- DGEList(counts = counts, group = pDat$ERCC_spike)
y$samples <- cbind(y$samples, pDat)

########################################################
# compute the spike in statistics of the experiment:

percent_reads <- data.frame(RL = colnames(counts) %>% str_extract(., pattern = "RL\\d+"),
                            ERCC_reads = colSums(ERCC_counts), gene_reads = colSums(gene_counts),
                            total_counts = colSums(counts) ) %>% 
  mutate(perc_ERCC = 100*ERCC_reads/total_counts) %>%
  left_join(., mdat) %>% dplyr::select(RL, perc_ERCC, sample, total_counts,Total_RNA_amount_.ng.) %>% 
  dplyr::rename(Total_RNA_amount = Total_RNA_amount_.ng.)

percent_reads$sample <- factor(percent_reads$sample, levels = c("WT_noDox","WT_Dox","WT_Dox_wd", "Mut_noDox","Mut_Dox"))


# Plot the % of ERCC reads per total reads, the total amount of RNA after extraction (Qubit RNA BR kit) and total read counts

ggsave(
  plot_grid(ggplot(percent_reads, aes(x = sample, y = perc_ERCC, color = sample)) + geom_point(size = 3, alpha = 0.8) + theme_minimal() +
              ylab("Percentage of ERCC reads / library") + xlab("") + theme(legend.position = "none"),
            ggplot(percent_reads, aes(x = sample, y = Total_RNA_amount, color = sample)) + geom_point(size = 3, alpha = 0.8) + theme_minimal() +
              ylab("Amount of RNA after 50,000 cell \n+ 2ul spike-in (ng)") + xlab("") + theme(legend.position = "none"),
            ggplot(percent_reads, aes(x = sample, y = total_counts, color = sample)) + geom_point(size = 3, alpha = 0.8) + theme_minimal() +
              ylab("Total read counts (genes + ERCC)") + xlab("") + theme(legend.position = "none"), cols = 1),
  filename = "~/Dropbox/0_ZF598/Plots/Percent_Spike_point.pdf", height = 8, width = 5) 

wilcox.test(percent_reads$perc_ERCC[percent_reads$sample == "WT_noDox"],percent_reads$perc_ERCC[percent_reads$sample == "WT_Dox"], paired = F)


########################################################
# Cluster counts by spearman

d.correlation <- as.dist(1 - cor(gene_counts,method=c("spearman")))
fit <- hclust(d.correlation, method="complete")
pdf("~/Dropbox/0_ZF598/Plots/RNAseq_clustering_byDEgenes.pdf", height = 20)
plot(fit) # display dendogram
dev.off()

########################################################
# Perform differentially expression tests with DEseq2

### check differential expression using DEseq2

rownames(gene_counts) %>% length()
rownames(gene_counts) %>% unique() %>% length()
#no duplicate genes

# delete genes that have 0 counts in all conditions 
allZero <- rowSums(gene_counts==0)==ncol(gene_counts)
gene_counts <- gene_counts[!allZero,]

# delimit metadata
coldata <- data.frame(sample=colnames(gene_counts)) %>%
  mutate(condition = ifelse(grepl("D3awt_P2A_GFP_noDox", sample), "Ctrl_NoDox_wt",
                            ifelse(grepl("D3awt_P2A_GFP_Dox1000_3d_mRNAseq", sample), "Trt_Dox_wt",
                                   ifelse(grepl("D3amut_P2A_GFP_noDo", sample), "Ctrl_NoDox_mut",
                                          ifelse(grepl("D3amut_P2A_GFP_Dox1000", sample), "Trt_Dox_mut","Withdrawal_wt"))))) 

dds_all <- DESeqDataSetFromMatrix(countData = gene_counts,
                                  colData = coldata,
                                  design = ~ condition)

# 1st comparison: noDox vs Dox
dds <- dds_all

dds <- dds[, dds$condition == "Ctrl_NoDox_wt" | 
             dds$condition == "Trt_Dox_wt" ]
dds$condition <- droplevels(dds$condition)

#filter low expressed genes (at least 10 counts in total, or more than 0 in at least 50% of the samples)
dds <- dds[rowSums(counts(dds)) >= 10,]
dds_good <- dds[rowSums(counts(dds)==0) <= ncol(counts(dds))*0.5,]

dds_good <- estimateSizeFactors(dds_good)

rowData(dds_good)$control_expr <- rowMeans(counts(dds_good, normalize = TRUE)[,dds_good$condition == "Ctrl_NoDox_wt"])
rowData(dds_good)$condition_expr <- rowMeans(counts(dds_good, normalize = TRUE)[,dds_good$condition == "Trt_Dox_wt"])

#######
#do the test
dds_good <- DESeq(dds_good)
res_good <- results(dds_good)

fdrlevel.de <- 0.05
sum(res_good$padj < fdrlevel.de, na.rm = TRUE)
#11350 DE genes
data.frame(res_good) %>% mutate(change_direction = ifelse(log2FoldChange < 0, "downregulated", "upregulated")) %>%
  filter(padj < fdrlevel.de) %>% .$change_direction %>% table()
#downregulated   upregulated
#5446          5904
res_good$control_expr <- rowData(dds_good)$control_expr
res_good$condition_expr <- rowData(dds_good)$condition_expr


dds_good <- dds_good[!is.na(res_good$padj)]
res_good <- res_good %>% na.omit()

saveRDS(res_good, file = "res_good.DEseq2.rds")

#############
### test mutant
dds_mut <- dds_all[, dds_all$condition == "Ctrl_NoDox_mut" | 
                     dds_all$condition == "Trt_Dox_mut" ]
dds_mut$condition <- droplevels(dds_mut$condition)

#filter low expressed
dds_mut <- dds_mut[rowSums(counts(dds_mut)) >= 10,]
dds_mut <- dds_mut[rowSums(counts(dds_mut)==0) <= ncol(counts(dds_mut))*0.5,]

dds_mut <- estimateSizeFactors(dds_mut)

rowData(dds_mut)$control_expr <- rowMeans(counts(dds_mut, normalize = TRUE)[,dds_mut$condition == "Ctrl_NoDox_mut"])
rowData(dds_mut)$condition_expr <- rowMeans(counts(dds_mut, normalize = TRUE)[,dds_mut$condition == "Trt_Dox_mut"])


#######
#do the test
dds_mut <- DESeq(dds_mut)
res_mut <- results(dds_mut)

sum(res_mut$padj < fdrlevel.de, na.rm = TRUE)
#5096 DE genes
data.frame(res_mut) %>% mutate(change_direction = ifelse(log2FoldChange < 0, "downregulated", "upregulated")) %>%
  filter(padj < fdrlevel.de) %>% .$change_direction %>% table()
#downregulated   upregulated
#2547          2549
res_mut$control_expr <- rowData(dds_mut)$control_expr
res_mut$condition_expr <- rowData(dds_mut)$condition_expr


dds_mut <- dds_mut[!is.na(res_mut$padj)]
res_mut <- res_mut %>% na.omit()

saveRDS(res_mut, file = "res_mut.DEseq2.rds")


#############
### test NoDoxes
dds_noDox <- dds_all[, dds_all$condition == "Ctrl_NoDox_wt" | 
                       dds_all$condition == "Ctrl_NoDox_mut" ]
dds_noDox$condition <- droplevels(dds_noDox$condition)

#filter low expressed
dds_noDox <- dds_noDox[rowSums(counts(dds_noDox)) >= 10,]
dds_noDox <- dds_noDox[rowSums(counts(dds_noDox)==0) <= ncol(counts(dds_noDox))*0.5,]

dds_noDox <- estimateSizeFactors(dds_noDox)

rowData(dds_noDox)$control_expr <- rowMeans(counts(dds_noDox, normalize = TRUE)[,dds_noDox$condition == "Ctrl_NoDox_wt"])


#######
#do the test corrected
dds_noDox <- DESeq(dds_noDox)
res_noDox <- results(dds_noDox)

sum(res_noDox$padj < fdrlevel.de, na.rm = TRUE)
#11 DE genes
data.frame(res_noDox) %>% mutate(change_direction = ifelse(log2FoldChange < 0, "downregulated", "upregulated")) %>%
  filter(padj < fdrlevel.de) %>% .$change_direction %>% table()
#downregulated   upregulated
#7          4
res_noDox$control_expr <- rowData(dds_noDox)$control_expr

dds_noDox <- dds_noDox[!is.na(res_noDox$padj)]
res_noDox <- res_noDox %>% na.omit()


#######
# Plotting the MA plots of differential expression (Supplementary Figure 4e)

genes_to_plot_MA <- c("SOX2","GAPDH","DNMT3A")

ggMAplot <- function(deseq_res, lab = "(WT noDox/Dox)"){
  deseq_res$gene_id <- rownames(deseq_res)
  ggplot(data.frame(deseq_res), aes(x = baseMean, y = log2FoldChange)) + 
    geom_point(size = 0.5, alpha = 0.5, mapping = aes(color = padj < 0.05)) +
    ylab(paste0("log Fold Change", lab)) +
    xlab("mean normalised expression") +
    geom_hline(yintercept=0, col="grey") +
    geom_hline(yintercept=1, col="grey", linetype="dashed") +
    geom_hline(yintercept=-1, col="grey", linetype="dashed") +
    geom_point(data =  data.frame(deseq_res[deseq_res$gene_id %in% genes_to_plot_MA,]), col="darkred") +
    geom_label(data =  data.frame(deseq_res[deseq_res$gene_id %in% genes_to_plot_MA,]), col="darkred",
               aes(label=deseq_res$gene_id[deseq_res$gene_id %in% genes_to_plot_MA]), hjust=0.5, vjust=1.25) +
    annotate(geom="text", x=1e+04, y=-6, label=paste("Downregulated genes =",nrow(deseq_res[deseq_res$log2FoldChange < 0 & deseq_res$padj < 0.05,])),color="black") +
    annotate(geom="text", x=1e+04, y=6, label=paste("Upregulated genes =",nrow(deseq_res[deseq_res$log2FoldChange > 0 & deseq_res$padj < 0.05,])),color="black") +
    theme_minimal() + scale_x_log10()
}

ggMA_wt_noDox_Dox <- ggMAplot(deseq_res = res_good, lab = "(WT Dox vs noDox)")
ggMA_mut_noDox_Dox <- ggMAplot(deseq_res = res_mut, lab = "(Mutant Dox vs noDox)")
ggMA_mut_noDoxes <- ggMAplot(deseq_res = res_noDox, lab = "(WT noDox vs Mutant noDox)")


ggsave(plot_grid(ggMA_wt_noDox_Dox, ggMA_mut_noDox_Dox, ggMA_noDoxes, cols = 1),
       filename = "MAplots.png", height = 10, width = 5)

##############################
## Upset DEseq2 results as in Figure 3b
library(UpSetR)
geneIDs <- c(rownames(res_good), rownames(res_mut), rownames(res_noDox), rownames(res_wd)) %>% unique()


overlap_DEgene_UpDown_sets <- function(res_df , genelist = geneIDs){
  logical_up <- genelist %in% rownames(res_df)[res_df$padj < fdrlevel.de & res_df$log2FoldChange > 0] %>% as.integer()
  logical_down <- genelist %in% rownames(res_df)[res_df$padj < fdrlevel.de & res_df$log2FoldChange < 0] %>% as.integer()
  logicals <- cbind(logical_up, logical_down)
  return(logicals)
}

wt_vs_mut_overlaps <- cbind(overlap_DEgene_UpDown_sets(res_df = res_good),
                            overlap_DEgene_UpDown_sets(res_df = res_mut)) %>% data.frame()
colnames(wt_vs_mut_overlaps) <- c("WT_Upregulated_Dox","WT_Downregulated_Dox",
                                  "Mut_Upregulated_Dox","Mut_Downregulated_Dox")

pdf("Upset_DE_genes_by_direction.pdf", useDingbats = FALSE)
upset(wt_vs_mut_overlaps, 
      show.numbers = FALSE,
      sets.x.label = "Total number of DE genes (fdr < 0.05)",
      mainbar.y.label = "Number of DE genes in intersection",
      #group.by = "sets",
      #sets.bar.color = "#56B4E9",
      order.by = "freq")
dev.off()

##############
# Are ZF-mut genes changing in the same direction than ZF-D3A genes? Supplementary Figure 4f


res_mut$gene_id <- rownames(res_mut)
res_good$gene_id <- rownames(res_good)

#join results mutant and wt, filter for DE in mutant
mutant_vs_wt <- left_join(data.frame(res_mut),data.frame(res_good), by = "gene_id" ) %>% 
  dplyr::select(gene_id, log2FoldChange.x, log2FoldChange.y, padj.x, padj.y) %>%
  dplyr::filter(padj.x < 0.05 & padj.y < 0.05)


ggsave(ggplot(mutant_vs_wt, aes(x = log2FoldChange.x, log2FoldChange.y)) + 
         geom_bin2d() + 
         theme_bw() +
         geom_smooth(method = "lm", span = 0.25) +
         annotate(geom="text", x=-1, y=10, label=paste("spearman r =",round(cor(mutant_vs_wt$log2FoldChange.x, mutant_vs_wt$log2FoldChange.y, method = "spearman"),digits = 2)),color="black") +
         xlab("log2(mRNA fold change)\nZF-D3A mutant Dox vs no Dox") + 
         ylab("log2(mRNA fold change)\nZF-D3A wt Dox vs no Dox"),
       filename = "Mutant_vs_WT_foldChanges_ggplot.pdf", height = 3, width = 4)



