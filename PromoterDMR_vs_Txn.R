library(stringr)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(tidyr)
library(dplyr)
library(cowplot)
library(ChIPpeakAnno)
library(GenomicFeatures)


# Import features from UCSC annotation
genePath <- "genes.gtf"
txdb <- makeTxDbFromGFF(file = genePath)
#promoters extracted by transcript, not gene, to take into account alternative promoters
promoters <- promoters(transcripts, upstream = 2000, downstream = 200)


# Import DMRs 
dmrs_filtered <- readRDS("NoDox_Dox.dmrs_filtered.rds")

# Import DEseq result from noDox vs Dox for both ZF-D3A-wt and ZF-D3A-mut
res_good <- readRDS("res_good.DEseq2.rds")
res_mut <- readRDS("res_mut.DEseq2.rds")

# Link DMRs to promoters
overlap_and_add_metadata_from_gr2 <- function (gr1, gr2){
  gr2 <- add_loci(gr2)
  hits <- findOverlaps(gr1, gr2, ignore.strand = TRUE)
  gr_out <- gr1[hits@from]
  gr_out@elementMetadata <- cbind(gr_out@elementMetadata,gr2[hits@to]@elementMetadata )
  gr_out <- unique(gr_out)
  return(gr_out)
}

promoters_2 <- unlist(promoters)
promoters_2$gene_id <- names(promoters_2) %>% as.character()

# get the overlaps and add the information from the DMRs
dmrs_in_promoters <- overlap_and_add_metadata_from_gr2(gr1 = promoters_2, gr2 = dmrs_filtered) %>% mcols() %>% data.frame() %>%
  dplyr::select(gene_id, delta, CpG, qval, stat,NoDox_mCG, Dox_mCG, loci )

# get the information of the DEseq2 noDox vs Dox
deseq_noDox_Dox <- res_good %>% data.frame() %>% dplyr::select(log2FoldChange, padj, gene_id, control_expr, condition_expr) %>% 
  dplyr::rename(padj_DEseq2 = padj )

# associate DEseq to promoter-DMRs
dmrsa2 <- left_join(dmrs_in_promoters, deseq_noDox_Dox) %>% 
  dplyr::filter(control_expr > 50 | condition_expr > 50) %>%
  mutate(sig = padj_DEseq2 < 0.05, FC = 2^log2FoldChange, stat = -stat) %>%
  mutate(ptocol2 = ifelse(sig, "Differentially Expressed", 
                          "Not Differentially Expressed"),
         ptcol3 = ifelse(qval > 0.01, "Not significant (DMR)", 
                         ifelse(sig, "DMR and DE", "Not significant (DE)"))) %>%
  dplyr::filter(!is.na(log2FoldChange) & !is.na(sig)) %>%
  unique()


genes_to_plot <- c("SOX2","GAPDH","PSME1")
fdrlevel.dmr.highconf <- 0.01

gg_scatter_stat <- ggplot(dmrsa2, aes(x = stat, y = log2FoldChange)) +
  geom_hline(yintercept=0, col="black") +
  geom_hline(yintercept=0, col="white", linetype="dashed") +
  geom_point(size=0.5, alpha=0.75, aes(color = ptcol3)) + 
  theme_bw() + 
  xlab("Region test statistic") +
  ylab("log2 fold change mRNA abundance") +
  scale_color_manual(values=c("red", "black", "grey")) +
  geom_smooth(method = "loess", span = 0.25) +
  labs(color="Significance")  +
  geom_point(data =  dmrsa2[dmrsa2$gene_id %in% genes_to_plot,], col="darkred") +
  geom_label(data =  dmrsa2[dmrsa2$gene_id %in% genes_to_plot,], col="darkred",
             aes(label=dmrsa2$gene_id[dmrsa2$gene_id %in% genes_to_plot]), hjust=0.5, vjust=1.25) +
  ggtitle("dmrseq DMRs on promoters") +
  ylim(-3,3) +
  geom_vline(xintercept=min(dmrsa2$stat[dmrsa2$qval<fdrlevel.dmr.highconf]), 
             linetype="dashed", color="grey20") +
  theme(legend.position="none")

gg_scatter_delta <- ggplot(dmrsa2, aes(x = delta, y = log2FoldChange)) +
  geom_hline(yintercept=0, col="black") +
  geom_hline(yintercept=0, col="white", linetype="dashed") +
  geom_point(size=0.5, alpha=0.75, aes(color = ptcol2)) + 
  theme_bw() + 
  xlab("mCG Dox - mCG noDox") +
  ylab("log2 fold change mRNA abundance") +
  scale_color_manual(values=c("red", "black", "grey")) +
  geom_smooth(method = "loess", span = 0.25) +
  labs(color="Significance")  +
  geom_point(data =  dmrsa2[dmrsa2$gene_id %in% genes_to_plot,], col="darkred") +
  geom_label(data =  dmrsa2[dmrsa2$gene_id %in% genes_to_plot,], col="darkred",
             aes(label=dmrsa2$gene_id[dmrsa2$gene_id %in% genes_to_plot]), hjust=0.5, vjust=1.25) +
  ggtitle("dmrseq DMRs on promoters") +
  ylim(-3,3) +
  #geom_vline(xintercept=min(dmrsa2$stat[dmrsa2$qval<fdrlevel.dmr.highconf]), 
  #          linetype="dashed", color="grey20") +
  theme(legend.position="none")

ggsave(gg_scatter_stat, filename = "Scatter_DMRs_promoter_stat.pdf")
ggsave(gg_scatter_delta, filename = "Scatter_DMRs_promoter_deltamCG.pdf")

##############################################
# Plot again discarding the ZF affected genes

dmrsa3 <- dmrsa2 %>% filter(!gene_id %in% res_mut$gene_id[res_mut$padj < 0.05] ) %>% unique()

gg_scatter_stat_filt <- ggplot(dmrsa3, aes(x = stat, y = log2FoldChange)) +
  geom_hline(yintercept=0, col="black") +
  geom_hline(yintercept=0, col="white", linetype="dashed") +
  geom_point(size=0.5, alpha=0.75, aes(color = ptcol3)) + 
  theme_bw() + 
  xlab("Region test statistic") +
  ylab("log2 fold change mRNA abundance") +
  scale_color_manual(values=c("red", "black", "grey")) +
  geom_smooth(method = "loess", span = 0.25) +
  labs(color="Significance")  +
  geom_point(data =  dmrsa3[dmrsa3$gene_id %in% genes_to_plot,], col="darkred") +
  geom_label(data =  dmrsa3[dmrsa3$gene_id %in% genes_to_plot,], col="darkred",
             aes(label=dmrsa3$gene_id[dmrsa3$gene_id %in% genes_to_plot]), hjust=0.5, vjust=1.25) +
  ggtitle("dmrseq DMRs on promoters") +
  ylim(-3,3) +
  geom_vline(xintercept=min(dmrsa3$stat[dmrsa3$qval<fdrlevel.dmr.highconf]), 
             linetype="dashed", color="grey20") +
  theme(legend.position="none")

gg_scatter_delta_filt <- ggplot(dmrsa3, aes(x = delta, y = log2FoldChange)) +
  geom_hline(yintercept=0, col="black") +
  geom_hline(yintercept=0, col="white", linetype="dashed") +
  geom_point(size=0.5, alpha=0.75, aes(color = ptcol3)) + 
  theme_bw() + 
  xlab("mCG Dox - mCG noDox") +
  ylab("log2 fold change mRNA abundance") +
  scale_color_manual(values=c("red", "black", "grey")) +
  geom_smooth(method = "loess", span = 0.25) +
  labs(color="Significance")  +
  geom_point(data =  dmrsa3[dmrsa3$gene_id %in% genes_to_plot,], col="darkred") +
  geom_label(data =  dmrsa3[dmrsa3$gene_id %in% genes_to_plot,], col="darkred",
             aes(label=dmrsa3$gene_id[dmrsa3$gene_id %in% genes_to_plot]), hjust=0.5, vjust=1.25) +
  ggtitle("dmrseq DMRs on promoters") +
  ylim(-3,3) +
  #geom_vline(xintercept=min(dmrsa2$stat[dmrsa2$qval<fdrlevel.dmr.highconf]), 
  #          linetype="dashed", color="grey20") +
  theme(legend.position="none")

ggsave(gg_scatter_stat_filt, filename = "~/Dropbox/0_ZF598/Plots/Scatter_DMRs_promoter_stat_NoDEG_mutant.pdf")
ggsave(gg_scatter_stat_filt, filename = "~/Dropbox/0_ZF598/Plots/Scatter_DMRs_promoter_stat_NoDEG_mutant.png")
ggsave(gg_scatter_delta_filt, filename = "~/Dropbox/0_ZF598/Plots/Scatter_DMRs_promoter_deltamCG_NoDEG_mutant.pdf")


######################################
### check proportion of genes that are upregulated, downregulated or not differentially expressed

compare_gene_set_DMR <- function(DMR_genes, normalisation, method_name, fdr_de = 0.05 ){
  df <- normalisation[rownames(normalisation) %in% DMR_genes,] %>% data.frame() %>%
    mutate(gene_id = rownames(.)) %>% 
    mutate(DE = ifelse(padj < fdr_de,"YES","NO")) %>%
    mutate(sign = ifelse(log2FoldChange < 0,"down","up")) %>%
    mutate(class = ifelse(DE == "NO", paste0("Not Diff. Expr (DEseq2 qval > ",fdr_de,")"),ifelse(sign == "up",paste0("Upregulated (DEseq2 qval < ",fdr_de,")"),paste0("Downregulated (DEseq2 qval < ",fdr_de,")")))) %>%
    dplyr::select(gene_id,class) %>% unique() %>% .$class %>% table() %>% data.frame() %>% 
    dplyr::rename(category = ".") %>%  
    mutate(perc = Freq/sum(Freq)*100) %>% 
    mutate(Method = method_name)
  return(df) 
  
}


A <- compare_gene_set_DMR(DMR_genes = unique(dmrsa2$gene_id[dmrsa2$qval < 0.01]), normalisation = res_good, method_name = "dmrseq DMRs (fdr < 0.01)\n All genes", fdr_de = 0.05)
B <- compare_gene_set_DMR(DMR_genes = unique(dmrsa3$gene_id[dmrsa3$qval < 0.01]), normalisation = res_good, method_name = "dmrseq DMRs (fdr < 0.01)\n Genes not differentially expressed\nin Dox ZF_D3A_mutant", fdr_de = 0.05)

final_comparison <- rbind(B, A)

cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

gg_proportions <- ggplot(data = final_comparison, 
                         aes(x = Method,y = as.numeric(perc), fill = category, label = paste0(Freq, "\n(",round(perc, digits = 1),"%)"))) + 
  geom_bar( stat = "identity" ) +
  geom_text(size = 3, color = "white", position = position_stack(vjust = 0.5)) +
  scale_y_continuous() + 
  
  # X axis label
  xlab(label = "") + 
  
  # Y axis label
  ylab(label = "percentage (%)") +
  
  
  #theme with white background
  theme_bw() + 
  #eliminates background, gridlines, and chart border
  theme( axis.text.x = element_text(angle = 45, hjust = 1)
         ,plot.background = element_blank()
         ,panel.grid.major = element_blank()
         ,panel.grid.minor = element_blank()
         ,panel.border = element_blank()
  ) +
  
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  
  # Edit legend
  scale_fill_manual(values = cbPalette) +
  coord_flip()

ggsave(gg_proportions, filename = "Proportion_genes_DMRpromoter.pdf",
       height = 3, width = 10)
