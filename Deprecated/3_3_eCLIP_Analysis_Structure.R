################################################################################
## eCLIP Analysis V5: 
## Written by Soon Yi
## Created: 2024-05-23
## Last Edited: 2024-05-24
## Figure 3
################################################################################

library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(corrplot)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(GenomicRanges)
library(IRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

## Custom Functions:
################################################################################
# Min-Max normalization:
min_max_norm = function(x, a = 0, b = 1) {
  a + (x - min(x, na.rm = TRUE)) * (b - a) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Specificity index calculation:
SI_Calculation = function(x) {
  x / (median(x, na.rm = FALSE))
}

# Scramble given sequence
ScrambleDNA = function(seq) {
  scrambled_seq = sample(unlist(strsplit(as.character(seq), "")))
  return(DNAString(paste(scrambled_seq, collapse = "")))
}

# Generate N random distances with a minimum of X nt away
generate_random_distances = function(N, X) {
  distances_negative = runif(N / 2, min = -1000, max = -X)
  distances_positive = runif(N / 2, min = X, max = 1000)
  distances = sample(c(distances_negative, distances_positive))
  return(distances)
}

################################################################################

## Set up basic parameters:
################################################################################
# baseDir = '~/Desktop/Genomics/Specificity/'

baseDir = '~/Desktop/Genomics/Specificity/Data_Final/'
sampleTable_eCLIP = read_csv(paste0(baseDir, 'sample_table_eCLIP.csv'), col_names = T, show_col_types = F)
sampleTable_eCLIP = data.frame(sampleTable_eCLIP)

cell_lines = c('K562', 'HepG2')
# K562_RBPs = c('HNRNPC', 'EIF4G2', 'RBFOX2')
# HepG2_RBPs = c('HNRNPC', 'PCBP2', 'RBFOX2')
K562_RBPs = as.character(unique((sampleTable_eCLIP %>% filter(Cell == 'K562'))[, 'RBP']))
HepG2_RBPs = as.character(unique((sampleTable_eCLIP %>% filter(Cell == 'HepG2'))[, 'RBP']))


DATE = '20240524_3'
extension = 25
################################################################################

## 0523:   50nt extension both ways, RNAplfold options default
## 0524:   50nt extension only 5'-end upstream, RNAplfold -W 25 (because sometimes len < 70nt)
## 0524_2: 25nt extension only 5'-end upstream, RNAplfold -W 25 (because sometimes len < 70nt)
## 0524_3: 25nt extension only 5'-end upstream, RNAplfold -W 25 -u 1 (and use lunp)

################################################################################
for (cell_line in cell_lines) {
  print(paste0('Working on ', cell_line, '....'))
  if (cell_line == 'K562') {RBPs = K562_RBPs} else {RBPs = HepG2_RBPs}
  for (RBP in RBPs) {
    # Read peaks and make GR object
    print(paste0('Working on ', RBP, '....'))
    print(paste0('     Reading and processing peak data....'))
    file_ID = sampleTable_eCLIP[sampleTable_eCLIP$RBP == RBP & sampleTable_eCLIP$Cell == cell_line, ]$eCLIP
    Peak = fread(paste0(baseDir, 'eCLIP_peak/raw/', file_ID,  '.bed'))
    Peak = Peak %>% mutate(V1 = case_when(V1 == 'chrMT' ~ 'chrM', TRUE ~ as.character(V1)))
    Peak = Peak %>% filter(V1 %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 
                                     'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
    Peak_GR = GRanges(seqnames = Rle(Peak$V1),
                      ranges = IRanges(start = Peak$V2 - extension, end = Peak$V3),
                      strand = Rle(Peak$V6))
    
    # Get sequence of peaks and count K-mers
    Peak_Seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, Peak_GR, as.character = TRUE)
    Peak_Seqs = DNAStringSet(Peak_Seqs)
    names(Peak_Seqs) = paste0(Peak$V4, "_", seq_along(Peak$V4))
    writeXStringSet(Peak_Seqs, paste0(baseDir, 'eCLIP_peak/peak_fasta/', DATE, '/', cell_line, '_', RBP, '_eCLIP_Peaks_', extension, 'nt.fa'))
  }
}
################################################################################

## Run bash /Users/soonyi/Repos/Specificity_BITs/Figure3/3_3_eCLIP_Analysis_Structure.sh

# DATE = '20240523'
# DATE = '20240524'
# DATE = '20240524_2'
# DATE = '20240524_3'
DATE = '20240524_4'

Peaks_BPPs = list()
Peaks_MFEs = list()

for (cell_line in cell_lines) {
  if (cell_line == 'K562') {RBPs = K562_RBPs} else {RBPs = HepG2_RBPs}
  for (RBP in RBPs) {
    print(paste0('Working on ', RBP, '....'))
    print(paste0('     Reading and processing peak data....'))
    files = list.files(path = paste0(baseDir, 'eCLIP_peak/peak_fasta/', DATE, '/', cell_line, '_', RBP, '_eCLIP'))
    # avg_bpp = c()
    mfes = c()
    for (file in files) {
      # upp = fread(paste0(baseDir, 'eCLIP_peak/peak_fasta/', DATE, '/', cell_line, '_', RBP, '_eCLIP', '/', file))
      # avg_bpp = c(avg_bpp, mean(1 - upp$`l=1`[c(25:length(upp$`l=1`))]))
      MFE = fread(paste0(baseDir, 'eCLIP_peak/peak_fasta/', DATE, '/', cell_line, '_', RBP, '_eCLIP', '/', file))
      mfes = c(mfes, as.numeric(gsub("[()]", "", MFE[1, 2])))
    }
    # Peaks_BPPs[paste0(cell_line, '_', RBP)] = list(avg_bpp)
    Peaks_MFEs[paste0(cell_line, '_', RBP)] = list(mfes)
  }
}


# plot_df = stack(Peaks_BPPs)
plot_df = stack(Peaks_MFEs)
plot_df = plot_df[complete.cases(plot_df), ]

order_list = c('K562_EIF4G2',
               'K562_EWSR1',
               'K562_HNRNPK_1',
               'K562_HNRNPK_2',
               'K562_PUM1',
               'K562_RBFOX2',
               'HepG2_RBFOX2',
               'K562_RBM22',
               'HepG2_RBM22',
               'HepG2_SRSF9',
               'K562_TAF15',
               'HepG2_TAF15',
               'K562_TRA2A',
               'HepG2_TRA2A',
               'HepG2_FUBP3',
               'K562_FUS',
               'HepG2_FUS',
               'K562_HNRNPC',
               'HepG2_HNRNPC',
               'K562_HNRNPL',
               'HepG2_HNRNPL',
               'K562_IGF2BP2',
               'K562_KHSRP',
               'HepG2_KHSRP',
               'K562_PCBP1',
               'HepG2_PCBP1',
               'HepG2_PCBP2',
               'K562_TIA1',
               'HepG2_TIA1',
               'K562_AKAP8L',
               'K562_APOBEC3C',
               'HepG2_EIF3D',
               'HepG2_IGF2BP3',
               'K562_LIN28B',
               'HepG2_LIN28B',
               'HepG2_RBM5',
               'K562_SAFB2',
               'K562_TROVE2',
               'HepG2_TROVE2',
               'K562_XRCC6',
               'HepG2_XRCC6',
               'K562_ZRANB2')
plot_df$ind = factor(plot_df$ind, levels = order_list)

colors = c('K562_EIF4G2' = 'antiquewhite2',
           'K562_EWSR1' = 'antiquewhite2',
           'K562_HNRNPK_1' = 'antiquewhite2',
           'K562_HNRNPK_2' = 'antiquewhite2',
           'K562_PUM1' = 'antiquewhite2',
           'K562_RBFOX2' = 'antiquewhite2',
           'HepG2_RBFOX2' = 'antiquewhite2',
           'K562_RBM22' = 'antiquewhite2',
           'HepG2_RBM22' = 'antiquewhite2',
           'HepG2_SRSF9' = 'antiquewhite2',
           'K562_TAF15' = 'antiquewhite2',
           'HepG2_TAF15' = 'antiquewhite2',
           'K562_TRA2A' = 'antiquewhite2',
           'HepG2_TRA2A' = 'antiquewhite2',
           'HepG2_FUBP3' = 'darkseagreen4',
           'K562_FUS' = 'darkseagreen4',
           'HepG2_FUS' = 'darkseagreen4',
           'K562_HNRNPC' = 'darkseagreen4',
           'HepG2_HNRNPC' = 'darkseagreen4',
           'K562_HNRNPL' = 'darkseagreen4',
           'HepG2_HNRNPL' = 'darkseagreen4',
           'K562_IGF2BP2' = 'darkseagreen4',
           'K562_KHSRP' = 'darkseagreen4',
           'HepG2_KHSRP' = 'darkseagreen4',
           'K562_PCBP1' = 'darkseagreen4',
           'HepG2_PCBP1' = 'darkseagreen4',
           'HepG2_PCBP2' = 'darkseagreen4',
           'K562_TIA1' = 'darkseagreen4',
           'HepG2_TIA1' = 'darkseagreen4',
           'K562_AKAP8L' = 'grey',
           'K562_APOBEC3C' = 'grey',
           'HepG2_EIF3D' = 'grey',
           'HepG2_IGF2BP3' = 'grey',
           'K562_LIN28B' = 'grey',
           'HepG2_LIN28B' = 'grey',
           'HepG2_RBM5' = 'grey',
           'K562_SAFB2' = 'grey',
           'K562_TROVE2' = 'grey',
           'HepG2_TROVE2' = 'grey',
           'K562_XRCC6' = 'grey',
           'HepG2_XRCC6' = 'grey',
           'K562_ZRANB2' = 'grey')


ggplot(plot_df, aes(x = ind, y = values, fill = ind)) +
  geom_boxplot(outliers = FALSE, notch = TRUE) +
  labs(x = "Group", y = "Value", title = "Boxplot of Values") +
  scale_fill_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")


wilcox.test(Peaks_BPPs$K562_EIF4G2, Peaks_BPPs$K562_HNRNPC)
ks.test(Peaks_BPPs$K562_EIF4G2, Peaks_BPPs$K562_HNRNPC)





t.test(MFE_means[c('K562_EIF4G2',
                        'K562_EWSR1',
                        'K562_HNRNPK_1',
                        'K562_HNRNPK_2',
                        'K562_PUM1',
                        'K562_RBFOX2',
                        'HepG2_RBFOX2',
                        'K562_RBM22',
                        'HepG2_RBM22',
                        'HepG2_SRSF9',
                        'K562_TAF15',
                        'HepG2_TAF15',
                        'K562_TRA2A',
                        'HepG2_TRA2A')], 
            MFE_means[c('HepG2_FUBP3',
                        'K562_FUS',
                        'HepG2_FUS',
                        'K562_HNRNPC',
                        'HepG2_HNRNPC',
                        'K562_HNRNPL',
                        'HepG2_HNRNPL',
                        'K562_IGF2BP2',
                        'K562_KHSRP',
                        'HepG2_KHSRP',
                        'K562_PCBP1',
                        'HepG2_PCBP1',
                        'HepG2_PCBP2',
                        'K562_TIA1',
                        'HepG2_TIA1')])




unlist()


wilcox.test(unlist(Peaks_BPPs[c('K562_EIF4G2',
                           'K562_EWSR1',
                           'K562_HNRNPK_1',
                           'K562_HNRNPK_2',
                           'K562_PUM1',
                           'K562_RBFOX2',
                           'HepG2_RBFOX2',
                           'K562_RBM22',
                           'HepG2_RBM22',
                           'HepG2_SRSF9',
                           'K562_TAF15',
                           'HepG2_TAF15',
                           'K562_TRA2A',
                           'HepG2_TRA2A')]), 
       unlist(Peaks_BPPs[c('HepG2_FUBP3',
                           'K562_FUS',
                           'HepG2_FUS',
                           'K562_HNRNPC',
                           'HepG2_HNRNPC',
                           'K562_HNRNPL',
                           'HepG2_HNRNPL',
                           'K562_IGF2BP2',
                           'K562_KHSRP',
                           'HepG2_KHSRP',
                           'K562_PCBP1',
                           'HepG2_PCBP1',
                           'HepG2_PCBP2',
                           'K562_TIA1',
                           'HepG2_TIA1')]))
