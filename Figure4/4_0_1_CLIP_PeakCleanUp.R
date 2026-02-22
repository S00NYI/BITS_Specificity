################################################################################
## Peak Cleaning Script for CLIP:
## Written by Soon Yi
## Created: 2025-08-09
## Last Edited: 2025-08-09
################################################################################

library(dplyr)
library(stringr)
library(data.table)

## Set up basic parameters:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Sequencing/COLLAPSED/CombinedPeaks/'
peakBaseName = 'Combined_peakCoverage.txt'
tagcountBaseName = 'tag_counts.txt'

sampleNames = c('hnRNPC_WT', 'hnRNPC_Mut', 'RBM25_WT', 'RBM25_WT_old', 'RBM25_WT_new', 'RBM25_Mut', 'RBM25_Mut_old', 'RBM25_Mut_new', 'hnRNPC_WT_inRBM25_WT', 'hnRNPC_WT_inRBM25_Mut')
################################################################################

## Filter peaks and generate isolated files for filtered peaks:
################################################################################
## Filters
## Peaks for eCLIP (combined), eCLIP K562, eCLIP HepG2, and HeLa hnRNPC CLIP all cleaned up already
# eCLIP: Total_BC == 4, HepG2_avg_RPM > 0.2, K562_avg_RPM > 0.2
# HepG2: HepG2_BC == 2, HepG2_avg_RPM > 0.3
# K562: K562_BC == 2, K562_avg_RPM > 0.3
# Hela: Hela >= 2, hela_avg_RPM > 2.0

## Peaks for our samples, filter with high stringency:
# BC: at least half of the samples

# BC_filter = c(1, 1, 1, 1, 1, 1)

for (sampleName in sampleNames) {
  peakTemp = data.frame(fread(paste0(baseDir, sampleName, '_', peakBaseName)))
  tagcountTemp = data.frame(fread(paste0(baseDir, sampleName, '_', tagcountBaseName)))
  rownames(tagcountTemp) = tagcountTemp$V1
  tagcountTemp$V1 = NULL
  colnames(tagcountTemp) = c('depth')
  
  sample_cols = colnames(peakTemp)[7:length(colnames(peakTemp))]
  
  peakTemp$BC = apply(peakTemp[sample_cols], 1, function(x) sum(x != 0))
  for (bed in sample_cols) {
    peakTemp[, paste0(bed, '_RPM')] = peakTemp[, bed] / tagcountTemp[bed, 'depth'] * 1e6
  }
  # peakTemp$avg_RPM = rowMeans(peakTemp[paste0(sample_cols, '_RPM')], na.rm = TRUE)
  peakTemp$avg_RPM = ifelse(peakTemp$BC > 0, rowSums(peakTemp[paste0(sample_cols, '_RPM')]) / peakTemp$BC, 0)

  
  # peakTemp = peakTemp %>% filter(BC >= BC_filter[sampleNames == sampleName])
  # peakTemp = peakTemp %>% filter(avg_RPM >= median(peakTemp$avg_RPM))
  
  peakTemp = peakTemp %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
  peakTemp = peakTemp %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                            'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
  
  write.table(peakTemp, str_replace(str_replace(paste0(baseDir, sampleName, '_', peakBaseName), '.txt', '_filtered.txt'), 'CombinedPeaks', 'FilteredPeaks'), sep = '\t', quote = F, row.names = FALSE)
  
}
################################################################################
