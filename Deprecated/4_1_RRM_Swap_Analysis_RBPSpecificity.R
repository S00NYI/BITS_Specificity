################################################################################
## Motif Analysis using RBPSpecificity Package:
## Written by Soon Yi
## Created: 2025-08-04
## Last Edited: 2025-08-04
## Figure 4
################################################################################

library(corrplot)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(GenomicRanges)
library(IRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)



library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(data.table)

library(RBPSpecificity)
library(Biostrings)
library(seqLogo)
library(VennDiagram)

## Set up basic parameters:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Sequencing/COLLAPSED/Final'
outputDir = '~/Desktop/Genomics/Specificity/AnalysisOutput/COLLAPSED_PEAKS/'
extensions = c(25)

K = 5
################################################################################

## Read files for filtered peaks:
################################################################################
## eCLIP
Peak_eCLIP = fread(paste0(baseDir, '/eCLIP_peaks_filtered.txt'))
Peak_eCLIP = Peak_eCLIP %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_eCLIP = Peak_eCLIP %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                              'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))


## K562
Peak_K562 = fread(paste0(baseDir, '/K562_peaks_filtered.txt'))
Peak_K562 = Peak_K562 %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_K562 = Peak_K562 %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                            'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))

## HepG2
Peak_HepG2 = fread(paste0(baseDir, '/HepG2_peaks_filtered.txt'))
Peak_HepG2 = Peak_HepG2 %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_HepG2 = Peak_HepG2 %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                              'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))

## Hela
Peak_HeLa = fread(paste0(baseDir, '/HeLa_peaks_filtered.txt'))
Peak_HeLa = Peak_HeLa %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_HeLa = Peak_HeLa %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                            'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))

## Our data
hnRNPC_BC = 3

Peak_hnRNPC_WT = fread(paste0(baseDir, '/hnRNPC_WT_Combined_peakCoverage_filtered_annotated.txt'))
colnames(Peak_hnRNPC_WT)[colnames(Peak_hnRNPC_WT) == "chrom"] = "chr"
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                                      'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
# sample_cols = c("hnRNPC_WT1.sorted.bed", "hnRNPC_WT2.sorted.bed", "hnRNPC_WT3.sorted.bed", "hnRNPC_WT4.sorted.bed", "hnRNPC_WT5.sorted.bed")
# rpm_cols = paste0(sample_cols, "_RPM")
# Peak_hnRNPC_WT$BC = apply(Peak_hnRNPC_WT[, ..sample_cols], 1, function(x) sum(x != 0))
# Peak_hnRNPC_WT$avg_RPM = rowMeans(Peak_hnRNPC_WT[, ..rpm_cols])
# Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% filter(hnRNPC_noRT.sorted.bed == 0)
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% filter(avg_RPM >= median(Peak_hnRNPC_WT$avg_RPM))
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% filter(BC >= hnRNPC_BC)

Peak_hnRNPC_Mut = fread(paste0(baseDir, '/hnRNPC_Mut_Combined_peakCoverage_filtered_annotated.txt'))
colnames(Peak_hnRNPC_Mut)[colnames(Peak_hnRNPC_Mut) == "chrom"] = "chr"
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                                        'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
# sample_cols = c("hnRNPC_Mut1.sorted.bed", "hnRNPC_Mut2.sorted.bed", "hnRNPC_Mut3.sorted.bed", "hnRNPC_Mut4.sorted.bed", "hnRNPC_Mut5.sorted.bed")
# rpm_cols = paste0(sample_cols, "_RPM")
# Peak_hnRNPC_Mut$BC = apply(Peak_hnRNPC_Mut[, ..sample_cols], 1, function(x) sum(x != 0))
# Peak_hnRNPC_Mut$avg_RPM = rowMeans(Peak_hnRNPC_Mut[, ..rpm_cols])
# Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% filter(hnRNPC_noRT.sorted.bed == 0)
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% filter(avg_RPM >= median(Peak_hnRNPC_Mut$avg_RPM))
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% filter(BC >= hnRNPC_BC)


RBM25_BC = 6

Peak_RBM25_WT = fread(paste0(baseDir, '/RBM25_WT_Combined_peakCoverage_filtered_annotated.txt'))
colnames(Peak_RBM25_WT)[colnames(Peak_RBM25_WT) == "chrom"] = "chr"
Peak_RBM25_WT = Peak_RBM25_WT %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_RBM25_WT = Peak_RBM25_WT %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                                    'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
# sample_cols = c("RBM25_WT1.bed", "RBM25_WT2.bed", "RBM25_WT3.bed", "RBM25_WT4.bed", "RBM25_WT5.bed")
# rpm_cols = paste0(sample_cols, "_RPM")
# Peak_RBM25_WT$BC = apply(Peak_RBM25_WT[, ..sample_cols], 1, function(x) sum(x != 0))
# Peak_RBM25_WT$avg_RPM = rowMeans(Peak_RBM25_WT[, ..rpm_cols])
# Peak_RBM25_WT = Peak_RBM25_WT %>% filter(RBM25_noRT.sorted.bed == 0)
Peak_RBM25_WT = Peak_RBM25_WT %>% filter(avg_RPM >= median(Peak_RBM25_WT$avg_RPM))
Peak_RBM25_WT = Peak_RBM25_WT %>% filter(BC >= RBM25_BC)


Peak_RBM25_Mut = fread(paste0(baseDir, '/RBM25_Mut_Combined_peakCoverage_filtered_annotated.txt'))
colnames(Peak_RBM25_Mut)[colnames(Peak_RBM25_Mut) == "chrom"] = "chr"
Peak_RBM25_Mut = Peak_RBM25_Mut %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_RBM25_Mut = Peak_RBM25_Mut %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                                      'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
# sample_cols = c("RBM25_Mut1.bed", "RBM25_Mut2.bed", "RBM25_Mut3.bed", "RBM25_Mut4.bed", "RBM25_Mut5.bed")
# rpm_cols = paste0(sample_cols, "_RPM")
# Peak_RBM25_Mut$BC = apply(Peak_RBM25_Mut[, ..sample_cols], 1, function(x) sum(x != 0))
# Peak_RBM25_Mut$avg_RPM = rowMeans(Peak_RBM25_Mut[, ..rpm_cols])
# Peak_RBM25_Mut = Peak_RBM25_Mut %>% filter(RBM25_noRT.sorted.bed == 0)
Peak_RBM25_Mut = Peak_RBM25_Mut %>% filter(avg_RPM >= median(Peak_RBM25_Mut$avg_RPM))
Peak_RBM25_Mut = Peak_RBM25_Mut %>% filter(BC >= RBM25_BC)


# RBM25_BC_old = 2
#
# Peak_RBM25_WT_old = fread(paste0(baseDir, '/RBM25_WT_old_Combined_peakCoverage_filtered_annotated.txt'))
# colnames(Peak_RBM25_WT_old)[colnames(Peak_RBM25_WT_old) == "chrom"] = "chr"
# Peak_RBM25_WT_old = Peak_RBM25_WT_old %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
# Peak_RBM25_WT_old = Peak_RBM25_WT_old %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
#                                                             'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
# # sample_cols = c("RBM25_WT1.bed", "RBM25_WT2.bed", "RBM25_WT3.bed", "RBM25_WT4.bed", "RBM25_WT5.bed")
# # rpm_cols = paste0(sample_cols, "_RPM")
# # Peak_RBM25_WT_old$BC = apply(Peak_RBM25_WT_old[, ..sample_cols], 1, function(x) sum(x != 0))
# # Peak_RBM25_WT_old$avg_RPM = rowMeans(Peak_RBM25_WT_old[, ..rpm_cols])
# # Peak_RBM25_WT_old = Peak_RBM25_WT_old %>% filter(RBM25_noRT.sorted.bed == 0)
# Peak_RBM25_WT_old = Peak_RBM25_WT_old %>% filter(avg_RPM >= median(Peak_RBM25_WT_old$avg_RPM))
# Peak_RBM25_WT_old = Peak_RBM25_WT_old %>% filter(BC >= RBM25_BC_old)
#
# Peak_RBM25_Mut_old = fread(paste0(baseDir, '/RBM25_Mut_old_Combined_peakCoverage_filtered_annotated.txt'))
# colnames(Peak_RBM25_Mut_old)[colnames(Peak_RBM25_Mut_old) == "chrom"] = "chr"
# Peak_RBM25_Mut_old = Peak_RBM25_Mut_old %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
# Peak_RBM25_Mut_old = Peak_RBM25_Mut_old %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
#                                                               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
# # sample_cols = c("RBM25_Mut1.bed", "RBM25_Mut2.bed", "RBM25_Mut3.bed", "RBM25_Mut4.bed", "RBM25_Mut5.bed")
# # rpm_cols = paste0(sample_cols, "_RPM")
# # Peak_RBM25_Mut_old$BC = apply(Peak_RBM25_Mut_old[, ..sample_cols], 1, function(x) sum(x != 0))
# # Peak_RBM25_Mut_old$avg_RPM = rowMeans(Peak_RBM25_Mut_old[, ..rpm_cols])
# # Peak_RBM25_Mut_old = Peak_RBM25_Mut_old %>% filter(RBM25_noRT.sorted.bed == 0)
# Peak_RBM25_Mut_old = Peak_RBM25_Mut_old %>% filter(avg_RPM >= median(Peak_RBM25_Mut_old$avg_RPM))
# Peak_RBM25_Mut_old = Peak_RBM25_Mut_old %>% filter(BC >= RBM25_BC_old)
#
#
# RBM25_BC_new = 2
#
# Peak_RBM25_WT_new = fread(paste0(baseDir, '/RBM25_WT_new_Combined_peakCoverage_filtered_annotated.txt'))
# colnames(Peak_RBM25_WT_new)[colnames(Peak_RBM25_WT_new) == "chrom"] = "chr"
# Peak_RBM25_WT_new = Peak_RBM25_WT_new %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
# Peak_RBM25_WT_new = Peak_RBM25_WT_new %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
#                                                             'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
# # sample_cols = c("RBM25_WT1.bed", "RBM25_WT2.bed", "RBM25_WT3.bed", "RBM25_WT4.bed", "RBM25_WT5.bed")
# # rpm_cols = paste0(sample_cols, "_RPM")
# # Peak_RBM25_WT_new$BC = apply(Peak_RBM25_WT_new[, ..sample_cols], 1, function(x) sum(x != 0))
# # Peak_RBM25_WT_new$avg_RPM = rowMeans(Peak_RBM25_WT_new[, ..rpm_cols])
# # Peak_RBM25_WT_new = Peak_RBM25_WT_new %>% filter(RBM25_noRT.sorted.bed == 0)
# Peak_RBM25_WT_new = Peak_RBM25_WT_new %>% filter(avg_RPM >= median(Peak_RBM25_WT_new$avg_RPM))
# Peak_RBM25_WT_new = Peak_RBM25_WT_new %>% filter(BC >= RBM25_BC_new)
#
# Peak_RBM25_Mut_new = fread(paste0(baseDir, '/RBM25_Mut_new_Combined_peakCoverage_filtered_annotated.txt'))
# colnames(Peak_RBM25_Mut_new)[colnames(Peak_RBM25_Mut_new) == "chrom"] = "chr"
# Peak_RBM25_Mut_new = Peak_RBM25_Mut_new %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
# Peak_RBM25_Mut_new = Peak_RBM25_Mut_new %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
#                                                               'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
# # sample_cols = c("RBM25_Mut1.bed", "RBM25_Mut2.bed", "RBM25_Mut3.bed", "RBM25_Mut4.bed", "RBM25_Mut5.bed")
# # rpm_cols = paste0(sample_cols, "_RPM")
# # Peak_RBM25_Mut_new$BC = apply(Peak_RBM25_Mut_new[, ..sample_cols], 1, function(x) sum(x != 0))
# # Peak_RBM25_Mut_new$avg_RPM = rowMeans(Peak_RBM25_Mut_new[, ..rpm_cols])
# # Peak_RBM25_Mut_new = Peak_RBM25_Mut_new %>% filter(RBM25_noRT.sorted.bed == 0)
# Peak_RBM25_Mut_new = Peak_RBM25_Mut_new %>% filter(avg_RPM >= median(Peak_RBM25_Mut_new$avg_RPM))
# Peak_RBM25_Mut_new = Peak_RBM25_Mut_new %>% filter(BC >= RBM25_BC_new)

################################################################################

## Stacked Bar
################################################################################
## Get Annotation counts:
countAnnotation = function(peak_matrix, annotation_column, new_column_name = NULL, annotation_to_skip = NULL, fraction = NULL) {
  temp = data.frame(table(peak_matrix[, ..annotation_column]), row.names = 1)
  if(!is.null(new_column_name)) {
    colnames(temp) = new_column_name
  }

  if(!is.null(annotation_to_skip)) {
    temp = temp[rownames(temp) != annotation_to_skip, , drop = FALSE]
  }

  if(!is.null(fraction)) {
    temp = temp/sum(temp)
  }

  return(temp)
}

## Fill Annotation counts if anything is mixing:
fillAnnotation = function(annotation_counts, annotation_list) {
  colnames((annotation_counts))

  temp = data.frame(Sample = numeric(length(annotation_list)))
  rownames(temp) = annotation_list
  temp2 = merge(temp, annotation_counts, by = "row.names", all = TRUE)
  temp2[is.na(temp2)] = 0
  temp2$Sample = NULL

  rownames(temp2) = temp2$Row.names
  temp2 = temp2[annotation_list, -1, drop = FALSE]

  return(temp2)
}

## Plot Stacked bar:
plotStackedBar = function(annotation_counts, sample_list, sample_label, title, y_lim = NULL, y_tick = NULL) {
  plot = ggplot(annotation_counts %>% filter(Source %in% sample_list), aes(fill = Annotation, y=Freq, x=Source)) +
    geom_bar(position='stack', stat='identity') +
    scale_x_discrete(labels = sample_label) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Set3") +
    theme_bw() +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14, face = 'bold'),
          legend.text = element_text(size=14))

  if(!is.null(y_lim)) {
    plot = plot + ylim(y_lim)
  }

  if (!is.null(y_tick)) {
    plot = plot + scale_y_continuous(breaks = seq(0, y_lim[2], by=y_tick), limits=c(0, y_lim[2]))
  }

  return(plot)
}


CLIP_List = c('eCLIP', 'HeLa', 'hnRNPC_WT', 'hnRNPC_Mut', 'RBM25_WT', 'RBM25_Mut')
All_Annotation_List = c("5'UTR", "CDS", "3'UTR", "intron", 'rRNA', 'tRNA', 'ncRNA', 'TE', "Other")

Peak_eCLIP = Peak_eCLIP %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_eCLIP = Peak_eCLIP %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_HeLa = Peak_HeLa %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_HeLa = Peak_HeLa %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_RBM25_WT = Peak_RBM25_WT %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_RBM25_WT = Peak_RBM25_WT %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_RBM25_Mut = Peak_RBM25_Mut %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_RBM25_Mut = Peak_RBM25_Mut %>% filter(finalized_annotation %in% All_Annotation_List)

# Peak_RBM25_WT_old = Peak_RBM25_WT_old %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
# Peak_RBM25_WT_old = Peak_RBM25_WT_old %>% filter(finalized_annotation %in% All_Annotation_List)
#
# Peak_RBM25_Mut_old = Peak_RBM25_Mut_old %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
# Peak_RBM25_Mut_old = Peak_RBM25_Mut_old %>% filter(finalized_annotation %in% All_Annotation_List)
#
# Peak_RBM25_WT_new = Peak_RBM25_WT_new %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
# Peak_RBM25_WT_new = Peak_RBM25_WT_new %>% filter(finalized_annotation %in% All_Annotation_List)
#
# Peak_RBM25_Mut_new = Peak_RBM25_Mut_new %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
# Peak_RBM25_Mut_new = Peak_RBM25_Mut_new %>% filter(finalized_annotation %in% All_Annotation_List)

Counts_Peak_eCLIP = countAnnotation(Peak_eCLIP, 'finalized_annotation', 'eCLIP', fraction = TRUE)
Counts_Peak_eCLIP['tRNA', 'eCLIP'] = 0
Counts_Peak_HeLa = countAnnotation(Peak_HeLa, 'finalized_annotation', 'HeLa', fraction = TRUE)
Counts_Peak_hnRNPC_WT = countAnnotation(Peak_hnRNPC_WT, 'finalized_annotation', 'hnRNPC_WT', fraction = TRUE)
Counts_Peak_hnRNPC_Mut = countAnnotation(Peak_hnRNPC_Mut, 'finalized_annotation', 'hnRNPC_Mut', fraction = TRUE)
Counts_Peak_RBM25_WT = countAnnotation(Peak_RBM25_WT, 'finalized_annotation', 'RBM25_WT', fraction = TRUE)
Counts_Peak_RBM25_Mut = countAnnotation(Peak_RBM25_Mut, 'finalized_annotation', 'RBM25_Mut', fraction = TRUE)
# Counts_Peak_RBM25_WT_old = countAnnotation(Peak_RBM25_WT_old, 'finalized_annotation', 'RBM25_WT_old', fraction = TRUE)
# Counts_Peak_RBM25_Mut_old = countAnnotation(Peak_RBM25_Mut_old, 'finalized_annotation', 'RBM25_Mut_old', fraction = TRUE)
# Counts_Peak_RBM25_WT_new = countAnnotation(Peak_RBM25_WT_new, 'finalized_annotation', 'RBM25_WT_new', fraction = TRUE)
# Counts_Peak_RBM25_Mut_new = countAnnotation(Peak_RBM25_Mut_new, 'finalized_annotation', 'RBM25_Mut_new', fraction = TRUE)

Counts_Peak_hnRNPC_WT = fillAnnotation(Counts_Peak_hnRNPC_WT, All_Annotation_List)
Counts_Peak_hnRNPC_Mut = fillAnnotation(Counts_Peak_hnRNPC_Mut, All_Annotation_List)
Counts_Peak_RBM25_WT = fillAnnotation(Counts_Peak_RBM25_WT, All_Annotation_List)
Counts_Peak_RBM25_Mut = fillAnnotation(Counts_Peak_RBM25_Mut, All_Annotation_List)
# Counts_Peak_RBM25_WT_old = fillAnnotation(Counts_Peak_RBM25_WT_old, All_Annotation_List)
# Counts_Peak_RBM25_Mut_old = fillAnnotation(Counts_Peak_RBM25_Mut_old, All_Annotation_List)
# Counts_Peak_RBM25_WT_new = fillAnnotation(Counts_Peak_RBM25_WT_new, All_Annotation_List)
# Counts_Peak_RBM25_Mut_new = fillAnnotation(Counts_Peak_RBM25_Mut_new, All_Annotation_List)

PeakDistribution_combined = cbind(Counts_Peak_eCLIP, Counts_Peak_HeLa,
                                  Counts_Peak_hnRNPC_WT, Counts_Peak_hnRNPC_Mut,
                                  Counts_Peak_RBM25_WT, Counts_Peak_RBM25_Mut)
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>% gather(key = "Source", value = "Freq", CLIP_List) %>% dplyr::select(Source, Freq, Annotation)
PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = CLIP_List)
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = All_Annotation_List)

plotStackedBar(PeakDistribution_combined, CLIP_List, CLIP_List, 'CLIP')
################################################################################

## Motif Analysis
################################################################################
DATE = format(Sys.Date(), "%Y%m%d")
if (!dir.exists(paste0(outputDir, DATE, '/'))) {
  dir.create(paste0(outputDir, DATE, '/'), recursive = TRUE, showWarnings = FALSE)
  print(paste("Created directory:", paste0(outputDir, DATE, '/')))
}

samplesList = c('eCLIP_hnRNPC_WT',
                'K562_hnRNP_WT', 'HepG2_hnRNPC_WT',
                'HeLa_hnRNPC_WT',
                'HEK293_hnRNPC_WT', 'HEK293_hnRNPC_Mut',
                'HEK293_RBM25_WT', 'HEK293_RBM25_Mut')

CLIP_NAME = c()
NT_EXT = c()

UUUUU_IS = c()
UUUUU_Rank = c()

GGGGG_IS = c()
GGGGG_Rank = c()

TOP_MOTIF = c()


for (extension in extensions) {
  print(paste0('     Working on ', extension, 'nt extension....'))
  for (sampleName in samplesList) {
    # Read peaks and make GR object
    print(paste0('Working on ', sampleName, '....'))
    print(paste0('     Reading and processing peak data....'))

    if (sampleName == 'eCLIP_hnRNPC_WT') {
      Peak = Peak_eCLIP
    } else if (sampleName == 'K562_hnRNPC_WT') {
      Peak = Peak_K562
    } else if (sampleName == 'HepG2_hnRNPC_WT') {
      Peak = Peak_HepG2
    } else if (sampleName == 'HeLa_hnRNPC_WT') {
      Peak = Peak_HeLa
    } else if (sampleName == 'HEK293_hnRNPC_WT') {
      Peak = Peak_hnRNPC_WT
    } else if (sampleName == 'HEK293_hnRNPC_Mut') {
      Peak = Peak_hnRNPC_Mut
    } else if (sampleName == 'HEK293_RBM25_WT') {
      Peak = Peak_RBM25_WT
    } else if (sampleName == 'HEK293_RBM25_Mut') {
      Peak = Peak_RBM25_Mut
    } else if (sampleName == 'HEK293_RBM25_WT_old') {
      Peak = Peak_RBM25_WT_old
    } else if (sampleName == 'HEK293_RBM25_Mut_old') {
      Peak = Peak_RBM25_Mut_old
    } else if (sampleName == 'HEK293_RBM25_WT_new') {
      Peak = Peak_RBM25_WT_new
    } else if (sampleName == 'HEK293_RBM25_Mut_new') {
      Peak = Peak_RBM25_Mut_new
    }

    CLIP_NAME = c(CLIP_NAME, sampleName)
    NT_EXT = c(NT_EXT, extension)

    enrichment = motifEnrichment(peak_data = Peak, 'hg38', K,
                                 extension = extension,
                                 Bkg_number = 100, Bkg_dist = 500, max_shift_dist = 1000,
                                 log_transform = FALSE)
    IS = returnIS(enrichment, return_type = 'all')
    MS = returnMS(enrichment, return_type = 'all')

    UUUUU_IS = c(UUUUU_IS, returnIS(enrichment, motif = 'TTTTT'))
    GGGGG_IS = c(GGGGG_IS, returnIS(enrichment, motif = 'GGGGG'))
    TOP_MOTIF = c(TOP_MOTIF, enrichment$MOTIF[enrichment$Score == max(enrichment$Score)])

    ordered_IS = IS[order(-IS$IS), ]
    row.names(ordered_IS) = NULL
    UUUUU_Rank = c(UUUUU_Rank, as.numeric(rownames(ordered_IS[ordered_IS$MOTIF=='TTTTT', ])))
    GGGGG_Rank = c(GGGGG_Rank, as.numeric(rownames(ordered_IS[ordered_IS$MOTIF=='GGGGG', ])))

    write.csv(enrichment,
              paste0(outputDir, DATE, '/', sampleName, '_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), quote = F)
    write.csv(IS,
              paste0(outputDir, DATE, '/', sampleName, '_ISperMotif_', extension, 'ntExt_', as.character(K), 'mer.csv'), quote = F)
    write.csv(MS,
              paste0(outputDir, DATE, '/', sampleName, '_MSperMotif_', extension, 'ntExt_', as.character(K), 'mer.csv'), quote = F)


  }
  TopMotif_Data_Compiled = data.frame(CLIP = CLIP_NAME,
                                      EXTN = NT_EXT,
                                      TOP = TOP_MOTIF,
                                      UUUUU_RANK = UUUUU_Rank,
                                      UUUUU_IS = UUUUU_IS,
                                      GGGGG_RANK = GGGGG_Rank,
                                      GGGGG_IS = GGGGG_IS)

  write.csv(TopMotif_Data_Compiled,
            paste0(outputDir, DATE, '/', K, 'mer_RBMSwap_Analysis_', extension, 'ntExt.csv'), quote = F)

}



################################################################################

## PWM
################################################################################
DATE = '20250821'
extension = 25


hnRNPC_WT = read_csv(paste0(outputDir, DATE, '/', 'HEK293_hnRNPC_WT_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)
hnRNPC_Mut = read_csv(paste0(outputDir, DATE, '/', 'HEK293_hnRNPC_Mut_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)
RBM25_WT = read_csv(paste0(outputDir, DATE, '/', 'HEK293_RBM25_WT_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)
RBM25_Mut = read_csv(paste0(outputDir, DATE, '/', 'HEK293_RBM25_Mut_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)

numMotif = 10

## For hnRNPC WT
data = hnRNPC_WT
data = data.frame(data[, c('MOTIF', 'Score')])
data = data[order(-data[, 'Score']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## For hnRNPC Mut
data = hnRNPC_Mut
data = data.frame(data[, c('MOTIF', 'Score')])
data = data[order(-data[, 'Score']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## For RBM25 WT
data = RBM25_WT
data = data.frame(data[, c('MOTIF', 'Score')])
data = data[order(-data[, 'Score']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## For RBM25 Mut
data = RBM25_Mut
data = data.frame(data[, c('MOTIF', 'Score')])
data = data[order(-data[, 'Score']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)
################################################################################

## Histogram
################################################################################
ggplot() +
  geom_histogram(data = hnRNPC_WT, aes(x = Score), binwidth = 0.01, fill = "darkgoldenrod1", alpha = 0.5) +
  geom_histogram(data = hnRNPC_Mut, aes(x = Score), binwidth = 0.01, fill = "cornflowerblue", alpha = 0.5) +
  labs(title = "hnRNPC", x = "Normalized Counts", y = "Frequency") +
  scale_y_continuous(limits = c(0, 350), breaks = seq(0, 550, by = 50)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))

ggplot() +
  geom_histogram(data = RBM25_WT, aes(x = Score), binwidth = 0.01, fill = "darkgoldenrod1", alpha = 0.5) +
  geom_histogram(data = RBM25_Mut, aes(x = Score), binwidth = 0.01, fill = "cornflowerblue", alpha = 0.5) +
  labs(title = "RBM25", x = "Normalized Counts", y = "Frequency") +
  scale_y_continuous(limits = c(0, 350), breaks = seq(0, 500, by = 50)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))
################################################################################

## MS Analysis
################################################################################

plotMS(hnRNPC_WT, output_type = 'matrix')
plotMS(hnRNPC_Mut, output_type = 'matrix')
plotMS(RBM25_WT, output_type = 'matrix')
plotMS(RBM25_Mut, output_type = 'matrix')


plotMS(RBM25_WT, motif = 'TTTTT', output_type = 'matrix')
plotMS(RBM25_Mut, motif = 'TTTTT', output_type = 'matrix')


returnMS(hnRNPC_WT, motif = 'TTTTT', output_type = 'number')    # 0.84
returnMS(hnRNPC_Mut, motif = 'TTTTT',output_type = 'number')    # 0.81
returnMS(RBM25_WT, motif = 'TTTTT',output_type = 'number')      # 0.23
returnMS(RBM25_Mut, motif = 'TTTTT',output_type = 'number')     # 0.46


returnIS(hnRNPC_WT, motif = 'TTTTT')     # 52.6
returnIS(hnRNPC_Mut, motif = 'TTTTT')    # 50.5
returnIS(RBM25_WT, motif = 'TTTTT')      # 1.82
returnIS(RBM25_Mut, motif = 'TTTTT')     # 2.46
returnIS(RBM25_Mut, motif = 'GGGGG')     # 2.02



returnIS(RBM25_WT)      # 2.35
returnIS(RBM25_Mut)     # 2.46

################################################################################

## Gene Level Analysis
################################################################################
hnRNPC_targetGenes_WT = unique(Peak_hnRNPC_WT$gene)
hnRNPC_targetGenes_Mut = unique(Peak_hnRNPC_Mut$gene)

hnRNPC_targetGenes_venn = venn.diagram(list(WT = hnRNPC_targetGenes_WT, Mut = hnRNPC_targetGenes_Mut),
                                       category.names = c("hnRNPC WT", "hnRNPC Mut"),
                                       filename = NULL,
                                       output = TRUE,
                                       imagetype = "png",
                                       fill = c("skyblue", "salmon"),
                                       col = "transparent",
                                       alpha = 0.5,
                                       cex = 1.5,
                                       cat.cex = 1.2,
                                       cat.pos = c(-20, 20),
                                       cat.dist = c(0.05, 0.05))
grid.newpage()
grid.draw(hnRNPC_targetGenes_venn)

RBM25_targetGenes_WT = unique(Peak_RBM25_WT$gene)
RBM25_targetGenes_Mut = unique(Peak_RBM25_Mut$gene)

RBM25_targetGenes_venn = venn.diagram(list(WT = RBM25_targetGenes_WT, Mut = RBM25_targetGenes_Mut),
                                       category.names = c("RBM25 WT", "RBM25 Mut"),
                                       filename = NULL,
                                       output = TRUE,
                                       imagetype = "png",
                                       fill = c("skyblue", "salmon"),
                                       col = "transparent",
                                       alpha = 0.5,
                                       cex = 1.5,
                                       cat.cex = 1.2,
                                       cat.pos = c(-20, 20),
                                       cat.dist = c(0.05, 0.05))
grid.newpage()
grid.draw(RBM25_targetGenes_venn)

hnRNPCWT_RBM25MWT_targetGenes_venn = venn.diagram(list(WT = hnRNPC_targetGenes_WT, Mut = RBM25_targetGenes_WT),
                                                  category.names = c("hnRNPC WT", "RBM25 WT"),
                                                  filename = NULL,
                                                  output = TRUE,
                                                  imagetype = "png",
                                                  fill = c("skyblue", "salmon"),
                                                  col = "transparent",
                                                  alpha = 0.5,
                                                  cex = 1.5,
                                                  cat.cex = 1.2,
                                                  cat.pos = c(-20, 20),
                                                  cat.dist = c(0.05, 0.05))
grid.newpage()
grid.draw(hnRNPCWT_RBM25MWT_targetGenes_venn)

hnRNPCWT_RBM25Mut_targetGenes_venn = venn.diagram(list(WT = hnRNPC_targetGenes_WT, Mut = RBM25_targetGenes_Mut),
                                       category.names = c("hnRNPC WT", "RBM25 Mut"),
                                       filename = NULL,
                                       output = TRUE,
                                       imagetype = "png",
                                       fill = c("skyblue", "salmon"),
                                       col = "transparent",
                                       alpha = 0.5,
                                       cex = 1.5,
                                       cat.cex = 1.2,
                                       cat.pos = c(-20, 20),
                                       cat.dist = c(0.05, 0.05))
grid.newpage()
grid.draw(hnRNPCWT_RBM25Mut_targetGenes_venn)

hnRNPCMut_RBM25WT_targetGenes_venn = venn.diagram(list(WT = hnRNPC_targetGenes_Mut, Mut = RBM25_targetGenes_WT),
                                                  category.names = c("hnRNPC Mut", "RBM25 WT"),
                                                  filename = NULL,
                                                  output = TRUE,
                                                  imagetype = "png",
                                                  fill = c("skyblue", "salmon"),
                                                  col = "transparent",
                                                  alpha = 0.5,
                                                  cex = 1.5,
                                                  cat.cex = 1.2,
                                                  cat.pos = c(-20, 20),
                                                  cat.dist = c(0.05, 0.05))
grid.newpage()
grid.draw(hnRNPCMut_RBM25WT_targetGenes_venn)

hnRNPC_3way_targetGenes_venn = venn.diagram(list(WT = hnRNPC_targetGenes_WT, Mut = hnRNPC_targetGenes_Mut, Mut2 = RBM25_targetGenes_Mut),
                                                  category.names = c("hnRNPC WT", "hnRNPC Mut", "RBM25 Mut"),
                                                  filename = NULL,
                                                  output = TRUE,
                                                  imagetype = "png",
                                                  fill = c("skyblue", "salmon", 'gold'))
grid.newpage()
grid.draw(hnRNPC_3way_targetGenes_venn)

hnRNPC_targetGenes_overlap = intersect(hnRNPC_targetGenes_WT, hnRNPC_targetGenes_Mut)

targetGenes_RPM_inhnRNPCWT = Peak_hnRNPC_WT %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(hnRNPCWT = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM_inhnRNPCMut = Peak_hnRNPC_Mut %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(hnRNPCMut = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM = inner_join(targetGenes_RPM_inhnRNPCWT, targetGenes_RPM_inhnRNPCMut, by = "gene")
targetGenes_RPM =  targetGenes_RPM[!(targetGenes_RPM$gene == ""), ]

library(biomaRt)
library(ggsignif)
ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids = unique(targetGenes_RPM$gene)

gene_coords = getBM(attributes = c('ensembl_gene_id', 'start_position', 'end_position'),
                    filters = 'ensembl_gene_id',
                    values = gene_ids,
                    mart = ensembl)

gene_lengths = gene_coords %>%
  mutate(length = end_position - start_position) %>% # Calculate length
  group_by(ensembl_gene_id) %>%                      # Group by gene ID
  summarise(gene_length = max(length))               # Find the max length for each group

targetGenes_RPM = left_join(targetGenes_RPM,
                            gene_lengths,
                            by = c("gene" = "ensembl_gene_id"))

targetGenes_RPM$hnRNPCWT_RPKM = targetGenes_RPM$hnRNPCWT/targetGenes_RPM$gene_length * 1000
targetGenes_RPM$hnRNPCMut_RPKM = targetGenes_RPM$hnRNPCMut/targetGenes_RPM$gene_length * 1000

p_val = wilcox.test(targetGenes_RPM$hnRNPCWT_RPKM, targetGenes_RPM$hnRNPCMut_RPKM)
p_val = scales::pvalue(p_val$p.value)

targetGenes_long = pivot_longer(targetGenes_RPM,
                                cols = c("hnRNPCWT_RPKM", "hnRNPCMut_RPKM"),
                                names_to = "condition",
                                values_to = "RPKM")

ggplot(targetGenes_long, aes(x = condition, y = RPKM, fill = condition)) +
  geom_boxplot(outliers = TRUE, notch = TRUE) +
  scale_y_log10() +
  geom_signif(comparisons = list(c('hnRNPCWT_RPKM', 'hnRNPCMut_RPKM')),
              annotations = p_val,
              y_position = log10(5000),
              tip_length = 0.01,
              map_signif_level = TRUE) +
  labs(title = "RPKM for hnRNPC target genes",
       x = "Cell Lines",
       y = "RPKM per Gene") +
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))

ggplot(targetGenes_RPM, aes(x = hnRNPCWT_RPKM, y = hnRNPCMut_RPKM)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "RPKM per hnRNPC target genes in WT vs. in Mut)",
       x = "Log(RPKM) in hnRNPC WT",
       y = "Log(RPKM) in hnRNPC Mut") +
  coord_fixed() + # Ensure a 1:1 aspect ratio for the axes
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))



## for RBM25
RBM25_targetGenes_overlap = intersect(RBM25_targetGenes_WT, RBM25_targetGenes_Mut)

targetGenes_RPM_inRBM25WT = Peak_RBM25_WT %>%
  filter(gene %in% RBM25_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(RBM25WT = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM_inRBM25Mut = Peak_RBM25_Mut %>%
  filter(gene %in% RBM25_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(RBM25Mut = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM = inner_join(targetGenes_RPM_inRBM25WT, targetGenes_RPM_inRBM25Mut, by = "gene")
targetGenes_RPM =  targetGenes_RPM[!(targetGenes_RPM$gene == ""), ]

ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids = unique(targetGenes_RPM$gene)

gene_coords = getBM(attributes = c('ensembl_gene_id', 'start_position', 'end_position'),
                    filters = 'ensembl_gene_id',
                    values = gene_ids,
                    mart = ensembl)

gene_lengths = gene_coords %>%
  mutate(length = end_position - start_position) %>% # Calculate length
  group_by(ensembl_gene_id) %>%                      # Group by gene ID
  summarise(gene_length = max(length))               # Find the max length for each group

targetGenes_RPM = left_join(targetGenes_RPM,
                            gene_lengths,
                            by = c("gene" = "ensembl_gene_id"))

targetGenes_RPM$RBM25WT_RPKM = targetGenes_RPM$RBM25WT/targetGenes_RPM$gene_length * 1000
targetGenes_RPM$RBM25Mut_RPKM = targetGenes_RPM$RBM25Mut/targetGenes_RPM$gene_length * 1000

p_val = wilcox.test(targetGenes_RPM$RBM25WT_RPKM, targetGenes_RPM$RBM25Mut_RPKM)
p_val = scales::pvalue(p_val$p.value)

targetGenes_long = pivot_longer(targetGenes_RPM,
                                cols = c("RBM25WT_RPKM", "RBM25Mut_RPKM"),
                                names_to = "condition",
                                values_to = "RPKM")

ggplot(targetGenes_long, aes(x = condition, y = RPKM, fill = condition)) +
  geom_boxplot(outliers = TRUE, notch = TRUE) +
  scale_y_log10() +
  geom_signif(comparisons = list(c('RBM25WT_RPKM', 'RBM25Mut_RPKM')),
              annotations = p_val,
              y_position = log10(5000),
              tip_length = 0.01,
              map_signif_level = TRUE) +
  labs(title = "RPKM for RBM25 target genes",
       x = "Cell Lines",
       y = "RPKM per Gene") +
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))

ggplot(targetGenes_RPM, aes(x = RBM25WT_RPKM, y = RBM25Mut_RPKM)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "RPKM per RBM25 target genes in WT vs. in Mut)",
       x = "Log(RPKM) in RBM25 WT",
       y = "Log(RPKM) in RBM25 Mut") +
  coord_fixed() + # Ensure a 1:1 aspect ratio for the axes
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))
################################################################################



## Peak Scatter Plot
################################################################################

peaks2GR = function(peaks_df) {
  peaks_df = peaks_df[, c('chr', 'start', 'end', 'strand')]
  peaks_df$start = as.integer(peaks_df$start)
  peaks_df$end = as.integer(peaks_df$end)
  peaks_GR = GRanges(peaks_df)

  return(peaks_GR)
}

GR_Peak_eCLIP = peaks2GR(Peak_eCLIP)
GR_Peak_HeLa = peaks2GR(Peak_HeLa)
GR_Peak_hnRNPC_WT = peaks2GR(Peak_hnRNPC_WT)
GR_Peak_hnRNPC_Mut = peaks2GR(Peak_hnRNPC_Mut)
GR_Peak_RBM25_WT = peaks2GR(Peak_RBM25_WT)
GR_Peak_RBM25_Mut = peaks2GR(Peak_RBM25_Mut)

# hnRNPC_WT vs hnRNPC_Mut
overlaps = as.data.frame(findOverlaps(GR_Peak_hnRNPC_WT, GR_Peak_hnRNPC_Mut, minoverlap = 15))

overlaps_hnRNPC_WT = Peak_hnRNPC_WT[overlaps$queryHits]
overlaps_hnRNPC_Mut = Peak_hnRNPC_Mut[overlaps$subjectHits]

overlaps_RPM = data.frame(hnRNPC_WT = overlaps_hnRNPC_WT$avg_RPM,
                          hnRNPC_Mut = overlaps_hnRNPC_Mut$avg_RPM)

RMSDist = sqrt(mean((overlaps_RPM$hnRNPC_Mut - overlaps_RPM$hnRNPC_WT)^2)) # for magnitude
MDist = mean(overlaps_RPM$hnRNPC_WT - overlaps_RPM$hnRNPC_Mut) # for direction
L2FC = log2(mean(overlaps_RPM$hnRNPC_Mut/overlaps_RPM$hnRNPC_WT))

# Calculate the *perpendicular* distance from the y=x line (as before)
custom_color_function = function(dist, pos) {
  if (dist < 0.001) { # Threshold for "near" (adjust as needed)
    return("black")
  } else {
    if (pos == "Above") {
      # Use pmax to avoid 0 index
      index <- pmax(1, pmin(100, round(dist * 100))) # Correct indexing
      return(colorRampPalette(c("black", "#ed6677"))(100)[index])
    } else {
      index <- pmax(1, pmin(100, round(dist * 100))) # Correct indexing
      return(colorRampPalette(c("black", "#258942"))(100)[index])
    }
  }
}

overlaps_RPM$distance = ((abs(overlaps_RPM$hnRNPC_Mut - overlaps_RPM$hnRNPC_WT) / sqrt(2)))
overlaps_RPM$position = ifelse(overlaps_RPM$hnRNPC_Mut > overlaps_RPM$hnRNPC_WT, "Above", "Below")
overlaps_RPM$color = mapply(custom_color_function, overlaps_RPM$distance, overlaps_RPM$position)

ggplot(overlaps_RPM,
       aes(x = hnRNPC_WT, y = hnRNPC_Mut, color = color)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_color_identity() +  # Use the pre-calculated colors
  scale_x_log10(limits = c(1, 1000)) +
  scale_y_log10(limits = c(1, 1000)) +
  labs(title = paste0("hnRNPC WT vs Mut RPM Comparison (N: ", (nrow(overlaps_RPM)), ", RMSD: ", (RMSDist), ", MD: ", (MDist), ")"),
       x = "hnRNPC WT Average RPM",
       y = "hnRNPC Mut Average RPM") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))

overlaps = as.data.frame(findOverlaps(GR_Peak_RBM25_WT, GR_Peak_RBM25_Mut))
overlaps_RBM25_WT = Peak_RBM25_WT[overlaps$queryHits]
overlaps_RBM25_Mut = Peak_RBM25_Mut[overlaps$subjectHits]

overlaps_RPM = data.frame(RBM25_WT = overlaps_RBM25_WT$avg_RPM,
                          RBM25_Mut = overlaps_RBM25_Mut$avg_RPM)

RMSDist = sqrt(mean((overlaps_RPM$RBM25_Mut - overlaps_RPM$RBM25_WT)^2)) # for magnitude
MDist = mean(overlaps_RPM$RBM25_WT - overlaps_RPM$RBM25_Mut) # for direction
L2FC = log2(mean(overlaps_RPM$RBM25_Mut/overlaps_RPM$RBM25_WT))

overlaps_RPM$distance = (abs(overlaps_RPM$RBM25_Mut - overlaps_RPM$RBM25_WT) / sqrt(2))
overlaps_RPM$position = ifelse(overlaps_RPM$RBM25_Mut > overlaps_RPM$RBM25_WT, "Above", "Below")
overlaps_RPM$color = mapply(custom_color_function, overlaps_RPM$distance, overlaps_RPM$position)

ggplot(overlaps_RPM,
       aes(x = RBM25_WT, y = RBM25_Mut, color = color)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_color_identity() +  # Use the pre-calculated colors
  scale_x_log10(limits = c(1, 1000)) +
  scale_y_log10(limits = c(1, 1000)) +
  labs(title = paste0("RBM25 WT vs Mut RPM Comparison (N: ", (nrow(overlaps_RPM)), ", RMSD: ", (RMSDist), ", MD: ", (MDist), ")"),
       x = "RBM25 WT Average RPM",
       y = "RBM25 Mut Average RPM") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))

################################################################################

## IS Ratio Analysis
################################################################################
hnRNPC_Mut2WT_ISRatio = data.frame(MOTIF = hnRNPC_WT$MOTIF,
                                   WT_IS = returnIS(hnRNPC_WT, return_type = 'all')$IS,
                                   Mut_IS = returnIS(hnRNPC_Mut, return_type = 'all')$IS,
                                   IS_Ratio = (returnIS(hnRNPC_Mut, return_type = 'all')$IS+0.001)/(returnIS(hnRNPC_WT, return_type = 'all')$IS+0.001))

RBM25_Mut2WT_ISRatio = data.frame(MOTIF = RBM25_WT$MOTIF,
                                  WT_IS = returnIS(RBM25_WT, return_type = 'all')$IS,
                                  Mut_IS = returnIS(RBM25_Mut, return_type = 'all')$IS,
                                  IS_Ratio = (returnIS(RBM25_Mut, return_type = 'all')$IS+0.001)/(returnIS(RBM25_WT, return_type = 'all')$IS+0.001))

numMotif = 20

## For hnRNPC, top 20 SI changes
data = data.frame(hnRNPC_Mut2WT_ISRatio[, c('MOTIF', 'IS_Ratio')])
data = data[order(-data[, 'IS_Ratio']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## For hnRNPC, bottom 20 SI changes
data = data[order(data[, 'IS_Ratio']), ]

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## For RBM25, top 20 SI changes
data = data.frame(RBM25_Mut2WT_ISRatio[, c('MOTIF', 'IS_Ratio')])
data = data[order(-data[, 'IS_Ratio']), ]
rownames(data) = NULL

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

## For RBM25, bottom 20 SI changes
data = data[order(data[, 'IS_Ratio']), ]

motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')
seqLogo(PWM, ic.scale = F)

# Generate all 5-mers
K = 5
NTs = c("A", "C", "G", "T")
KMer_polyT = expand.grid(replicate(K, NTs, simplify = FALSE))
KMer_polyT = apply(KMer_polyT, 1, paste, collapse = "")
KMer_polyT = KMer_polyT[grepl("TTT", KMer_polyT)]

RBM25_polyT = RBM25_Mut2WT_ISRatio %>% filter(MOTIF %in% KMer_polyT) %>% filter(MOTIF != "TTTTT")

################################################################################


## PCA
################################################################################
All_peaks = fread('~/Desktop/Genomics/Specificity/Sequencing/ALL/All_Combined_peakCoverage.txt')
count_matrix = as.matrix(All_peaks[, 7:ncol(All_peaks)])
count_matrix = round(count_matrix)

# Create sample metadata (colData) for grouping and coloring
sample_names = colnames(count_matrix)
sample_info = data.frame(
  row.names = sample_names,
  condition = gsub("(\\d+)$", "", sample_names)
)

# Create a DESeqDataSet object
dds = DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ 1) # Use design = ~ 1 for QC

# Apply the variance stabilizing transformation (VST)
vsd = vst(dds, blind = TRUE)

pca_data = plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar = round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples") +
  theme_bw(base_size = 14) +
  # Use ggrepel to add sample labels without them overlapping
  # geom_text_repel(box.padding = 0.5, max.overlaps = Inf) +
  # coord_fixed()  +
  xlim(-50, 50) + ylim(-50, 50) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))



################################################################################


