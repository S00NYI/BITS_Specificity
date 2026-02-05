################################################################################
## Competition CLIP Analysis using RBPSpecificity Package:
## Written by Soon Yi
## Created: 2025-08-04
## Last Edited: 2025-09-03
## Figure 7
################################################################################

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(data.table)
library(rtracklayer)
library(BSgenome)
library(Biostrings)

library(seqLogo)
library(RBPSpecificity)
library(ggsignif)
library(biomaRt)
library(ggrepel)

## Set up basic parameters:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Sequencing/COLLAPSED/Final'
outputDir = '~/Desktop/Genomics/Specificity/AnalysisOutput/COLLAPSED_PEAKS/'
extensions = c(50)

K = 5
################################################################################

## Read files for filtered peaks:
################################################################################
chr_list = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
             'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM')

## eCLIP
Peak_eCLIP = fread(paste0(baseDir, '/eCLIP_peaks_filtered.txt'))
Peak_eCLIP = Peak_eCLIP %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_eCLIP = Peak_eCLIP %>% filter(chr %in% chr_list)

## Hela
Peak_HeLa = fread(paste0(baseDir, '/HeLa_peaks_filtered.txt'))
Peak_HeLa = Peak_HeLa %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_HeLa = Peak_HeLa %>% filter(chr %in% chr_list)

## Our data
hnRNPC_BC = 3

Peak_hnRNPC_inRBM25WT = fread(paste0(baseDir, '/hnRNPC_WT_inRBM25_WT_Combined_peakCoverage_filtered_annotated.txt'))
colnames(Peak_hnRNPC_inRBM25WT)[colnames(Peak_hnRNPC_inRBM25WT) == "chrom"] = "chr"
Peak_hnRNPC_inRBM25WT = Peak_hnRNPC_inRBM25WT %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_hnRNPC_inRBM25WT = Peak_hnRNPC_inRBM25WT %>% filter(chr %in% chr_list)
Peak_hnRNPC_inRBM25WT = Peak_hnRNPC_inRBM25WT %>% filter(avg_RPM >= median(Peak_hnRNPC_inRBM25WT$avg_RPM))
Peak_hnRNPC_inRBM25WT = Peak_hnRNPC_inRBM25WT %>% filter(BC >= hnRNPC_BC)

Peak_hnRNPC_inRBM25Mut = fread(paste0(baseDir, '/hnRNPC_WT_inRBM25_Mut_Combined_peakCoverage_filtered_annotated.txt'))
colnames(Peak_hnRNPC_inRBM25Mut)[colnames(Peak_hnRNPC_inRBM25Mut) == "chrom"] = "chr"
Peak_hnRNPC_inRBM25Mut = Peak_hnRNPC_inRBM25Mut %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_hnRNPC_inRBM25Mut = Peak_hnRNPC_inRBM25Mut %>% filter(chr %in% chr_list)
Peak_hnRNPC_inRBM25Mut = Peak_hnRNPC_inRBM25Mut %>% filter(avg_RPM >= median(Peak_hnRNPC_inRBM25Mut$avg_RPM))
Peak_hnRNPC_inRBM25Mut = Peak_hnRNPC_inRBM25Mut %>% filter(BC >= hnRNPC_BC)

Peak_hnRNPC_WT = fread(paste0(baseDir, '/hnRNPC_WT_Combined_peakCoverage_filtered_annotated.txt'))
colnames(Peak_hnRNPC_WT)[colnames(Peak_hnRNPC_WT) == "chrom"] = "chr"
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% filter(chr %in% chr_list)
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% filter(avg_RPM >= median(Peak_hnRNPC_WT$avg_RPM))
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% filter(BC >= hnRNPC_BC)

Peak_hnRNPC_Mut = fread(paste0(baseDir, '/hnRNPC_Mut_Combined_peakCoverage_filtered_annotated.txt'))
colnames(Peak_hnRNPC_Mut)[colnames(Peak_hnRNPC_Mut) == "chrom"] = "chr"
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% filter(chr %in% chr_list)
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% filter(avg_RPM >= median(Peak_hnRNPC_Mut$avg_RPM))
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% filter(BC >= hnRNPC_BC)

RBM25_BC = 6

Peak_RBM25_WT = fread(paste0(baseDir, '/RBM25_WT_Combined_peakCoverage_filtered_annotated.txt'))
colnames(Peak_RBM25_WT)[colnames(Peak_RBM25_WT) == "chrom"] = "chr"
Peak_RBM25_WT = Peak_RBM25_WT %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_RBM25_WT = Peak_RBM25_WT %>% filter(chr %in% chr_list)
Peak_RBM25_WT = Peak_RBM25_WT %>% filter(avg_RPM >= median(Peak_RBM25_WT$avg_RPM))
Peak_RBM25_WT = Peak_RBM25_WT %>% filter(BC >= RBM25_BC)

Peak_RBM25_Mut = fread(paste0(baseDir, '/RBM25_Mut_Combined_peakCoverage_filtered_annotated.txt'))
colnames(Peak_RBM25_Mut)[colnames(Peak_RBM25_Mut) == "chrom"] = "chr"
Peak_RBM25_Mut = Peak_RBM25_Mut %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_RBM25_Mut = Peak_RBM25_Mut %>% filter(chr %in% chr_list)
Peak_RBM25_Mut = Peak_RBM25_Mut %>% filter(avg_RPM >= median(Peak_RBM25_Mut$avg_RPM))
Peak_RBM25_Mut = Peak_RBM25_Mut %>% filter(BC >= RBM25_BC)
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

CLIP_List = c('eCLIP', 'HeLa', 'hnRNPC_WT', 'hnRNPC_Mut', 'RBM25_WT', 'RBM25_Mut', 'hnRNPC_in_RBM25WT', 'hnRNPC_in_RBM25Mut')
All_Annotation_List = c("5'UTR", "CDS", "3'UTR", "intron", 'rRNA', 'tRNA', 'ncRNA', 'TE', "Other")

Peak_eCLIP = Peak_eCLIP %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_eCLIP = Peak_eCLIP %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_HeLa = Peak_HeLa %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_HeLa = Peak_HeLa %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_hnRNPC_inRBM25WT = Peak_hnRNPC_inRBM25WT %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_hnRNPC_inRBM25WT = Peak_hnRNPC_inRBM25WT %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_hnRNPC_inRBM25Mut = Peak_hnRNPC_inRBM25Mut %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_hnRNPC_inRBM25Mut = Peak_hnRNPC_inRBM25Mut %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_hnRNPC_WT = Peak_hnRNPC_WT %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_hnRNPC_Mut = Peak_hnRNPC_Mut %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_RBM25_WT = Peak_RBM25_WT %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_RBM25_WT = Peak_RBM25_WT %>% filter(finalized_annotation %in% All_Annotation_List)

Peak_RBM25_Mut = Peak_RBM25_Mut %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'promoter', 'Other', finalized_annotation))
Peak_RBM25_Mut = Peak_RBM25_Mut %>% filter(finalized_annotation %in% All_Annotation_List)

Counts_Peak_eCLIP = countAnnotation(Peak_eCLIP, 'finalized_annotation', 'eCLIP', fraction = TRUE)
Counts_Peak_eCLIP['tRNA', 'eCLIP'] = 0
Counts_Peak_HeLa = countAnnotation(Peak_HeLa, 'finalized_annotation', 'HeLa', fraction = TRUE)
Counts_Peak_Peak_hnRNPC_inRBM25WT = countAnnotation(Peak_hnRNPC_inRBM25WT, 'finalized_annotation', 'hnRNPC_in_RBM25WT', fraction = TRUE)
Counts_Peak_Peak_hnRNPC_inRBM25WT = fillAnnotation(Counts_Peak_Peak_hnRNPC_inRBM25WT, All_Annotation_List)
Counts_Peak_Peak_hnRNPC_inRBM25Mut = countAnnotation(Peak_hnRNPC_inRBM25Mut, 'finalized_annotation', 'hnRNPC_in_RBM25Mut', fraction = TRUE)
Counts_Peak_Peak_hnRNPC_inRBM25Mut = fillAnnotation(Counts_Peak_Peak_hnRNPC_inRBM25Mut, All_Annotation_List)
Counts_Peak_hnRNPC_WT = countAnnotation(Peak_hnRNPC_WT, 'finalized_annotation', 'hnRNPC_WT', fraction = TRUE)
Counts_Peak_hnRNPC_WT = fillAnnotation(Counts_Peak_hnRNPC_WT, All_Annotation_List)
Counts_Peak_hnRNPC_Mut = countAnnotation(Peak_hnRNPC_Mut, 'finalized_annotation', 'hnRNPC_Mut', fraction = TRUE)
Counts_Peak_hnRNPC_Mut = fillAnnotation(Counts_Peak_hnRNPC_Mut, All_Annotation_List)
Counts_Peak_RBM25_WT = countAnnotation(Peak_RBM25_WT, 'finalized_annotation', 'RBM25_WT', fraction = TRUE)
Counts_Peak_RBM25_WT = fillAnnotation(Counts_Peak_RBM25_WT, All_Annotation_List)
Counts_Peak_RBM25_Mut = countAnnotation(Peak_RBM25_Mut, 'finalized_annotation', 'RBM25_Mut', fraction = TRUE)
Counts_Peak_RBM25_Mut = fillAnnotation(Counts_Peak_RBM25_Mut, All_Annotation_List)

PeakDistribution_combined = cbind(Counts_Peak_eCLIP, Counts_Peak_HeLa,
                                  Counts_Peak_hnRNPC_WT, Counts_Peak_hnRNPC_Mut,
                                  Counts_Peak_RBM25_WT, Counts_Peak_RBM25_Mut,
                                  Counts_Peak_Peak_hnRNPC_inRBM25WT,
                                  Counts_Peak_Peak_hnRNPC_inRBM25Mut)

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

# samplesList = c('eCLIP_hnRNPC_WT', 'HeLa_hnRNPC_WT',
#                 'HEK293_hnRNPC_inRBM25WT',
#                 'HEK293_hnRNPC_inRBM25Mut')

samplesList = c('HEK293_RBM25_WT', 'HEK293_RBM25_Mut')

# samplesList = c('HEK293_hnRNPC_inRBM25WT',
#                 'HEK293_hnRNPC_inRBM25Mut')

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
    } else if (sampleName == 'HeLa_hnRNPC_WT') {
      Peak = Peak_HeLa
    } else if (sampleName == 'HEK293_hnRNPC_inRBM25WT') {
      Peak = Peak_hnRNPC_inRBM25WT
    } else if (sampleName == 'HEK293_hnRNPC_inRBM25Mut') {
      Peak = Peak_hnRNPC_inRBM25Mut
    } else if (sampleName == 'HEK293_RBM25_WT') {
      Peak = Peak_RBM25_WT
    } else if (sampleName == 'HEK293_RBM25_Mut') {
      Peak = Peak_RBM25_Mut
    }

    CLIP_NAME = c(CLIP_NAME, sampleName)
    NT_EXT = c(NT_EXT, extension)

    enrichment = motifEnrichment(peak_data = Peak, 'hg38', K, extension = extension,
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
            paste0(outputDir, DATE, '/', K, 'mer_hnRNPC_Competition_', extension, 'ntExt.csv'), quote = F)

}



################################################################################

## PWM
################################################################################
# DATE = '20250824'
# extension = 25

# eCLIP = read_csv(paste0(outputDir, DATE, '/', 'eCLIP_hnRNPC_WT_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)
# HeLa = read_csv(paste0(outputDir, DATE, '/', 'HeLa_hnRNPC_WT_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)
# hnRNPC_inRBM25WT = read_csv(paste0(outputDir, DATE, '/', 'HEK293_hnRNPC_inRBM25WT_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)
# hnRNPC_inRBM25Mut = read_csv(paste0(outputDir, DATE, '/', 'HEK293_hnRNPC_inRBM25Mut_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)

RBM25_WT = read_csv(paste0(outputDir, DATE, '/', 'HEK293_RBM25_WT_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)
RBM25_Mut = read_csv(paste0(outputDir, DATE, '/', 'HEK293_RBM25_Mut_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)

numMotif = 10

# ## For eCLIP
# data = eCLIP
# data = data.frame(data[, c('MOTIF', 'Score')])
# data = data[order(-data[, 'Score']), ]
# rownames(data) = NULL
#
# motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
# pwm = PWM(consensusMatrix(motif))
# pwm = pwm - min(pwm)
#
# pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
# pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
# pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
# pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
# pwm[, 5] = pwm[, 5]/sum(pwm[, 5])
#
# PWM = makePWM(pwm, alphabet = 'RNA')
# seqLogo(PWM, ic.scale = F)
#
# ## For HeLa
# data = HeLa
# data = data.frame(data[, c('MOTIF', 'Score')])
# data = data[order(-data[, 'Score']), ]
# rownames(data) = NULL
#
# motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
# pwm = PWM(consensusMatrix(motif))
# pwm = pwm - min(pwm)
#
# pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
# pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
# pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
# pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
# pwm[, 5] = pwm[, 5]/sum(pwm[, 5])
#
# PWM = makePWM(pwm, alphabet = 'RNA')
# seqLogo(PWM, ic.scale = F)
#
# ## For hnRNPC in RBM25 WT
# data = hnRNPC_inRBM25WT
# data = data.frame(data[, c('MOTIF', 'Score')])
# data = data[order(-data[, 'Score']), ]
# rownames(data) = NULL
#
# motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
# pwm = PWM(consensusMatrix(motif))
# pwm = pwm - min(pwm)
#
# pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
# pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
# pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
# pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
# pwm[, 5] = pwm[, 5]/sum(pwm[, 5])
#
# PWM = makePWM(pwm, alphabet = 'RNA')
# seqLogo(PWM, ic.scale = F)
#
# ## For hnRNPC in RBM25 Mut
# data = hnRNPC_inRBM25Mut
# data = data.frame(data[, c('MOTIF', 'Score')])
# data = data[order(-data[, 'Score']), ]
# rownames(data) = NULL
#
# motif = DNAStringSet(gsub("U", "T", data$MOTIF[1:numMotif]))
# pwm = PWM(consensusMatrix(motif))
# pwm = pwm - min(pwm)
#
# pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
# pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
# pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
# pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
# pwm[, 5] = pwm[, 5]/sum(pwm[, 5])
#
# PWM = makePWM(pwm, alphabet = 'RNA')
# seqLogo(PWM, ic.scale = F)

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
# ggplot() +
#   geom_histogram(data = eCLIP, aes(x = Score), binwidth = 0.01, fill = "darkgoldenrod1", alpha = 0.5) +
#   geom_histogram(data = HeLa, aes(x = Score), binwidth = 0.01, fill = "cornflowerblue", alpha = 0.5) +
#   labs(title = "hnRNPC", x = "Normalized Counts", y = "Frequency") +
#   # scale_y_continuous(limits = c(0, 400), breaks = seq(0, 600, by = 50)) +
#   theme_minimal() +
#   theme_bw() +
#   theme(axis.text = element_text(size=14),
#         axis.title = element_text(size=14, face = 'bold'),
#         legend.text = element_text(size=14))
#
# ggplot() +
#   # geom_histogram(data = eCLIP, aes(x = Score), binwidth = 0.01, fill = "black", alpha = 0.5) +
#   # geom_histogram(data = HeLa, aes(x = Score), binwidth = 0.01, fill = "salmon", alpha = 0.5) +
#   geom_histogram(data = hnRNPC_inRBM25WT, aes(x = Score), binwidth = 0.01, fill = "darkgoldenrod1", alpha = 0.5) +
#   geom_histogram(data = hnRNPC_inRBM25Mut, aes(x = Score), binwidth = 0.01, fill = "cornflowerblue", alpha = 0.5) +
#   labs(title = "RBM25", x = "Normalized Counts", y = "Frequency") +
#   # scale_y_continuous(limits = c(0, 400), breaks = seq(0, 500, by = 50)) +
#   theme_minimal() +
#   theme_bw() +
#   theme(axis.text = element_text(size=14),
#         axis.title = element_text(size=14, face = 'bold'),
#         legend.text = element_text(size=14))

ggplot() +
  # geom_histogram(data = eCLIP, aes(x = Score), binwidth = 0.01, fill = "black", alpha = 0.5) +
  # geom_histogram(data = HeLa, aes(x = Score), binwidth = 0.01, fill = "salmon", alpha = 0.5) +
  geom_histogram(data = RBM25_WT, aes(x = Score), binwidth = 0.01, fill = "darkgoldenrod1", alpha = 0.5) +
  geom_histogram(data = RBM25_Mut, aes(x = Score), binwidth = 0.01, fill = "cornflowerblue", alpha = 0.5) +
  labs(title = "RBM25", x = "Normalized Counts", y = "Frequency") +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 500, by = 50)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))


################################################################################

## MS Analysis
################################################################################
plotMS(eCLIP, output_type = 'matrix')
plotMS(HeLa, output_type = 'matrix')
plotMS(hnRNPC_inRBM25WT, output_type = 'matrix')
plotMS(hnRNPC_inRBM25Mut, output_type = 'matrix')

returnMS(eCLIP, output_type = 'number')                 # 0.89
returnMS(HeLa, output_type = 'number')                  # 0.74
returnMS(hnRNPC_inRBM25WT, output_type = 'number')      # 0.77
returnMS(hnRNPC_inRBM25Mut, output_type = 'number')     # 0.77

returnIS(eCLIP)                                         # 41.2
returnIS(HeLa)                                          # 24.4
returnIS(hnRNPC_inRBM25WT)                              # 42.7
returnIS(hnRNPC_inRBM25Mut)                             # 44.2

################################################################################

## Gene Level Analysis hnRNPC WT vs hnRNPC WT in RBM25 WT OE
################################################################################
hnRNPC_targetGenes = unique(Peak_hnRNPC_WT$gene)
hnRNPC_targetGenes_inRBM25WT = unique(Peak_hnRNPC_inRBM25WT$gene)
hnRNPC_targetGenes_overlap = intersect(hnRNPC_targetGenes, hnRNPC_targetGenes_inRBM25WT)

targetGenes_RPM_inEndo = Peak_hnRNPC_WT %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inEndo = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM_inRBM25WT = Peak_hnRNPC_inRBM25WT %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inRBM25WT = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM = inner_join(targetGenes_RPM_inEndo, targetGenes_RPM_inRBM25WT, by = "gene")
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

gene_names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = targetGenes_RPM$gene,
                   mart = ensembl)
targetGenes_RPM = left_join(targetGenes_RPM, gene_names, by = c("gene" = "ensembl_gene_id"))

targetGenes_RPM$inEndo_RPKM = targetGenes_RPM$inEndo/targetGenes_RPM$gene_length * 1000
targetGenes_RPM$inRBM25WT_RPKM = targetGenes_RPM$inRBM25WT/targetGenes_RPM$gene_length * 1000

targetGenes_RPM$L2FC_RPKM = log2(targetGenes_RPM$inRBM25WT_RPKM/targetGenes_RPM$inEndo_RPKM)
targetGenes_RPM = targetGenes_RPM[!is.na(targetGenes_RPM$L2FC_RPKM), ]
targetGenes_RPM = targetGenes_RPM[!(targetGenes_RPM$external_gene_name == ''), ]

p_val = wilcox.test(targetGenes_RPM$inEndo_RPKM, targetGenes_RPM$inRBM25WT_RPKM)
p_val = scales::pvalue(p_val$p.value)

targetGenes_long = pivot_longer(targetGenes_RPM,
                                 cols = c("inEndo_RPKM", "inRBM25WT_RPKM"),
                                 names_to = "condition",
                                 values_to = "RPKM")

ggplot(targetGenes_long, aes(x = condition, y = RPKM, fill = condition)) +
  geom_boxplot(outliers = TRUE, notch = TRUE) +
  scale_y_log10() +
  geom_signif(comparisons = list(c('inEndo_RPKM', 'inRBM25WT_RPKM')),
              annotations = p_val,
              y_position = log10(10000),
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

genes_to_highlight = c("ENSG00000169760", "ENSG00000101911", "ENSG00000172318", "ENSG00000101596", "ENSG00000184903", "ENSG00000184047")
highlight_df = targetGenes_RPM %>% filter(gene %in% genes_to_highlight)

ggplot(targetGenes_RPM, aes(x = inEndo_RPKM, y = inRBM25WT_RPKM)) +
  geom_hex(bins = 25, color = 'black') + scale_fill_viridis_c() +
  # geom_point(alpha = 0.1) +
  # geom_point(data = highlight_df, color = "red", size = 4, shape = 18) +
  # geom_text_repel(data = highlight_df, aes(label = external_gene_name),
  #                 color = "black",
  #                 box.padding = 0.5,
  #                 max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, color = "salmon", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "RPKM per hnRNPC target genes",
       x = "Log(RPKM) in hnRNPC WT CLIP",
       y = "Log(RPKM) in RBM25 WT") +
  coord_fixed() + # Ensure a 1:1 aspect ratio for the axes
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))
################################################################################

## Gene Level Analysis hnRNPC WT vs hnRNPC WT in RBM25 Mut OE
################################################################################
hnRNPC_targetGenes = unique(Peak_hnRNPC_WT$gene)
hnRNPC_targetGenes_inRBM25Mut = unique(Peak_hnRNPC_inRBM25Mut$gene)
hnRNPC_targetGenes_overlap = intersect(hnRNPC_targetGenes, hnRNPC_targetGenes_inRBM25Mut)

targetGenes_RPM_inEndo = Peak_hnRNPC_WT %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inEndo = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM_inRBM25Mut = Peak_hnRNPC_inRBM25Mut %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inRBM25Mut = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM = inner_join(targetGenes_RPM_inEndo, targetGenes_RPM_inRBM25Mut, by = "gene")
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

gene_names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = targetGenes_RPM$gene,
                   mart = ensembl)
targetGenes_RPM = left_join(targetGenes_RPM, gene_names, by = c("gene" = "ensembl_gene_id"))

targetGenes_RPM$inEndo_RPKM = targetGenes_RPM$inEndo/targetGenes_RPM$gene_length * 1000
targetGenes_RPM$inRBM25Mut_RPKM = targetGenes_RPM$inRBM25Mut/targetGenes_RPM$gene_length * 1000

targetGenes_RPM$L2FC_RPKM = log2(targetGenes_RPM$inRBM25Mut_RPKM/targetGenes_RPM$inEndo_RPKM)
targetGenes_RPM = targetGenes_RPM[!is.na(targetGenes_RPM$L2FC_RPKM), ]
targetGenes_RPM = targetGenes_RPM[!(targetGenes_RPM$external_gene_name == ''), ]

p_val = wilcox.test(targetGenes_RPM$inEndo_RPKM, targetGenes_RPM$inRBM25Mut_RPKM)
p_val = scales::pvalue(p_val$p.value)

targetGenes_long = pivot_longer(targetGenes_RPM,
                                cols = c("inEndo_RPKM", "inRBM25Mut_RPKM"),
                                names_to = "condition",
                                values_to = "RPKM")

ggplot(targetGenes_long, aes(x = condition, y = RPKM, fill = condition)) +
  geom_boxplot(outliers = TRUE, notch = TRUE) +
  scale_y_log10() +
  geom_signif(comparisons = list(c('inEndo_RPKM', 'inRBM25Mut_RPKM')),
              annotations = p_val,
              y_position = log10(10000),
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

genes_to_highlight = c("ENSG00000169760", "ENSG00000101911", "ENSG00000172318", "ENSG00000101596", "ENSG00000184903", "ENSG00000184047")
highlight_df = targetGenes_RPM %>% filter(gene %in% genes_to_highlight)

ggplot(targetGenes_RPM, aes(x = inEndo_RPKM, y = inRBM25Mut_RPKM)) +
  geom_hex(bins = 25, color = 'black') + scale_fill_viridis_c() +
  # geom_point(alpha = 0.1) +
  # geom_point(data = highlight_df, color = "red", size = 4, shape = 18) +
  # geom_text_repel(data = highlight_df, aes(label = external_gene_name),
  #                 color = "black",
  #                 box.padding = 0.5,
  #                 max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, color = "salmon", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "RPKM per hnRNPC target genes",
       x = "Log(RPKM) in hnRNPC WT CLIP",
       y = "Log(RPKM) in RBM25 Mut") +
  coord_fixed() + # Ensure a 1:1 aspect ratio for the axes
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))
################################################################################

## Gene Level Analysis hnRNPC WT in RBM25 WT OE vs in RBM25 Mut OE
################################################################################
hnRNPC_targetGenes_inRBM25WT = unique(Peak_hnRNPC_inRBM25WT$gene)
hnRNPC_targetGenes_inRBM25Mut = unique(Peak_hnRNPC_inRBM25Mut$gene)
hnRNPC_targetGenes_overlap = intersect(hnRNPC_targetGenes_inRBM25WT, hnRNPC_targetGenes_inRBM25Mut)

targetGenes_RPM_inRBM25WT = Peak_hnRNPC_inRBM25WT %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inRBM25WT = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM_inRBM25Mut = Peak_hnRNPC_inRBM25Mut %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inRBM25Mut = sum(avg_RPM, na.rm = TRUE))

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

gene_names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = targetGenes_RPM$gene,
                   mart = ensembl)
targetGenes_RPM = left_join(targetGenes_RPM, gene_names, by = c("gene" = "ensembl_gene_id"))

targetGenes_RPM$inRBM25WT_RPKM = targetGenes_RPM$inRBM25WT/targetGenes_RPM$gene_length * 1000
targetGenes_RPM$inRBM25Mut_RPKM = targetGenes_RPM$inRBM25Mut/targetGenes_RPM$gene_length * 1000

targetGenes_RPM$L2FC_RPKM = log2(targetGenes_RPM$inRBM25Mut_RPKM/targetGenes_RPM$inRBM25WT_RPKM)
targetGenes_RPM = targetGenes_RPM[!is.na(targetGenes_RPM$L2FC_RPKM), ]
targetGenes_RPM = targetGenes_RPM[!(targetGenes_RPM$external_gene_name == ''), ]

p_val = wilcox.test(targetGenes_RPM$inRBM25WT_RPKM, targetGenes_RPM$inRBM25Mut_RPKM)
p_val = scales::pvalue(p_val$p.value)

targetGenes_long = pivot_longer(targetGenes_RPM,
                                cols = c("inRBM25WT_RPKM", "inRBM25Mut_RPKM"),
                                names_to = "condition",
                                values_to = "RPKM")

ggplot(targetGenes_long, aes(x = condition, y = RPKM, fill = condition)) +
  geom_boxplot(outliers = TRUE, notch = TRUE) +
  scale_y_log10() +
  geom_signif(comparisons = list(c('inRBM25WT_RPKM', 'inRBM25Mut_RPKM')),
              annotations = p_val,
              y_position = log10(1000),
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

genes_to_highlight = c("ENSG00000169760", "ENSG00000101911", "ENSG00000172318", "ENSG00000101596", "ENSG00000184903", "ENSG00000184047")
highlight_df = targetGenes_RPM %>% filter(gene %in% genes_to_highlight)

ggplot(targetGenes_RPM, aes(x = inRBM25WT_RPKM, y = inRBM25Mut_RPKM)) +
  geom_hex(bins = 25, color = 'black') + scale_fill_viridis_c() +
  # geom_point(alpha = 0.1) +
  geom_point(data = highlight_df, color = "red", size = 4, shape = 18) +
  geom_text_repel(data = highlight_df, aes(label = external_gene_name),
                  color = "black",
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, color = "salmon", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "RPKM per hnRNPC target genes in RBM25 WT vs. in RBM25 Mut)",
       x = "Log(RPKM) in RBM25 WT",
       y = "Log(RPKM) in RBM25 Mut") +
  coord_fixed() + # Ensure a 1:1 aspect ratio for the axes
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))
################################################################################

## Gene Level Analysis #3 Endo vs RBM25 WT OE vs RBM25 Mut OE
################################################################################
hnRNPC_targetGenes_inWT = unique(Peak_hnRNPC_WT$gene)
hnRNPC_targetGenes_inRBM25WT = unique(Peak_hnRNPC_inRBM25WT$gene)
hnRNPC_targetGenes_inRBM25Mut = unique(Peak_hnRNPC_inRBM25Mut$gene)

hnRNPC_targetGenes_overlap = intersect(hnRNPC_targetGenes_inRBM25WT, hnRNPC_targetGenes_inRBM25Mut)
hnRNPC_targetGenes_overlap = intersect(hnRNPC_targetGenes_inWT, hnRNPC_targetGenes_overlap)

targetGenes_RPM_inWT= Peak_hnRNPC_WT %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inWT = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM_inRBM25WT = Peak_hnRNPC_inRBM25WT %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inRBM25WT = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM_inRBM25Mut = Peak_hnRNPC_inRBM25Mut %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inRBM25Mut = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM = inner_join(targetGenes_RPM_inRBM25WT, targetGenes_RPM_inRBM25Mut, by = "gene")
targetGenes_RPM = inner_join(targetGenes_RPM, targetGenes_RPM_inWT, by = "gene")

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

gene_names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = targetGenes_RPM$gene,
                   mart = ensembl)
targetGenes_RPM = left_join(targetGenes_RPM, gene_names, by = c("gene" = "ensembl_gene_id"))

targetGenes_RPM$inWT_RPKM = targetGenes_RPM$inWT/targetGenes_RPM$gene_length * 1000
targetGenes_RPM$inRBM25WT_RPKM = targetGenes_RPM$inRBM25WT/targetGenes_RPM$gene_length * 1000
targetGenes_RPM$inRBM25Mut_RPKM = targetGenes_RPM$inRBM25Mut/targetGenes_RPM$gene_length * 1000

# targetGenes_RPM$L2FC_RPKM = log2(targetGenes_RPM$inRBM25Mut_RPKM/targetGenes_RPM$inRBM25WT_RPKM)
# targetGenes_RPM = targetGenes_RPM[!is.na(targetGenes_RPM$L2FC_RPKM), ]
targetGenes_RPM = targetGenes_RPM[!(targetGenes_RPM$external_gene_name == ''), ]

p_val_WT_RWT = wilcox.test(targetGenes_RPM$inWT_RPKM, targetGenes_RPM$inRBM25WT_RPKM)
p_val_WT_RWT = scales::pvalue(p_val_WT_RWT$p.value)

p_val_WT_RMut = wilcox.test(targetGenes_RPM$inWT_RPKM, targetGenes_RPM$inRBM25Mut_RPKM)
p_val_WT_RMut = scales::pvalue(p_val_WT_RMut$p.value)

p_val_RWT_RMut = wilcox.test(targetGenes_RPM$inRBM25WT_RPKM, targetGenes_RPM$inRBM25Mut_RPKM)
p_val_RWT_RMut = scales::pvalue(p_val_RWT_RMut$p.value)

targetGenes_long = pivot_longer(targetGenes_RPM,
                                cols = c("inWT_RPKM", "inRBM25WT_RPKM", "inRBM25Mut_RPKM"),
                                names_to = "condition",
                                values_to = "RPKM")

ggplot(targetGenes_long, aes(x = condition, y = RPKM, fill = condition)) +
  geom_boxplot(outliers = TRUE, notch = TRUE) +
  scale_y_log10() +
  # geom_signif(comparisons = list(c('inWT_RPKM', 'inRBM25WT_RPKM'), c('inWT_RPKM', 'inRBM25Mut_RPKM'), c('inRBM25WT_RPKM', 'inRBM25Mut_RPKM')),
  #             annotations = c(p_val_WT_RWT, p_val_WT_RMut, p_val_RWT_RMut),
  #             y_position = c(log10(10000), log10(50000), log10(10000)),
  #             tip_length = 0.01,
  #             map_signif_level = TRUE) +
  labs(title = "RPKM for hnRNPC target genes",
       x = "Cell Lines",
       y = "RPKM per Gene") +
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))

genes_to_highlight = c("ENSG00000169760", "ENSG00000101911", "ENSG00000172318", "ENSG00000101596", "ENSG00000184903", "ENSG00000184047")
highlight_df = targetGenes_RPM %>% filter(gene %in% genes_to_highlight)

ggplot(targetGenes_RPM, aes(x = inRBM25WT_RPKM, y = inRBM25Mut_RPKM)) +
  geom_hex(bins = 25, color = 'black') + scale_fill_viridis_c() +
  # geom_point(alpha = 0.1) +
  geom_point(data = highlight_df, color = "red", size = 4, shape = 18) +
  geom_text_repel(data = highlight_df, aes(label = external_gene_name),
                  color = "black",
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, color = "salmon", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "RPKM per hnRNPC target genes in RBM25 WT vs. in RBM25 Mut)",
       x = "Log(RPKM) in RBM25 WT",
       y = "Log(RPKM) in RBM25 Mut") +
  coord_fixed() + # Ensure a 1:1 aspect ratio for the axes
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))
################################################################################

## Filter hnRNPC in RBM25 OE
################################################################################
GR_Peak_hnRNPC_WT = makeGRangesFromDataFrame(Peak_hnRNPC_WT, keep.extra.columns = TRUE)
GR_Peak_hnRNPC_WT_inRBM25WT = makeGRangesFromDataFrame(Peak_hnRNPC_inRBM25WT, keep.extra.columns = TRUE)
GR_Peak_hnRNPC_WT_inRBM25Mut = makeGRangesFromDataFrame(Peak_hnRNPC_inRBM25Mut, keep.extra.columns = TRUE)

GR_Peak_hnRNPC_WT_inRBM25WT = subsetByOverlaps(GR_Peak_hnRNPC_WT_inRBM25WT, GR_Peak_hnRNPC_WT)
GR_Peak_hnRNPC_WT_inRBM25Mut = subsetByOverlaps(GR_Peak_hnRNPC_WT_inRBM25Mut, GR_Peak_hnRNPC_WT)

Peak_hnRNPC_WT_inRBM25WT_filtered = as.data.frame(GR_Peak_hnRNPC_WT_inRBM25WT)
Peak_hnRNPC_WT_inRBM25WT_filtered = Peak_hnRNPC_WT_inRBM25WT_filtered[, c('seqnames', 'start', 'end', 'peak_names', 'score', 'strand')]

Peak_hnRNPC_WT_inRBM25Mut_filtered = as.data.frame(GR_Peak_hnRNPC_WT_inRBM25Mut)
Peak_hnRNPC_WT_inRBM25Mut_filtered = Peak_hnRNPC_WT_inRBM25Mut_filtered[, c('seqnames', 'start', 'end', 'peak_names', 'score', 'strand')]

DATE = format(Sys.Date(), "%Y%m%d")
if (!dir.exists(paste0(outputDir, DATE, '/'))) {
  dir.create(paste0(outputDir, DATE, '/'), recursive = TRUE, showWarnings = FALSE)
  print(paste("Created directory:", paste0(outputDir, DATE, '/')))
}

write.table(Peak_hnRNPC_WT_inRBM25WT_filtered,
            paste0(outputDir, DATE, '/Peak_hnRNPC_WT_inRBM25_WT.bed'), sep = '\t', col.names = FALSE, row.names = FALSE, quote = F)

write.table(Peak_hnRNPC_WT_inRBM25Mut_filtered,
            paste0(outputDir, DATE, '/Peak_hnRNPC_WT_inRBM25_Mut.bed'), sep = '\t', col.names = FALSE, row.names = FALSE, quote = F)

################################################################################


## Subgroup Analysis

## Peaks: + in Endo, + in WT OE, - in Mut OE
## Peaks: + in Endo, - in WT OE, + in Mut OE
## Peaks: + in Endo, + in WT OE, + in Mut OE
## Peaks: + in Endo, - in WT OE, - in Mut OE
## Peaks: - in Endo, + in WT OE, + in Mut OE
## Peaks: - in Endo, + in WT OE, - in Mut OE
################################################################################
GR_Peak_hnRNPC_WT = makeGRangesFromDataFrame(Peak_hnRNPC_WT, keep.extra.columns = TRUE)
GR_Peak_RBM25_WT = makeGRangesFromDataFrame(Peak_RBM25_WT, keep.extra.columns = TRUE)
GR_Peak_RBM25_Mut = makeGRangesFromDataFrame(Peak_RBM25_Mut, keep.extra.columns = TRUE)

Peak_hnRNPC_WT_Unique = subsetByOverlaps(GR_Peak_hnRNPC_WT, GR_Peak_RBM25_WT, invert = TRUE)
Peak_hnRNPC_WT_Unique = subsetByOverlaps(Peak_hnRNPC_WT_Unique, GR_Peak_RBM25_Mut, invert = TRUE)

Peak_RBM25_WT_Unique = subsetByOverlaps(GR_Peak_RBM25_WT, GR_Peak_hnRNPC_WT, invert = TRUE)
Peak_RBM25_WT_Unique = subsetByOverlaps(Peak_RBM25_WT_Unique, GR_Peak_RBM25_Mut, invert = TRUE)

Peak_RBM25_Mut_Unique = subsetByOverlaps(GR_Peak_RBM25_Mut, GR_Peak_hnRNPC_WT, invert = TRUE)
Peak_RBM25_Mut_Unique = subsetByOverlaps(Peak_RBM25_Mut_Unique, GR_Peak_RBM25_WT, invert = TRUE)




Gene_RBM25_WT_Unique = setdiff(Peak_RBM25_WT_Unique$gene, Peak_RBM25_Mut_Unique$gene)
Gene_RBM25_Mut_Unique = setdiff(Peak_RBM25_Mut_Unique$gene, Peak_RBM25_WT_Unique$gene)

DATE = format(Sys.Date(), "%Y%m%d")
if (!dir.exists(paste0(outputDir, DATE, '/'))) {
  dir.create(paste0(outputDir, DATE, '/'), recursive = TRUE, showWarnings = FALSE)
  print(paste("Created directory:", paste0(outputDir, DATE, '/')))
}

samplesList = c('Peak_hnRNPC_WT_Unique',
                'Peak_RBM25_WT_Unique',
                'Peak_RBM25_Mut_Unique')

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

    if (sampleName == 'Peak_RBM25_WT_Unique') {
      Peak = Peak_RBM25_WT_Unique
    } else if (sampleName == 'Peak_RBM25_Mut_Unique') {
      Peak = Peak_RBM25_Mut_Unique
    }

    CLIP_NAME = c(CLIP_NAME, sampleName)
    NT_EXT = c(NT_EXT, extension)

    enrichment = motifEnrichment(peak_data = Peak, 'hg38', K, extension = extension, Bkg_number = 100, max_shift_dist = 250, Bkg_dist = 50)
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
            paste0(outputDir, DATE, '/', K, 'mer_Unique_Peaks_', extension, 'ntExt.csv'), quote = F)

}

write.table(as.data.frame(Peak_hnRNPC_WT_Unique),
          paste0(outputDir, DATE, 'Peak_hnRNPC_WT_Unique.bed'), sep = '\t', col.names = FALSE, row.names = FALSE, quote = F)

write.table(as.data.frame(Peak_RBM25_WT_Unique),
            paste0(outputDir, DATE, 'Peak_hnRNPC_WT_Unique.bed'), sep = '\t', col.names = FALSE, row.names = FALSE, quote = F)

write.table(as.data.frame(Peak_RBM25_Mut_Unique),
            paste0(outputDir, DATE, 'Peak_hnRNPC_WT_Unique.bed'), sep = '\t', col.names = FALSE, row.names = FALSE, quote = F)


hnRNPC_WT_Unique = read_csv(paste0(outputDir, DATE, '/', 'Peak_hnRNPC_WT_Unique_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)
RBM25_WT_Unique = read_csv(paste0(outputDir, DATE, '/', 'Peak_RBM25_WT_Unique_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)
RBM25_Mut_Unique = read_csv(paste0(outputDir, DATE, '/', 'Peak_RBM25_Mut_Unique_NormalizedEnrichment_', extension, 'ntExt_', as.character(K), 'mer.csv'), col_names = T, show_col_types = F)

numMotif = 25

## For RBM25_WT_Unique
data = RBM25_WT_Unique
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

## For RBM25_Mut_Unique
data = RBM25_Mut_Unique
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

## Gene Level Analysis hnRNPC WT vs hnRNPC WT in RBM25 WT OE
################################################################################
hnRNPC_targetGenes = unique(Peak_hnRNPC_WT$gene)
hnRNPC_targetGenes_inRBM25WT = unique(Peak_hnRNPC_inRBM25WT$gene)
hnRNPC_targetGenes_overlap = intersect(hnRNPC_targetGenes, hnRNPC_targetGenes_inRBM25WT)
hnRNPC_targetGenes_overlap = intersect(hnRNPC_targetGenes_overlap, Gene_RBM25_WT_Unique)

targetGenes_RPM_inEndo = Peak_hnRNPC_WT %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inEndo = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM_inRBM25WT = Peak_hnRNPC_inRBM25WT %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inRBM25WT = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM = inner_join(targetGenes_RPM_inEndo, targetGenes_RPM_inRBM25WT, by = "gene")
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

gene_names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = targetGenes_RPM$gene,
                   mart = ensembl)
targetGenes_RPM = left_join(targetGenes_RPM, gene_names, by = c("gene" = "ensembl_gene_id"))

targetGenes_RPM$inEndo_RPKM = targetGenes_RPM$inEndo/targetGenes_RPM$gene_length * 1000
targetGenes_RPM$inRBM25WT_RPKM = targetGenes_RPM$inRBM25WT/targetGenes_RPM$gene_length * 1000

targetGenes_RPM$L2FC_RPKM = log2(targetGenes_RPM$inRBM25WT_RPKM/targetGenes_RPM$inEndo_RPKM)
targetGenes_RPM = targetGenes_RPM[!is.na(targetGenes_RPM$L2FC_RPKM), ]
targetGenes_RPM = targetGenes_RPM[!(targetGenes_RPM$external_gene_name == ''), ]

p_val = wilcox.test(targetGenes_RPM$inEndo_RPKM, targetGenes_RPM$inRBM25WT_RPKM)
p_val = scales::pvalue(p_val$p.value)

targetGenes_long = pivot_longer(targetGenes_RPM,
                                cols = c("inEndo_RPKM", "inRBM25WT_RPKM"),
                                names_to = "condition",
                                values_to = "RPKM")

ggplot(targetGenes_long, aes(x = condition, y = RPKM, fill = condition)) +
  geom_boxplot(outliers = TRUE, notch = TRUE) +
  scale_y_log10() +
  geom_signif(comparisons = list(c('inEndo_RPKM', 'inRBM25WT_RPKM')),
              annotations = p_val,
              y_position = log10(10000),
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

genes_to_highlight = c("ENSG00000169760", "ENSG00000101911", "ENSG00000172318", "ENSG00000101596", "ENSG00000184903", "ENSG00000184047")
highlight_df = targetGenes_RPM %>% filter(gene %in% genes_to_highlight)

ggplot(targetGenes_RPM, aes(x = inEndo_RPKM, y = inRBM25WT_RPKM)) +
  # geom_hex(bins = 25, color = 'black') + scale_fill_viridis_c() +
  geom_point(alpha = 0.1) +
  # geom_point(data = highlight_df, color = "red", size = 4, shape = 18) +
  # geom_text_repel(data = highlight_df, aes(label = external_gene_name),
  #                 color = "black",
  #                 box.padding = 0.5,
  #                 max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, color = "salmon", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "RPKM per hnRNPC target genes",
       x = "Log(RPKM) in hnRNPC WT CLIP",
       y = "Log(RPKM) in RBM25 WT") +
  coord_fixed() + # Ensure a 1:1 aspect ratio for the axes
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))
################################################################################

## Gene Level Analysis hnRNPC WT vs hnRNPC WT in RBM25 Mut OE
################################################################################
hnRNPC_targetGenes = unique(Peak_hnRNPC_WT$gene)
hnRNPC_targetGenes_inRBM25Mut = unique(Peak_hnRNPC_inRBM25Mut$gene)
hnRNPC_targetGenes_overlap = intersect(hnRNPC_targetGenes, hnRNPC_targetGenes_inRBM25Mut)
hnRNPC_targetGenes_overlap = intersect(hnRNPC_targetGenes_overlap, Gene_RBM25_Mut_Unique)

targetGenes_RPM_inEndo = Peak_hnRNPC_WT %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inEndo = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM_inRBM25Mut = Peak_hnRNPC_inRBM25Mut %>%
  filter(gene %in% hnRNPC_targetGenes_overlap) %>%
  group_by(gene) %>%
  summarise(inRBM25Mut = sum(avg_RPM, na.rm = TRUE))

targetGenes_RPM = inner_join(targetGenes_RPM_inEndo, targetGenes_RPM_inRBM25Mut, by = "gene")
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

gene_names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = targetGenes_RPM$gene,
                   mart = ensembl)
targetGenes_RPM = left_join(targetGenes_RPM, gene_names, by = c("gene" = "ensembl_gene_id"))

targetGenes_RPM$inEndo_RPKM = targetGenes_RPM$inEndo/targetGenes_RPM$gene_length * 1000
targetGenes_RPM$inRBM25Mut_RPKM = targetGenes_RPM$inRBM25Mut/targetGenes_RPM$gene_length * 1000

targetGenes_RPM$L2FC_RPKM = log2(targetGenes_RPM$inRBM25Mut_RPKM/targetGenes_RPM$inEndo_RPKM)
targetGenes_RPM = targetGenes_RPM[!is.na(targetGenes_RPM$L2FC_RPKM), ]
targetGenes_RPM = targetGenes_RPM[!(targetGenes_RPM$external_gene_name == ''), ]

p_val = wilcox.test(targetGenes_RPM$inEndo_RPKM, targetGenes_RPM$inRBM25Mut_RPKM)
p_val = scales::pvalue(p_val$p.value)

targetGenes_long = pivot_longer(targetGenes_RPM,
                                cols = c("inEndo_RPKM", "inRBM25Mut_RPKM"),
                                names_to = "condition",
                                values_to = "RPKM")

ggplot(targetGenes_long, aes(x = condition, y = RPKM, fill = condition)) +
  geom_boxplot(outliers = TRUE, notch = TRUE) +
  scale_y_log10() +
  geom_signif(comparisons = list(c('inEndo_RPKM', 'inRBM25Mut_RPKM')),
              annotations = p_val,
              y_position = log10(10000),
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

genes_to_highlight = c("ENSG00000169760", "ENSG00000101911", "ENSG00000172318", "ENSG00000101596", "ENSG00000184903", "ENSG00000184047")
highlight_df = targetGenes_RPM %>% filter(gene %in% genes_to_highlight)

ggplot(targetGenes_RPM, aes(x = inEndo_RPKM, y = inRBM25Mut_RPKM)) +
  # geom_hex(bins = 25, color = 'black') + scale_fill_viridis_c() +
  geom_point(alpha = 0.1) +
  # geom_point(data = highlight_df, color = "red", size = 4, shape = 18) +
  # geom_text_repel(data = highlight_df, aes(label = external_gene_name),
  #                 color = "black",
  #                 box.padding = 0.5,
  #                 max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, color = "salmon", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "RPKM per hnRNPC target genes",
       x = "Log(RPKM) in hnRNPC WT CLIP",
       y = "Log(RPKM) in RBM25 Mut") +
  coord_fixed() + # Ensure a 1:1 aspect ratio for the axes
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))
################################################################################









## PCA
################################################################################
library(DESeq2)
peakPath = '/Users/soonyi/Desktop/Genomics/Specificity/Sequencing/ALL/All_Combined_peakCoverage.txt'
All_peaks = fread(peakPath)
count_matrix = as.matrix(All_peaks[, 7:ncol(All_peaks)])
count_matrix = round(count_matrix)
colnames(count_matrix) = c('hnRNPC_inRBM25_MUT1', 'hnRNPC_inRBM25_MUT2', 'hnRNPC_inRBM25_MUT3',
                           'hnRNPC_inRBM25_WT1', 'hnRNPC_inRBM25_WT2', 'hnRNPC_inRBM25_WT3',
                           'RBM25_MUT1', 'RBM25_MUT2', 'RBM25_MUT3',
                           'RBM25_WT1', 'RBM25_WT2', 'RBM25_WT3',
                           'hnRNPC_MUT1', 'hnRNPC_MUT2', 'hnRNPC_MUT3',
                           'hnRNPC_MUT4', 'hnRNPC_MUT5', 'hnRNPC_WT1',
                           'hnRNPC_WT2', 'hnRNPC_WT3', 'hnRNPC_WT4',
                           'hnRNPC_WT5', 'RBM25_MUT4', 'RBM25_MUT5',
                           'RBM25_MUT6', 'RBM25_MUT7', 'RBM25_MUT8',
                           'RBM25_WT4', 'RBM25_WT5', 'RBM25_WT6',
                           'RBM25_WT7', 'RBM25_WT8')

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
vsd = varianceStabilizingTransformation(dds, blind = TRUE)

pca_data = plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar = round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples") +
  theme_bw(base_size = 14) +
  # geom_text_repel(box.padding = 0.5, max.overlaps = Inf) +
  # coord_fixed()  +
  xlim(-120, 120) + ylim(-120, 120) +
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
GR_Peak_hnRNPC_inRBM25WT = peaks2GR(Peak_hnRNPC_inRBM25WT)
GR_Peak_hnRNPC_inRBM25Mut = peaks2GR(Peak_hnRNPC_inRBM25Mut)

# Overlap between hnRNPC in RBM25WT and RBM25MUT
hnRNPC_Peaks_overlap = as.data.frame(findOverlaps(GR_Peak_hnRNPC_inRBM25WT,
                                                  GR_Peak_hnRNPC_inRBM25Mut, minoverlap = 15))

hnRNPC_overlapPeaks_inRBM25WT = Peak_hnRNPC_inRBM25WT[hnRNPC_Peaks_overlap$queryHits]
hnRNPC_overlapPeaks_inRBM25Mut = Peak_hnRNPC_inRBM25Mut[hnRNPC_Peaks_overlap$subjectHits]

hnRNPC_overlapPeaks_RPM = data.frame(inRBM25WT = hnRNPC_overlapPeaks_inRBM25WT$avg_RPM,
                                     inRBM25Mut = hnRNPC_overlapPeaks_inRBM25Mut$avg_RPM)

hnRNPC_overlapPeaks_RPM_long = pivot_longer(hnRNPC_overlapPeaks_RPM,
                                            cols = c("inRBM25WT", "inRBM25Mut"),
                                            names_to = "condition",
                                            values_to = "avg_RPM")

p_val = t.test(hnRNPC_overlapPeaks_RPM$inRBM25WT, hnRNPC_overlapPeaks_RPM$inRBM25Mut)
p_val = scales::pvalue(p_val$p.value)

ggplot(hnRNPC_overlapPeaks_RPM_long, aes(x = condition, y = avg_RPM, fill = condition)) +
  geom_boxplot(outliers = FALSE, notch = TRUE) +
  scale_y_log10() +
  geom_signif(comparisons = list(c('inRBM25WT', 'inRBM25Mut')),
              annotations = p_val,
              y_position = log10(2.5),
              tip_length = 0.01,
              map_signif_level = TRUE) +
  labs(title = "RPM per hnRNPC Peak",
       x = "Cell Lines",
       y = "RPM per Peak") +
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))

plot_data = hnRNPC_overlapPeaks_RPM
plot_data$inRBM25WT = plot_data$inRBM25WT + 0.001
plot_data$inRBM25Mut = plot_data$inRBM25Mut + 0.001

# Create the 2D binned heatmap
ggplot(plot_data, aes(x = inRBM25WT, y = inRBM25Mut)) +
  geom_bin2d(bins = 100) +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c(trans = "log10", option = "C") +
  geom_abline(slope = 1, intercept = 0, color = "white", linetype = "dashed") +
  labs(title = "Heatmap of hnRNPC Peak RPM (RBM25 WT vs. RBM25 Mut)",
       x = "Log(RPM) in RBM25 WT",
       y = "Log(RPM) in RBM25 Mut",
       fill = "Frequency") +
  coord_fixed() +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))

ggplot(hnRNPC_overlapPeaks_RPM, aes(x = inRBM25WT, y = inRBM25Mut)) +
  geom_point(alpha = 0.4) +
  # Add a diagonal y=x line for reference. Points above the line are higher in Mut.
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  # Use a log scale for both axes to spread out the dense cluster at the origin
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "hnRNPC Peak RPM in RBM25 WT vs. in RBM25 Mut)",
       x = "Log(RPM) in RBM25 WT",
       y = "Log(RPM) in RBM25 Mut") +
  coord_fixed() + # Ensure a 1:1 aspect ratio for the axes
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))


################################################################################

## Peak Overlap with RBM25 Mutant
################################################################################
# RBM25_BC_new = 1

Peak_RBM25_WT_new = fread(paste0(baseDir, '/RBM25_WT_new_Combined_peakCoverage_filtered_annotated.txt'))
colnames(Peak_RBM25_WT_new)[colnames(Peak_RBM25_WT_new) == "chrom"] = "chr"
Peak_RBM25_WT_new = Peak_RBM25_WT_new %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_RBM25_WT_new = Peak_RBM25_WT_new %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                                            'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
# Peak_RBM25_WT_new = Peak_RBM25_WT_new %>% filter(avg_RPM >= median(Peak_RBM25_WT_new$avg_RPM))
# Peak_RBM25_WT_new = Peak_RBM25_WT_new %>% filter(BC >= RBM25_BC_new)

Peak_RBM25_Mut_new = fread(paste0(baseDir, '/RBM25_Mut_new_Combined_peakCoverage_filtered_annotated.txt'))
colnames(Peak_RBM25_Mut_new)[colnames(Peak_RBM25_Mut_new) == "chrom"] = "chr"
Peak_RBM25_Mut_new = Peak_RBM25_Mut_new %>% mutate(chr = case_when(chr == 'chrMT' ~ 'chrM', TRUE ~ as.character(chr)))
Peak_RBM25_Mut_new = Peak_RBM25_Mut_new %>% filter(chr %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                                              'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))
# Peak_RBM25_Mut_new = Peak_RBM25_Mut_new %>% filter(avg_RPM >= median(Peak_RBM25_Mut_new$avg_RPM))
# Peak_RBM25_Mut_new = Peak_RBM25_Mut_new %>% filter(BC >= RBM25_BC_new)

GR_Peak_RBM25_WT = peaks2GR(Peak_RBM25_WT_new)
GR_Peak_RBM25_Mut = peaks2GR(Peak_RBM25_Mut_new)

final_peaks = Reduce(intersect, list(GR_Peak_hnRNPC_inRBM25WT, GR_Peak_hnRNPC_inRBM25Mut, GR_Peak_RBM25_Mut))
export(final_peaks, paste0(outputDir, 'competition_CLIP_vs_RBM25Mut.bed'), format = "bed")

overlaps_hnRNPC_inRBM25WT = findOverlaps(GR_Peak_hnRNPC_inRBM25WT, final_peaks)
overlaps_hnRNPC_inRBM25Mut = findOverlaps(GR_Peak_hnRNPC_inRBM25Mut, final_peaks)

# Use the queryHits to get the row indices from the original data frames
overlaps_hnRNPC_inRBM25WT_RPM = Peak_hnRNPC_inRBM25WT$avg_RPM[queryHits(overlaps_hnRNPC_inRBM25WT)]
overlaps_hnRNPC_inRBM25Mut_RPM = Peak_hnRNPC_inRBM25Mut$avg_RPM[queryHits(overlaps_hnRNPC_inRBM25Mut)]

# Create a new data frame specifically for plotting the RPM values
plot_df = data.frame(
  avg_RPM = c(overlaps_hnRNPC_inRBM25WT_RPM, overlaps_hnRNPC_inRBM25Mut_RPM),
  condition = c(rep("WT", length(overlaps_hnRNPC_inRBM25WT_RPM)),
                rep("Mutant", length(overlaps_hnRNPC_inRBM25Mut_RPM)))
)

ggplot(plot_df, aes(x = condition, y = avg_RPM, fill = condition)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) + # Hides outlier points for clarity
  # Use a log scale if the RPM values have a wide range
  scale_y_log10() +
  # geom_jitter(width = 0.1, alpha = 0.1) + # Optional: shows individual points
  labs(
    title = "Comparison of RPM for Subsetted Peaks",
    x = "Condition",
    y = "Average RPM (Log Scale)"
  ) +
  coord_fixed() + # Ensure a 1:1 aspect ratio for the axes
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))

ggplot(plot_df, aes(x = condition, y = avg_RPM, fill = condition)) +
  geom_violin(trim = FALSE) +
  # Optional: Add a small boxplot inside the violin
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # Use a log scale to better visualize data with a wide range
  scale_y_log10() +
  labs(
    title = "Distribution of RPM for Subsetted Peaks",
    x = "Condition",
    y = "Average RPM (Log Scale)"
  ) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))


overlaps = findOverlaps(GR_Peak_hnRNPC_inRBM25WT, GR_Peak_hnRNPC_inRBM25Mut, minoverlap = 15)
scatter_df = data.frame(
  WT_RPM = Peak_hnRNPC_inRBM25WT$avg_RPM[queryHits(overlaps)],
  Mut_RPM = Peak_hnRNPC_inRBM25Mut$avg_RPM[subjectHits(overlaps)]
)

pseudocount = 1
scatter_df_colored = scatter_df %>%
  mutate(change = case_when(
    (Mut_RPM + pseudocount) / (WT_RPM + pseudocount) >= 1.5 ~ "Upregulated in Mutant",
    (WT_RPM + pseudocount) / (Mut_RPM + pseudocount) >= 1.5 ~ "Upregulated in WT",
    TRUE ~ "No significant change"
  )
  )

# common_peaks_gr = Reduce(intersect, list(GR_Peak_hnRNPC_inRBM25WT, GR_Peak_hnRNPC_inRBM25Mut, GR_Peak_RBM25_Mut))
# overlaps_hnRNPC_inRBM25WT = findOverlaps(GR_Peak_hnRNPC_inRBM25WT, common_peaks_gr, minoverlap = 15)
# overlaps_hnRNPC_inRBM25Mut = findOverlaps(GR_Peak_hnRNPC_inRBM25Mut, common_peaks_gr, minoverlap = 15)
#
# overlaps_hnRNPC_inRBM25WT = data.frame(
#   common_idx = subjectHits(overlaps_hnRNPC_inRBM25WT),
#   WT_RPM = Peak_hnRNPC_inRBM25WT$avg_RPM[queryHits(overlaps_hnRNPC_inRBM25WT)])
#
# overlaps_hnRNPC_inRBM25Mut = data.frame(
#   common_idx = subjectHits(overlaps_hnRNPC_inRBM25Mut),
#   Mut_RPM = Peak_hnRNPC_inRBM25Mut$avg_RPM[queryHits(overlaps_hnRNPC_inRBM25Mut)])
#
# overlaps_hnRNPC_inRBM25WT = aggregate(WT_RPM ~ common_idx, data = overlaps_hnRNPC_inRBM25WT, FUN = mean)
# overlaps_hnRNPC_inRBM25Mut = aggregate(Mut_RPM ~ common_idx, data = overlaps_hnRNPC_inRBM25Mut, FUN = mean)
# scatter_df = merge(overlaps_hnRNPC_inRBM25WT, overlaps_hnRNPC_inRBM25Mut, by = "common_idx")
#
# pseudocount = 1
# scatter_df_colored = scatter_df %>%
#   mutate(change = case_when(
#       (Mut_RPM + pseudocount) / (WT_RPM + pseudocount) >= 1.5 ~ "Upregulated in Mutant",
#       (WT_RPM + pseudocount) / (Mut_RPM + pseudocount) >= 1.5 ~ "Upregulated in WT",
#       TRUE                                                  ~ "No significant change"
#     )
#   )

# 3. Create the scatter plot
ggplot(scatter_df_colored, aes(x = WT_RPM, y = Mut_RPM, color = change)) +
  geom_point(alpha = 0.6) +
  # Add a y=x line for reference
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  # Use log scales to spread out points
  scale_x_log10() +
  scale_y_log10() +
  # Manually set the colors for each category
  scale_color_manual(
    name = "Fold Change",
    values = c(
      "Upregulated in Mutant" = "red",
      "Upregulated in WT" = "blue",
      "No significant change" = "grey"
    )
  ) +
  labs(
    title = "hnRNPC Peak RPM Comparison (RBM25 WT vs. RBM25 Mutant)",
    x = "hnRNPC log(RPM) per Peak in RBM25 WT",
    y = "hnRNPC log(RPM) per Peak in RBM25 Mut"
  ) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'))

# Find the row indices from your scatter plot data for each category
mutant_up_indices = which(scatter_df_colored$change == "Upregulated in Mutant")
wt_up_indices = which(scatter_df_colored$change == "Upregulated in WT")

# Use these indices to get the corresponding indices from the original GRanges object
gr_upregulated_in_mutant = GR_Peak_hnRNPC_inRBM25WT[queryHits(overlaps)[mutant_up_indices]]
gr_upregulated_in_wt = GR_Peak_hnRNPC_inRBM25WT[queryHits(overlaps)[wt_up_indices]]

export(gr_upregulated_in_mutant, paste0(outputDir, "upregulated_in_mutant.bed"), format = "bed")
export(gr_upregulated_in_wt, paste0(outputDir, "upregulated_in_wt.bed"), format = "bed")

# Get the original row indices from Peak_hnRNPC_inRBM25WT for each category
wt_indices_for_wt_up = queryHits(overlaps)[wt_up_indices]
wt_indices_for_mut_up = queryHits(overlaps)[mutant_up_indices]

# Subset the original data frame to get the gene names and coordinates
wt_up_peaks_with_genes = Peak_hnRNPC_inRBM25WT[
  wt_indices_for_wt_up,
  c("chr", "start", "end", "score", "strand", "gene", "external_gene_name", "finalized_annotation", "BC",
    "Luna_4460_Pool_hnRNPC_inRBM_WT1", "Luna_4460_Pool_hnRNPC_inRBM_WT2", "Luna_4460_Pool_hnRNPC_inRBM_WT3",
    "Luna_4460_Pool_hnRNPC_inRBM_WT1_RPM", "Luna_4460_Pool_hnRNPC_inRBM_WT2_RPM", "Luna_4460_Pool_hnRNPC_inRBM_WT3_RPM", "avg_RPM")
]

mut_up_peaks_with_genes = Peak_hnRNPC_inRBM25Mut[
  wt_indices_for_mut_up,
  c("chr", "start", "end", "score", "strand", "gene", "external_gene_name", "finalized_annotation", "BC",
    "Luna_4460_Pool_hnRNPC_inRBM_MUT1", "Luna_4460_Pool_hnRNPC_inRBM_MUT2", "Luna_4460_Pool_hnRNPC_inRBM_MUT3",
    "Luna_4460_Pool_hnRNPC_inRBM_MUT1_RPM", "Luna_4460_Pool_hnRNPC_inRBM_MUT2_RPM", "Luna_4460_Pool_hnRNPC_inRBM_MUT3_RPM", "avg_RPM")
]

################################################################################
