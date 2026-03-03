################################################################################
## CLIP Peak Filtering
## Author: Soon Yi with Antigravity
## Date: February 2026
################################################################################

# === Load Libraries ===
library(data.table)
library(tidyverse)
library(tibble)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggsignif)

## 1. Basic Setup:
################################################################################
# BASE_DIR = "F:/Specificity/CLIP/Analysis/"
# INPUT_DIR = paste0(BASE_DIR, "output/")
# OUTPUT_DIR = paste0(BASE_DIR, "output/")

BASE_DIR = "/Volumes/1TB_Data/Specificity/CLIP/Analysis/"
INPUT_DIR = paste0(BASE_DIR, "output2/")
OUTPUT_DIR = paste0(BASE_DIR, "output2/")

PEAK_MATRIX_FILE = paste0(INPUT_DIR, "meta_peakMatrix_normalized_annotated_hnRNPC_eCLIP.txt")

################################################################################

## 2. Load and filter annotated peak matrix:
################################################################################
peaksMatrix = fread(PEAK_MATRIX_FILE)
# peaksMatrix = peaksMatrix %>% filter(grouped_annotation != "unannotated")

peaks_hnRNPC = peaksMatrix %>% filter(BC_hnRNPC_eCLIP == 2 & BC_Control == 0)
peaks_hnRNPC = peaks_hnRNPC %>% filter(nTC_hnRNPC_eCLIP >= quantile(peaks_hnRNPC$nTC_hnRNPC_eCLIP)[3])

fwrite(peaks_hnRNPC, paste0(OUTPUT_DIR, "peaks_hnRNPC_eCLIP", ".txt"), sep = "\t")

################################################################################

## 3. Load peak matrix:
################################################################################
peaks_hnRNPC_eCLIP = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_eCLIP.txt"))

################################################################################

## 4. Visualization - Stacked Bar Graphs:
################################################################################
## Custom Functions 
# Get Annotation counts:
countAnnotation = function(peak_matrix, annotation_column, new_column_name = NULL, annotation_to_skip = NULL, fraction = NULL) {
  temp = data.frame(table(peak_matrix[[annotation_column]]), row.names = 1)
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

fillAnnotation = function(annotation_counts, annotation_list) {
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

## Global Variables:
CLIP_List = c('hnRNPC_eCLIP')
All_Annotation_List = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", 'ncRNA', "TE", "Other", "retained_intron", 'unannotated')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_hnRNPC_eCLIP = countAnnotation(peaks_hnRNPC_eCLIP, 'grouped_annotation', 'hnRNPC_eCLIP', 'unannotated')

PeakDistribution_combined = PeakDistribution_hnRNPC_eCLIP
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>% gather(key = "Source", value = "Freq", all_of(CLIP_List)) %>% dplyr::select(Source, Freq, Annotation)
PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = CLIP_List)
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = All_Annotation_List)

peakDistribution_hnRNPC = plotStackedBar(PeakDistribution_combined, c('hnRNPC_eCLIP'), c('HNRNPC eCLIP'), 'HNRNPC eCLIP Peak Distribution') 

## Add proportional plot for all samples:
PeakDistribution_prop = PeakDistribution_combined %>% 
  group_by(Source) %>% 
  mutate(Freq = Freq / sum(Freq)) %>% 
  ungroup()

peakDistribution_group = plotStackedBar(PeakDistribution_prop, CLIP_List, 
                                        c('hnRNPC_eCLIP'), 
                                        'All Samples Peak Distribution (Proportion)')

# Display plots in RStudio
print(peakDistribution_hnRNPC)
print(peakDistribution_group)

################################################################################









## 1. Basic Setup:
################################################################################
# BASE_DIR = "F:/Specificity/CLIP/Analysis/"
# INPUT_DIR = paste0(BASE_DIR, "output/")
# OUTPUT_DIR = paste0(BASE_DIR, "output/")

BASE_DIR = "/Volumes/1TB_Data/Specificity/CLIP/Analysis/"
INPUT_DIR = paste0(BASE_DIR, "output2/")
OUTPUT_DIR = paste0(BASE_DIR, "output2/")

PEAK_MATRIX_FILE = paste0(INPUT_DIR, "meta_peakMatrix_normalized_annotated_specificity_CLIP.txt")

################################################################################

## 2. Load and filter annotated peak matrix:
################################################################################
peaksMatrix = fread(PEAK_MATRIX_FILE)
# peaksMatrix = peaksMatrix %>% filter(grouped_annotation != "unannotated")

peaks_hnRNPC_WT = peaksMatrix %>% filter(BC_hnRNPC_WT >= 4)
peaks_hnRNPC_WT = peaks_hnRNPC_WT %>% filter(nTC_hnRNPC_WT >= quantile(peaks_hnRNPC_WT$nTC_hnRNPC_WT)[3])

peaks_hnRNPC_Mut = peaksMatrix %>% filter(BC_hnRNPC_Mut >= 4)
peaks_hnRNPC_Mut = peaks_hnRNPC_Mut %>% filter(nTC_hnRNPC_Mut >= quantile(peaks_hnRNPC_Mut$nTC_hnRNPC_Mut)[3])

peaks_hnRNPC_WT_inRBM25_WT = peaksMatrix %>% filter(BC_hnRNPC_WT_inRBM25_WT >= 3)
peaks_hnRNPC_WT_inRBM25_WT = peaks_hnRNPC_WT_inRBM25_WT %>% filter(nTC_hnRNPC_WT_inRBM25_WT >= 3*quantile(peaks_hnRNPC_WT_inRBM25_WT$nTC_hnRNPC_WT_inRBM25_WT)[4])

peaks_hnRNPC_WT_inRBM25_Mut = peaksMatrix %>% filter(BC_hnRNPC_WT_inRBM25_Mut >= 3)
peaks_hnRNPC_WT_inRBM25_Mut = peaks_hnRNPC_WT_inRBM25_Mut %>% filter(nTC_hnRNPC_WT_inRBM25_Mut >= 3*quantile(peaks_hnRNPC_WT_inRBM25_Mut$nTC_hnRNPC_WT_inRBM25_Mut)[4])

peaks_hnRNPC_WT_inKD = peaksMatrix %>% filter(BC_hnRNPC_WT_inKD >= 3)
peaks_hnRNPC_WT_inKD = peaks_hnRNPC_WT_inKD %>% filter(nTC_hnRNPC_WT_inKD >= quantile(peaks_hnRNPC_WT_inKD$nTC_hnRNPC_WT_inKD)[3])

peaks_hnRNPC_Mut_inKD = peaksMatrix %>% filter(BC_hnRNPC_Mut_inKD >= 3)
peaks_hnRNPC_Mut_inKD = peaks_hnRNPC_Mut_inKD %>% filter(nTC_hnRNPC_Mut_inKD >= quantile(peaks_hnRNPC_Mut_inKD$nTC_hnRNPC_Mut_inKD)[3])

peaks_RBM25_WT = peaksMatrix %>% filter(BC_RBM25_WT >= 3)
peaks_RBM25_WT = peaks_RBM25_WT %>% filter(nTC_RBM25_WT >= quantile(peaks_RBM25_WT$nTC_RBM25_WT)[3])

peaks_RBM25_Mut = peaksMatrix %>% filter(BC_RBM25_Mut >= 3)
peaks_RBM25_Mut = peaks_RBM25_Mut %>% filter(nTC_RBM25_Mut >= quantile(peaks_RBM25_Mut$nTC_RBM25_Mut)[3])

fwrite(peaks_hnRNPC_WT, paste0(OUTPUT_DIR, "peaks_hnRNPC_WT", ".txt"), sep = "\t")
fwrite(peaks_hnRNPC_Mut, paste0(OUTPUT_DIR, "peaks_hnRNPC_Mut", ".txt"), sep = "\t")
fwrite(peaks_hnRNPC_WT_inRBM25_WT, paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inRBM25_WT", ".txt"), sep = "\t")
fwrite(peaks_hnRNPC_WT_inRBM25_Mut, paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inRBM25_Mut", ".txt"), sep = "\t")
fwrite(peaks_hnRNPC_WT_inKD, paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inKD", ".txt"), sep = "\t")
fwrite(peaks_hnRNPC_Mut_inKD, paste0(OUTPUT_DIR, "peaks_hnRNPC_Mut_inKD", ".txt"), sep = "\t")
fwrite(peaks_RBM25_WT, paste0(OUTPUT_DIR, "peaks_RBM25_WT", ".txt"), sep = "\t")
fwrite(peaks_RBM25_Mut, paste0(OUTPUT_DIR, "peaks_RBM25_Mut", ".txt"), sep = "\t")

################################################################################

## 3. Load peak matrix:
################################################################################
peaks_hnRNPC_WT = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_WT.txt"))
peaks_hnRNPC_Mut = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_Mut.txt"))
peaks_hnRNPC_WT_inRBM25_WT = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inRBM25_WT.txt"))
peaks_hnRNPC_WT_inRBM25_Mut = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inRBM25_Mut.txt"))
peaks_hnRNPC_WT_inKD = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inKD.txt"))
peaks_hnRNPC_Mut_inKD = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_Mut_inKD.txt"))
peaks_RBM25_WT = fread(paste0(OUTPUT_DIR, "peaks_RBM25_WT.txt"))
peaks_RBM25_Mut = fread(paste0(OUTPUT_DIR, "peaks_RBM25_Mut.txt"))
################################################################################

## 4. Visualization - Stacked Bar Graphs:
################################################################################
## Custom Functions 
# Get Annotation counts:
countAnnotation = function(peak_matrix, annotation_column, new_column_name = NULL, annotation_to_skip = NULL, fraction = NULL) {
  temp = data.frame(table(peak_matrix[[annotation_column]]), row.names = 1)
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

fillAnnotation = function(annotation_counts, annotation_list) {
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

## Global Variables:
CLIP_List = c('hnRNPC_WT', 'hnRNPC_Mut', 'hnRNPC_WT_inRBM25_WT', 'hnRNPC_WT_inRBM25_Mut', 'hnRNPC_WT_inKD', 'hnRNPC_Mut_inKD', 'RBM25_WT', 'RBM25_Mut')
All_Annotation_List = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", 'ncRNA', "TE", "Other", "retained_intron", 'unannotated')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_hnRNPC_WT = countAnnotation(peaks_hnRNPC_WT, 'grouped_annotation', 'hnRNPC_WT', 'unannotated')
PeakDistribution_hnRNPC_Mut = countAnnotation(peaks_hnRNPC_Mut, 'grouped_annotation', 'hnRNPC_Mut', 'unannotated')
PeakDistribution_RBM25_WT = countAnnotation(peaks_RBM25_WT, 'grouped_annotation', 'RBM25_WT', 'unannotated')
PeakDistribution_RBM25_Mut = countAnnotation(peaks_RBM25_Mut, 'grouped_annotation', 'RBM25_Mut', 'unannotated')
PeakDistribution_hnRNPC_WT_inRBM25_WT = countAnnotation(peaks_hnRNPC_WT_inRBM25_WT, 'grouped_annotation', 'hnRNPC_WT_inRBM25_WT', 'unannotated')
PeakDistribution_hnRNPC_WT_inRBM25_Mut = countAnnotation(peaks_hnRNPC_WT_inRBM25_Mut, 'grouped_annotation', 'hnRNPC_WT_inRBM25_Mut', 'unannotated')
PeakDistribution_hnRNPC_WT_inKD = countAnnotation(peaks_hnRNPC_WT_inKD, 'grouped_annotation', 'hnRNPC_WT_inKD', 'unannotated')
PeakDistribution_hnRNPC_Mut_inKD = countAnnotation(peaks_hnRNPC_Mut_inKD, 'grouped_annotation', 'hnRNPC_Mut_inKD', 'unannotated')

PeakDistribution_hnRNPC_Mut_inKD = fillAnnotation(PeakDistribution_hnRNPC_Mut_inKD, rownames(PeakDistribution_hnRNPC_WT))
PeakDistribution_RBM25_WT = fillAnnotation(PeakDistribution_RBM25_WT, rownames(PeakDistribution_hnRNPC_WT))
PeakDistribution_RBM25_Mut = fillAnnotation(PeakDistribution_RBM25_Mut, rownames(PeakDistribution_hnRNPC_WT))

PeakDistribution_combined = cbind(PeakDistribution_hnRNPC_WT, PeakDistribution_hnRNPC_Mut, 
                                  PeakDistribution_hnRNPC_WT_inRBM25_WT, PeakDistribution_hnRNPC_WT_inRBM25_Mut,
                                  PeakDistribution_hnRNPC_WT_inKD, PeakDistribution_hnRNPC_Mut_inKD,
                                  PeakDistribution_RBM25_WT, PeakDistribution_RBM25_Mut)
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>% gather(key = "Source", value = "Freq", all_of(CLIP_List)) %>% dplyr::select(Source, Freq, Annotation)
PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = CLIP_List)
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = All_Annotation_List)

peakDistribution_hnRNPC = plotStackedBar(PeakDistribution_combined, c('hnRNPC_WT', 'hnRNPC_Mut'), c('HNRNPC WT', 'HNRNPC Mut'), 'HNRNPC WT vs Mut Peak Distribution') 
peakDistribution_hnRNPC_inRBM25 = plotStackedBar(PeakDistribution_combined, c('hnRNPC_WT_inRBM25_WT', 'hnRNPC_WT_inRBM25_Mut'), c('HNRNPC in RBM25 WT', 'HNRNPC in RBM25 Mut'), 'HNRNPC in RBM25 WT vs Mut Peak Distribution') 
peakDistribution_hnRNPC_inKD = plotStackedBar(PeakDistribution_combined, c('hnRNPC_WT_inKD', 'hnRNPC_Mut_inKD'), c('HNRNPC WT in KD', 'HNRNPC Mut in KD'), 'HNRNPC WT vs Mut in KD Peak Distribution')
peakDistribution_RBM25 = plotStackedBar(PeakDistribution_combined, c('RBM25_WT', 'RBM25_Mut'), c('RBM25 WT', 'RBM25 Mut'), 'RBM25 WT vs Mut Peak Distribution')

## Add proportional plot for all samples:
PeakDistribution_prop = PeakDistribution_combined %>% 
  group_by(Source) %>% 
  mutate(Freq = Freq / sum(Freq)) %>% 
  ungroup()

peakDistribution_group = plotStackedBar(PeakDistribution_prop, CLIP_List, 
                                        c('HNRNPC WT', 'HNRNPC Mut', 'HNRNPC WT in RBM25 WT', 'HNRNPC WT in RBM25 Mut', 'HNRNPC WT in KD', 'HNRNPC Mut in KD', 'RBM25 WT', 'RBM25 Mut'), 
                                        'All Samples Peak Distribution (Proportion)')

# Display plots in RStudio
print(peakDistribution_hnRNPC)
print(peakDistribution_hnRNPC_inRBM25)
print(peakDistribution_hnRNPC_inKD)
print(peakDistribution_RBM25)
print(peakDistribution_group)

################################################################################
