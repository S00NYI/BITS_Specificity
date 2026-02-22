## CoCLIP Analysis: 
## Specificity Analysis
## Written by Soon Yi
## Use Homer output for subsequent analysis.
## Last Edit: 2023-11-15

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

## Custom Functions
####################################################################################################################
## Filter peaks based on the designated criteria:
filterPeakMatrix = function(peak_matrix, sample_list, info_columns, BC_criteria, rowSum_criteria = NULL) {
  temp = peak_matrix[, c(info_columns, sample_list)]
  
  if (!is.null(rowSum_criteria)) {
    temp$rowAvg = rowSums(temp[, sample_list]) / length(sample_list)
    temp = temp[(rowSums(temp[, paste0(unique(sub("_[0-9]+$", "", sample_list)), '_BC')]) >= BC_criteria) & (temp$rowAvg > (median(temp$rowAvg) * rowSum_criteria)), ]
    temp$rowAvg = NULL
  } else {
    temp = temp[(rowSums(temp[, paste0(unique(sub("_[0-9]+$", "", sample_list)), '_BC')]) >= BC_criteria), ]
  }
  
  return(temp)
}

## Count motif occurrence:
motifCounts = function(peak_matrix, motifs) {
  result_container = list()
  for (motif in motifs) {
    org_motif = motif
    motif = gsub("U", "T", motif)
    motif_positions = sapply(motif, function(motif) gregexpr(paste0("(?=", motif, ")"), peak_matrix$sequence, perl=TRUE))
    motif_positions = lapply(motif_positions, function(x) { attributes(x) <- NULL; x })
    counts = numeric(length(motif_positions))
    for (i in 1:length(motif_positions)) {
      sublist = motif_positions[[i]]
      if (length(sublist) == 1 && sublist[[1]] == -1) {
        counts[i] = 0
      } else {
        counts[i] = length(sublist)
      }
    }
    result_container[[org_motif]] = counts
  }
  result_df = data.frame(result_container)
  colnames(result_df) = motifs
  return(result_df)
}

## Min-Max normalization:
min_max_norm = function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

## Specificity index calculation:
SI_Calculation = function(x) {
  x / (median(x, na.rm = FALSE))
}
####################################################################################################################

## Load peak matrix and clean up:
####################################################################################################################
peaksMatrix_PATH = '/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/'
peaksMatrix_FILE = 'Combined_peakCoverage_groomed_normalized_annotated.txt'

peaksMatrix = read_delim(paste0(peaksMatrix_PATH, peaksMatrix_FILE), show_col_types = FALSE)
peaksMatrix = peaksMatrix %>% mutate_at('TOTAL_BC', as.numeric)
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'unannotated', 'UnAn', finalized_annotation))
peaksMatrix = peaksMatrix %>% mutate(grouped_annotation = ifelse(grouped_annotation == 'unannotated', 'UnAn', grouped_annotation))
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'downstream 10K', 'DS10K', finalized_annotation))
peaksMatrix = peaksMatrix %>% mutate(grouped_annotation = ifelse(grouped_annotation == 'downstream 10K', 'DS10K', grouped_annotation))
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'ncRNA_Retained_intron', 'nC_RI', finalized_annotation))
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'CDS_Retained_intron', 'CDS_RI', finalized_annotation))

## Column organization:
inert_columns = c('chrom', 'start', 'end', 'peak_names', 'score', 'strand', 
                  "gene", "external_gene_name", "annotation", "finalized_annotation", "grouped_annotation", "annotation_count", "TOTAL_TagCount")
BC_columns = c("Nuc_F_M_BC", "Nuc_F_S_BC", "Cyto_F_M_BC", "Cyto_F_S_BC", 
               "NLS_I_M_BC", "NLS_I_S_BC", "NES_I_M_BC", "NES_I_S_BC", "G3BP_I_M_BC", "G3BP_I_S_BC",
               "NLS_E_M_BC", "NLS_E_S_BC", "NES_E_M_BC", "NES_E_S_BC", "G3BP_E_M_BC", "G3BP_E_S_BC",
               "TOTAL_BC")

Nuc_F_M = c('Nuc_F_M_1', 'Nuc_F_M_2', 'Nuc_F_M_3')
Nuc_F_S = c('Nuc_F_S_1', 'Nuc_F_S_2', 'Nuc_F_S_3')
Cyto_F_M = c('Cyto_F_M_1', 'Cyto_F_M_2', 'Cyto_F_M_3')
Cyto_F_S = c('Cyto_F_S_1', 'Cyto_F_S_2', 'Cyto_F_S_3')

NLS_I_M = c('NLS_I_M_1', 'NLS_I_M_2')
NES_I_M = c('NES_I_M_1', 'NES_I_M_2')
G3BP_I_M = c('G3BP_I_M_1', 'G3BP_I_M_2', 'G3BP_I_M_3', 'G3BP_I_M_4')

NLS_I_S = c('NLS_I_S_1', 'NLS_I_S_2')
NES_I_S = c('NES_I_S_1', 'NES_I_S_2')
G3BP_I_S = c('G3BP_I_S_1', 'G3BP_I_S_2', 'G3BP_I_S_3', 'G3BP_I_S_4', 'G3BP_I_S_5')

NLS_E_M = c('NLS_E_M_1', 'NLS_E_M_2', 'NLS_E_M_3', 'NLS_E_M_4')
NLS_E_S = c('NLS_E_S_1', 'NLS_E_S_2', 'NLS_E_S_3', 'NLS_E_S_4')
NES_E_M = c('NES_E_M_1', 'NES_E_M_2', 'NES_E_M_3', 'NES_E_M_4')
NES_E_S = c('NES_E_S_1', 'NES_E_S_2', 'NES_E_S_3', 'NES_E_S_4')
G3BP_E_M = c('G3BP_E_M_1', 'G3BP_E_M_2', 'G3BP_E_M_3', 'G3BP_E_M_4', 'G3BP_E_M_5', 'G3BP_E_M_6')
G3BP_E_S = c('G3BP_E_S_1', 'G3BP_E_S_2', 'G3BP_E_S_3', 'G3BP_E_S_4', 'G3BP_E_S_5', 'G3BP_E_S_6', 'G3BP_E_S_7')

## Add row sum columns for further filtering:
peaksMatrix$F_rowSum = rowSums(peaksMatrix[, c(Nuc_F_M, Nuc_F_S, Cyto_F_M, Cyto_F_S)])
peaksMatrix$I_rowSum = rowSums(peaksMatrix[, c(NLS_I_M, NLS_I_S, NES_I_M, NES_I_S, G3BP_I_M, G3BP_I_S)])
peaksMatrix$E_rowSum = rowSums(peaksMatrix[, c(NLS_E_M, NLS_E_S, NES_E_M, NES_E_S, G3BP_E_M, G3BP_E_S)])
rowSum_columns = c('F_rowSum', 'I_rowSum', 'E_rowSum')

## Add pseudocount:
pseudoCount = min(peaksMatrix[, colnames(peaksMatrix)[6:63]][peaksMatrix[, colnames(peaksMatrix)[6:63]] != 0], na.rm = TRUE)
peaksMatrix[, colnames(peaksMatrix)[6:63]] = peaksMatrix[, colnames(peaksMatrix)[6:63]] + pseudoCount

####################################################################################################################

## Filter Criteria:
####################################################################################################################
BC_Threshold_F = 2
BC_Threshold_I = 4
BC_Threshold_E = 2
BC_Threshold_E_SG = 3

rowSum_Multiplier_F = 4
rowSum_Multiplier_I = 2
rowSum_Multiplier_E = 2
####################################################################################################################

## Subset of Peaks for downstream analysis:
####################################################################################################################
## Fractionation CLIP
Peak_F_Nuc_M = filterPeakMatrix(peaksMatrix, Nuc_F_M, c(inert_columns, BC_columns), BC_Threshold_F, rowSum_Multiplier_F)
Peak_F_Nuc_S = filterPeakMatrix(peaksMatrix, Nuc_F_S, c(inert_columns, BC_columns), BC_Threshold_F, rowSum_Multiplier_F)
Peak_F_Cyt_M = filterPeakMatrix(peaksMatrix, Cyto_F_M, c(inert_columns, BC_columns), BC_Threshold_F, rowSum_Multiplier_F)
Peak_F_Cyt_S = filterPeakMatrix(peaksMatrix, Cyto_F_S, c(inert_columns, BC_columns), BC_Threshold_F, rowSum_Multiplier_F)

## CoCLIP
Peak_Co_Input_M = filterPeakMatrix(peaksMatrix, c(NLS_I_M, NES_I_M, G3BP_I_M), c(inert_columns, BC_columns), BC_Threshold_I, rowSum_Multiplier_I)
Peak_Co_Input_S = filterPeakMatrix(peaksMatrix, c(NLS_I_S, NES_I_S, G3BP_I_S), c(inert_columns, BC_columns), BC_Threshold_I, rowSum_Multiplier_I)
Peak_Co_NLS_M = filterPeakMatrix(peaksMatrix, NLS_E_M, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_NLS_S = filterPeakMatrix(peaksMatrix, NLS_E_S, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_NES_M = filterPeakMatrix(peaksMatrix, NES_E_M, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_NES_S = filterPeakMatrix(peaksMatrix, NES_E_S, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_G3BP_M = filterPeakMatrix(peaksMatrix, G3BP_E_M, c(inert_columns, BC_columns), BC_Threshold_E_SG, rowSum_Multiplier_E)
Peak_Co_G3BP_S = filterPeakMatrix(peaksMatrix, G3BP_E_S, c(inert_columns, BC_columns), BC_Threshold_E_SG, rowSum_Multiplier_E)
####################################################################################################################

## Build peak sequence table and motif enrichment ranks based on peaksMatrix:
####################################################################################################################
baseDir = '/Users/soonyi/Desktop/Genomics/Specificity/Data/'
extensions = c(0, 15, 40, 65, 90)
Ks = c(5)
# Ks = c(4, 6, 7)
NTs = c("A", "C", "G", "T")

for (extension in extensions) {
  print(paste0('Working on ', extension, 'nt extension for coCLIP....'))
  peaksGR = peaksMatrix[, c('chrom', 'start', 'end', 'strand', 'peak_names')]
  peaksGR = peaksGR %>% mutate(chrom = ifelse(chrom == "chrMT", "chrM", chrom))
  peaksGR$start = as.integer(peaksGR$start) + 1 - extension
  peaksGR$end = as.integer(peaksGR$end) + extension
  peaksGR = GRanges(peaksGR)
  
  peaksGR_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, peaksGR, as.character = TRUE)
  peaksGR_seqs = cbind(peaksMatrix[, c('chrom', 'start', 'end', 'strand', 'peak_names', BC_columns, 'grouped_annotation', 'finalized_annotation')], data.frame(sequence = peaksGR_seqs))
  
  Samples = c('I_M', 'NLS_M', 'NES_M', 'G3BP_M', 'I_S','NLS_S', 'NES_S', 'G3BP_S')
  I_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_Input_M$peak_names))
  NLS_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NLS_M$peak_names))
  NES_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NES_M$peak_names))
  G3BP_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_G3BP_M$peak_names))
  I_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_Input_S$peak_names))
  NLS_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NLS_S$peak_names))
  NES_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NES_S$peak_names))
  G3BP_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_G3BP_S$peak_names))
  
  for (K in Ks) {
    print(paste0('Working on ', K, 'mer....'))
    # Make all K-mers
    combinations = expand.grid(replicate(K, NTs, simplify = FALSE))
    kmers = apply(combinations, 1, paste, collapse = "")
    
    raw_count = data.frame(MOTIF = kmers)
    nor_count = data.frame(MOTIF = kmers)
    
    for (sample in Samples) {
      INPUT = get(sample)
      raw_count[[sample]] = rep(0, 4^K)
      nor_count[[sample]] = rep(0, 4^K)
      num_peak = nrow(INPUT)
      for (kmer in kmers) {
        kmer_counts = colSums(motifCounts(INPUT, kmer))
        raw_count[which(raw_count$MOTIF == kmer), sample] = kmer_counts
        nor_count[which(nor_count$MOTIF == kmer), sample] = kmer_counts/num_peak
      }
    }
    write.csv(raw_count, paste0(baseDir, 'CLIP_Analysis/', K ,'mer/coCLIP_Raw_', K, 'mer_', extension ,'nt_PeakExtension_Counts.csv'), quote = F)
    write.csv(nor_count, paste0(baseDir, 'CLIP_Analysis/', K ,'mer/coCLIP_Normalized_', K, 'mer_', extension ,'nt_PeakExtension_Counts.csv'), quote = F)
    normalized_dataframe = nor_count
    normalized_dataframe = normalized_dataframe %>% mutate_if(is.numeric, min_max_norm) %>% mutate_if(is.numeric, round, 10)
    write.csv(normalized_dataframe, paste0(baseDir, 'CLIP_Analysis/', K ,'mer/coCLIP_Scaled_', K, 'mer_', extension ,'nt_PeakExtension_Counts_Normalized.csv'), quote = F)
    specificity_dataframe = normalized_dataframe %>% mutate_if(is.numeric, SI_Calculation) %>% mutate_if(is.numeric, round, 10)
    write.csv(specificity_dataframe, paste0(baseDir, 'CLIP_Analysis/', K ,'mer/coCLIP_', K, 'mer_', extension ,'nt_PeakExtension_SpecificityIndex.csv'), quote = F)
  }
}


eCLIP = read.delim('/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/motifEnrichment/eCLIP_bed/ELAVL1_eCLIP_ENCFF566LNK.bed', header = F)
eCLIP = eCLIP[, c('V1', 'V2', 'V3', 'V4', 'V6')]
colnames(eCLIP) = c('chrom', 'start', 'end', 'name', 'strand')
eCLIP = eCLIP %>% filter(chrom %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'))

for (extension in extensions) {
  print(paste0('Working on ', extension, 'nt extension for eCLIP....'))
  peaksGR = eCLIP[, c('chrom', 'start', 'end', 'strand', 'name')]
  peaksGR$start = as.integer(peaksGR$start) + 1 - extension
  peaksGR$end = as.integer(peaksGR$end) + extension
  peaksGR = GRanges(peaksGR)
  peaksGR_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, peaksGR, as.character = TRUE)
  peaksGR_seqs = cbind(eCLIP, data.frame(sequence = peaksGR_seqs))
  for (K in Ks) {
    print(paste0('Working on ', K, 'mer....'))
    # Make all K-mers
    combinations = expand.grid(replicate(K, NTs, simplify = FALSE))
    kmers = apply(combinations, 1, paste, collapse = "")
    
    raw_count_eCLIP = data.frame(MOTIF = kmers)
    nor_count_eCLIP = data.frame(MOTIF = kmers)
    raw_count_eCLIP$eCLIP = rep(0, 4^K)
    nor_count_eCLIP$eCLIP = rep(0, 4^K)
    num_peak = nrow(peaksGR_seqs)
    
    for (kmer in kmers) {
      kmer_counts = colSums(motifCounts(peaksGR_seqs, kmer))
      raw_count_eCLIP[which(raw_count_eCLIP$MOTIF == kmer), 'eCLIP'] = kmer_counts
      nor_count_eCLIP[which(nor_count_eCLIP$MOTIF == kmer), 'eCLIP'] = kmer_counts/num_peak
    }

    write.csv(raw_count_eCLIP, paste0(baseDir, 'CLIP_Analysis/', K ,'mer/eCLIP_Raw_', K, 'mer_', extension ,'nt_PeakExtension_Counts.csv'), quote = F)
    write.csv(nor_count_eCLIP, paste0(baseDir, 'CLIP_Analysis/', K ,'mer/eCLIP_Normalized_', K, 'mer_', extension ,'nt_PeakExtension_Counts.csv'), quote = F)
    normalized_dataframe_eCLIP = nor_count_eCLIP
    normalized_dataframe_eCLIP = normalized_dataframe_eCLIP %>% mutate_if(is.numeric, min_max_norm) %>% mutate_if(is.numeric, round, 10)
    write.csv(normalized_dataframe_eCLIP, paste0(baseDir, 'CLIP_Analysis/', K ,'mer/eCLIP_Scaled_', K, 'mer_', extension ,'nt_PeakExtension_Counts_Normalized.csv'), quote = F)
    specificity_dataframe_eCLIP = normalized_dataframe_eCLIP %>% mutate_if(is.numeric, SI_Calculation) %>% mutate_if(is.numeric, round, 10)
    write.csv(specificity_dataframe_eCLIP, paste0(baseDir, 'CLIP_Analysis/', K ,'mer/eCLIP_', K, 'mer_', extension ,'nt_PeakExtension_SpecificityIndex.csv'), quote = F)
  }
}

####################################################################################################################



## Plotting
####################################################################################################################
baseDir = '/Users/soonyi/Desktop/Genomics/Specificity/Data/'
extensions = c(0, 15, 40, 65, 90)
extension = extensions[3]

eCLIP_Scaled_5mer = read.delim(paste0(baseDir, 'CLIP_Analysis/5mer/eCLIP_Scaled_5mer_', extension, 'nt_PeakExtension_Counts_Normalized.csv'), sep = ',')
coCLIP_Scaled_5mer = read.delim(paste0(baseDir, 'CLIP_Analysis/5mer/coCLIP_Scaled_5mer_', extension, 'nt_PeakExtension_Counts_Normalized.csv'), sep = ',')

Scaled_5mer = data.frame(cbind(data.frame(MOTIF = eCLIP_Scaled_5mer$MOTIF, eCLIP = eCLIP_Scaled_5mer$eCLIP), coCLIP_Scaled_5mer[, c(-1, -2)]))

mock_data = Scaled_5mer %>%
  select(I_M, eCLIP, NLS_M, NES_M, G3BP_M) %>%
  pivot_longer(cols = c("eCLIP", "NLS_M", "NES_M", "G3BP_M"), names_to = "Condition", values_to = "Value")

mock_colors = c("eCLIP" = "ivory4", "NLS_M" = "skyblue", "NES_M" = "darkseagreen2", "G3BP_M" = "salmon")

ggplot(mock_data, aes(x = Value, y = I_M, color = Condition)) +
  geom_point(alpha = 1, shape = 16) +
  scale_color_manual(values = mock_colors) + 
  theme_minimal() +
  labs(x = "Condition Value", y = "I_M Value", title = "Scatterplot of I_M vs Conditions") +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(mock_data, aes(x = Value, fill = Condition)) +
  geom_histogram(binwidth = 0.005, alpha = 0.75, position = 'identity') +
  scale_fill_manual(values = mock_colors) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 400), 
                     breaks = seq(0, 400, by = 100))  +
  labs(x = "Value", y = "Count", title = "Overlapping Histograms") +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(Scaled_5mer, aes(I_M)) +
  geom_histogram(binwidth = 0.005, alpha = 0.75) +
  scale_y_continuous(limits = c(0, 400), 
                     breaks = seq(0, 400, by = 100))  +
  labs(x = "Value", y = "Count", title = "Input")  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

stress_data = Scaled_5mer %>%
  select(I_S, NLS_S, NES_S, G3BP_S) %>%
  pivot_longer(cols = c("NLS_S", "NES_S", "G3BP_S"), names_to = "Condition", values_to = "Value")

stress_colors = c("NLS_S" = "skyblue", "NES_S" = "darkseagreen2", "G3BP_S" = "salmon")

ggplot(stress_data, aes(x = Value, y = I_S, color = Condition)) +
  geom_point(alpha = 1, shape = 16) +
  scale_color_manual(values = stress_colors) + 
  theme_minimal() +
  labs(x = "Condition Value", y = "I_S Value", title = "Scatterplot of I_S vs Conditions") +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(stress_data, aes(x = Value, fill = Condition)) +
  geom_histogram(binwidth = 0.005, alpha = 0.75, position = 'identity') +
  scale_fill_manual(values = stress_colors) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, by = 100))  +
  labs(x = "Value", y = "Count", title = "Overlapping Histograms") +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

ggplot(Scaled_5mer, aes(I_S)) +
  geom_histogram(binwidth = 0.005, alpha = 0.75) +
  scale_y_continuous(limits = c(0, 400), 
                     breaks = seq(0, 400, by = 100))  +
  labs(x = "Value", y = "Count", title = "Input")  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

for (extension in extensions) {
  eCLIP_SI_5mer = read.delim(paste0(baseDir, 'CLIP_Analysis/5mer/eCLIP_5mer_', extension, 'nt_PeakExtension_SpecificityIndex.csv'), sep = ',')
  coCLIP_SI_5mer = read.delim(paste0(baseDir, 'CLIP_Analysis/5mer/coCLIP_5mer_', extension, 'nt_PeakExtension_SpecificityIndex.csv'), sep = ',')
  
  combined_df = data.frame(MOTIF = eCLIP_SI_5mer$MOTIF, eCLIP = eCLIP_SI_5mer$eCLIP, coCLIP_SI_5mer[, -c(1, 2)])
  df_name = paste0("SI_5mer_", (20 + extension*2), "nt")
  
  assign(df_name, combined_df, envir = .GlobalEnv)
}

dataframes = list(SI_5mer_20nt, SI_5mer_50nt, SI_5mer_100nt, SI_5mer_150nt, SI_5mer_200nt)
extract_values = function(df, motif) {
  return(df[df$MOTIF == motif, c('eCLIP', 'I_M', 'NLS_M', 'NES_M', 'G3BP_M', 'I_S', 'NLS_S', 'NES_S', 'G3BP_S')])
}
values_list = lapply(dataframes, extract_values, motif = 'AAAAA')
SI_5mer_AAAAA = do.call(rbind, values_list)
row.names(SI_5mer_AAAAA) = c('20nt', '50nt', '100nt', '150nt', '200nt')

values_list = lapply(dataframes, extract_values, motif = 'TTTTT')
SI_5mer_UUUUU = do.call(rbind, values_list)
row.names(SI_5mer_UUUUU) = c('20nt', '50nt', '100nt', '150nt', '200nt')

SI_5mer_AAAAA$ext = rownames(SI_5mer_AAAAA)
SI_5mer_AAAAA$ext = factor(SI_5mer_AAAAA$ext, levels = c("20nt", "50nt", "100nt", "150nt", "200nt"))
AAAAA_plotData = SI_5mer_AAAAA %>%
  pivot_longer(cols = -ext, names_to = "Column", values_to = "Value") %>%
  mutate(Column = factor(Column, levels = c('eCLIP', 'I_M', 'NLS_M', 'NES_M', 'G3BP_M', 'I_S', 'NLS_S', 'NES_S', 'G3BP_S')))

colors = c("grey", "black", "skyblue", "darkseagreen2", "salmon", "black", "skyblue", "darkseagreen2", "salmon")
line_types = c("solid", "solid", "solid", "solid", "solid", "dotted", "dotted", "dotted", "dotted")

# Create the line graph
ggplot(AAAAA_plotData, aes(x = ext, y = Value, group = Column, color = Column, linetype = Column)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = line_types) +
  theme_minimal() +
  labs(x = "Extensions (nt)", y = "SI", title = "AAAAA SI per Extension") +
  scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, by = 100))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

SI_5mer_UUUUU$ext = rownames(SI_5mer_UUUUU)
SI_5mer_UUUUU$ext = factor(SI_5mer_UUUUU$ext, levels = c("20nt", "50nt", "100nt", "150nt", "200nt"))
UUUUU_plotData = SI_5mer_UUUUU %>%
  pivot_longer(cols = -ext, names_to = "Column", values_to = "Value") %>%
  mutate(Column = factor(Column, levels = c('eCLIP', 'I_M', 'NLS_M', 'NES_M', 'G3BP_M', 'I_S', 'NLS_S', 'NES_S', 'G3BP_S')))

colors = c("grey", "black", "skyblue", "darkseagreen2", "salmon", "black", "skyblue", "darkseagreen2", "salmon")
line_types = c("solid", "solid", "solid", "solid", "solid", "dotted", "dotted", "dotted", "dotted")

# Create the line graph
ggplot(UUUUU_plotData, aes(x = ext, y = Value, group = Column, color = Column, linetype = Column)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = line_types) +
  theme_minimal() +
  labs(x = "Extensions (nt)", y = "SI", title = "UUUUU SI per Extension") +
  scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, by = 100))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

SI_selection = which((SI_5mer_50nt$I_M >= 1) &
                       (SI_5mer_50nt$NLS_M >= 1) &
                       (SI_5mer_50nt$NES_M >= 1) &
                       (SI_5mer_50nt$G3BP_M >= 1) &
                       (SI_5mer_50nt$I_S >= 1) &
                       (SI_5mer_50nt$NLS_S >= 1) &
                       (SI_5mer_50nt$NES_S >= 1) &
                       (SI_5mer_50nt$G3BP_S >= 1))

SI_FC = data.frame(MOTIF = SI_5mer_50nt$MOTIF,
                       Input = (SI_5mer_50nt$I_S + 1e-10)/(SI_5mer_50nt$I_M + 1e-10),
                       NLS = (SI_5mer_50nt$NLS_S + 1e-10)/(SI_5mer_50nt$NLS_M + 1e-10),
                       NES = (SI_5mer_50nt$NES_S + 1e-10)/(SI_5mer_50nt$NES_M + 1e-10),
                       G3BP = (SI_5mer_50nt$G3BP_S + 1e-10)/(SI_5mer_50nt$G3BP_M + 1e-10))

SI_FC = SI_FC[SI_selection, ]
SI_L2FC = SI_FC %>% mutate_if(is.numeric, log2)
SI_5mer_50nt_Selection = SI_5mer_50nt[SI_selection, ]
SI_5mer_50nt_Selection_L10FC = SI_5mer_50nt_Selection %>% mutate_if(is.numeric, log10)
####################################################################################################################

## Locus specific analysis:
####################################################################################################################
baseDir = '/Users/soonyi/Desktop/Genomics/Specificity/Data/'
extensions = c(0, 15, 40, 65, 90)
LOCI = c("5'UTR", "CDS", "3'UTR", "intron")

for (extension in extensions) {
  print(paste0('Working on ', extension, 'nt extension for coCLIP....'))
  for (locus in LOCI) {
    filtered_peaksMatrix = peaksMatrix %>% filter(grouped_annotation == locus)
    print(paste0('Looking at ', locus, ' specific peaks (', nrow(filtered_peaksMatrix) , ' peaks) for coCLIP....'))
    
    peaksGR = filtered_peaksMatrix[, c('chrom', 'start', 'end', 'strand', 'peak_names')]
    peaksGR = peaksGR %>% mutate(chrom = ifelse(chrom == "chrMT", "chrM", chrom))
    peaksGR$start = as.integer(peaksGR$start) + 1 - extension
    peaksGR$end = as.integer(peaksGR$end) + extension
    peaksGR = GRanges(peaksGR)
    
    peaksGR_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, peaksGR, as.character = TRUE)
    peaksGR_seqs = cbind(filtered_peaksMatrix[, c('chrom', 'start', 'end', 'strand', 'peak_names', BC_columns, 'grouped_annotation', 'finalized_annotation')], data.frame(sequence = peaksGR_seqs))
    
    Samples = c('I_M', 'NLS_M', 'NES_M', 'G3BP_M', 'I_S','NLS_S', 'NES_S', 'G3BP_S')
    I_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_Input_M$peak_names))
    NLS_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NLS_M$peak_names))
    NES_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NES_M$peak_names))
    G3BP_M = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_G3BP_M$peak_names))
    I_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_Input_S$peak_names))
    NLS_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NLS_S$peak_names))
    NES_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_NES_S$peak_names))
    G3BP_S = peaksGR_seqs %>% filter(peak_names %in% c(Peak_Co_G3BP_S$peak_names))
    
    for (K in Ks) {
      print(paste0('Working on ', K, 'mer....'))
      # Make all K-mers
      combinations = expand.grid(replicate(K, NTs, simplify = FALSE))
      kmers = apply(combinations, 1, paste, collapse = "")
      
      raw_count = data.frame(MOTIF = kmers)
      nor_count = data.frame(MOTIF = kmers)
      
      for (sample in Samples) {
        INPUT = get(sample)
        raw_count[[sample]] = rep(0, 4^K)
        nor_count[[sample]] = rep(0, 4^K)
        num_peak = nrow(INPUT)
        for (kmer in kmers) {
          kmer_counts = colSums(motifCounts(INPUT, kmer))
          raw_count[which(raw_count$MOTIF == kmer), sample] = kmer_counts
          nor_count[which(nor_count$MOTIF == kmer), sample] = kmer_counts/num_peak
        }
      }
      if (locus == "5'UTR") {
        locus_fn = "UTR5"
      } else if (locus == "3'UTR") {
        locus_fn = "UTR3"
      } else if (locus == 'intron') {
        locus_fn = "INTRON"
      } else {
        locus_fn = locus
      }
      write.csv(raw_count, paste0(baseDir, 'CLIP_Analysis/', K ,'mer/Locus_Specific/', locus_fn, '_coCLIP_Raw_', K, 'mer_', extension ,'nt_PeakExtension_Counts.csv'), quote = F)
      write.csv(nor_count, paste0(baseDir, 'CLIP_Analysis/', K ,'mer/Locus_Specific/', locus_fn, '_coCLIP_Normalized_', K, 'mer_', extension ,'nt_PeakExtension_Counts.csv'), quote = F)
      normalized_dataframe = nor_count
      normalized_dataframe = normalized_dataframe %>% mutate_if(is.numeric, min_max_norm) %>% mutate_if(is.numeric, round, 10)
      write.csv(normalized_dataframe, paste0(baseDir, 'CLIP_Analysis/', K ,'mer/Locus_Specific/', locus_fn, '_coCLIP_Scaled_', K, 'mer_', extension ,'nt_PeakExtension_Counts_Normalized.csv'), quote = F)
      specificity_dataframe = normalized_dataframe %>% mutate_if(is.numeric, SI_Calculation) %>% mutate_if(is.numeric, round, 10)
      write.csv(specificity_dataframe, paste0(baseDir, 'CLIP_Analysis/', K ,'mer/Locus_Specific/', locus_fn, '_coCLIP_', K, 'mer_', extension ,'nt_PeakExtension_SpecificityIndex.csv'), quote = F)
    } 
  }
}
####################################################################################################################

## Locus specific analysis plotting:
####################################################################################################################
for (locus in LOCI) {
  for (extension in extensions) {
    if (locus == "5'UTR") {
      locus_fn = "UTR5"
    } else if (locus == "3'UTR") {
      locus_fn = "UTR3"
    } else if (locus == 'intron') {
      locus_fn = "INTRON"
    } else {
      locus_fn = locus
    }
    temp = read.delim(paste0(baseDir, 'CLIP_Analysis/5mer/Locus_Specific/', locus_fn ,'_coCLIP_5mer_', extension, 'nt_PeakExtension_SpecificityIndex.csv'), sep = ',')
    temp = temp[, -1]
    df_name = paste0(locus_fn, "_SI_5mer_", (20 + extension*2), "nt")
    assign(df_name, temp, envir = .GlobalEnv)
  }
}

## UTR5
UTR5_SI_selection = which((UTR5_SI_5mer_50nt$I_M >= 1) &
                          (UTR5_SI_5mer_50nt$NLS_M >= 1) &
                          (UTR5_SI_5mer_50nt$NES_M >= 1) &
                          (UTR5_SI_5mer_50nt$G3BP_M >= 1) &
                          (UTR5_SI_5mer_50nt$I_S >= 1) &
                          (UTR5_SI_5mer_50nt$NLS_S >= 1) &
                          (UTR5_SI_5mer_50nt$NES_S >= 1) &
                          (UTR5_SI_5mer_50nt$G3BP_S >= 1))

UTR5_SI_FC = data.frame(MOTIF = UTR5_SI_5mer_50nt$MOTIF,
                        Input = (UTR5_SI_5mer_50nt$I_S + 1e-10)/(UTR5_SI_5mer_50nt$I_M + 1e-10),
                        NLS = (UTR5_SI_5mer_50nt$NLS_S + 1e-10)/(UTR5_SI_5mer_50nt$NLS_M + 1e-10),
                        NES = (UTR5_SI_5mer_50nt$NES_S + 1e-10)/(UTR5_SI_5mer_50nt$NES_M + 1e-10),
                        G3BP = (UTR5_SI_5mer_50nt$G3BP_S + 1e-10)/(UTR5_SI_5mer_50nt$G3BP_M + 1e-10))

UTR5_SI_FC = UTR5_SI_FC[UTR5_SI_selection, ]
UTR5_SI_L2FC = UTR5_SI_FC %>% mutate_if(is.numeric, log2)
UTR5_SI_5mer_50nt_Selection = UTR5_SI_5mer_50nt[UTR5_SI_selection, ]
UTR5_SI_5mer_50nt_Selection_L10FC = UTR5_SI_5mer_50nt_Selection %>% mutate_if(is.numeric, log10)

## CDS
CDS_SI_selection = which((CDS_SI_5mer_50nt$I_M >= 1) &
                            (CDS_SI_5mer_50nt$NLS_M >= 1) &
                            (CDS_SI_5mer_50nt$NES_M >= 1) &
                            (CDS_SI_5mer_50nt$G3BP_M >= 1) &
                            (CDS_SI_5mer_50nt$I_S >= 1) &
                            (CDS_SI_5mer_50nt$NLS_S >= 1) &
                            (CDS_SI_5mer_50nt$NES_S >= 1) &
                            (CDS_SI_5mer_50nt$G3BP_S >= 1))

CDS_SI_FC = data.frame(MOTIF = CDS_SI_5mer_50nt$MOTIF,
                        Input = (CDS_SI_5mer_50nt$I_S + 1e-10)/(CDS_SI_5mer_50nt$I_M + 1e-10),
                        NLS = (CDS_SI_5mer_50nt$NLS_S + 1e-10)/(CDS_SI_5mer_50nt$NLS_M + 1e-10),
                        NES = (CDS_SI_5mer_50nt$NES_S + 1e-10)/(CDS_SI_5mer_50nt$NES_M + 1e-10),
                        G3BP = (CDS_SI_5mer_50nt$G3BP_S + 1e-10)/(CDS_SI_5mer_50nt$G3BP_M + 1e-10))

CDS_SI_FC = CDS_SI_FC[CDS_SI_selection, ]
CDS_SI_L2FC = CDS_SI_FC %>% mutate_if(is.numeric, log2)
CDS_SI_5mer_50nt_Selection = CDS_SI_5mer_50nt[CDS_SI_selection, ]
CDS_SI_5mer_50nt_Selection_L10FC = CDS_SI_5mer_50nt_Selection %>% mutate_if(is.numeric, log10)

## UTR3
UTR3_SI_selection = which((UTR3_SI_5mer_50nt$I_M >= 1) &
                            (UTR3_SI_5mer_50nt$NLS_M >= 1) &
                            (UTR3_SI_5mer_50nt$NES_M >= 1) &
                            (UTR3_SI_5mer_50nt$G3BP_M >= 1) &
                            (UTR3_SI_5mer_50nt$I_S >= 1) &
                            (UTR3_SI_5mer_50nt$NLS_S >= 1) &
                            (UTR3_SI_5mer_50nt$NES_S >= 1) &
                            (UTR3_SI_5mer_50nt$G3BP_S >= 1))

UTR3_SI_FC = data.frame(MOTIF = UTR3_SI_5mer_50nt$MOTIF,
                        Input = (UTR3_SI_5mer_50nt$I_S + 1e-10)/(UTR3_SI_5mer_50nt$I_M + 1e-10),
                        NLS = (UTR3_SI_5mer_50nt$NLS_S + 1e-10)/(UTR3_SI_5mer_50nt$NLS_M + 1e-10),
                        NES = (UTR3_SI_5mer_50nt$NES_S + 1e-10)/(UTR3_SI_5mer_50nt$NES_M + 1e-10),
                        G3BP = (UTR3_SI_5mer_50nt$G3BP_S + 1e-10)/(UTR3_SI_5mer_50nt$G3BP_M + 1e-10))

UTR3_SI_FC = UTR3_SI_FC[UTR3_SI_selection, ]
UTR3_SI_L2FC = UTR3_SI_FC %>% mutate_if(is.numeric, log2)
UTR3_SI_5mer_50nt_Selection = UTR3_SI_5mer_50nt[UTR3_SI_selection, ]
UTR3_SI_5mer_50nt_Selection_L10FC = UTR3_SI_5mer_50nt_Selection %>% mutate_if(is.numeric, log10)

## Intron
INTRON_SI_selection = which((INTRON_SI_5mer_50nt$I_M >= 1) &
                            (INTRON_SI_5mer_50nt$NLS_M >= 1) &
                            (INTRON_SI_5mer_50nt$NES_M >= 1) &
                            (INTRON_SI_5mer_50nt$G3BP_M >= 1) &
                            (INTRON_SI_5mer_50nt$I_S >= 1) &
                            (INTRON_SI_5mer_50nt$NLS_S >= 1) &
                            (INTRON_SI_5mer_50nt$NES_S >= 1) &
                            (INTRON_SI_5mer_50nt$G3BP_S >= 1))

INTRON_SI_FC = data.frame(MOTIF = INTRON_SI_5mer_50nt$MOTIF,
                        Input = (INTRON_SI_5mer_50nt$I_S + 1e-10)/(INTRON_SI_5mer_50nt$I_M + 1e-10),
                        NLS = (INTRON_SI_5mer_50nt$NLS_S + 1e-10)/(INTRON_SI_5mer_50nt$NLS_M + 1e-10),
                        NES = (INTRON_SI_5mer_50nt$NES_S + 1e-10)/(INTRON_SI_5mer_50nt$NES_M + 1e-10),
                        G3BP = (INTRON_SI_5mer_50nt$G3BP_S + 1e-10)/(INTRON_SI_5mer_50nt$G3BP_M + 1e-10))

INTRON_SI_FC = INTRON_SI_FC[INTRON_SI_selection, ]
INTRON_SI_L2FC = INTRON_SI_FC %>% mutate_if(is.numeric, log2)
INTRON_SI_5mer_50nt_Selection = INTRON_SI_5mer_50nt[INTRON_SI_selection, ]
INTRON_SI_5mer_50nt_Selection_L10FC = INTRON_SI_5mer_50nt_Selection %>% mutate_if(is.numeric, log10)

##
SI_FC_long = SI_L2FC %>%
  select(MOTIF, Input, NLS, NES, G3BP) %>%
  pivot_longer(cols = c("Input", "NLS", "NES", "G3BP"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("Input", "NLS", "NES", "G3BP")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index FC", title = "Total") +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


##
UTR5_SI_FC_long = UTR5_SI_L2FC %>%
  select(MOTIF, Input, NLS, NES, G3BP) %>%
  pivot_longer(cols = c("Input", "NLS", "NES", "G3BP"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("Input", "NLS", "NES", "G3BP")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(UTR5_SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index FC", title = "5'UTR") +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


##
CDS_SI_FC_long = CDS_SI_L2FC %>%
  select(MOTIF, Input, NLS, NES, G3BP) %>%
  pivot_longer(cols = c("Input", "NLS", "NES", "G3BP"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("Input", "NLS", "NES", "G3BP")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(CDS_SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index FC", title = "CDS") +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


##
UTR3_SI_FC_long = UTR3_SI_L2FC %>%
  select(MOTIF, Input, NLS, NES, G3BP) %>%
  pivot_longer(cols = c("Input", "NLS", "NES", "G3BP"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("Input", "NLS", "NES", "G3BP")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(UTR3_SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index FC", title = "3'UTR") +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


##
INTRON_SI_FC_long = INTRON_SI_L2FC %>%
  select(MOTIF, Input, NLS, NES, G3BP) %>%
  pivot_longer(cols = c("Input", "NLS", "NES", "G3BP"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("Input", "NLS", "NES", "G3BP")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(INTRON_SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index FC", title = "Intron") +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

####################################################################################################################


####################################################################################################################
##
SI_FC_long = SI_5mer_50nt_Selection_L10FC %>%
  select(MOTIF, I_M, NLS_M, NES_M, G3BP_M) %>%
  pivot_longer(cols = c("I_M", "NLS_M", "NES_M", "G3BP_M"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("I_M", "NLS_M", "NES_M", "G3BP_M")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index Log10", title = "Total") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


##
UTR5_SI_FC_long = UTR5_SI_5mer_50nt_Selection_L10FC %>%
  select(MOTIF, I_M, NLS_M, NES_M, G3BP_M) %>%
  pivot_longer(cols = c("I_M", "NLS_M", "NES_M", "G3BP_M"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("I_M", "NLS_M", "NES_M", "G3BP_M")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(UTR5_SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index FC", title = "5'UTR") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


##
CDS_SI_FC_long = CDS_SI_5mer_50nt_Selection_L10FC %>%
  select(MOTIF, I_M, NLS_M, NES_M, G3BP_M) %>%
  pivot_longer(cols = c("I_M", "NLS_M", "NES_M", "G3BP_M"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("I_M", "NLS_M", "NES_M", "G3BP_M")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(CDS_SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index Log10", title = "CDS") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


##
UTR3_SI_FC_long = UTR3_SI_5mer_50nt_Selection_L10FC %>%
  select(MOTIF, I_M, NLS_M, NES_M, G3BP_M) %>%
  pivot_longer(cols = c("I_M", "NLS_M", "NES_M", "G3BP_M"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("I_M", "NLS_M", "NES_M", "G3BP_M")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(UTR3_SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index Log10", title = "3'UTR") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


##
INTRON_SI_FC_long = INTRON_SI_5mer_50nt_Selection_L10FC %>%
  select(MOTIF, I_M, NLS_M, NES_M, G3BP_M) %>%
  pivot_longer(cols = c("I_M", "NLS_M", "NES_M", "G3BP_M"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("I_M", "NLS_M", "NES_M", "G3BP_M")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(INTRON_SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index Log10", title = "Intron") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
####################################################################################################################

####################################################################################################################
##
SI_FC_long = SI_5mer_50nt_Selection_L10FC %>%
  select(MOTIF, I_S, NLS_S, NES_S, G3BP_S) %>%
  pivot_longer(cols = c("I_S", "NLS_S", "NES_S", "G3BP_S"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("I_S", "NLS_S", "NES_S", "G3BP_S")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index Log10", title = "Total") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


##
UTR5_SI_FC_long = UTR5_SI_5mer_50nt_Selection_L10FC %>%
  select(MOTIF, I_S, NLS_S, NES_S, G3BP_S) %>%
  pivot_longer(cols = c("I_S", "NLS_S", "NES_S", "G3BP_S"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("I_S", "NLS_S", "NES_S", "G3BP_S")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(UTR5_SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index FC", title = "5'UTR") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


##
CDS_SI_FC_long = CDS_SI_5mer_50nt_Selection_L10FC %>%
  select(MOTIF, I_S, NLS_S, NES_S, G3BP_S) %>%
  pivot_longer(cols = c("I_S", "NLS_S", "NES_S", "G3BP_S"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("I_S", "NLS_S", "NES_S", "G3BP_S")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(CDS_SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index Log10", title = "CDS") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


##
UTR3_SI_FC_long = UTR3_SI_5mer_50nt_Selection_L10FC %>%
  select(MOTIF, I_S, NLS_S, NES_S, G3BP_S) %>%
  pivot_longer(cols = c("I_S", "NLS_S", "NES_S", "G3BP_S"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("I_S", "NLS_S", "NES_S", "G3BP_S")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(UTR3_SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index Log10", title = "3'UTR") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


##
INTRON_SI_FC_long = INTRON_SI_5mer_50nt_Selection_L10FC %>%
  select(MOTIF, I_S, NLS_S, NES_S, G3BP_S) %>%
  pivot_longer(cols = c("I_S", "NLS_S", "NES_S", "G3BP_S"), names_to = "Category", values_to = "Value") %>%
  mutate(Category = factor(Category, levels = c("I_S", "NLS_S", "NES_S", "G3BP_S")))

colors = c("black", "skyblue", "darkseagreen2", "salmon")

ggplot(INTRON_SI_FC_long, aes(x = Category, y = Value, color = Category)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 1, shape = 16) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Category", y = "Specificity Index Log10", title = "Intron") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 5, by = 1))  +
  theme(legend.title = element_blank()) +
  theme_minimal() +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
####################################################################################################################





max(normalized_dataframe[, 'I_M'])/median(normalized_dataframe[, 'I_M'])
max(normalized_dataframe[, 'NES_M'])/median(normalized_dataframe[, 'NES_M'])
max(normalized_dataframe[, 'NLS_M'])/median(normalized_dataframe[, 'NLS_M'])
max(normalized_dataframe[, 'G3BP_M'])/median(normalized_dataframe[, 'G3BP_M'])
max(normalized_dataframe[, 'I_S'])/median(normalized_dataframe[, 'I_S'])
max(normalized_dataframe[, 'NES_S'])/median(normalized_dataframe[, 'NES_S'])
max(normalized_dataframe[, 'NLS_S'])/median(normalized_dataframe[, 'NLS_S'])
max(normalized_dataframe[, 'G3BP_S'])/median(normalized_dataframe[, 'G3BP_S'])

normalized_dataframe[which(normalized_dataframe$MOTIF == 'AAAAA'), 'I_M']/median(normalized_dataframe[, 'I_M'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'AAAAA'), 'NES_M']/median(normalized_dataframe[, 'NES_M'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'AAAAA'), 'NLS_M']/median(normalized_dataframe[, 'NLS_M'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'AAAAA'), 'G3BP_M']/median(normalized_dataframe[, 'G3BP_M'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'AAAAA'), 'I_S']/median(normalized_dataframe[, 'I_S'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'AAAAA'), 'NES_S']/median(normalized_dataframe[, 'NES_S'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'AAAAA'), 'NLS_S']/median(normalized_dataframe[, 'NLS_S'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'AAAAA'), 'G3BP_S']/median(normalized_dataframe[, 'G3BP_S'])

normalized_dataframe[which(normalized_dataframe$MOTIF == 'TTTTT'), 'I_M']/median(normalized_dataframe[, 'I_M'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'TTTTT'), 'NES_M']/median(normalized_dataframe[, 'NES_M'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'TTTTT'), 'NLS_M']/median(normalized_dataframe[, 'NLS_M'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'TTTTT'), 'G3BP_M']/median(normalized_dataframe[, 'G3BP_M'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'TTTTT'), 'I_S']/median(normalized_dataframe[, 'I_S'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'TTTTT'), 'NES_S']/median(normalized_dataframe[, 'NES_S'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'TTTTT'), 'NLS_S']/median(normalized_dataframe[, 'NLS_S'])
normalized_dataframe[which(normalized_dataframe$MOTIF == 'TTTTT'), 'G3BP_S']/median(normalized_dataframe[, 'G3BP_S'])


specificity_fc = data.frame(MOTIF = normalized_dataframe$MOTIF)
specificity_fc$Input = (normalized_dataframe$I_S + 1e-12) / (normalized_dataframe$I_M + 1e-12)
specificity_fc$NES = (normalized_dataframe$NES_S + 1e-12) / (normalized_dataframe$NES_M + 1e-12)
specificity_fc$NLS = (normalized_dataframe$NLS_S + 1e-12) / (normalized_dataframe$NLS_M + 1e-12)
specificity_fc$G3BP = (normalized_dataframe$G3BP_S + 1e-12) / (normalized_dataframe$G3BP_M + 1e-12)

fc_u = mean(as.matrix(specificity_fc[, c(2, 3, 4, 5)]))
fc_sd = sd(as.matrix(specificity_fc[, c(2, 3, 4, 5)]))

specificity_fc_z= data.frame(MOTIF = normalized_dataframe$MOTIF)
specificity_fc_z$Input = (specificity_fc$Input - fc_u)/fc_sd
specificity_fc_z$NES = (specificity_fc$NES - fc_u)/fc_sd
specificity_fc_z$NLS = (specificity_fc$NLS - fc_u)/fc_sd
specificity_fc_z$G3BP = (specificity_fc$G3BP - fc_u)/fc_sd

####################################################################################################################

## :
####################################################################################################################

histogram_plot = ggplot(normalized_dataframe, aes(x = I_M)) + geom_histogram(bins = 500)
hist_data = ggplot_build(histogram_plot)$data[[1]]
I_M_hist = as.data.frame(hist_data)[, c('x', 'y')]

histogram_plot = ggplot(normalized_dataframe, aes(x = I_S)) + geom_histogram(bins = 500)
hist_data = ggplot_build(histogram_plot)$data[[1]]
I_S_hist = as.data.frame(hist_data)[, c('x', 'y')]

histogram_plot = ggplot(normalized_dataframe, aes(x = NES_M)) + geom_histogram(bins = 500)
hist_data = ggplot_build(histogram_plot)$data[[1]]
NES_M_hist = as.data.frame(hist_data)[, c('x', 'y')]

histogram_plot = ggplot(normalized_dataframe, aes(x = NES_S)) + geom_histogram(bins = 500)
hist_data = ggplot_build(histogram_plot)$data[[1]]
NES_S_hist = as.data.frame(hist_data)[, c('x', 'y')]

histogram_plot = ggplot(normalized_dataframe, aes(x = NLS_M)) + geom_histogram(bins = 500)
hist_data = ggplot_build(histogram_plot)$data[[1]]
NLS_M_hist = as.data.frame(hist_data)[, c('x', 'y')]

histogram_plot = ggplot(normalized_dataframe, aes(x = NLS_S)) + geom_histogram(bins = 500)
hist_data = ggplot_build(histogram_plot)$data[[1]]
NLS_S_hist = as.data.frame(hist_data)[, c('x', 'y')]

histogram_plot = ggplot(normalized_dataframe, aes(x = G3BP_M)) + geom_histogram(bins = 500)
hist_data = ggplot_build(histogram_plot)$data[[1]]
G3BP_M_hist = as.data.frame(hist_data)[, c('x', 'y')]

histogram_plot = ggplot(normalized_dataframe, aes(x = G3BP_S)) + geom_histogram(bins = 500)
hist_data = ggplot_build(histogram_plot)$data[[1]]
G3BP_S_hist = as.data.frame(hist_data)[, c('x', 'y')]

t.test()

I_comparison = rbind(I_M_hist$y, I_S_hist$y)
chisq.test(I_comparison)

ks.test(normalized_dataframe$I_M, normalized_dataframe$I_S)
ks.test(normalized_dataframe$NES_M, normalized_dataframe$NES_S)
ks.test(normalized_dataframe$NLS_M, normalized_dataframe$NLS_S)
ks.test(normalized_dataframe$G3BP_M, normalized_dataframe$G3BP_S)

chisq.test(I_M_hist$y, I_S_hist$y, )
chisq.test(NES_M_hist$y, NES_S_hist$y)
chisq.test(NLS_M_hist$y, NLS_S_hist$y)
chisq.test(G3BP_M_hist$y, G3BP_S_hist$y)

t.test(normalized_dataframe$I_M, normalized_dataframe$I_S)
t.test(normalized_dataframe$NES_M, normalized_dataframe$NES_S)
t.test(normalized_dataframe$NLS_M, normalized_dataframe$NLS_S)
t.test(normalized_dataframe$G3BP_M, normalized_dataframe$G3BP_S)

wilcox.test(I_M_hist$y, I_S_hist$y, )
wilcox.test(NES_M_hist$y, NES_S_hist$y)
wilcox.test(NLS_M_hist$y, NLS_S_hist$y)
wilcox.test(G3BP_M_hist$y, G3BP_S_hist$y)

####################################################################################################################

####################################################################################################################

ggplot(normalized_dataframe_eCLIP, aes(eCLIP)) +
  geom_histogram(bins = 500) +
  ylim(0, 150)
  

max(normalized_dataframe_eCLIP[, 'eCLIP'])/median(normalized_dataframe_eCLIP[, 'eCLIP'])
normalized_dataframe_eCLIP[which(normalized_dataframe_eCLIP$MOTIF == 'AAAAA'), 'eCLIP']/median(normalized_dataframe_eCLIP[, 'eCLIP'])
normalized_dataframe_eCLIP[which(normalized_dataframe_eCLIP$MOTIF == 'TTTTT'), 'eCLIP']/median(normalized_dataframe_eCLIP[, 'eCLIP'])

####################################################################################################################








