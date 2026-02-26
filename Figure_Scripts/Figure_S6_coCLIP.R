################################################################################
## CoCLIP Analysis
## Written by Soon Yi
## Created: 2026-02-24
## Last Edited: 2026-02-24
## Figure S6
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
library(RBPSpecificity)
library(seqLogo)
library(ggseqlogo)

## 1. Basic Setup:
################################################################################
BASE_DIR = "/Volumes/1TB_Data/Specificity/CoCLIP/"
OUTPUT_DIR = paste0(BASE_DIR, "output/")

PEAK_MATRIX_FILE = paste0(BASE_DIR, "Combined_peakCoverage_groomed_normalized_annotated.txt")
################################################################################

## 2. Custom Functions:
################################################################################

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

# Compute per-position nucleotide frequency around peak centers (vectorized):
nucleotideFrequencyAroundPeaks = function(peaks_df, genome, window = 200, extension = c(0, 0)) {
  
  peaks_df_sense = peaks_df[peaks_df$strand == '+', ]
  peaks_df_antisense = peaks_df[peaks_df$strand == '-', ]
  
  ## Apply strand-aware extension: ext[1] = 5' end, ext[2] = 3' end
  if (nrow(peaks_df_sense) > 0) {
    peaks_df_sense$start = peaks_df_sense$start - extension[1]
    peaks_df_sense$end   = peaks_df_sense$end   + extension[2]
  }
  if (nrow(peaks_df_antisense) > 0) {
    peaks_df_antisense$start = peaks_df_antisense$start - extension[2]
    peaks_df_antisense$end   = peaks_df_antisense$end   + extension[1]
  }
  
  all_seqs = DNAStringSet()
  
  ## Sense strand: extract window around peak center
  if (nrow(peaks_df_sense) > 0) {
    chroms = as.character(peaks_df_sense$chr)
    centers = floor((peaks_df_sense$start + peaks_df_sense$end) / 2)
    
    starts = centers - window
    ends = centers + window - 1
    
    ## Filter valid ranges
    valid_chr = chroms %in% seqnames(genome)
    chr_lens = rep(NA_integer_, length(chroms))
    chr_lens[valid_chr] = seqlengths(genome)[chroms[valid_chr]]
    valid = valid_chr & starts >= 1 & ends <= chr_lens
    
    if (sum(valid) > 0) {
      gr = GRanges(seqnames = chroms[valid],
                   ranges = IRanges(start = starts[valid], end = ends[valid]),
                   strand = '+')
      seqs = getSeq(genome, gr)
      all_seqs = c(all_seqs, seqs)
    }
  }
  
  ## Antisense strand: extract and reverse complement
  if (nrow(peaks_df_antisense) > 0) {
    chroms = as.character(peaks_df_antisense$chr)
    centers = floor((peaks_df_antisense$start + peaks_df_antisense$end) / 2)
    
    starts = centers - window
    ends = centers + window - 1
    
    valid_chr = chroms %in% seqnames(genome)
    chr_lens = rep(NA_integer_, length(chroms))
    chr_lens[valid_chr] = seqlengths(genome)[chroms[valid_chr]]
    valid = valid_chr & starts >= 1 & ends <= chr_lens
    
    if (sum(valid) > 0) {
      gr = GRanges(seqnames = chroms[valid],
                   ranges = IRanges(start = starts[valid], end = ends[valid]),
                   strand = '-')
      seqs = getSeq(genome, gr)
      ## getSeq(strand='-') automatically reverse complements to orient 5'->3' related to RNA
      all_seqs = c(all_seqs, seqs)
    }
  }
  
  if (length(all_seqs) == 0) {
    warning("No valid sequences extracted")
    return(NULL)
  }
  
  ## consensusMatrix: rows = A/C/G/T, columns = positions
  cm = consensusMatrix(all_seqs, as.prob = TRUE)
  
  ## Extract only A, C, G, T rows (ignore N, etc.)
  positions = seq(-window, window - 1)
  
  result = data.frame(
    position = positions,
    A = cm["A", ],
    C = cm["C", ],
    G = cm["G", ],
    U = cm["T", ]   ## T on DNA = U on RNA
  )
  
  return(result)
}

# Plot stacked barplot of nucleotide composition around peaks:
plotNucleotideBar = function(freq_df, sampleName = NULL, window_lim = NULL, nuc_colors = NULL, x_ticks = NULL) {
  
  if (is.null(nuc_colors)) {
    nuc_colors = c("A" = "#4DAF4A", "C" = "#377EB8", "G" = "#FF7F00", "U" = "#E41A1C")
  }
  
  plot_data = freq_df %>%
    pivot_longer(cols = c(A, C, G, U), names_to = "Nucleotide", values_to = "Fraction")
  
  ## Set nucleotide order (bottom to top of stack)
  plot_data$Nucleotide = factor(plot_data$Nucleotide, levels = c("A", "C", "G", "U"))
  
  plot = ggplot(plot_data, aes(x = position, y = Fraction, fill = Nucleotide)) +
    geom_col(position = 'stack', width = 1) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = 0.25, color = "gray40", linetype = "dotted", linewidth = 0.3) +
    scale_fill_manual(values = nuc_colors) +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = 'bold'),
          legend.text = element_text(size = 14),
          panel.grid.minor = element_blank()) +
    labs(x = "Distance from peak center (nt)",
         y = "Nucleotide Fraction")
  
  if (!is.null(sampleName)) {
    plot = plot + ggtitle(paste0("Sequence Context around Peaks: ", sampleName))
  }
  
  if (!is.null(x_ticks) && length(x_ticks) == 2) {
    minor_val = x_ticks[1]
    major_val = x_ticks[2]
    max_pos = max(abs(freq_df$position))
    plot = plot + scale_x_continuous(breaks = seq(-max_pos, max_pos, by = major_val),
                                     minor_breaks = seq(-max_pos, max_pos, by = minor_val))
  }
  
  if (!is.null(window_lim)) {
    plot = plot + xlim(c(-window_lim, window_lim))
  }
  
  return(plot)
}

# Plot sequence logo using ggseqlogo (for narrow window):
plotSeqLogo = function(peaks_df, genome, window = 15, sampleName = NULL, extension = c(0, 0), method = 'bits', x_ticks = NULL) {
  
  peaks_df_sense = peaks_df[peaks_df$strand == '+', ]
  peaks_df_antisense = peaks_df[peaks_df$strand == '-', ]
  
  ## Apply strand-aware extension: ext[1] = 5' end, ext[2] = 3' end
  if (nrow(peaks_df_sense) > 0) {
    peaks_df_sense$start = peaks_df_sense$start - extension[1]
    peaks_df_sense$end   = peaks_df_sense$end   + extension[2]
  }
  if (nrow(peaks_df_antisense) > 0) {
    peaks_df_antisense$start = peaks_df_antisense$start - extension[2]
    peaks_df_antisense$end   = peaks_df_antisense$end   + extension[1]
  }
  
  all_seqs = c()
  
  ## Sense strand:
  if (nrow(peaks_df_sense) > 0) {
    chroms = as.character(peaks_df_sense$chr)
    centers = floor((peaks_df_sense$start + peaks_df_sense$end) / 2)
    starts = centers - window
    ends = centers + window
    
    valid_chr = chroms %in% seqnames(genome)
    chr_lens = rep(NA_integer_, length(chroms))
    chr_lens[valid_chr] = seqlengths(genome)[chroms[valid_chr]]
    valid = valid_chr & starts >= 1 & ends <= chr_lens
    
    if (sum(valid) > 0) {
      gr = GRanges(seqnames = chroms[valid],
                   ranges = IRanges(start = starts[valid], end = ends[valid]),
                   strand = '+')
      seqs = getSeq(genome, gr)
      all_seqs = c(all_seqs, as.character(seqs))
    }
  }
  
  ## Antisense strand:
  if (nrow(peaks_df_antisense) > 0) {
    chroms = as.character(peaks_df_antisense$chr)
    centers = floor((peaks_df_antisense$start + peaks_df_antisense$end) / 2)
    starts = centers - window
    ends = centers + window
    
    valid_chr = chroms %in% seqnames(genome)
    chr_lens = rep(NA_integer_, length(chroms))
    chr_lens[valid_chr] = seqlengths(genome)[chroms[valid_chr]]
    valid = valid_chr & starts >= 1 & ends <= chr_lens
    
    if (sum(valid) > 0) {
      gr = GRanges(seqnames = chroms[valid],
                   ranges = IRanges(start = starts[valid], end = ends[valid]),
                   strand = '-')
      seqs = getSeq(genome, gr)
      ## getSeq(strand='-') automatically reverse complements to orient 5'->3' related to RNA
      all_seqs = c(all_seqs, as.character(seqs))
    }
  }
  
  if (length(all_seqs) == 0) {
    warning("No valid sequences for logo")
    return(NULL)
  }
  
  ## Replace T with U for RNA display
  all_seqs = gsub("T", "U", all_seqs)
  
  title = "Sequence Logo around Peak Center"
  if (!is.null(sampleName)) {
    title = paste0(title, ": ", sampleName)
  }
  
  plot = ggseqlogo(all_seqs, method = method, seq_type = 'rna') +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = 'bold'),
          plot.title = element_text(size = 14, face = 'bold')) +
    ggtitle(title) +
    xlab(paste0("Position relative to peak center (±", window, " nt)"))
  
  if (!is.null(x_ticks) && length(x_ticks) == 2) {
    minor_val = x_ticks[1]
    major_val = x_ticks[2]
    # In ggseqlogo, positions are 1-indexed by default, but we can set scale_x_continuous explicitly
    # To align with -window to +window, ggseqlogo usually plots 1:(2*window).
    # We override the x-axis scale completely to match the actual distance:
    plot = plot + scale_x_continuous(breaks = seq(1, 2*window + 1, by = major_val),
                                     minor_breaks = seq(1, 2*window + 1, by = minor_val),
                                     labels = seq(-window, window, by = major_val))
  }
  
  return(plot)
}

################################################################################

## 3. Load and filter annotated peak matrix:
################################################################################
peaksMatrix = read_delim(PEAK_MATRIX_FILE)
peaksMatrix = peaksMatrix %>% mutate_at('TOTAL_BC', as.numeric)
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'unannotated', 'UnAn', finalized_annotation))
peaksMatrix = peaksMatrix %>% mutate(grouped_annotation = ifelse(grouped_annotation == 'unannotated', 'UnAn', grouped_annotation))
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'downstream 10K', 'DS10K', finalized_annotation))
peaksMatrix = peaksMatrix %>% mutate(grouped_annotation = ifelse(grouped_annotation == 'downstream 10K', 'DS10K', grouped_annotation))
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'ncRNA_Retained_intron', 'nC_RI', finalized_annotation))
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'CDS_Retained_intron', 'CDS_RI', finalized_annotation))

## further consolidate ncRNA:
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'miRNA', 'Other', finalized_annotation))
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'scaRNA', 'Other', finalized_annotation))
peaksMatrix = peaksMatrix %>% mutate(finalized_annotation = ifelse(finalized_annotation == 'nC_RI', 'Other', finalized_annotation))

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
################################################################################

## Filter Criteria:
################################################################################
BC_Threshold_F = 2
BC_Threshold_I = 4
BC_Threshold_E = 2
BC_Threshold_E_SG = 3

rowSum_Multiplier_F = 4
rowSum_Multiplier_I = 2
rowSum_Multiplier_E = 2
################################################################################

## Subset of Peaks for downstream analysis:
################################################################################
## CoCLIP
Peak_Co_Input_M = filterPeakMatrix(peaksMatrix, c(NLS_I_M, NES_I_M, G3BP_I_M), c(inert_columns, BC_columns), BC_Threshold_I, rowSum_Multiplier_I)
Peak_Co_Input_S = filterPeakMatrix(peaksMatrix, c(NLS_I_S, NES_I_S, G3BP_I_S), c(inert_columns, BC_columns), BC_Threshold_I, rowSum_Multiplier_I)
Peak_Co_NLS_M = filterPeakMatrix(peaksMatrix, NLS_E_M, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_NLS_S = filterPeakMatrix(peaksMatrix, NLS_E_S, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_NES_M = filterPeakMatrix(peaksMatrix, NES_E_M, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_NES_S = filterPeakMatrix(peaksMatrix, NES_E_S, c(inert_columns, BC_columns), BC_Threshold_E, rowSum_Multiplier_E)
Peak_Co_G3BP_M = filterPeakMatrix(peaksMatrix, G3BP_E_M, c(inert_columns, BC_columns), BC_Threshold_E_SG, rowSum_Multiplier_E)
Peak_Co_G3BP_S = filterPeakMatrix(peaksMatrix, G3BP_E_S, c(inert_columns, BC_columns), BC_Threshold_E_SG, rowSum_Multiplier_E)

colnames(Peak_Co_Input_M)[c(1, 4)] = c('chr', 'name')
colnames(Peak_Co_Input_S)[c(1, 4)] = c('chr', 'name')
colnames(Peak_Co_NLS_M)[c(1, 4)] = c('chr', 'name')
colnames(Peak_Co_NLS_S)[c(1, 4)] = c('chr', 'name')
colnames(Peak_Co_NES_M)[c(1, 4)] = c('chr', 'name')
colnames(Peak_Co_NES_S)[c(1, 4)] = c('chr', 'name')
colnames(Peak_Co_G3BP_M)[c(1, 4)] = c('chr', 'name')
colnames(Peak_Co_G3BP_S)[c(1, 4)] = c('chr', 'name')


fwrite(Peak_Co_Input_M, paste0(OUTPUT_DIR, "Peak_Co_Input_M", ".txt"), sep = "\t")
fwrite(Peak_Co_Input_S, paste0(OUTPUT_DIR, "Peak_Co_Input_S", ".txt"), sep = "\t")
fwrite(Peak_Co_NLS_M, paste0(OUTPUT_DIR, "Peak_Co_NLS_M", ".txt"), sep = "\t")
fwrite(Peak_Co_NLS_S, paste0(OUTPUT_DIR, "Peak_Co_NLS_S", ".txt"), sep = "\t")
fwrite(Peak_Co_NES_M, paste0(OUTPUT_DIR, "Peak_Co_NES_M", ".txt"), sep = "\t")
fwrite(Peak_Co_NES_S, paste0(OUTPUT_DIR, "Peak_Co_NES_S", ".txt"), sep = "\t")
fwrite(Peak_Co_G3BP_M, paste0(OUTPUT_DIR, "Peak_Co_G3BP_M", ".txt"), sep = "\t")
fwrite(Peak_Co_G3BP_S, paste0(OUTPUT_DIR, "Peak_Co_G3BP_S", ".txt"), sep = "\t")

################################################################################

## 4. Load peak matrix:
################################################################################
Peak_Co_Input_M = fread(paste0(OUTPUT_DIR, "Peak_Co_Input_M.txt"))
Peak_Co_Input_S = fread(paste0(OUTPUT_DIR, "Peak_Co_Input_S.txt"))
Peak_Co_NLS_M = fread(paste0(OUTPUT_DIR, "Peak_Co_NLS_M.txt"))
Peak_Co_NLS_S = fread(paste0(OUTPUT_DIR, "Peak_Co_NLS_S.txt"))
Peak_Co_NES_M = fread(paste0(OUTPUT_DIR, "Peak_Co_NES_M.txt"))
Peak_Co_NES_S = fread(paste0(OUTPUT_DIR, "Peak_Co_NES_S.txt"))
Peak_Co_G3BP_M = fread(paste0(OUTPUT_DIR, "Peak_Co_G3BP_M.txt"))
Peak_Co_G3BP_S = fread(paste0(OUTPUT_DIR, "Peak_Co_G3BP_S.txt"))
################################################################################

################################################################################
SAMPLES = c("Peak_Co_Input_M", "Peak_Co_Input_S",
            "Peak_Co_NLS_M", "Peak_Co_NLS_S",
            "Peak_Co_NES_M", "Peak_Co_NES_S",
            "Peak_Co_G3BP_M", "Peak_Co_G3BP_S")

K_MER = 5
METHOD = 'subtract'

# Output PNG dimensions (easily adjustable)
AFFINITY_WIDTH  = 16
AFFINITY_HEIGHT = 10
AFFINITY_DPI    = 300
PWM_WIDTH       = 1200
PWM_HEIGHT      = 1600
PWM_RES         = 150
PWM_SD          = 2

# Define test configurations
LOCAL_EXTENSIONS = list(
  # "ext0_0" = c(0, 0),
  # "ext0_n5" = c(0, -5),
  # "ext0_n10" = c(0, -10),
  # "ext0_n15" = c(0, -15),
  # "ext5_n5" = c(5, -5),
  # "ext5_n10" = c(5, -10),
  # "ext5_n15" = c(5, -15),
  "ext10_n5" = c(10, -5)
  # "ext10_n10" = c(10, -10),
  # "ext10_n15" = c(10, -15),
  # "ext15_n5" = c(15, -5),
  # "ext15_n10" = c(15, -10),
  # "ext15_n15" = c(15, -15),
  # "ext25_n5" = c(25, -5),
  # "ext25_n10" = c(25, -10),
  # "ext25_n15" = c(25, -15)
)
################################################################################

## 5. Motif Enrichment:
################################################################################

SCRAMBLE_OPTIONS = c(FALSE)

for (ext_name in names(LOCAL_EXTENSIONS)) {
  ext = LOCAL_EXTENSIONS[[ext_name]]
  
  for (scramble in SCRAMBLE_OPTIONS) {
    scramble_suffix = ifelse(scramble, "scrambleON", "scrambleOFF")
    CONDITION_DIR = paste0(OUTPUT_DIR, "local_", ext_name, "_", scramble_suffix, "/")
    if (!dir.exists(CONDITION_DIR)) dir.create(CONDITION_DIR)
    
    message(paste0("Running local background: ", ext_name, " (scramble=", scramble, ")"))
    
    for (SAMPLE in SAMPLES) {
      FILE_PATH = paste0(OUTPUT_DIR, SAMPLE, ".txt")
      peaks = fread(FILE_PATH)
      peak_data = peaks[, .(chr, start, end, name, score, strand)]
      
      motif_enrich = motifEnrichment(peak_data, 'hg38',
                                     K = K_MER,
                                     extension = ext,
                                     nucleic_acid_type = 'RNA',
                                     enrichment_method = METHOD,
                                     background_type = "local",
                                     scramble_bkg = scramble)
      
      fwrite(motif_enrich, paste0(CONDITION_DIR, "motifEnrichment_", SAMPLE, ".txt"), sep = "\t")
    }
  }
}

################################################################################

## 6. Affinity Distribution & PWM (per condition):
################################################################################
# Build all condition directory paths
local_dirs = c()
for (ext_name in names(LOCAL_EXTENSIONS)) {
  local_dirs = c(local_dirs,
                 paste0(OUTPUT_DIR, "local_", ext_name, "_scrambleON/"),
                 paste0(OUTPUT_DIR, "local_", ext_name, "_scrambleOFF/"))
}
all_condition_dirs = c(local_dirs)

returnPWM = function(data, sd_multipler, plot = TRUE) {
  data = data.frame(data[, c('MOTIF', 'Score')])
  data = data[order(-data[, 'Score']), ]
  rownames(data) = NULL
  
  data = data %>% filter(Score > (mean(data$Score) + sd_multipler*sd(data$Score)))
  
  if (nrow(data) == 0) {
    print("No motifs passed the filter.")
    return(NULL)
  }
  
  motif = DNAStringSet(gsub("U", "T", data$MOTIF))
  pfm_raw = consensusMatrix(motif)
  
  bases = c("A", "C", "G", "T")
  pfm = matrix(0, nrow = 4, ncol = ncol(pfm_raw))
  rownames(pfm) = bases
  colnames(pfm) = colnames(pfm_raw)
  
  common_rows = intersect(bases, rownames(pfm_raw))
  if (length(common_rows) > 0) {
    pfm[common_rows, ] = pfm_raw[common_rows, ]
  }
  
  ppm = prop.table(pfm, margin = 2)
  PWM = makePWM(ppm, alphabet = 'RNA')
  if (plot) seqLogo(PWM, ic.scale = F)
  
  ppm_log_ppm = ppm * log2(ppm)
  ppm_log_ppm[is.nan(ppm_log_ppm)] = 0
  H_l = -colSums(ppm_log_ppm)
  IC_l = 2 - H_l
  
  invisible(list(H_motif = sum(H_l), IC_motif = sum(IC_l), ppm = ppm, PWM = PWM, TopMotif = data$MOTIF[1]))
}

capturePWMGrob = function(data, sd_multiplier, title) {
  result = returnPWM(data, sd_multiplier, plot = FALSE)
  if (is.null(result)) {
    grob = grid::textGrob(paste0(title, "\n(No motifs passed filter)"),
                          gp = grid::gpar(fontsize = 11))
  } else {
    logo_grob = grid::grid.grabExpr(seqLogo(result$PWM, ic.scale = FALSE), warn = 0)
    title_grob = grid::textGrob(title, gp = grid::gpar(fontsize = 10, fontface = "bold"))
    grob = gridExtra::arrangeGrob(title_grob, logo_grob, heights = c(0.08, 0.92))
  }
  return(grob)
}

for (CONDITION_DIR in all_condition_dirs) {
  if (!dir.exists(CONDITION_DIR)) next
  
  condition_name = basename(dirname(paste0(CONDITION_DIR, "/")))
  
  motifEnrichment_Co_Input_M = fread(paste0(CONDITION_DIR, "motifEnrichment_Peak_Co_Input_M.txt"))
  motifEnrichment_Co_Input_S = fread(paste0(CONDITION_DIR, "motifEnrichment_Peak_Co_Input_S.txt"))
  motifEnrichment_Co_NLS_M = fread(paste0(CONDITION_DIR, "motifEnrichment_Peak_Co_NLS_M.txt"))
  motifEnrichment_Co_NLS_S = fread(paste0(CONDITION_DIR, "motifEnrichment_Peak_Co_NLS_S.txt"))
  motifEnrichment_Co_NES_M = fread(paste0(CONDITION_DIR, "motifEnrichment_Peak_Co_NES_M.txt"))
  motifEnrichment_Co_NES_S = fread(paste0(CONDITION_DIR, "motifEnrichment_Peak_Co_NES_S.txt"))
  motifEnrichment_Co_G3BP_M = fread(paste0(CONDITION_DIR, "motifEnrichment_Peak_Co_G3BP_M.txt"))
  motifEnrichment_Co_G3BP_S = fread(paste0(CONDITION_DIR, "motifEnrichment_Peak_Co_G3BP_S.txt"))
  
  affinityDistribution_hnRNPC = ggplot() +
    geom_histogram(data = motifEnrichment_Co_Input_M, aes(x = Score), binwidth = 0.005, fill = "blue", alpha = 0.5) +
    geom_histogram(data = motifEnrichment_Co_Input_S, aes(x = Score), binwidth = 0.005, fill = "red", alpha = 0.5) +
    labs(title = paste0("HuR Input Mock vs Stress (", condition_name, ")"), x = "Motif Enrichment", y = "Frequency") +
    scale_y_continuous(limits = c(0, 750), breaks = seq(0, 750, by = 50)) +
    theme_bw() +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14, face = 'bold'))
  
  affinityDistribution_hnRNPC_inRBM25 = ggplot() +
    geom_histogram(data = motifEnrichment_Co_NLS_M, aes(x = Score), binwidth = 0.005, fill = "blue", alpha = 0.5) +
    geom_histogram(data = motifEnrichment_Co_NLS_S, aes(x = Score), binwidth = 0.005, fill = "red", alpha = 0.5) +
    labs(title = paste0("HuR NLS Mock vs Stress (", condition_name, ")"), x = "Motif Enrichment", y = "Frequency") +
    scale_y_continuous(limits = c(0, 750), breaks = seq(0, 750, by = 50)) +
    theme_bw() +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14, face = 'bold'))
  
  affinityDistribution_hnRNPC_inKD = ggplot() +
    geom_histogram(data = motifEnrichment_Co_NES_M, aes(x = Score), binwidth = 0.005, fill = "blue", alpha = 0.5) +
    geom_histogram(data = motifEnrichment_Co_NES_S, aes(x = Score), binwidth = 0.005, fill = "red", alpha = 0.5) +
    labs(title = paste0("HuR NES Mock vs Stress (", condition_name, ")"), x = "Motif Enrichment", y = "Frequency") +
    scale_y_continuous(limits = c(0, 750), breaks = seq(0, 750, by = 50)) +
    theme_bw() +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14, face = 'bold'))
  
  affinityDistribution_RBM25 = ggplot() +
    geom_histogram(data = motifEnrichment_Co_G3BP_M, aes(x = Score), binwidth = 0.005, fill = "blue", alpha = 0.5) +
    geom_histogram(data = motifEnrichment_Co_G3BP_S, aes(x = Score), binwidth = 0.005, fill = "red", alpha = 0.5) +
    labs(title = paste0("HuR G3BP Mock vs Stress (", condition_name, ")"), x = "Motif Enrichment", y = "Frequency") +
    scale_y_continuous(limits = c(0, 750), breaks = seq(0, 750, by = 50)) +
    theme_bw() +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14, face = 'bold'))
  
  # Save 4-panel affinity distribution (PNG + PDF)
  combined_affinity = cowplot::plot_grid(
    affinityDistribution_hnRNPC, affinityDistribution_RBM25,
    affinityDistribution_hnRNPC_inRBM25, affinityDistribution_hnRNPC_inKD,
    nrow = 2, ncol = 2
  )
  ggsave(paste0(CONDITION_DIR, "affinity_distributions.png"),
         combined_affinity, width = AFFINITY_WIDTH, height = AFFINITY_HEIGHT, dpi = AFFINITY_DPI)
  ggsave(paste0(CONDITION_DIR, "affinity_distributions.pdf"),
         combined_affinity, width = AFFINITY_WIDTH, height = AFFINITY_HEIGHT)
  message(paste0("Saved: ", CONDITION_DIR, "affinity_distributions.png/.pdf"))
  
  # Save 8-panel PWM logos (PNG + PDF)
  pwm_grobs = list(
    capturePWMGrob(motifEnrichment_Co_Input_M, PWM_SD, "HuR Input Mock"),
    capturePWMGrob(motifEnrichment_Co_Input_S, PWM_SD, "HuR Input Stress"),
    capturePWMGrob(motifEnrichment_Co_NLS_M, PWM_SD, "HuR NLS Mock"),
    capturePWMGrob(motifEnrichment_Co_NLS_S, PWM_SD, "HuR NLS Stress"),
    capturePWMGrob(motifEnrichment_Co_NES_M, PWM_SD, "HuR NES Mock"),
    capturePWMGrob(motifEnrichment_Co_NES_S, PWM_SD, "HuR NES Stress"),
    capturePWMGrob(motifEnrichment_Co_G3BP_M, PWM_SD, "HuR G3BP Mock"),
    capturePWMGrob(motifEnrichment_Co_G3BP_S, PWM_SD, "HuR G3BP Stress")
  )
  png(paste0(CONDITION_DIR, "pwm_logos.png"), width = PWM_WIDTH, height = PWM_HEIGHT, res = PWM_RES)
  gridExtra::grid.arrange(grobs = pwm_grobs, nrow = 4, ncol = 2,
                          top = grid::textGrob(condition_name, gp = grid::gpar(fontsize = 14, fontface = "bold")))
  dev.off()
  pdf(paste0(CONDITION_DIR, "pwm_logos.pdf"), width = PWM_WIDTH/PWM_RES, height = PWM_HEIGHT/PWM_RES)
  gridExtra::grid.arrange(grobs = pwm_grobs, nrow = 4, ncol = 2,
                          top = grid::textGrob(condition_name, gp = grid::gpar(fontsize = 14, fontface = "bold")))
  dev.off()
  message(paste0("Saved: ", CONDITION_DIR, "pwm_logos.png/.pdf"))
}
################################################################################

## 7. IS/MS Summary Tables:
################################################################################
for (CONDITION_DIR in all_condition_dirs) {
  if (!dir.exists(CONDITION_DIR)) next
  
  summary_data = data.frame(
    Sample = character(),
    TopMotif = character(),
    IS = numeric(),
    MS = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (SAMPLE in SAMPLES) {
    enrich_file = paste0(CONDITION_DIR, "motifEnrichment_", SAMPLE, ".txt")
    if (!file.exists(enrich_file)) next
    
    enrich_data = fread(enrich_file)
    
    # Get Top Motif
    top_motif = enrich_data[which.max(Score), MOTIF]
    
    # Get IS and MS
    is_val = returnIS(enrich_data)
    ms_val = returnMS(enrich_data, output_type = "number")
    
    summary_data = rbind(summary_data, data.frame(
      Sample = SAMPLE,
      TopMotif = top_motif,
      IS = is_val,
      MS = ms_val
    ))
  }
  
  fwrite(summary_data, paste0(CONDITION_DIR, "summary_IS_MS.txt"), sep = "\t")
  message(paste0("Saved summary: ", CONDITION_DIR, "summary_IS_MS.txt"))
}
################################################################################

## 8. Sequence Logos (Narrow View):
################################################################################
genome = BSgenome.Hsapiens.UCSC.hg38
narrow_window = 25
extension = c(0, 0)
x_ticks = c(5, 10)
logo_method = 'prob'
# logo_method = 'bits'
print(plotSeqLogo(Peak_Co_Input_M, genome, narrow_window, "HuR Input Mock", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(Peak_Co_Input_S, genome, narrow_window, "HuR Input Stress", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(Peak_Co_NLS_M, genome, narrow_window, "HuR NLS Mock", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(Peak_Co_NLS_S, genome, narrow_window, "HuR NLS Stress", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(Peak_Co_NES_M, genome, narrow_window, "HuR NES Mock", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(Peak_Co_NES_S, genome, narrow_window, "HuR NES Stress", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(Peak_Co_G3BP_M, genome, narrow_window, "HuR G3BP Stress", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(Peak_Co_G3BP_S, genome, narrow_window, "HuR G3BP Stress", extension = extension, method = logo_method, x_ticks = x_ticks))
################################################################################


################################################################################
t.test(motifEnrichment_Co_Input_M$Score, motifEnrichment_Co_Input_S$Score)
t.test(motifEnrichment_Co_NLS_M$Score, motifEnrichment_Co_NLS_S$Score)
t.test(motifEnrichment_Co_NES_M$Score, motifEnrichment_Co_NES_S$Score)
t.test(motifEnrichment_Co_G3BP_M$Score, motifEnrichment_Co_G3BP_S$Score)
################################################################################
