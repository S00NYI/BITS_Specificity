################################################################################
## CLIP Peak Sequence Context Analysis
## Author: Soon Yi with Antigravity
## Date: February 2026
################################################################################

# === Load Libraries ===
library(dplyr)
library(tidyr)
library(data.table)
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(ggseqlogo)

## 1. Basic Setup:
################################################################################
# BASE_DIR = "F:/Specificity/CLIP/Analysis/"
# INPUT_DIR = paste0(BASE_DIR, "output/")
# OUTPUT_DIR = paste0(BASE_DIR, "output/sequence_context/")

BASE_DIR = "/Volumes/1TB_Data/Specificity/CLIP/Analysis/"
INPUT_DIR = paste0(BASE_DIR, "output/")
OUTPUT_DIR = paste0(BASE_DIR, "output/sequence_context/")

genome = BSgenome.Hsapiens.UCSC.hg38

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

CLIP_List = c('hnRNPC_WT', 'hnRNPC_Mut', 'hnRNPC_WT_inRBM25_WT', 
              'hnRNPC_WT_inRBM25_Mut', 'hnRNPC_WT_inKD', 'hnRNPC_Mut_inKD',
              'RBM25_WT', 'RBM25_Mut')
################################################################################

## 2. Custom Functions:
################################################################################

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
plotNucleotideBar = function(freq_df, sampleName = NULL, window_lim = NULL, 
                             nuc_colors = NULL, x_ticks = NULL) {
  
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

## 3. Load Pre-filtered Peaks:
################################################################################
peaks_hnRNPC_WT = fread(paste0(INPUT_DIR, "peaks_hnRNPC_WT.txt"))
peaks_hnRNPC_Mut = fread(paste0(INPUT_DIR, "peaks_hnRNPC_Mut.txt"))
peaks_hnRNPC_WT_inRBM25_WT = fread(paste0(INPUT_DIR, "peaks_hnRNPC_WT_inRBM25_WT.txt"))
peaks_hnRNPC_WT_inRBM25_Mut = fread(paste0(INPUT_DIR, "peaks_hnRNPC_WT_inRBM25_Mut.txt"))
peaks_hnRNPC_WT_inKD = fread(paste0(INPUT_DIR, "peaks_hnRNPC_WT_inKD.txt"))
peaks_hnRNPC_Mut_inKD = fread(paste0(INPUT_DIR, "peaks_hnRNPC_Mut_inKD.txt"))
peaks_RBM25_WT = fread(paste0(INPUT_DIR, "peaks_RBM25_WT.txt"))
peaks_RBM25_Mut = fread(paste0(INPUT_DIR, "peaks_RBM25_Mut.txt"))
################################################################################

## 4. Compute Nucleotide Frequencies (Broad View):
################################################################################
broad_window = 500

freq_hnRNPC_WT = nucleotideFrequencyAroundPeaks(peaks_hnRNPC_WT, genome, broad_window)
freq_hnRNPC_Mut = nucleotideFrequencyAroundPeaks(peaks_hnRNPC_Mut, genome, broad_window)
freq_hnRNPC_WT_inRBM25_WT = nucleotideFrequencyAroundPeaks(peaks_hnRNPC_WT_inRBM25_WT, genome, broad_window)
freq_hnRNPC_WT_inRBM25_Mut = nucleotideFrequencyAroundPeaks(peaks_hnRNPC_WT_inRBM25_Mut, genome, broad_window)
freq_hnRNPC_WT_inKD = nucleotideFrequencyAroundPeaks(peaks_hnRNPC_WT_inKD, genome, broad_window)
freq_hnRNPC_Mut_inKD = nucleotideFrequencyAroundPeaks(peaks_hnRNPC_Mut_inKD, genome, broad_window)
freq_RBM25_WT = nucleotideFrequencyAroundPeaks(peaks_RBM25_WT, genome, broad_window)
freq_RBM25_Mut = nucleotideFrequencyAroundPeaks(peaks_RBM25_Mut, genome, broad_window)
################################################################################

## 5. Stacked Barplots (Broad View):
################################################################################
print(plotNucleotideBar(freq_hnRNPC_WT, "hnRNPC WT"))
print(plotNucleotideBar(freq_hnRNPC_Mut, "hnRNPC Mut"))
print(plotNucleotideBar(freq_hnRNPC_WT_inRBM25_WT, "hnRNPC WT in RBM25 WT"))
print(plotNucleotideBar(freq_hnRNPC_WT_inRBM25_Mut, "hnRNPC WT in RBM25 Mut"))
print(plotNucleotideBar(freq_hnRNPC_WT_inKD, "hnRNPC WT in KD"))
print(plotNucleotideBar(freq_hnRNPC_Mut_inKD, "hnRNPC Mut in KD"))
print(plotNucleotideBar(freq_RBM25_WT, "RBM25 WT"))
print(plotNucleotideBar(freq_RBM25_Mut, "RBM25 Mut"))
################################################################################

## 6. Sequence Logos (Narrow View):
################################################################################
narrow_window = 50
extension = c(0, 0)
x_ticks = c(5, 10)
# logo_method = 'prob'
logo_method = 'bits'
print(plotSeqLogo(peaks_hnRNPC_WT, genome, narrow_window, "hnRNPC WT", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(peaks_hnRNPC_Mut, genome, narrow_window, "hnRNPC Mut", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(peaks_hnRNPC_WT_inRBM25_WT, genome, narrow_window, "hnRNPC WT in RBM25 WT", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(peaks_hnRNPC_WT_inRBM25_Mut, genome, narrow_window, "hnRNPC WT in RBM25 Mut", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(peaks_hnRNPC_WT_inKD, genome, narrow_window, "hnRNPC WT in KD", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(peaks_hnRNPC_Mut_inKD, genome, narrow_window, "hnRNPC Mut in KD", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(peaks_RBM25_WT, genome, narrow_window, "RBM25 WT", extension = extension, method = logo_method, x_ticks = x_ticks))
print(plotSeqLogo(peaks_RBM25_Mut, genome, narrow_window, "RBM25 Mut", extension = extension, method = logo_method, x_ticks = x_ticks))
################################################################################

## 7. Difference Plots (WT vs Mut):
################################################################################
## Compute per-nucleotide difference: Mut - WT
## Positive = enriched in Mut; Negative = enriched in WT

plotNucleotideDiff = function(freq_a, freq_b, label_a = "A", label_b = "B", nuc_colors = NULL) {
  
  if (is.null(nuc_colors)) {
    nuc_colors = c("A" = "#4DAF4A", "C" = "#377EB8", "G" = "#FF7F00", "U" = "#E41A1C")
  }
  
  diff_df = data.frame(
    position = freq_a$position,
    A = freq_b$A - freq_a$A,
    C = freq_b$C - freq_a$C,
    G = freq_b$G - freq_a$G,
    U = freq_b$U - freq_a$U
  )
  
  diff_long = diff_df %>%
    pivot_longer(cols = c(A, C, G, U), names_to = "Nucleotide", values_to = "Difference")
  diff_long$Nucleotide = factor(diff_long$Nucleotide, levels = c("A", "C", "G", "U"))
  
  plot = ggplot(diff_long, aes(x = position, y = Difference, color = Nucleotide)) +
    geom_hline(yintercept = 0, color = "gray40", linetype = "solid") +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_smooth(span = 0.3, se = FALSE, linewidth = 1) +
    scale_color_manual(values = nuc_colors) +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = 'bold'),
          legend.text = element_text(size = 14)) +
    labs(title = paste0("Nucleotide Difference: ", label_b, " - ", label_a),
         x = "Distance from peak center (nt)",
         y = "Fraction Difference")
  
  return(plot)
}

## hnRNPC WT vs Mut:
print(plotNucleotideDiff(freq_hnRNPC_WT, freq_hnRNPC_Mut, "hnRNPC WT", "hnRNPC Mut"))

## hnRNPC WT in RBM25 WT vs Mut:
print(plotNucleotideDiff(freq_hnRNPC_WT_inRBM25_WT, freq_hnRNPC_WT_inRBM25_Mut, 
                         "hnRNPC WT in RBM25 WT", "hnRNPC WT in RBM25 Mut"))

## hnRNPC WT in KD vs hnRNPC Mut in KD:
print(plotNucleotideDiff(freq_hnRNPC_WT_inKD, freq_hnRNPC_Mut_inKD, 
                         "hnRNPC WT in KD", "hnRNPC Mut in KD"))

## RBM25 WT vs Mut:
print(plotNucleotideDiff(freq_RBM25_WT, freq_RBM25_Mut, "RBM25 WT", "RBM25 Mut"))

## hnRNPC WT vs hnRNPC WT in RBM25 WT:
print(plotNucleotideDiff(freq_hnRNPC_WT, freq_hnRNPC_WT_inRBM25_WT, 
                         "hnRNPC WT", "hnRNPC WT in RBM25 WT"))

## hnRNPC WT vs hnRNPC WT in KD:
print(plotNucleotideDiff(freq_hnRNPC_WT, freq_hnRNPC_WT_inKD, 
                         "hnRNPC WT", "hnRNPC WT in KD"))
################################################################################
