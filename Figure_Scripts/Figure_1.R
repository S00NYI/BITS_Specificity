################################################################################
## Generate Distribution from RBNS using RBPSpecificity Package:
## Written by Soon Yi
## Created: 2025-08-08
## Last Edited: 2026-02-21
## Figure 1
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(seqLogo)
library(Biostrings)
library(RBPSpecificity)
library(tibble)
library(ggrepel) 
library(pheatmap)
library(RColorBrewer)

## Set some basic parameters:
################################################################################
K = 5
baseDir = '~/Repos/BITS_Specificity/Dataset/Analysis/'
RBNS = read_csv(paste0(baseDir, 'RBNS/RBNS_normalized_', K, 'mer.csv'), col_names = T, show_col_types = F)
colnames(RBNS) = c('Motif', colnames(RBNS)[2:ncol(RBNS)])
RBPs = colnames(RBNS)[2:ncol(RBNS)]
################################################################################

## Set custom functions:
################################################################################
returnPWM = function(data, sd_multipler) {
  
  # --- 1. Filter Data ---
  data = data.frame(data[, c('MOTIF', 'Score')])
  data = data[order(-data[, 'Score']), ]
  rownames(data) = NULL
  
  data = data %>% filter(Score > (mean(data$Score) + sd_multipler*sd(data$Score)))
  
  if (nrow(data) == 0) {
    print("No motifs passed the filter.")
    return(NULL)
  }
  
  # --- 2. Get PFM (Position Frequency Matrix) ---
  motif = DNAStringSet(gsub("U", "T", data$MOTIF))
  pfm_raw = consensusMatrix(motif)
  
  # --- 3. NEW: Create a clean 4-row PFM (A, C, G, T) ---
  # This handles cases where 'N' (or other) is present,
  # or where a base is missing entirely.
  
  # Define the 4 required bases
  bases = c("A", "C", "G", "T")
  
  # Create an empty 4-row matrix with the correct dimensions
  pfm = matrix(0, nrow = 4, ncol = ncol(pfm_raw))
  rownames(pfm) = bases
  colnames(pfm) = colnames(pfm_raw)
  
  # Find the common rows between raw PFM and our 4 bases
  common_rows = intersect(bases, rownames(pfm_raw))
  
  # Copy the counts from the raw pfm to the new one
  if (length(common_rows) > 0) {
    pfm[common_rows, ] = pfm_raw[common_rows, ]
  }
  # --- End of NEW section ---
  
  # --- 4. Create PPM (Position Probability Matrix) ---
  # Use the *clean* pfm
  ppm = prop.table(pfm, margin = 2)
  
  # Handle 0 probabilities to avoid NaN in log2(0)
  # (prop.table already makes this, but good to be safe)
  if(any(rowSums(ppm) == 0)) {
    print("Warning: One or more bases (A,C,G,T) are not present in any filtered motif.")
  }
  
  # --- 5. Plot SeqLogo ---
  PWM = makePWM(ppm, alphabet = 'RNA')
  seqLogo(PWM, ic.scale = F) # ic.scale=F plots entropy height
  
  print(paste0("Used ", as.character(nrow(data)), " Top Motifs...."))
  
  # --- 6. Calculate Entropy and IC ---
  
  # We must handle P_b,l = 0, because 0 * log2(0) = 0.
  # R calculates 0 * log2(0) as NaN (0 * -Inf), so we replace NaNs with 0.
  
  # Matrix of P_b,l * log2(P_b,l)
  ppm_log_ppm = ppm * log2(ppm)
  ppm_log_ppm[is.nan(ppm_log_ppm)] = 0
  
  # H(l) = Shannon entropy per position (column)
  H_l = -colSums(ppm_log_ppm)
  
  # H_max = max possible entropy (log2(4) for DNA/RNA = 2 bits)
  H_max = 2
  
  # IC(l) = Information Content per position
  IC_l = H_max - H_l
  
  # --- 7. Calculate Total Motif Entropy and IC ---
  H_motif = sum(H_l)
  IC_motif = sum(IC_l)
  
  # --- 8. Print Results ---
  print(paste0("Total Motif Entropy (H_motif): ", round(H_motif, 4), " bits"))
  print(paste0("Total Information Content (IC_motif): ", round(IC_motif, 4), " bits"))
  
  # Return the results
  invisible(list(H_motif = H_motif, IC_motif = IC_motif, ppm = ppm))
}

Motif_Variants = function(motif, nucleotides) {
  variants = c()
  # Loop over each position in the motif
  for (i in 1:nchar(motif)) {
    # Loop over each nucleotide
    for (nuc in nucleotides) {
      if (substr(motif, i, i) != nuc) {
        # Construct the new variant
        variant = paste0(substr(motif, 1, i - 1), nuc, substr(motif, i + 1, nchar(motif)))
        variants = c(variants, variant)
      }
    }
  }
  return(variants)
}
################################################################################

## Panel B
################################################################################
RBM25 = 'RBM25'
RBM25 = RBNS[, c('Motif', 'RBM25')]
colnames(RBM25) = c('MOTIF', 'Score')

returnPWM(RBM25, 2)

ggplot() +
  geom_histogram(data = RBM25, aes(x = Score), binwidth = 0.001, fill = "black", alpha = 1.0) +
  labs(title = "RBM25", x = "RBNS Score", y = "Frequency") +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 25)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))

plotMS(RBM25)
returnMS(RBM25, output_type = 'number')      # 0.22
returnIS(RBM25)                              # 2.45
################################################################################

## Panel C
################################################################################
hnRNPC = 'hnRNPC'
hnRNPC = RBNS[, c('Motif', 'HNRNPC')]
colnames(hnRNPC) = c('MOTIF', 'Score')

returnPWM(hnRNPC, 2)

ggplot() +
  geom_histogram(data = hnRNPC, aes(x = Score), binwidth = 0.001, fill = "black", alpha = 1.0) +
  labs(title = "hnRNPC", x = "RBNS Score", y = "Frequency") +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 25)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))

plotMS(hnRNPC)
returnMS(hnRNPC, output_type = 'number')     # 0.86
returnIS(hnRNPC)                             # 50.8
################################################################################


## Load and process RBNS data for panel D and E:
################################################################################
# Read RBNS file:
RBNS_4mer = read_csv(paste0(baseDir, 'RBNS/RBNS_singleRBD_normalized_4mer.csv'), col_names = T, show_col_types = F)
RBNS_5mer = read_csv(paste0(baseDir, 'RBNS/RBNS_singleRBD_normalized_5mer.csv'), col_names = T, show_col_types = F)
RBNS_6mer = read_csv(paste0(baseDir, 'RBNS/RBNS_singleRBD_normalized_6mer.csv'), col_names = T, show_col_types = F)
RBNS_7mer = read_csv(paste0(baseDir, 'RBNS/RBNS_singleRBD_normalized_7mer.csv'), col_names = T, show_col_types = F)

RBPs_4mer = colnames(RBNS_4mer)[2:ncol(RBNS_4mer)]
RBPs_5mer = colnames(RBNS_5mer)[2:ncol(RBNS_5mer)]
RBPs_6mer = colnames(RBNS_6mer)[2:ncol(RBNS_6mer)]
RBPs_7mer = colnames(RBNS_7mer)[2:ncol(RBNS_7mer)]

# Calculate Specificity Indexes for RBPs:
InherentSpecificity = as.data.frame(matrix(NA, nrow = 26, ncol = 5))
colnames(InherentSpecificity) = c('RBP', '4mer', '5mer', '6mer', '7mer')
InherentSpecificity$RBP = RBPs_5mer

Ks = c(4, 5, 6, 7)

for (K in Ks) {
  RBNS_Data = get(paste0('RBNS_', K, 'mer'))
  RBP_Data = get(paste0('RBPs_', K, 'mer'))
  for (RBP in RBP_Data) {
    temp = data.frame(MOTIF = unlist(RBNS_Data[, 'Motif']),
                      Score = unlist(RBNS_Data[, RBP]))
    rownames(temp) = NULL
    InherentSpecificity[which(InherentSpecificity$RBP == RBP), paste0(K, 'mer')] = returnIS(temp)
  }
}

# Calculate Mutational Sensitivity for RBPs:
MutationalSensitivity = as.data.frame(matrix(NA, nrow = 26, ncol = 5))
colnames(MutationalSensitivity) = c('RBP', '4mer', '5mer', '6mer', '7mer')
MutationalSensitivity$RBP = RBPs_5mer

Ks = c(4, 5, 6, 7)

for (K in Ks) {
  RBNS_Data = get(paste0('RBNS_', K, 'mer'))
  RBP_Data = get(paste0('RBPs_', K, 'mer'))
  for (RBP in RBP_Data) {
    temp = data.frame(MOTIF = unlist(RBNS_Data[, 'Motif']),
                      Score = unlist(RBNS_Data[, RBP]))
    rownames(temp) = NULL
    MutationalSensitivity[which(MutationalSensitivity$RBP == RBP), paste0(K, 'mer')] = returnMS(temp, output_type = 'number')
  }
}

################################################################################

## Panel D: IS vs MS Between Kmers:
################################################################################
IS_MS_5 = data.frame(RBP = InherentSpecificity$RBP, IS = InherentSpecificity$'5mer', MS = MutationalSensitivity$'5mer')

r = cor(IS_MS_5$IS, IS_MS_5$MS, use = 'complete.obs')
ggplot(IS_MS_5, aes(x = IS, y = MS, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "5mer Specificity Index", y = "5mer Mutational Sensitivity", title = "5-mer: SI vs MS") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 100, by = 10)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
################################################################################

## Panel E: Heatmap of correlations:
################################################################################
IS_Data = InherentSpecificity[, c('4mer', '5mer', '6mer', '7mer')]
colnames(IS_Data) = c('IS_4', 'IS_5', 'IS_6', 'IS_7')
MS_Data = MutationalSensitivity[, c('4mer', '5mer', '6mer', '7mer')]
colnames(MS_Data) = c('MS_4', 'MS_5', 'MS_6', 'MS_7')

# Generate SI heatmap
KMer_Data = IS_Data
CorrMatrix = cor(KMer_Data,  use = "pairwise.complete.obs")
CorrMatrix = matrix(round(CorrMatrix,2), nrow = ncol(KMer_Data))
rownames(CorrMatrix) = colnames(KMer_Data)
colnames(CorrMatrix) = colnames(KMer_Data)

color_palette = (colorRampPalette(brewer.pal(9, "GnBu"))(100))
# color_palette = (colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
breaks = seq(0.5, 1, length.out = length(color_palette) + 1)
pheatmap(CorrMatrix, cluster_rows = F, cluster_cols = F, display_numbers = T, color = color_palette, breaks = breaks)

# Generate MS heatmap
KMer_Data = MS_Data
CorrMatrix = cor(KMer_Data,  use = "pairwise.complete.obs")
CorrMatrix = matrix(round(CorrMatrix,2), nrow = ncol(KMer_Data))
rownames(CorrMatrix) = colnames(KMer_Data)
colnames(CorrMatrix) = colnames(KMer_Data)

color_palette = (colorRampPalette(brewer.pal(9, "GnBu"))(100))
# color_palette = (colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
breaks = seq(0.5, 1, length.out = length(color_palette) + 1)
pheatmap(CorrMatrix, cluster_rows = F, cluster_cols = F, display_numbers = T, color = color_palette, breaks = breaks)
################################################################################
