################################################################################
## Generate PWM/Logo from RBNS
## Written by Soon Yi
## Created: 2023-10-20
## Last Edited: 2025-10-30
## Figure 1
################################################################################

library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(Biostrings)
library(seqLogo)


## Set some basic parameters:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/Data_Final/'

K = 5
RBNS = read_csv(paste0(baseDir, 'RBNS/RBNS_normalized_', K, 'mer.csv'), col_names = T, show_col_types = F)
colnames(RBNS) = c('Motif', colnames(RBNS)[2:ncol(RBNS)])
RBPs = colnames(RBNS)[2:ncol(RBNS)]
################################################################################

## For hnRNPC and EIF4G2:
################################################################################

# returnPWM = function(data, sd_multipler) {
#   data = data.frame(data[, c('MOTIF', 'Score')])
#   data = data[order(-data[, 'Score']), ]
#   rownames(data) = NULL
#
#   data = data %>% filter(Score > (mean(data$Score) + sd_multipler*sd(data$Score)))
#
#   motif = DNAStringSet(gsub("U", "T", data$MOTIF))
#   pwm = PWM(consensusMatrix(motif))
#   pwm = pwm - min(pwm)
#
#   pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
#   pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
#   pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
#   pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
#   pwm[, 5] = pwm[, 5]/sum(pwm[, 5])
#
#   PWM = makePWM(pwm, alphabet = 'RNA')
#   seqLogo(PWM, ic.scale = F)
#   print(paste0("Used ", as.character(nrow(data)), " Top Motifs...."))
# }


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

HNRNPC = data.frame(RBNS[, c('Motif', 'HNRNPC')])
colnames(HNRNPC) = c('MOTIF', 'Score')

RBM25 = data.frame(RBNS[, c('Motif', 'RBM25')])
colnames(RBM25) = c('MOTIF', 'Score')

# EIF4G2 = data.frame(RBNS[, c('Motif', 'EIF4G2')])
# colnames(EIF4G2) = c('MOTIF', 'Score')

## 2 standard deviation above mean --> 95 percentile
returnPWM(HNRNPC, 2)
returnPWM(RBM25, 2)
# returnPWM(EIF4G2, 3)



################################################################################





