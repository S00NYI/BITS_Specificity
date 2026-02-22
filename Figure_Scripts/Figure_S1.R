################################################################################
## Generate Distribution from RBNS using RBPSpecificity Package:
## Written by Soon Yi
## Created: 2025-08-08
## Last Edited: 2026-02-21
## Figure S1
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(seqLogo)
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

## Load and Process Data
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

# Calculate Inherent Specificityes for RBPs:
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

## Plot Inherent Specificity rank across different K-mers:
################################################################################
long_data_rank = InherentSpecificity %>%
  gather(key = "mer", value = "value", -RBP) %>%
  group_by(mer) %>%
  mutate(rank = rank(-value)) %>%
  ungroup()

long_data_rank[which(is.na(long_data_rank$value)), ]$rank = NA

# In Heatmap Form:
ranked_df = long_data_rank %>%
  select(RBP, mer, rank) %>%
  pivot_wider(names_from = mer, values_from = rank, values_fill = NA) %>%
  arrange(RBP)

heatmap_data = ranked_df %>% column_to_rownames(var = 'RBP')
avg_ranks = rowMeans(heatmap_data, na.rm = TRUE)
ordered_heatmap_data = heatmap_data[order(avg_ranks), ]

SI_data_for_heatmap = InherentSpecificity %>% column_to_rownames(var = 'RBP') %>% as.matrix()
SI_data_for_heatmap = round(SI_data_for_heatmap, 1)
SI_data_for_heatmap = SI_data_for_heatmap[rownames(ordered_heatmap_data), ]

pheatmap(ordered_heatmap_data,
         display_numbers = SI_data_for_heatmap,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)),
         na_col = "grey",
         number_color = 'black')

################################################################################

## Plot Mutational Sensitivity rank across different K-mers:
################################################################################
long_data_rank = MutationalSensitivity %>%
  gather(key = "mer", value = "value", -RBP) %>%
  group_by(mer) %>%
  mutate(rank = rank(-value)) %>%
  ungroup()

long_data_rank[which(is.na(long_data_rank$value)), ]$rank = NA

# In Heatmap Form:
ranked_df = long_data_rank %>%
  select(RBP, mer, rank) %>%
  pivot_wider(names_from = mer, values_from = rank, values_fill = NA) %>%
  arrange(RBP)

heatmap_data = ranked_df %>% column_to_rownames(var = 'RBP')
avg_ranks = rowMeans(heatmap_data, na.rm = TRUE)
ordered_heatmap_data = heatmap_data[order(avg_ranks), ]

MS_data_for_heatmap = MutationalSensitivity %>% column_to_rownames(var = 'RBP') %>% as.matrix()
MS_data_for_heatmap = round(MS_data_for_heatmap, 2)
MS_data_for_heatmap = MS_data_for_heatmap[rownames(ordered_heatmap_data), ]

pheatmap(ordered_heatmap_data,
         display_numbers = MS_data_for_heatmap,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)),
         na_col = "grey",
         number_color = 'black')

################################################################################

## Plot IS vs MS Between K-mers:
################################################################################
IS_MS_4 = data.frame(RBP = InherentSpecificity$RBP, IS = InherentSpecificity$'4mer', MS = MutationalSensitivity$'4mer')
IS_MS_6 = data.frame(RBP = InherentSpecificity$RBP, IS = InherentSpecificity$'6mer', MS = MutationalSensitivity$'6mer')
IS_MS_7 = data.frame(RBP = InherentSpecificity$RBP, IS = InherentSpecificity$'7mer', MS = MutationalSensitivity$'7mer')

r = cor(IS_MS_4$IS, IS_MS_4$MS, use = 'complete.obs')
ggplot(IS_MS_4, aes(x = IS, y = MS, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Specificity Index", y = "4mer Mutational Sensitivity", title = "4-mer: SI vs MS") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 100, by = 5)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


r = cor(IS_MS_6$IS, IS_MS_6$MS, use = 'complete.obs')
ggplot(IS_MS_6, aes(x = IS, y = MS, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "6mer Specificity Index", y = "6mer Mutational Sensitivity", title = "6-mer: SI vs MS") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 160), breaks = seq(0, 160, by = 20)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) 

r = cor(IS_MS_7$IS, IS_MS_7$MS, use = 'complete.obs')
ggplot(IS_MS_7, aes(x = IS, y = MS, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "7mer Specificity Index", y = "7mer Mutational Sensitivity", title = "7-mer: SI vs MS") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 320), breaks = seq(0, 320, by = 40)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
################################################################################



## NOT INCLUDED: Plot Inherent Specificity across different K-mers:
################################################################################
long_data = InherentSpecificity %>%
  pivot_longer(cols = c('4mer', '5mer', '6mer', '7mer'),
               names_to = 'mer',
               values_to = 'value')

ggplot(long_data, aes(x = mer, y = value, group = RBP, color = RBP)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "K-Mer", y = "Inherent Specificity", title = "SI per RBP across Different K-mers") +
  scale_y_continuous(limits = c(0, 275),
                     breaks = seq(0, 300, by = 25))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))
################################################################################

## NOT INCLUDED: Plot Inherent Specificity Comparison Between Kmers: 
################################################################################
# 4 v 5
r = cor(InherentSpecificity$'4mer', InherentSpecificity$'5mer', use = 'complete.obs')
ggplot(InherentSpecificity, aes(x = `4mer`, y = `5mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Inherent Specificity", y = "5mer Inherent Specificity", title = "SI Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 100, by = 20)) +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 100, by = 20)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# 4 v 6
r = cor(InherentSpecificity$'4mer', InherentSpecificity$'6mer', use = 'complete.obs')
ggplot(InherentSpecificity, aes(x = `4mer`, y = `6mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Inherent Specificity", y = "6mer Inherent Specificity", title = "SI Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 25)) +
  scale_x_continuous(limits = c(0, 125), breaks = seq(0, 125, by = 25)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# 4 v 7
r = cor(InherentSpecificity$'4mer', InherentSpecificity$'7mer', use = 'complete.obs')
ggplot(InherentSpecificity, aes(x = `4mer`, y = `7mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Inherent Specificity", y = "7mer Inherent Specificity", title = "SI Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 275), breaks = seq(0, 300, by = 25)) +
  scale_x_continuous(limits = c(0, 275), breaks = seq(0, 300, by = 25)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# 5 v 6
r = cor(InherentSpecificity$'5mer', InherentSpecificity$'6mer', use = 'complete.obs')
ggplot(InherentSpecificity, aes(x = `5mer`, y = `6mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "5mer Inherent Specificity", y = "6mer Inherent Specificity", title = "SI Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 125), breaks = seq(0, 150, by = 25)) +
  scale_x_continuous(limits = c(0, 125), breaks = seq(0, 150, by = 25)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# 5 v 7
r = cor(InherentSpecificity$'5mer', InherentSpecificity$'7mer', use = 'complete.obs')
ggplot(InherentSpecificity, aes(x = `5mer`, y = `7mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "5mer Inherent Specificity", y = "7mer Inherent Specificity", title = "SI Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 275), breaks = seq(0, 275, by = 25)) +
  scale_x_continuous(limits = c(0, 275), breaks = seq(0, 275, by = 25)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# 6 v 7
r = cor(InherentSpecificity$'6mer', InherentSpecificity$'7mer', use = 'complete.obs')
ggplot(InherentSpecificity, aes(x = `6mer`, y = `7mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "6mer Inherent Specificity", y = "7mer Inherent Specificity", title = "SI Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 275), breaks = seq(0, 275, by = 25)) +
  scale_x_continuous(limits = c(0, 275), breaks = seq(0, 275, by = 25)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
################################################################################

## NOT INCLUDED: Plot Mutational Sensitivity across different K-mers:
################################################################################
long_data = MutationalSensitivity %>%
  pivot_longer(cols = c('4mer', '5mer', '6mer', '7mer'),
               names_to = 'mer',
               values_to = 'value')

ggplot(long_data, aes(x = mer, y = value, group = RBP, color = RBP)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "K-Mer", y = "Mutational Sensitivity", title = "MS per RBP across Different K-mers") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.1))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))
################################################################################

## NOT INCLUDED: Mutational Sensitivity Comparison Between Kmers:
################################################################################
# 4 v 5
r = cor(MutationalSensitivity$'4mer', MutationalSensitivity$'5mer', use = 'complete.obs')
ggplot(MutationalSensitivity, aes(x = `4mer`, y = `5mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Mutational Sensitivity", y = "5mer Mutational Sensitivity", title = "MS Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# 4 v 6
r = cor(MutationalSensitivity$'4mer', MutationalSensitivity$'6mer', use = 'complete.obs')
ggplot(MutationalSensitivity, aes(x = `4mer`, y = `6mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Mutational Sensitivity", y = "6mer Mutational Sensitivity", title = "MS Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) 

# 4 v 7
r = cor(MutationalSensitivity$'4mer', MutationalSensitivity$'7mer', use = 'complete.obs')
ggplot(MutationalSensitivity, aes(x = `4mer`, y = `7mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "4mer Mutational Sensitivity", y = "7mer Mutational Sensitivity", title = "MS Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# 5 v 6
r = cor(MutationalSensitivity$'5mer', MutationalSensitivity$'6mer', use = 'complete.obs')
ggplot(MutationalSensitivity, aes(x = `5mer`, y = `6mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "5mer Mutational Sensitivity", y = "6mer Mutational Sensitivity", title = "MS Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# 5 v 7
r = cor(MutationalSensitivity$'5mer', MutationalSensitivity$'7mer', use = 'complete.obs')
ggplot(MutationalSensitivity, aes(x = `5mer`, y = `7mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "5mer Mutational Sensitivity", y = "7mer Mutational Sensitivity", title = "MS Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

# 6 v 7
r = cor(MutationalSensitivity$'6mer', MutationalSensitivity$'7mer', use = 'complete.obs')
ggplot(MutationalSensitivity, aes(x = `6mer`, y = `7mer`, label = RBP)) +
  geom_point() +
  geom_text_repel() +
  labs(x = "6mer Mutational Sensitivity", y = "7mer Mutational Sensitivity", title = "MS Comparison between K-mers") +
  annotate("text", x = Inf, y = Inf, label = paste0("R == ", format(r, digits = 2)), parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme_bw() +
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))
################################################################################
