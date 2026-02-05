################################################################################
## Generate Distribution from RBNS using RBPSpecificity Package:
## Written by Soon Yi
## Created: 2025-08-08
## Last Edited: 2025-08-25
## Figure 1
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RBPSpecificity)

## Set some basic parameters:
################################################################################
K = 5
baseDir = '~/Desktop/Genomics/Specificity/Data_Final/'
RBNS = read_csv(paste0(baseDir, 'RBNS/RBNS_normalized_', K, 'mer.csv'), col_names = T, show_col_types = F)
colnames(RBNS) = c('Motif', colnames(RBNS)[2:ncol(RBNS)])
RBPs = colnames(RBNS)[2:ncol(RBNS)]
################################################################################

# Read RBNS file:
hnRNPC = 'hnRNPC'
hnRNPC = RBNS[, c('Motif', 'HNRNPC')]
colnames(hnRNPC) = c('MOTIF', 'Score')

RBM25 = 'RBM25'
RBM25 = RBNS[, c('Motif', 'RBM25')]
colnames(RBM25) = c('MOTIF', 'Score')


ggplot() +
  geom_histogram(data = hnRNPC, aes(x = Score), binwidth = 0.001, fill = "black", alpha = 1.0) +
  labs(title = "hnRNPC", x = "RBNS Score", y = "Frequency") +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 25)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))

ggplot() +
  geom_histogram(data = RBM25, aes(x = Score), binwidth = 0.001, fill = "black", alpha = 1.0) +
  labs(title = "RBM25", x = "RBNS Score", y = "Frequency") +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 25)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))

hnRNPC_MS = plotMS(hnRNPC)
RBM25_MS = plotMS(RBM25)

returnMS(hnRNPC, output_type = 'number')     # 0.86
returnMS(RBM25, output_type = 'number')      # 0.22


returnIS(hnRNPC)                             # 50.8
returnIS(RBM25)                              # 2.45














# Histogram Comparisons
################################################################################
# Parameters for the first distribution
num_elements1 = 1024
mean_score1 = 20
sd_score1 = 5

# Parameters for the second distribution
num_elements2 = 1024
mean_score2 = 30
sd_score2 = 20

# General plot parameters
num_bins = 150
wid_bins = 0.5

scores1 = rnorm(n = num_elements1, mean = mean_score1, sd = sd_score1)
scores2 = rnorm(n = num_elements2, mean = mean_score2, sd = sd_score2)

# Create two separate data frames, adding a column to identify the source distribution.
data1 = data.frame(Score = scores1, Distribution = "Group 1")
data2 = data.frame(Score = scores2, Distribution = "Group 2")
combined_data = rbind(data1, data2)

ggplot(combined_data, aes(x = Score, fill = Distribution)) +
  geom_histogram(bins = num_bins,
                 # binwidth = wid_bins,
                 alpha = 1,
                 position = "identity",
                 color = "white") +
  scale_fill_manual(values = c("Group 1" = "#bf616a", "Group 2" = "#5e81ac")) +
  labs(title = "Comparison of Two Normal Distributions",
       x = "Score",
       y = "Frequency (Count)",
       fill = "Distribution Group")  +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  scale_x_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))

