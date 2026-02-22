################################################################################
## PWM_to_Score 
## Written by Soon Yi
## Created: 2023-11-08
## Last Edited: 2023-11-22
## Figure 1
################################################################################

library(dplyr)
library(Biostrings)

# Set some basic parameters:
################################################################################
baseDir = '~/Desktop/Genomics/Specificity/'
K = 5
RBNS = read_csv(paste0(baseDir, 'Data/RBNS/RBNS_normalized_', K, 'mer.csv'), col_names = T, show_col_types = F)
colnames(RBNS) = c('Motif', colnames(RBNS)[2:ncol(RBNS)])
RBPs = colnames(RBNS)[2:ncol(RBNS)]

# All combinations of K-mers:
nts_id = c('A', 'C', 'G', 'U')
Kmers = expand.grid(rep(list(nts_id), K))
colnames(Kmers) = c(1, 2, 3, 4, 5)
################################################################################

################################################################################
RBP = 'HNRNPC'
RBP = 'EIF4G2'

data = data.frame(RBNS[, c('Motif', RBP)])
data = data[order(-data[, RBP]), ]
rownames(data) = NULL

numMotif = 10

motif = DNAStringSet(gsub("U", "T", data$Motif[1:numMotif]))
pwm = PWM(consensusMatrix(motif))
pwm = pwm - min(pwm)

pwm[, 1] = pwm[, 1]/sum(pwm[, 1])
pwm[, 2] = pwm[, 2]/sum(pwm[, 2])
pwm[, 3] = pwm[, 3]/sum(pwm[, 3])
pwm[, 4] = pwm[, 4]/sum(pwm[, 4])
pwm[, 5] = pwm[, 5]/sum(pwm[, 5])

PWM = makePWM(pwm, alphabet = 'RNA')

# Container for scores:
Scores = matrix(NA, nrow = 1024, ncol = 5)
colnames(Scores) = c(1, 2, 3, 4, 5)

# Calculate scores per motif:
for (id in 1:nrow(Kmers)) {
  Kmer = Kmers[id, ]
  for (pos in 1:ncol(Kmer)) {
    Scores[id, pos] = PWM@pwm[Kmer[[pos]], pos]
  }
}

PWM_Motif_Scores = data.frame(Motif = apply(Kmers, 1, paste, collapse = ""),
                              Score = rowSums(Scores))

PWM_Motif_Scores_Normalized = PWM_Motif_Scores
PWM_Motif_Scores_Normalized$Score = (PWM_Motif_Scores_Normalized$Score - min(PWM_Motif_Scores_Normalized$Score)) / (max(PWM_Motif_Scores_Normalized$Score) - min(PWM_Motif_Scores_Normalized$Score))

baseDir = '~/Desktop/Genomics/Specificity/Data/Specific_RBP/'
write.csv(PWM_Motif_Scores_Normalized, paste0(baseDir, RBP, '_Counts_PWM_Motif_ScoresNormed_new.csv'), quote = F)
################################################################################
