################################################################################
## CLIP Peak Motif Analysis
## Author: Soon Yi with Antigravity
## Date: February 2026
################################################################################

# === Load Libraries ===
library(RBPSpecificity)
library(ggplot2)
library(dplyr)
library(data.table)
library(seqLogo)
library(Biostrings)
library(cowplot)
library(grid)
library(gridExtra)

## 1. Basic Setup:
################################################################################
# BASE_DIR = "F:/Specificity/CLIP/Analysis/"
# INPUT_DIR = paste0(BASE_DIR, "output/")
# OUTPUT_BASE = paste0(BASE_DIR, "output/")

BASE_DIR = "/Volumes/1TB_Data/Specificity/CLIP/Analysis/"
INPUT_DIR = paste0(BASE_DIR, "output/")
OUTPUT_BASE = paste0(BASE_DIR, "output/")


SAMPLES = c("hnRNPC_WT", "hnRNPC_Mut",
           "hnRNPC_WT_inRBM25_WT", "hnRNPC_WT_inRBM25_Mut",
           "hnRNPC_WT_inKD", "hnRNPC_Mut_inKD",
           "RBM25_WT", "RBM25_Mut")

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

  "ext15_n10" = c(15, -10),
  "ext15_n15" = c(15, -15),
  "ext20_n10" = c(20, -10),
  "ext20_n15" = c(20, -15),
  "ext25_n10" = c(25, -10),
  "ext25_n15" = c(25, -15),
  "ext30_n10" = c(30, -10),
  "ext30_n15" = c(30, -15)
)

## Bootstrap ext15_n10 or ext15_n15 multiple times to get statistically significant result for comparison.


# GLOBAL_REGIONS = c("WholeGenome")
################################################################################

## 2. Motif Enrichment:
################################################################################

# A. Local Background Tests (with scramble ON and OFF)
SCRAMBLE_OPTIONS = c(FALSE)
# SCRAMBLE_OPTIONS = c(FALSE)

for (ext_name in names(LOCAL_EXTENSIONS)) {
    ext = LOCAL_EXTENSIONS[[ext_name]]

    for (scramble in SCRAMBLE_OPTIONS) {
        scramble_suffix = ifelse(scramble, "scrambleON", "scrambleOFF")
        CONDITION_DIR = paste0(OUTPUT_BASE, "local_", ext_name, "_", scramble_suffix, "/")
        if (!dir.exists(CONDITION_DIR)) dir.create(CONDITION_DIR)

        message(paste0("Running local background: ", ext_name, " (scramble=", scramble, ")"))

        for (SAMPLE in SAMPLES) {
            FILE_PATH = paste0(INPUT_DIR, "peaks_", SAMPLE, ".txt")
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

# # B. Global Background Tests
# for (region in GLOBAL_REGIONS) {
#     CONDITION_DIR = paste0(OUTPUT_BASE, "global_", region, "/")
#     if (!dir.exists(CONDITION_DIR)) dir.create(CONDITION_DIR)
# 
#     message(paste0("Running global background: ", region))
# 
#     for (SAMPLE in SAMPLES) {
#         FILE_PATH = paste0(INPUT_DIR, "peaks_", SAMPLE, ".txt")
#         peaks = fread(FILE_PATH)
#         peak_data = peaks[, .(chr, start, end, name, score, strand)]
# 
#         motif_enrich = motifEnrichment(peak_data, 'hg38',
#                                         K = K_MER,
#                                         extension = c(0, 0),
#                                         nucleic_acid_type = 'RNA',
#                                         enrichment_method = METHOD,
#                                         background_type = "global",
#                                         global_background_region = region)
# 
#         fwrite(motif_enrich, paste0(CONDITION_DIR, "motifEnrichment_", SAMPLE, ".txt"), sep = "\t")
#     }
# }
################################################################################

## 3. Affinity Distribution & PWM (per condition):
################################################################################
# Build all condition directory paths
local_dirs = c()
for (ext_name in names(LOCAL_EXTENSIONS)) {
    local_dirs = c(local_dirs,
                   paste0(OUTPUT_BASE, "local_", ext_name, "_scrambleON/"),
                   paste0(OUTPUT_BASE, "local_", ext_name, "_scrambleOFF/"))
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

    motifEnrichment_hnRNPC_WT = fread(paste0(CONDITION_DIR, "motifEnrichment_hnRNPC_WT.txt"))
    motifEnrichment_hnRNPC_Mut = fread(paste0(CONDITION_DIR, "motifEnrichment_hnRNPC_Mut.txt"))
    motifEnrichment_hnRNPC_inRBM25_WT = fread(paste0(CONDITION_DIR, "motifEnrichment_hnRNPC_WT_inRBM25_WT.txt"))
    motifEnrichment_hnRNPC_inRBM25_Mut = fread(paste0(CONDITION_DIR, "motifEnrichment_hnRNPC_WT_inRBM25_Mut.txt"))
    motifEnrichment_hnRNPC_WT_inKD = fread(paste0(CONDITION_DIR, "motifEnrichment_hnRNPC_WT_inKD.txt"))
    motifEnrichment_hnRNPC_Mut_inKD = fread(paste0(CONDITION_DIR, "motifEnrichment_hnRNPC_Mut_inKD.txt"))
    motifEnrichment_RBM25_WT = fread(paste0(CONDITION_DIR, "motifEnrichment_RBM25_WT.txt"))
    motifEnrichment_RBM25_Mut = fread(paste0(CONDITION_DIR, "motifEnrichment_RBM25_Mut.txt"))

    affinityDistribution_hnRNPC = ggplot() +
      geom_histogram(data = motifEnrichment_hnRNPC_WT, aes(x = Score), binwidth = 0.005, fill = "blue", alpha = 0.5) +
      geom_histogram(data = motifEnrichment_hnRNPC_Mut, aes(x = Score), binwidth = 0.005, fill = "red", alpha = 0.5) +
      labs(title = paste0("hnRNPC WT vs Mut (", condition_name, ")"), x = "Motif Enrichment", y = "Frequency") +
      scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by = 50)) +
      theme_bw() +
      theme(axis.text = element_text(size=14),
            axis.title = element_text(size=14, face = 'bold'))

    affinityDistribution_hnRNPC_inRBM25 = ggplot() +
      geom_histogram(data = motifEnrichment_hnRNPC_WT, aes(x = Score), binwidth = 0.005, fill = "black", alpha = 0.25) +
      geom_histogram(data = motifEnrichment_hnRNPC_inRBM25_WT, aes(x = Score), binwidth = 0.005, fill = "blue", alpha = 0.5) +
      geom_histogram(data = motifEnrichment_hnRNPC_inRBM25_Mut, aes(x = Score), binwidth = 0.005, fill = "red", alpha = 0.5) +
      labs(title = paste0("hnRNPC WT in RBM25 (", condition_name, ")"), x = "Motif Enrichment", y = "Frequency") +
      scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by = 50)) +
      theme_bw() +
      theme(axis.text = element_text(size=14),
            axis.title = element_text(size=14, face = 'bold'))
    
    affinityDistribution_hnRNPC_inKD = ggplot() +
      geom_histogram(data = motifEnrichment_hnRNPC_WT, aes(x = Score), binwidth = 0.005, fill = "black", alpha = 0.25) +
      geom_histogram(data = motifEnrichment_hnRNPC_WT_inKD, aes(x = Score), binwidth = 0.005, fill = "blue", alpha = 0.5) +
      geom_histogram(data = motifEnrichment_hnRNPC_Mut_inKD, aes(x = Score), binwidth = 0.005, fill = "red", alpha = 0.5) +
      labs(title = paste0("hnRNPC WT vs Mut in KD (", condition_name, ")"), x = "Motif Enrichment", y = "Frequency") +
      scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by = 50)) +
      theme_bw() +
      theme(axis.text = element_text(size=14),
            axis.title = element_text(size=14, face = 'bold'))

    affinityDistribution_RBM25 = ggplot() +
      geom_histogram(data = motifEnrichment_RBM25_WT, aes(x = Score), binwidth = 0.005, fill = "blue", alpha = 0.5) +
      geom_histogram(data = motifEnrichment_RBM25_Mut, aes(x = Score), binwidth = 0.005, fill = "red", alpha = 0.5) +
      labs(title = paste0("RBM25 WT vs Mut (", condition_name, ")"), x = "Motif Enrichment", y = "Frequency") +
      scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by = 50)) +
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
        capturePWMGrob(motifEnrichment_hnRNPC_WT, PWM_SD, "hnRNPC WT"),
        capturePWMGrob(motifEnrichment_hnRNPC_Mut, PWM_SD, "hnRNPC Mut"),
        capturePWMGrob(motifEnrichment_RBM25_WT, PWM_SD, "RBM25 WT"),
        capturePWMGrob(motifEnrichment_RBM25_Mut, PWM_SD, "RBM25 Mut"),
        capturePWMGrob(motifEnrichment_hnRNPC_inRBM25_WT, PWM_SD, "hnRNPC WT in RBM25 WT"),
        capturePWMGrob(motifEnrichment_hnRNPC_inRBM25_Mut, PWM_SD, "hnRNPC WT in RBM25 Mut"),
        capturePWMGrob(motifEnrichment_hnRNPC_WT_inKD, PWM_SD, "hnRNPC WT in KD"),
        capturePWMGrob(motifEnrichment_hnRNPC_Mut_inKD, PWM_SD, "hnRNPC Mut in KD")
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


## 5. IS/MS Summary Tables:
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

## 6. Focused Analysis: ext10_n15 scrambleOFF
################################################################################
FOCUS_DIR = paste0(OUTPUT_BASE, "local_ext15_n15_scrambleOFF/")

motifEnrichment_hnRNPC_WT = fread(paste0(FOCUS_DIR, "motifEnrichment_hnRNPC_WT.txt"))
motifEnrichment_hnRNPC_Mut = fread(paste0(FOCUS_DIR, "motifEnrichment_hnRNPC_Mut.txt"))
motifEnrichment_hnRNPC_WT_inRBM25_WT = fread(paste0(FOCUS_DIR, "motifEnrichment_hnRNPC_WT_inRBM25_WT.txt"))
motifEnrichment_hnRNPC_WT_inRBM25_Mut = fread(paste0(FOCUS_DIR, "motifEnrichment_hnRNPC_WT_inRBM25_Mut.txt"))
motifEnrichment_hnRNPC_WT_inKD = fread(paste0(FOCUS_DIR, "motifEnrichment_hnRNPC_WT_inKD.txt"))
motifEnrichment_hnRNPC_Mut_inKD = fread(paste0(FOCUS_DIR, "motifEnrichment_hnRNPC_Mut_inKD.txt"))
motifEnrichment_RBM25_WT = fread(paste0(FOCUS_DIR, "motifEnrichment_RBM25_WT.txt"))
motifEnrichment_RBM25_Mut = fread(paste0(FOCUS_DIR, "motifEnrichment_RBM25_Mut.txt"))

# 6A. Affinity Distributions
affinityDistribution_hnRNPC = ggplot() +
  geom_histogram(data = motifEnrichment_hnRNPC_WT, aes(x = Score), binwidth = 0.005, fill = "blue", alpha = 0.5) +
  geom_histogram(data = motifEnrichment_hnRNPC_Mut, aes(x = Score), binwidth = 0.005, fill = "red", alpha = 0.5) +
  labs(title = "hnRNPC WT vs Mut", x = "Motif Enrichment", y = "Frequency") +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by = 50)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))

affinityDistribution_RBM25 = ggplot() +
  geom_histogram(data = motifEnrichment_RBM25_WT, aes(x = Score), binwidth = 0.005, fill = "blue", alpha = 0.5) +
  geom_histogram(data = motifEnrichment_RBM25_Mut, aes(x = Score), binwidth = 0.005, fill = "red", alpha = 0.5) +
  labs(title = "RBM25 WT vs Mut", x = "Motif Enrichment", y = "Frequency") +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by = 50)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))

affinityDistribution_hnRNPC_inRBM25 = ggplot() +
  # geom_histogram(data = motifEnrichment_hnRNPC_WT, aes(x = Score), binwidth = 0.005, fill = "black", alpha = 0.25) +
  geom_histogram(data = motifEnrichment_hnRNPC_WT_inRBM25_WT, aes(x = Score), binwidth = 0.005, fill = "blue", alpha = 0.5) +
  geom_histogram(data = motifEnrichment_hnRNPC_WT_inRBM25_Mut, aes(x = Score), binwidth = 0.005, fill = "red", alpha = 0.5) +
  labs(title = "hnRNPC WT in RBM25", x = "Motif Enrichment", y = "Frequency") +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by = 50)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))

affinityDistribution_hnRNPC_inKD = ggplot() +
  # geom_histogram(data = motifEnrichment_hnRNPC_WT, aes(x = Score), binwidth = 0.005, fill = "black", alpha = 0.25) +
  geom_histogram(data = motifEnrichment_hnRNPC_WT_inKD, aes(x = Score), binwidth = 0.005, fill = "blue", alpha = 0.5) +
  geom_histogram(data = motifEnrichment_hnRNPC_Mut_inKD, aes(x = Score), binwidth = 0.005, fill = "red", alpha = 0.5) +
  labs(title = "hnRNPC WT in KD", x = "Motif Enrichment", y = "Frequency") +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by = 50)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14, face = 'bold'),
        legend.text = element_text(size=14))

print(affinityDistribution_hnRNPC)
print(affinityDistribution_RBM25)
print(affinityDistribution_hnRNPC_inRBM25)
print(affinityDistribution_hnRNPC_inKD)

# 6B. Mutational Sensitivity Plots
plotMS(motifEnrichment_hnRNPC_WT)
plotMS(motifEnrichment_hnRNPC_Mut)
plotMS(motifEnrichment_hnRNPC_WT_inRBM25_WT)
plotMS(motifEnrichment_hnRNPC_WT_inRBM25_Mut)
plotMS(motifEnrichment_hnRNPC_WT_inKD)
plotMS(motifEnrichment_hnRNPC_Mut_inKD)
plotMS(motifEnrichment_RBM25_WT)
plotMS(motifEnrichment_RBM25_Mut)


returnIS(motifEnrichment_hnRNPC_WT)
returnIS(motifEnrichment_hnRNPC_Mut)
returnIS(motifEnrichment_hnRNPC_WT_inRBM25_WT)
returnIS(motifEnrichment_hnRNPC_WT_inRBM25_Mut)
returnIS(motifEnrichment_hnRNPC_WT_inKD)
returnIS(motifEnrichment_hnRNPC_Mut_inKD)
returnIS(motifEnrichment_RBM25_WT)
returnIS(motifEnrichment_RBM25_Mut)

# 6C. Position Weight Matrices
returnPWM(motifEnrichment_hnRNPC_WT, 2)
returnPWM(motifEnrichment_hnRNPC_Mut, 2)
returnPWM(motifEnrichment_hnRNPC_WT_inRBM25_WT, 2)
returnPWM(motifEnrichment_hnRNPC_WT_inRBM25_Mut, 2)
returnPWM(motifEnrichment_hnRNPC_WT_inKD, 2)
returnPWM(motifEnrichment_hnRNPC_Mut_inKD, 2)
returnPWM(motifEnrichment_RBM25_WT, 2)
returnPWM(motifEnrichment_RBM25_Mut, 2)

t.test(motifEnrichment_hnRNPC_WT$Score, motifEnrichment_hnRNPC_Mut$Score)
t.test(motifEnrichment_hnRNPC_WT_inRBM25_WT$Score, motifEnrichment_hnRNPC_WT_inRBM25_Mut$Score)
t.test(motifEnrichment_hnRNPC_WT_inKD$Score, motifEnrichment_hnRNPC_Mut_inKD$Score)
t.test(motifEnrichment_RBM25_WT$Score, motifEnrichment_RBM25_Mut$Score)

# Boxplot comparison:

boxplot_data = data.frame(
  Score = c(motifEnrichment_hnRNPC_WT_inKD$Score, motifEnrichment_hnRNPC_Mut_inKD$Score),
  Sample = factor(c(rep("hnRNPC WT in KD", nrow(motifEnrichment_hnRNPC_WT_inKD)),
                    rep("hnRNPC Mut in KD", nrow(motifEnrichment_hnRNPC_Mut_inKD))),
                  levels = c("hnRNPC WT in KD", "hnRNPC Mut in KD"))
)
ggplot(boxplot_data, aes(x = Sample, y = Score, fill = Sample)) +
  geom_boxplot(alpha = 0.7, outliers = F, notch = T) +
  scale_fill_manual(values = c("hnRNPC WT in KD" = "blue", "hnRNPC Mut in KD" = "red")) +
  labs(title = "Motif Enrichment Score: hnRNPC WT vs Mut in KD",
       x = NULL, y = "Motif Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.position = "none")



t.test(motifEnrichment_hnRNPC_WT$Score, motifEnrichment_hnRNPC_Mut$Score)
t.test(motifEnrichment_hnRNPC_WT_inKD$Score, motifEnrichment_hnRNPC_Mut_inKD$Score)
t.test(motifEnrichment_hnRNPC_WT_inRBM25_WT$Score, motifEnrichment_hnRNPC_WT_inRBM25_Mut$Score)
t.test(motifEnrichment_RBM25_WT$Score, motifEnrichment_RBM25_Mut$Score)


wilcox.test(motifEnrichment_RBM25_WT$Score, 
            motifEnrichment_RBM25_Mut$Score, 
            paired = TRUE, 
            conf.int = TRUE)


################################################################################
