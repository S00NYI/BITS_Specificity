################################################################################
## U2AF2 Cofactor Competition on PTBP2 Transcript
## Written by Soon Yi
## Created: 2026-02-22
## Last Edited: 2026-02-22
## Figure 3G
################################################################################

library(data.table)
library(tidyr)
library(RBPEqBind)

library(readr)
library(dplyr)
library(ggplot2)

## Set up paths and parameters:
################################################################################
simDir = "~/Repos/BITS_Specificity/Dataset/Analysis/RBPEqBind_Simulation"
modelFile = file.path(simDir, "DATA_PROCESSED", "rnacompete_affinity_scores.csv")
fastaFile = file.path(simDir, "DATA_PROCESSED", "FASTA", "PTBP2.fa")
bgBedDir = file.path(simDir, "RESULTS_ANALYSIS", "BEDGRAPH_PROCESSED")
bindingSitesFile = file.path(simDir, "DATA_PROCESSED", "binding_sites.csv")

## Transcript/ROI info:
ptbp2_genomic_start = 97269727
roi_genomic_start = 97271580
roi_genomic_end = 97271780
roi_tx_start = roi_genomic_start - ptbp2_genomic_start + 1 # 1854
roi_tx_end = roi_genomic_end - ptbp2_genomic_start + 1 # 2054
xlim = c(roi_tx_start, roi_tx_end)

## Simulation parameters (validated from prior analysis):
u2af2_kd = 50 # nM
cofac_kd = 50 # nM
u2af2_conc = 500 # nM
hnrnpc_conc = 200 # nM
ptbp1_conc = 200 # nM
rna_conc = 6.75 # nM
k = 7

## Transcript name (from FASTA header):
tx_name = "PTBP2|chr1:97269727-97272451|+|len=2725"
################################################################################

## Load RNACompete Model
################################################################################
model_raw = read_csv(modelFile, show_col_types = FALSE)
colnames(model_raw)[1] = "Motif"
temp_model = file.path(simDir, "DATA_PROCESSED", "temp_model.csv")
write_csv(model_raw, temp_model)
raw_models = loadModel(temp_model)
################################################################################

## Condition 1: U2AF2 Alone
################################################################################
max_aff_1 = c(U2AF2 = 1 / u2af2_kd)
rbp_models_1 = setModel(raw_models, max_affinity = max_aff_1, min_affinity = 0.00001)

simResults_u2af2 = simulateBindingF(
    fasta_file = fastaFile,
    rbp_models = rbp_models_1["U2AF2"],
    protein_concs = c(U2AF2 = u2af2_conc),
    rna_conc = rna_conc,
    k = k
)
################################################################################

## Condition 2: U2AF2 + HNRNPC
################################################################################
max_aff_2 = c(U2AF2 = 1 / u2af2_kd, HNRNPC = 1 / cofac_kd)
rbp_models_2 = setModel(raw_models, max_affinity = max_aff_2, min_affinity = 0.00001)

simResults_hnrnpc = simulateBindingF(
    fasta_file = fastaFile,
    rbp_models = rbp_models_2[c("U2AF2", "HNRNPC")],
    protein_concs = c(U2AF2 = u2af2_conc, HNRNPC = hnrnpc_conc),
    rna_conc = rna_conc,
    k = k
)
################################################################################

## Condition 3: U2AF2 + PTBP1
################################################################################
max_aff_3 = c(U2AF2 = 1 / u2af2_kd, PTBP1 = 1 / cofac_kd)
rbp_models_3 = setModel(raw_models, max_affinity = max_aff_3, min_affinity = 0.00001)

simResults_ptbp1 = simulateBindingF(
    fasta_file = fastaFile,
    rbp_models = rbp_models_3[c("U2AF2", "PTBP1")],
    protein_concs = c(U2AF2 = u2af2_conc, PTBP1 = ptbp1_conc),
    rna_conc = rna_conc,
    k = k
)
################################################################################

## Simulation Heatmaps (U2AF2 density at ROI)
################################################################################
zlim = c(0, 0.00075)

plotHeatmap(simResults_u2af2,
    transcript = tx_name,
    rbps = "U2AF2",
    metric = "density",
    xlim = xlim,
    zlim = zlim
)

plotHeatmap(simResults_hnrnpc,
    transcript = tx_name,
    rbps = "U2AF2",
    metric = "density",
    xlim = xlim,
    zlim = zlim
)

plotHeatmap(simResults_ptbp1,
    transcript = tx_name,
    rbps = "U2AF2",
    metric = "density",
    xlim = xlim,
    zlim = zlim
)
################################################################################

## Load Experimental iCLIP Bedgraphs
################################################################################
load_bedgraph_roi = function(bg_file, target_chr, roi_start, roi_end) {
    d = read_tsv(bg_file, col_names = c("chr", "start", "end", "val"), show_col_types = FALSE)
    d_sub = d %>% filter(chr == target_chr, end > roi_start, start < roi_end)

    if (nrow(d_sub) == 0) {
        return(data.frame(pos = integer(), val = numeric()))
    }

    track = data.frame()
    for (i in seq_len(nrow(d_sub))) {
        positions = seq(max(d_sub$start[i], roi_start), min(d_sub$end[i], roi_end))
        track = rbind(track, data.frame(pos = positions, val = d_sub$val[i]))
    }
    return(track)
}

bg_u2af2 = load_bedgraph_roi(
    file.path(bgBedDir, "bg_u2af2_500nM.bedgraph"),
    "chr1", roi_genomic_start, roi_genomic_end
)
bg_hnrnpc = load_bedgraph_roi(
    file.path(bgBedDir, "bg_u2af2_500nM_hnrnpc_200nM.bedgraph"),
    "chr1", roi_genomic_start, roi_genomic_end
)
bg_ptbp1 = load_bedgraph_roi(
    file.path(bgBedDir, "bg_u2af2_500nM_ptbp1_200nM.bedgraph"),
    "chr1", roi_genomic_start, roi_genomic_end
)

bg_u2af2$condition = "U2AF2 only"
bg_hnrnpc$condition = "U2AF2 + HNRNPC"
bg_ptbp1$condition = "U2AF2 + PTBP1"

bg_combined = rbind(bg_u2af2, bg_hnrnpc, bg_ptbp1)
bg_combined$condition = factor(bg_combined$condition,
    levels = c("U2AF2 only", "U2AF2 + HNRNPC", "U2AF2 + PTBP1")
)
################################################################################

## Plot Experimental iCLIP Signal at ROI
################################################################################
## Load binding sites for annotation
binding_sites = read_csv(bindingSitesFile, show_col_types = FALSE) %>%
    filter(transcript == "PTBP2") %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    filter(end > roi_genomic_start, start < roi_genomic_end)

p_iclip = ggplot(bg_combined, aes(x = pos, y = val, color = condition)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~condition, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c(
        "U2AF2 only" = "#34495E",
        "U2AF2 + HNRNPC" = "#E91E63",
        "U2AF2 + PTBP1" = "#F39C12"
    )) +
    geom_rect(
        data = binding_sites,
        aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE, fill = "gold", alpha = 0.2
    ) +
    theme_minimal() +
    theme(
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "none"
    ) +
    labs(
        title = "Experimental iCLIP Signal: PTBP2 ROI",
        subtitle = "Gold regions = Experimentally validated binding sites",
        x = "Genomic Position (chr1)",
        y = "Signal (Normalized)"
    ) +
    scale_x_continuous(labels = scales::comma) + 
  ylim(0, 0.000025)
  # ylim(0, 0.0025)

print(p_iclip)
################################################################################


## Correlation: iCLIP vs Simulation at ROI
################################################################################

## Helper: extract per-position vectors aligned to the ROI
get_sim_density = function(simResults, tx_name, roi_start, roi_end) {
  d = simResults[transcript == tx_name & pos >= roi_start & pos <= roi_end]
  setkey(d, pos)
  return(d[, .(pos, density = U2AF2_density)])
}

get_bg_signal = function(bg_track, ptbp2_genomic_start) {
  bg_track$tx_pos = bg_track$pos - ptbp2_genomic_start + 1
  data.frame(pos = bg_track$tx_pos, signal = bg_track$val)
}


## Build aligned data for each condition
conditions = list(
  list(name = "U2AF2 only",       sim = simResults_u2af2,  bg = bg_u2af2),
  list(name = "U2AF2 + HNRNPC",   sim = simResults_hnrnpc, bg = bg_hnrnpc),
  list(name = "U2AF2 + PTBP1",    sim = simResults_ptbp1,  bg = bg_ptbp1)
)

cor_results = data.frame()
for (cond in conditions) {
  sim_d = get_sim_density(cond$sim, tx_name, roi_tx_start, roi_tx_end)
  bg_d = get_bg_signal(cond$bg, ptbp2_genomic_start)
  
  # Merge on transcript position
  merged = merge(sim_d, bg_d, by = "pos", all = FALSE)
  
  if (nrow(merged) > 10) {
    r_pearson  = cor(merged$density, merged$signal, method = "pearson")
    r_spearman = cor(merged$density, merged$signal, method = "spearman")
  } else {
    r_pearson  = NA
    r_spearman = NA
  }
  
  cor_results = rbind(cor_results, data.frame(
    Condition = cond$name,
    N_positions = nrow(merged),
    Pearson = round(r_pearson, 4),
    Spearman = round(r_spearman, 4)
  ))
}

print(cor_results)
################################################################################
