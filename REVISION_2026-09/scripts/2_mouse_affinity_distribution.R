################################################################################
## Task J: affinity (motif-enrichment) distribution for the mouse eIF4G2 HITS-CLIP
## 5-mer vector, with the top motif labelled.
##
## Matches the affinity-distribution idiom already used in the manuscript:
##   - Figure_1.R panels B/C  (geom_histogram of Score, theme_bw, size-14 axes)
##   - Figure4_5_CLIP/20260204_CLIP_Peak_Motif_Analysis.R  (x = "Motif Enrichment",
##     binwidth 0.005) -- same units the eCLIP panels use.
##
## Input : REVISION_2026-09/2_mouse_5mer_enrichment.csv  (1024 5-mers, U alphabet)
##         Regenerate it with scripts/2_mouse_5mer_enrichment.R if absent
##         (needs BSgenome.Mmusculus.UCSC.mm39; identical procedure to 2_mouse_CS_CVS).
## Output: a ggplot object `affinityDistribution_mouse`. No file is written here --
##         save it however you save your other panels (kept plot-only per your workflow).
################################################################################

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

REV_DIR = '/Users/soonyi/Repos/BITS_Specificity/REVISION_2026-09/'
enr = read_csv(paste0(REV_DIR, '2_mouse_5mer_enrichment.csv'), show_col_types = FALSE)
# columns: MOTIF (U alphabet), Score (log-scaled motif enrichment, the same quantity
# plotted for the eCLIP panels). CS = max(Score)/median(Score).

## --- affinity-distribution plot, top motif labelled ---------------------------
# Returns a ggplot; call it for any enrichment table with MOTIF + Score.
plotAffinityDistribution = function(data, title, binwidth = 0.005,
                                    bar_fill = 'black', label_col = 'firebrick', CVS = NA) {
  data = as.data.frame(data[, c('MOTIF', 'Score')])
  top_i     = which.max(data$Score)
  top_motif = data$MOTIF[top_i]
  top_score = data$Score[top_i]
  CS  = round(max(data$Score) / median(data$Score), 2)
  subtitle = paste0('top 5-mer ', top_motif, '  •  CS = ', CS,
                    if (!is.na(CVS)) paste0('  •  CVS = ', round(CVS, 2)) else '')

  ggplot(data, aes(x = Score)) +
    geom_histogram(binwidth = binwidth, fill = bar_fill, alpha = 1.0) +
    # mark + name the top 5-mer at its enrichment value
    geom_vline(xintercept = top_score, colour = label_col, linetype = 'dashed', linewidth = 0.5) +
    annotate('text', x = top_score, y = Inf, label = top_motif,
             colour = label_col, fontface = 'bold', size = 5, hjust = 1.15, vjust = 1.8) +
    labs(title = title, subtitle = subtitle, x = 'Motif Enrichment', y = 'Frequency') +
    theme_bw() +
    theme(plot.title    = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.text     = element_text(size = 14),
          axis.title    = element_text(size = 14, face = 'bold'),
          legend.text   = element_text(size = 14))
}

# CVS = 0.6404 comes from 2_mouse_CS_CVS (needs the per-variant computation, not derivable
# from the enrichment vector alone); pass it in so the subtitle can show it.
affinityDistribution_mouse = plotAffinityDistribution(
  enr, title = 'Mouse eIF4G2 HITS-CLIP', CVS = 0.6404)

dir.create(file.path(REV_DIR, 'figures'), showWarnings = FALSE)
ggsave(file.path(REV_DIR, 'figures', '2_mouse_affinity_distribution.pdf'),
       affinityDistribution_mouse, width = 5.2, height = 4)

## sanity echo
cat('top 5-mer :', enr$MOTIF[which.max(enr$Score)], '\n')
cat('CS        :', round(max(enr$Score) / median(enr$Score), 3), '\n')
