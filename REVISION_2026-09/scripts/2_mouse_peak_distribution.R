################################################################################
## Task K: peak (region) distribution for the mouse eIF4G2 HITS-CLIP peaks,
## shown next to the human EIF4G2 K562 eCLIP peaks for comparison.
##
## Matches the peak-distribution idiom already used in the manuscript:
##   Figure4_5_CLIP/20260204_CLIP_Peak_Analysis.R  ->  plotStackedBar()
##   (geom_bar position='stack', scale_fill_brewer "Set3", theme_bw, size-14 axes).
##
## Input : REVISION_2026-09/2_mouse_region_distribution.csv  (already computed in Task E)
## Output: ggplot objects `peakDistribution_peaks` (primary) and
##         `peakDistribution_tags`. No file written -- ggsave() as you prefer.
################################################################################

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

REV_DIR = '/Users/soonyi/Repos/BITS_Specificity/REVISION_2026-09/'
reg = read_csv(paste0(REV_DIR, '2_mouse_region_distribution.csv'), show_col_types = FALSE)
# columns: dataset, region, n, fraction, n_peaks_total

## --- cosmetic mappings (region order + human-readable region / source labels) --
region_levels = c('5UTR', 'CDS', '3UTR', 'intron', 'upstream_10kb', 'downstream_10kb', 'intergenic')
region_labels = c("5'UTR", 'CDS', "3'UTR", 'intron', 'upstream 10kb', 'downstream 10kb', 'intergenic')

source_labels = c(
  human_EIF4G2_K562_eCLIP        = 'eCLIP\n(human, K562)',
  mouse_eIF4G2_HITSCLIP_resting  = 'HITS-CLIP\n(mouse)'
)

## --- stacked-bar plotter (fraction, stacked to 1), matching plotStackedBar() ----
plotPeakDistribution = function(df, sources, title) {
  d = df %>%
    filter(dataset %in% names(sources)) %>%
    mutate(Annotation = factor(region, levels = region_levels, labels = region_labels),
           Source     = factor(dataset, levels = names(sources), labels = sources))
  ggplot(d, aes(fill = Annotation, y = fraction, x = Source)) +
    geom_bar(position = 'stack', stat = 'identity') +
    ggtitle(title) +
    scale_fill_brewer(palette = 'Set3') +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    coord_cartesian(ylim = c(0, 1)) +   # zoom, don't clip (fractions sum to ~1)
    labs(x = NULL, y = 'Fraction of peaks') +
    theme_bw() +
    theme(plot.title  = element_text(hjust = 0.5),
          axis.text   = element_text(size = 14),
          axis.title  = element_text(size = 14, face = 'bold'),
          legend.text = element_text(size = 14),
          legend.title = element_blank())
}

## primary: peak-level distribution (mouse resting vs human eCLIP)
peakDistribution_peaks = plotPeakDistribution(
  reg, source_labels, title = 'eIF4G2 peak region distribution')
peakDistribution_peaks

## secondary: tag-level distribution (mouse), if you want the read-level view too.
## The tag rows carry the "_TAGS" suffix in the same CSV.
tag_sources = c(mouse_eIF4G2_HITSCLIP_resting_TAGS = 'eIF4G2 HITS-CLIP\n(mouse, tags)')
peakDistribution_tags = reg %>%
  filter(dataset %in% names(tag_sources)) %>%
  mutate(Annotation = factor(region, levels = region_levels, labels = region_labels),
         Source     = factor(dataset, levels = names(tag_sources), labels = tag_sources)) %>%
  ggplot(aes(fill = Annotation, y = fraction, x = Source)) +
  geom_bar(position = 'stack', stat = 'identity') +
  ggtitle('eIF4G2 HITS-CLIP tag region distribution (mouse)') +
  scale_fill_brewer(palette = 'Set3') +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +   # zoom, don't clip (fractions sum to ~1)
  labs(x = NULL, y = 'Fraction of tags') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14), legend.title = element_blank())

dir.create(file.path(REV_DIR, 'figures'), showWarnings = FALSE)
ggsave(file.path(REV_DIR, 'figures', '2_mouse_peak_distribution.pdf'), peakDistribution_peaks, width = 5.5, height = 4.5)
ggsave(file.path(REV_DIR, 'figures', '2_mouse_tag_distribution.pdf'),  peakDistribution_tags,  width = 4.5, height = 4.5)

## sanity echo
cat('peak-level fractions:\n')
print(reg %>% filter(dataset == 'mouse_eIF4G2_HITSCLIP_resting') %>% select(region, fraction))
