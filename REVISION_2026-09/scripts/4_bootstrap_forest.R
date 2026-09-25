## Forest-plot strip of bootstrap 95% CIs for the reported correlations.
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
revDir = '/Users/soonyi/Repos/BITS_Specificity/REVISION_2026-09/'
d = read_csv(paste0(revDir, '4_bootstrap_all_CI.csv'), show_col_types = FALSE)
d = d[!grepl(' All', d$panel), ]                     # drop the pooled (All-18) rows

fig_of = function(p) sub(' .*', '', p)
d = d %>% mutate(
  figure = recode(fig_of(panel), Fig1D='Fig 1D', Fig2C='Fig 2C', Fig2D='Fig 2D',
                  Fig6A='Fig 6A', FigS6A='Fig S6A'),
  figure = factor(figure, levels = c('Fig 1D','Fig 2C','Fig 2D','Fig 6A','Fig S6A')),
  row = dplyr::case_when(
    grepl('Low',  panel) ~ 'Low context',
    grepl('High', panel) ~ 'High context',
    grepl(' All', panel) ~ 'All (pooled)',
    grepl('Fig1D', panel) ~ 'IS vs VS',
    grepl('Fig6A', panel) ~ 'CS K562 vs HepG2',
    grepl('FigS6A', panel) ~ 'CVS K562 vs HepG2'),
  excl0 = boot_lo > 0 | boot_hi < 0,
  lab = sprintf('%.2f [%.2f, %.2f]', r, boot_lo, boot_hi))
d$rowf = factor(d$panel, levels = rev(d$panel))      # unique key; relabel below

p = ggplot(d, aes(r, rowf, colour = excl0)) +
  scale_y_discrete(labels = setNames(d$row, d$panel)) +
  geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey55') +
  geom_errorbarh(aes(xmin = boot_lo, xmax = boot_hi), height = 0.25, linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(x = 1.02, label = lab), hjust = 0, size = 3, colour = 'grey25') +
  facet_grid(figure ~ ., scales = 'free_y', space = 'free_y', switch = 'y') +
  scale_colour_manual(values = c('TRUE' = '#F39C12', 'FALSE' = '#8492A6'),
                      labels = c('TRUE' = 'CI excludes 0', 'FALSE' = 'CI crosses 0'), name = NULL) +
  scale_x_continuous(limits = c(-1, 1.75), breaks = seq(-1, 1, 0.5)) +
  labs(title = 'Bootstrap 95% CIs (percentile) for the reported correlations',
       x = 'Pearson r', y = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        strip.text.y.left = element_text(angle = 0, face = 'bold'), strip.placement = 'outside',
        legend.position = 'top', axis.text = element_text(size = 10),
        panel.grid.major.y = element_blank())
dir.create(paste0(revDir, 'figures'), showWarnings = FALSE)
ggsave(paste0(revDir, 'figures/4_bootstrap_forest.pdf'), p, width = 8, height = 4.4)
cat('rows:', nrow(d), '\n')
