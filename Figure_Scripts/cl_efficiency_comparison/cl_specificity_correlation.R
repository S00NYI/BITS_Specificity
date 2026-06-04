################################################################################
## Crosslinking Efficiency vs RBP Specificity Metrics:
## Written by Soon Yi
## Created: 2026-06-03
## Correlate IS/VS/CS/CVS against %CL (abundance-controlled check included)
################################################################################

library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

## Set some basic parameters:
################################################################################
baseDir   = '~/Repos/BITS_Specificity/Figure_Scripts/cl_efficiency_comparison/'
xlsx_file = paste0(baseDir, 'Supplemental_Table_S1.xlsx')
csv_file  = paste0(baseDir, 'RBP_metrics.csv')
sheet     = 'Supplemental Table S1'
################################################################################

## Load and process data:
################################################################################
# xlsx: real header on row 10 (9 note rows above); keep per-protein block A:W.
# Reference columns by position, not name: the sheet has duplicate header names
# (a second block at col AB) so readxl appends ...1/...3/etc. and clean names
# no longer exist. Col C = gene; L,M,N = (%TP)input R1-3; U = Average %CL.
# Abundance (mean (%TP)input) is kept only to control for it in the partial
# correlation below; it is not a correlation target.
xl = read_excel(xlsx_file, sheet = sheet, skip = 9)
xl = xl[, 1:23]
xl = tibble(RBP       = xl[[3]],
            avg_CL    = xl[[21]],
            abundance = rowMeans(cbind(xl[[12]], xl[[13]], xl[[14]]), na.rm = TRUE))

# csv: BOM + trailing empty columns.
csv = read.csv(csv_file, fileEncoding = 'UTF-8-BOM', check.names = FALSE)
csv = csv[, c('RBP', 'IS', 'VS', 'CS', 'CVS')]

# join.
df = inner_join(csv, xl, by = 'RBP')

missing = setdiff(csv$RBP, df$RBP)
if (length(missing) > 0) warning('No xlsx match for: ', paste(missing, collapse = ', '))
message('Matched ', nrow(df), ' RBPs')
################################################################################

## Compute correlations + p-values: metric vs %CL (Pearson + Spearman):
################################################################################
metrics = c('IS', 'VS', 'CS', 'CVS')

results = data.frame(metric = metrics, n = nrow(df),
                     pearson_r = NA_real_, pearson_p = NA_real_,
                     spearman_rho = NA_real_, spearman_p = NA_real_)
for (i in seq_along(metrics)) {
  pc = cor.test(df[[metrics[i]]], df$avg_CL, method = 'pearson')
  sc = suppressWarnings(cor.test(df[[metrics[i]]], df$avg_CL, method = 'spearman'))
  results$pearson_r[i]    = pc$estimate
  results$pearson_p[i]    = pc$p.value
  results$spearman_rho[i] = sc$estimate
  results$spearman_p[i]   = sc$p.value
}
cat('Correlation: metric vs %CL\n')
print(results, digits = 3)
################################################################################

## Partial correlation: metric vs %CL, controlling for abundance:
################################################################################
# Asks "does specificity relate to crosslinking EFFICIENCY independent of
# abundance?". %CL is already an efficiency (abundance divided out by definition),
# but %CL still trends with abundance across proteins, so we residualize both
# log10(%CL) and the metric on log10(abundance), then correlate the residuals.
partial = data.frame(metric = metrics, n = nrow(df), partial_r = NA_real_, partial_p = NA_real_)
lz = log10(df$abundance)
for (i in seq_along(metrics)) {
  res_CL = resid(lm(log10(df$avg_CL) ~ lz))
  res_M  = resid(lm(df[[metrics[i]]] ~ lz))
  ct     = cor.test(res_CL, res_M, method = 'pearson')
  partial$partial_r[i] = ct$estimate
  partial$partial_p[i] = ct$p.value
}
cat('\nPartial correlation: metric vs log10(%CL) | log10(abundance)\n')
print(partial, digits = 3)
################################################################################

## Scatterplots: metric vs %CL, log-log axes, linear Pearson r annotated:
################################################################################
# NOTE: axes are log10 for readability, but the annotated R is the Pearson
# correlation on LINEAR values (as requested), so the implied fit is not the
# straight line through these log-log points.
# Axis limits are in DATA units (not log units) and must be strictly positive
# (a 0 lower bound -> NaN under log10). Edit the limits = c(...) per block.

## IS vs %CL
r   = cor(df$avg_CL, df$IS, use = 'complete.obs')
pv  = cor.test(df$avg_CL, df$IS)$p.value
sc  = suppressWarnings(cor.test(df$avg_CL, df$IS, method = 'spearman'))
ggplot(df, aes(x = avg_CL, y = IS, label = RBP)) +
  geom_point(size = 2) +
  geom_text_repel(size = 3) +
  scale_x_continuous(transform = "log10", limits = c(0.1, 100)) +
  scale_y_continuous(transform = "log10", limits = c(0.1, 100)) +
  labs(x = "%CL", y = "IS", title = "IS vs %CL") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("R == ", format(r, digits = 2), "~~~p == ", format(pv, digits = 2)),
           parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("rho == ", format(sc$estimate, digits = 2), "~~~p == ", format(sc$p.value, digits = 2)),
           parse = TRUE, hjust = 1.1, vjust = 3.6, size = 5, color = "blue") +
  theme_bw() +
  theme(axis.text   = element_text(size = 14),
        axis.title  = element_text(size = 14, face = 'bold'),
        plot.title  = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14))

## VS vs %CL
r   = cor(df$avg_CL, df$VS, use = 'complete.obs')
pv  = cor.test(df$avg_CL, df$VS)$p.value
sc  = suppressWarnings(cor.test(df$avg_CL, df$VS, method = 'spearman'))
ggplot(df, aes(x = avg_CL, y = VS, label = RBP)) +
  geom_point(size = 2) +
  geom_text_repel(size = 3) +
  scale_x_continuous(transform = "log10", limits = c(0.1, 100)) +
  scale_y_continuous(transform = "log10", limits = c(0.01, 1)) +
  labs(x = "%CL", y = "VS", title = "VS vs %CL") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("R == ", format(r, digits = 2), "~~~p == ", format(pv, digits = 2)),
           parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("rho == ", format(sc$estimate, digits = 2), "~~~p == ", format(sc$p.value, digits = 2)),
           parse = TRUE, hjust = 1.1, vjust = 3.6, size = 5, color = "blue") +
  theme_bw() +
  theme(axis.text   = element_text(size = 14),
        axis.title  = element_text(size = 14, face = 'bold'),
        plot.title  = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14))

## CS vs %CL
r   = cor(df$avg_CL, df$CS, use = 'complete.obs')
pv  = cor.test(df$avg_CL, df$CS)$p.value
sc  = suppressWarnings(cor.test(df$avg_CL, df$CS, method = 'spearman'))
ggplot(df, aes(x = avg_CL, y = CS, label = RBP)) +
  geom_point(size = 2) +
  geom_text_repel(size = 3) +
  scale_x_continuous(transform = "log10", limits = c(0.1, 100)) +
  scale_y_continuous(transform = "log10", limits = c(0.1, 100)) +
  labs(x = "%CL", y = "CS", title = "CS vs %CL") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("R == ", format(r, digits = 2), "~~~p == ", format(pv, digits = 2)),
           parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("rho == ", format(sc$estimate, digits = 2), "~~~p == ", format(sc$p.value, digits = 2)),
           parse = TRUE, hjust = 1.1, vjust = 3.6, size = 5, color = "blue") +
  theme_bw() +
  theme(axis.text   = element_text(size = 14),
        axis.title  = element_text(size = 14, face = 'bold'),
        plot.title  = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14))

## CVS vs %CL
r   = cor(df$avg_CL, df$CVS, use = 'complete.obs')
pv  = cor.test(df$avg_CL, df$CVS)$p.value
sc  = suppressWarnings(cor.test(df$avg_CL, df$CVS, method = 'spearman'))
ggplot(df, aes(x = avg_CL, y = CVS, label = RBP)) +
  geom_point(size = 2) +
  geom_text_repel(size = 3) +
  scale_x_continuous(transform = "log10", limits = c(0.1, 100)) +
  scale_y_continuous(transform = "log10", limits = c(0.1, 1)) +
  labs(x = "%CL", y = "CVS", title = "CVS vs %CL") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("R == ", format(r, digits = 2), "~~~p == ", format(pv, digits = 2)),
           parse = TRUE, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("rho == ", format(sc$estimate, digits = 2), "~~~p == ", format(sc$p.value, digits = 2)),
           parse = TRUE, hjust = 1.1, vjust = 3.6, size = 5, color = "blue") +
  theme_bw() +
  theme(axis.text   = element_text(size = 14),
        axis.title  = element_text(size = 14, face = 'bold'),
        plot.title  = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14))
################################################################################

# write.csv(results, paste0(baseDir, 'correlation_results.csv'), row.names = FALSE)
# write.csv(df,      paste0(baseDir, 'merged_data.csv'),         row.names = FALSE)
