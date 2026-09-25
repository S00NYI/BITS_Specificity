################################################################################
## Task A: K562 vs HepG2 reproducibility of CS and CVS (Reviewer 2, comment 2)
## Genome Biology revision, 2026-09
##
## Uses the stored per-cell-line 5-mer enrichment vectors that the manuscript
## figures were built from. No enrichment is recomputed here.
##   Fig. 2 set: Dataset/Analysis/eCLIP_peak/output/20260317/ (Figure_2.R)
##   Fig. 6 set: Dataset/Analysis/eCLIP_all/output/20260316/  (Figure_6.R)
## CS and CVS use RBPSpecificity 0.99.0 returnIS / returnMS(output_type = "number"),
## the installed build whose function bodies match commit 8b121b3 (HEAD on 2026-03-17).
################################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(RBPSpecificity)
})

baseDir = '/Users/soonyi/Repos/BITS_Specificity/Dataset/Analysis/'
outDir = '/Users/soonyi/Repos/BITS_Specificity/REVISION_2026-09/'
K = 5
stopifnot(as.character(packageVersion('RBPSpecificity')) == '0.99.0')

## Helpers (same steps as Figure_2.R)
################################################################################
min_max_norm = function(x, a = 0, b = 1) {
  a + (x - min(x, na.rm = TRUE)) * (b - a) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

read_vec = function(path) {
  x = read_csv(path, col_names = TRUE, show_col_types = FALSE)
  x = x[, c('MOTIF', 'Score')]
  x$MOTIF = gsub('T', 'U', x$MOTIF)
  data.frame(x %>% arrange(MOTIF))
}

metrics = function(v) {
  top = v$MOTIF[which.max(v$Score)]
  c(top = top,
    CS = suppressMessages(returnIS(v, top)),
    CVS = suppressMessages(returnMS(v, top, output_type = 'number')))
}

pair_stats = function(x, y, label, set) {
  ok = is.finite(x) & is.finite(y)
  p = cor.test(x[ok], y[ok], method = 'pearson')
  s = suppressWarnings(cor.test(x[ok], y[ok], method = 'spearman', exact = FALSE))
  data.frame(record_type = 'summary', set = set,
             statistic = paste0(label, c('_pearson_r', '_spearman_rho')),
             estimate = c(unname(p$estimate), unname(s$estimate)),
             p_value = c(p$p.value, s$p.value),
             n = sum(ok))
}

dist_stats = function(r, label, set) {
  o = order(r$value)
  data.frame(record_type = 'summary', set = set,
             statistic = paste0(label, c('_median', '_Q1', '_Q3', '_IQR', '_min', '_max', '_lowest3')),
             estimate = c(median(r$value), quantile(r$value, 0.25), quantile(r$value, 0.75),
                          IQR(r$value), min(r$value), max(r$value), NA),
             p_value = NA, n = nrow(r),
             note = c(rep(NA, 6), paste(sprintf('%s (%.3f)', r$RBP[o][1:3], r$value[o][1:3]), collapse = '; ')))
}

per_rbp = function(set, rbps, path_fun, acc_fun) {
  rows = lapply(rbps, function(RBP) {
    k = read_vec(path_fun('K562', RBP))
    h = read_vec(path_fun('HepG2', RBP))
    stopifnot(identical(k$MOTIF, h$MOTIF), nrow(k) == 4^K)
    mk = metrics(k)
    mh = metrics(h)
    merged = data.frame(MOTIF = k$MOTIF, Score = min_max_norm((k$Score + h$Score) / 2))
    mm = metrics(merged)
    cs = as.numeric(c(mk['CS'], mh['CS']))
    data.frame(record_type = 'per_RBP', set = set, RBP = RBP,
               accession_K562 = acc_fun('K562', RBP), accession_HepG2 = acc_fun('HepG2', RBP),
               top5mer_K562 = mk['top'], top5mer_HepG2 = mh['top'],
               top5mer_identical = mk['top'] == mh['top'],
               CS_K562 = cs[1], CS_HepG2 = cs[2],
               CS_fold_diff = max(cs) / min(cs),
               CS_gt_2fold = max(cs) / min(cs) > 2,
               CVS_K562 = as.numeric(mk['CVS']), CVS_HepG2 = as.numeric(mh['CVS']),
               vector_pearson_r = cor(k$Score, h$Score, method = 'pearson'),
               vector_spearman_rho = cor(k$Score, h$Score, method = 'spearman'),
               CS_mean_of_lines = mean(cs),
               CVS_mean_of_lines = mean(as.numeric(c(mk['CVS'], mh['CVS']))),
               CS_merged_vector = as.numeric(mm['CS']),
               CVS_merged_vector = as.numeric(mm['CVS']),
               top5mer_merged_vector = mm['top'],
               row.names = NULL)
  })
  bind_rows(rows)
}

summarise_set = function(tab, set) {
  bind_rows(
    pair_stats(tab$CS_K562, tab$CS_HepG2, 'CS', set),
    pair_stats(log2(tab$CS_K562), log2(tab$CS_HepG2), 'log2CS', set),
    pair_stats(tab$CVS_K562, tab$CVS_HepG2, 'CVS', set),
    dist_stats(data.frame(RBP = tab$RBP, value = tab$vector_pearson_r), 'vector_pearson_r', set),
    dist_stats(data.frame(RBP = tab$RBP, value = tab$vector_spearman_rho), 'vector_spearman_rho', set),
    data.frame(record_type = 'summary', set = set, statistic = 'n_CS_gt_2fold',
               estimate = sum(tab$CS_gt_2fold), n = nrow(tab),
               note = paste(tab$RBP[tab$CS_gt_2fold], collapse = '; ')),
    data.frame(record_type = 'summary', set = set, statistic = 'n_top5mer_identical',
               estimate = sum(tab$top5mer_identical), n = nrow(tab))
  )
}

## Fig. 2 set (41 datasets, sample_table_eCLIP.csv; enrichment run 20260317)
################################################################################
st = read.csv(paste0(baseDir, 'sample_table_eCLIP.csv'), fileEncoding = 'UTF-8-BOM')
fig2_both = sort(intersect(st$RBP[st$Cell == 'K562'], st$RBP[st$Cell == 'HepG2']))
fig2_path = function(cell, RBP) paste0(baseDir, 'eCLIP_peak/output/20260317/', cell,
                                       '_eCLIP_NormalizedEnrichment_25ntExt_', RBP, '_5mer.csv')
fig2_acc = function(cell, RBP) st$eCLIP[st$RBP == RBP & st$Cell == cell]
fig2 = per_rbp('Fig2_set', fig2_both, fig2_path, fig2_acc)

# Validation 1: merged-vector CS/CVS must equal the published Fig. 2 summary table.
fig2_pub = read_csv(paste0(baseDir, 'eCLIP_peak/output/20260317_5mer_RBNS_eCLIP_Analysis_25ntExt.csv'),
                    show_col_types = FALSE)
fig2 = fig2 %>%
  left_join(fig2_pub %>% select(RBP, published_CS = eCLIP_Specificity, published_CVS = eCLIP_Sensitivity),
            by = 'RBP')
fig2$published_merge_rule = ifelse(abs(fig2$CS_merged_vector - fig2$published_CS) < 1e-9, 'merged_vector',
                            ifelse(abs(fig2$CS_K562 - fig2$published_CS) < 1e-9, 'K562_only',
                            ifelse(abs(fig2$CS_HepG2 - fig2$published_CS) < 1e-9, 'HepG2_only', 'NO_MATCH')))

# Validation 2: per-line CS must equal the stored ISperMotif file at the top motif.
for (i in seq_len(nrow(fig2))) {
  for (cell in c('K562', 'HepG2')) {
    is_file = read_csv(paste0(baseDir, 'eCLIP_peak/output/20260317/', cell, '_eCLIP_ISperMotif_25ntExt_',
                              fig2$RBP[i], '_5mer.csv'), show_col_types = FALSE)
    top = gsub('U', 'T', fig2[i, paste0('top5mer_', cell)])
    stopifnot(abs(is_file[[3]][is_file[[2]] == top] - fig2[i, paste0('CS_', cell)]) < 1e-9)
  }
}

# Rerun noise floor: the same peak files were run again on 20260316 (Fig. 6 run).
for (cell in c('K562', 'HepG2')) {
  fig2[[paste0('rerun_vector_pearson_r_', cell)]] = sapply(fig2$RBP, function(RBP) {
    cor(read_vec(fig2_path(cell, RBP))$Score,
        read_vec(paste0(baseDir, 'eCLIP_all/output/20260316/', cell, '_eCLIP_Enrichment_', RBP, '_5mer.csv'))$Score)
  })
  fig2[[paste0('rerun_CS_', cell)]] = sapply(fig2$RBP, function(RBP) {
    as.numeric(metrics(read_vec(paste0(baseDir, 'eCLIP_all/output/20260316/', cell,
                                       '_eCLIP_Enrichment_', RBP, '_5mer.csv')))['CS'])
  })
}

## Fig. 6 set (eCLIP_all/eCLIP_list.csv; enrichment run 20260316)
################################################################################
el = read.csv(paste0(baseDir, 'eCLIP_all/eCLIP_list.csv'), fileEncoding = 'UTF-8-BOM')
fig6_both = sort(intersect(el$RBP[el$CELL == 'K562'], el$RBP[el$CELL == 'HepG2']))
fig6_path = function(cell, RBP) paste0(baseDir, 'eCLIP_all/output/20260316/', cell,
                                       '_eCLIP_Enrichment_', RBP, '_5mer.csv')
fig6_acc = function(cell, RBP) paste(el$ACCESSION[el$RBP == RBP & el$CELL == cell], collapse = '|')
fig6 = per_rbp('Fig6_set', fig6_both, fig6_path, fig6_acc)

# Validation 3: per-line CS/CVS must equal eCLIP_5mer_IS_VS_summary.csv; identify how IS_avg was built.
fig6_pub = read_csv(paste0(baseDir, 'eCLIP_all/output/20260316/eCLIP_5mer_IS_VS_summary.csv'),
                    show_col_types = FALSE)
fig6 = fig6 %>% left_join(fig6_pub, by = 'RBP')
fig6$published_CS = fig6$IS_avg
fig6$published_CVS = fig6$VS_avg
fig6$published_merge_rule = ifelse(abs(fig6$CS_mean_of_lines - fig6$IS_avg) < 1e-9, 'mean_of_line_values',
                            ifelse(abs(fig6$CS_merged_vector - fig6$IS_avg) < 1e-9, 'merged_vector', 'NO_MATCH'))
stopifnot(all(abs(fig6$CS_K562 - fig6$IS_K562) < 1e-9), all(abs(fig6$CS_HepG2 - fig6$IS_HepG2) < 1e-9),
          all(abs(fig6$CVS_K562 - fig6$VS_K562) < 1e-9), all(abs(fig6$CVS_HepG2 - fig6$VS_HepG2) < 1e-9))
fig6 = fig6 %>% select(-IS_K562, -IS_HepG2, -IS_avg, -VS_K562, -VS_HepG2, -VS_avg)

## Write one CSV
################################################################################
res = bind_rows(fig2, fig6, summarise_set(fig2, 'Fig2_set'), summarise_set(fig6, 'Fig6_set'))
write_csv(res, paste0(outDir, '4_bootstrap_cellline_reproducibility.csv'), na = '')

## Console report
################################################################################
options(width = 200)
cat('\nFig2 merge rule check:\n'); print(table(fig2$published_merge_rule)); print(fig2$RBP[fig2$published_merge_rule != 'merged_vector'])
cat('\nFig6 merge rule check:\n'); print(table(fig6$published_merge_rule))
cat('\nRerun noise (Fig2 datasets, 20260316 vs 20260317 vectors):\n')
print(summary(c(fig2$rerun_vector_pearson_r_K562, fig2$rerun_vector_pearson_r_HepG2)))
print(summary(abs(log2(c(fig2$rerun_CS_K562 / fig2$CS_K562, fig2$rerun_CS_HepG2 / fig2$CS_HepG2)))))
for (s in c('Fig2_set', 'Fig6_set')) {
  cat('\n====', s, '\n')
  print(res %>% filter(record_type == 'summary', set == s) %>% select(statistic, estimate, p_value, n, note))
}
cat('\nFig2 per-RBP:\n')
print(fig2 %>% select(RBP, top5mer_K562, top5mer_HepG2, CS_K562, CS_HepG2, CS_fold_diff, CVS_K562, CVS_HepG2,
                      vector_pearson_r, vector_spearman_rho, published_CS, published_merge_rule) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))))
cat('\nFig6 per-RBP (CS > 2-fold):\n')
print(fig6 %>% filter(CS_gt_2fold) %>% select(RBP, top5mer_K562, top5mer_HepG2, CS_K562, CS_HepG2, CS_fold_diff,
                                              CVS_K562, CVS_HepG2, vector_pearson_r) %>%
        arrange(desc(CS_fold_diff)) %>% mutate(across(where(is.numeric), ~ round(.x, 3))))
