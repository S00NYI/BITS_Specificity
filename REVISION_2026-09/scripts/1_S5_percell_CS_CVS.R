################################################################################
## Task B2: augment the Table S1 audit with per-cell-line and merged CS/CVS.
## (Reviewer 2 -- transparency on which cell line each value comes from.)
##
## The Task B CSV only carried the single "reported" CS (one cell line OR the merge).
## This adds, for every eCLIP RBP:
##   CS_K562, CVS_K562        (blank if no K562 dataset)
##   CS_HepG2, CVS_HepG2      (blank if no HepG2 dataset)
##   CS_merged, CVS_merged    (the value the paper reports: merged-vector when both
##                             cell lines exist, else the single cell line)
##   merge_status             merged | K562_only | HepG2_only
## RBNS-only RBPs (no eCLIP) get blank CS/CVS.
##
## Same procedure as Task A / Figure_2.R: read the stored per-cell 5-mer enrichment
## vectors, RBPSpecificity 0.99.0 returnIS/returnMS; merged = min-max( mean(K562,HepG2) ).
## Every value is validated against eCLIP_5mer_IS_VS_summary.csv before writing.
################################################################################

suppressPackageStartupMessages({ library(readr); library(dplyr); library(RBPSpecificity) })

baseDir = '/Users/soonyi/Repos/BITS_Specificity/Dataset/Analysis/'
revDir  = '/Users/soonyi/Repos/BITS_Specificity/REVISION_2026-09/'
K = 5
stopifnot(as.character(packageVersion('RBPSpecificity')) == '0.99.0')

## --- helpers (identical to Task A) --------------------------------------------
min_max_norm = function(x, a = 0, b = 1)
  a + (x - min(x, na.rm = TRUE)) * (b - a) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

read_vec = function(path) {
  x = read_csv(path, col_names = TRUE, show_col_types = FALSE)[, c('MOTIF', 'Score')]
  x$MOTIF = gsub('T', 'U', x$MOTIF)
  data.frame(x %>% arrange(MOTIF))
}
metrics = function(v) {
  top = v$MOTIF[which.max(v$Score)]
  list(top = top,
       CS  = as.numeric(suppressMessages(returnIS(v, top))),
       CVS = as.numeric(suppressMessages(returnMS(v, top, output_type = 'number'))))
}
enr_path = function(cell, RBP)
  paste0(baseDir, 'eCLIP_all/output/20260316/', cell, '_eCLIP_Enrichment_', RBP, '_5mer.csv')

## --- which RBPs have which cell lines -----------------------------------------
el = read.csv(paste0(baseDir, 'eCLIP_all/eCLIP_list.csv'), fileEncoding = 'UTF-8-BOM')
has = function(RBP, cell) RBP %in% el$RBP[el$CELL == cell]
rbps = sort(unique(el$RBP))

## --- compute per-cell + merged CS/CVS for every eCLIP RBP ---------------------
rows = lapply(rbps, function(RBP) {
  hasK = has(RBP, 'K562'); hasH = has(RBP, 'HepG2')
  mk = if (hasK) metrics(read_vec(enr_path('K562',  RBP))) else NULL
  mh = if (hasH) metrics(read_vec(enr_path('HepG2', RBP))) else NULL

  if (hasK && hasH) {
    k = read_vec(enr_path('K562', RBP)); h = read_vec(enr_path('HepG2', RBP))
    stopifnot(identical(k$MOTIF, h$MOTIF))
    mm = metrics(data.frame(MOTIF = k$MOTIF, Score = min_max_norm((k$Score + h$Score) / 2)))
    status = 'merged'
  } else if (hasK) { mm = mk; status = 'K562_only'
  } else            { mm = mh; status = 'HepG2_only' }

  data.frame(
    RBP = RBP, merge_status = status,
    CS_K562   = if (hasK) mk$CS  else NA, CVS_K562  = if (hasK) mk$CVS else NA,
    CS_HepG2  = if (hasH) mh$CS  else NA, CVS_HepG2 = if (hasH) mh$CVS else NA,
    CS_merged = mm$CS, CVS_merged = mm$CVS,
    top5mer_K562   = if (hasK) mk$top else NA,
    top5mer_HepG2  = if (hasH) mh$top else NA,
    top5mer_merged = mm$top, row.names = NULL)
})
cs = bind_rows(rows)

## --- validate against the published summary -----------------------------------
pub = read_csv(paste0(baseDir, 'eCLIP_all/output/20260316/eCLIP_5mer_IS_VS_summary.csv'),
               show_col_types = FALSE)
chk = cs %>% left_join(pub, by = 'RBP')
tol = 1e-6
ok_k = with(chk, is.na(CS_K562)  | abs(CS_K562  - IS_K562)  < tol)
ok_h = with(chk, is.na(CS_HepG2) | abs(CS_HepG2 - IS_HepG2) < tol)
ok_vk = with(chk, is.na(CVS_K562)  | abs(CVS_K562  - VS_K562)  < tol)
ok_vh = with(chk, is.na(CVS_HepG2) | abs(CVS_HepG2 - VS_HepG2) < tol)
merged_rows = chk %>% filter(merge_status == 'merged')
ok_a  = abs(merged_rows$CS_merged  - merged_rows$IS_avg) < tol
ok_va = abs(merged_rows$CVS_merged - merged_rows$VS_avg) < tol
cat('validation vs eCLIP_5mer_IS_VS_summary.csv:\n')
cat('  CS_K562 == IS_K562  :', all(ok_k),  '(', sum(!ok_k),  'mismatches )\n')
cat('  CS_HepG2 == IS_HepG2:', all(ok_h),  '(', sum(!ok_h),  'mismatches )\n')
cat('  CVS_K562 == VS_K562 :', all(ok_vk), '(', sum(!ok_vk), 'mismatches )\n')
cat('  CVS_HepG2 == VS_HepG2:', all(ok_vh),'(', sum(!ok_vh), 'mismatches )\n')
cat('  CS_merged == IS_avg  (both-cell RBPs):', all(ok_a),  '(', sum(!ok_a),  'mismatches )\n')
cat('  CVS_merged == VS_avg (both-cell RBPs):', all(ok_va), '(', sum(!ok_va), 'mismatches )\n')
stopifnot(all(ok_k), all(ok_h), all(ok_vk), all(ok_vh), all(ok_a), all(ok_va))

## --- join onto the Task B Table S1 audit (keeps RBNS-only rows, blank CS/CVS) --
tb = read_csv(paste0(revDir, '1_S1_cellline_audit.csv'), show_col_types = FALSE)
out = tb %>%
  left_join(cs %>% select(RBP, merge_status,
                          CS_K562, CVS_K562, CS_HepG2, CVS_HepG2, CS_merged, CVS_merged,
                          top5mer_K562, top5mer_HepG2, top5mer_merged), by = 'RBP')

# round the new numeric columns for the table (leave audit columns untouched)
num_new = c('CS_K562','CVS_K562','CS_HepG2','CVS_HepG2','CS_merged','CVS_merged')
out[num_new] = lapply(out[num_new], function(x) round(x, 4))

write_csv(out, paste0(revDir, '1_S5_percell_CS_CVS.csv'), na = '')

## --- report -------------------------------------------------------------------
cat('\nRBPs: total', nrow(out),
    '| eCLIP', sum(!is.na(out$merge_status)),
    '| merged', sum(out$merge_status == 'merged', na.rm = TRUE),
    '| K562_only', sum(out$merge_status == 'K562_only', na.rm = TRUE),
    '| HepG2_only', sum(out$merge_status == 'HepG2_only', na.rm = TRUE),
    '| RBNS-only (blank CS)', sum(is.na(out$merge_status)), '\n\n')
cat('Sample rows:\n')
print(as.data.frame(out %>% filter(RBP %in% c('EIF4G2','RBFOX2','HNRNPC','AARS','A1CF','GRSF1')) %>%
      select(RBP, cell_line, merge_status, CS_K562, CS_HepG2, CS_merged, CVS_merged, published_CS)))
