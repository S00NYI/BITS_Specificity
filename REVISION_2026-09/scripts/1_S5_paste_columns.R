################################################################################
## Task B2b: build the per-cell-line CS/CVS columns to paste into Table S5.
## Output rows are in EXACTLY Table S5's row order, with RBP names + the existing
## S5 CS/CVS as verification columns, so the paste can be checked row-for-row.
##   New columns to paste: CS_K562, CVS_K562, CS_HepG2, CVS_HepG2
##   Verification only    : S5_CS_existing, S5_CVS_existing, K562_flag, HepG2_flag
################################################################################

suppressPackageStartupMessages({ library(readxl); library(readr); library(dplyr) })

revDir = '/Users/soonyi/Repos/BITS_Specificity/REVISION_2026-09/'
s5 = read_excel('/Users/soonyi/Repos/BITS_Specificity/Dataset/Table_S5.xlsx')
b2 = read_csv(paste0(revDir, '1_S5_percell_CS_CVS.csv'), show_col_types = FALSE)

num = function(x) suppressWarnings(as.numeric(ifelse(x %in% c('.', '', 'NA'), NA, x)))

s5o = s5 %>% transmute(
  RBP,
  S5_CS_existing  = num(CS),
  S5_CVS_existing = num(CVS),
  K562_flag  = K562,
  HepG2_flag = HepG2,
  .s5_row = row_number())               # preserve S5 order

out = s5o %>%
  left_join(b2 %>% select(RBP, CS_K562, CVS_K562, CS_HepG2, CVS_HepG2,
                          CS_merged, merge_status), by = 'RBP') %>%
  arrange(.s5_row)

## --- alignment check: percell merged CS must equal the value already in S5 -------
chk = out %>% filter(!is.na(S5_CS_existing), !is.na(CS_merged))
gap = abs(chk$S5_CS_existing - chk$CS_merged)
cat('rows in S5:', nrow(out),
    '| eCLIP rows w/ CS in S5:', nrow(chk),
    '| max |S5_CS - percell_CS_merged|:', signif(max(gap), 3), '\n')
unmatched = out %>% filter(!is.na(S5_CS_existing), is.na(CS_merged))
if (nrow(unmatched)) { cat('!! S5 RBPs with a CS but no percell match:\n'); print(unmatched$RBP) }
# S5 stores CS at ~5-decimal precision, so allow rounding slack; a real row
# misalignment would differ by whole CS units, far above this.
stopifnot(max(gap) < 1e-3, nrow(unmatched) == 0)

## --- write paste-ready file (S5 order) ----------------------------------------
paste_ready = out %>%
  transmute(RBP,
            CS_K562  = round(CS_K562, 4),  CVS_K562  = round(CVS_K562, 4),
            CS_HepG2 = round(CS_HepG2, 4), CVS_HepG2 = round(CVS_HepG2, 4),
            S5_CS_existing = round(S5_CS_existing, 4),
            S5_CVS_existing = round(S5_CVS_existing, 4),
            K562_flag, HepG2_flag, merge_status)
write_csv(paste_ready, paste0(revDir, '1_S5_paste_columns.csv'), na = '')

cat('\nwrote 1_S5_paste_columns.csv (', nrow(paste_ready), 'rows, S5 order )\n')
cat('preview (first eCLIP rows of each kind):\n')
print(as.data.frame(paste_ready %>%
      filter(RBP %in% c('AARS','EIF4G2','HNRNPC','RBFOX2','GRSF1','MBNL1','A1CF')) ))
