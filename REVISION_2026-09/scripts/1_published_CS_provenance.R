################################################################################
## Pre-work check: which eCLIP run produced the published CS/CVS?
## Compares Table S5 (Dataset/Table_S5.xlsx) and the Fig. 2C/2D text values
## against every stored eCLIP 5-mer run, package-era and pre-package.
## Console output only.
################################################################################

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(RBPSpecificity)
})
options(width = 220)

repo = '/Users/soonyi/Repos/BITS_Specificity/Dataset/Analysis/'
old = '/Users/soonyi/Desktop/Genomics/Specificity/'
num = function(x) suppressWarnings(as.numeric(x))

## Published values
################################################################################
s5 = read_excel('/Users/soonyi/Repos/BITS_Specificity/Dataset/Table_S5.xlsx')
s5 = data.frame(RBP = s5$RBP, IS = num(s5$IS), VS = num(s5$VS), CS = num(s5$CS), CVS = num(s5$CVS))

st = read.csv(paste0(repo, 'sample_table_eCLIP.csv'), fileEncoding = 'UTF-8-BOM')
st$RBP0 = sub('_[0-9]+$', '', st$RBP)
structure = st %>% distinct(RBP0, structure_context) %>% rename(RBP = RBP0)

## Helpers (Figure_2.R steps)
################################################################################
min_max_norm = function(x, a = 0, b = 1) a + (x - min(x)) * (b - a) / (max(x) - min(x))

read_vec = function(path) {
  x = read_csv(path, show_col_types = FALSE, name_repair = 'minimal')
  x = x[, c('MOTIF', 'Score')]
  x$MOTIF = gsub('T', 'U', x$MOTIF)
  data.frame(x %>% arrange(MOTIF))
}

metrics = function(v) {
  top = v$MOTIF[which.max(v$Score)]
  c(CS = suppressMessages(returnIS(v, top)),
    CVS = suppressMessages(returnMS(v, top, output_type = 'number')))
}

# Rebuild Fig. 2 CS/CVS from per-cell-line vectors with the Figure_2.R merge rules.
fig2_from_vectors = function(dir) {
  path = function(cell, RBP) paste0(dir, cell, '_eCLIP_NormalizedEnrichment_25ntExt_', RBP, '_5mer.csv')
  bind_rows(lapply(unique(st$RBP0), function(R) {
    cells = st$Cell[st$RBP0 == R]
    if (R == 'HNRNPK') {
      v1 = read_vec(path('K562', 'HNRNPK_1'))
      v2 = read_vec(path('K562', 'HNRNPK_2'))
      v = data.frame(MOTIF = v1$MOTIF, Score = min_max_norm((v1$Score + v2$Score) / 2))
    } else if (length(cells) == 2 && R != 'XRCC6') {
      k = read_vec(path('K562', R))
      h = read_vec(path('HepG2', R))
      v = data.frame(MOTIF = k$MOTIF, Score = min_max_norm((h$Score + k$Score) / 2))
    } else {
      v = read_vec(path(if ('K562' %in% cells) 'K562' else 'HepG2', R))
    }
    m = metrics(v)
    data.frame(RBP = R, CS = m[['CS']], CVS = m[['CVS']])
  }))
}

read_summary = function(p) {
  x = read_csv(p, show_col_types = FALSE)
  data.frame(RBP = x$RBP, CS = num(x$eCLIP_Specificity), CVS = num(x$eCLIP_Sensitivity))
}

## Candidate runs
################################################################################
cands = list()
f6 = read_csv(paste0(repo, 'eCLIP_all/output/20260316/eCLIP_5mer_IS_VS_summary.csv'), show_col_types = FALSE)
cands[['2026-03-16 pkg 0.99.0 (Fig6 summary)']] =
  data.frame(RBP = f6$RBP, CS = coalesce(f6$IS_avg, f6$IS_K562, f6$IS_HepG2),
             CVS = coalesce(f6$VS_avg, f6$VS_K562, f6$VS_HepG2))
cands[['2026-03-17 pkg 0.99.0 (Fig2 summary)']] =
  read_summary(paste0(repo, 'eCLIP_peak/output/20260317_5mer_RBNS_eCLIP_Analysis_25ntExt.csv'))
cands[['2026-02-21 pkg (rebuilt, Fig2 merge)']] =
  tryCatch(fig2_from_vectors(paste0(old, 'AnalysisOutput/20260221/')),
           error = function(e) { message('20260221 skipped: ', conditionMessage(e)); NULL })
cands[['2025-09-13 pkg early (25nt)']] =
  read_summary(paste0(old, 'AnalysisOutput/20250913_5mer_RBNS_eCLIP_Analysis_25ntExt.csv'))
cands[['2025-08-26 pkg early (25nt)']] =
  read_summary(paste0(old, 'AnalysisOutput/20250826_5mer_RBNS_eCLIP_Analysis_25ntExt.csv'))
cands[['2025-04-27 pre-package script']] =
  read_summary(paste0(old, 'Data_Final/eCLIP_peak/Output/20250427_5_mer_RBNS_eCLIP_Analysis.csv'))
cands = cands[!vapply(cands, is.null, logical(1))]

# Sanity: the Fig. 2 rebuild reproduces the 20260317 summary from its own vectors.
chk = inner_join(fig2_from_vectors(paste0(repo, 'eCLIP_peak/output/20260317/')),
                 cands[['2026-03-17 pkg 0.99.0 (Fig2 summary)']], by = 'RBP')
cat('Rebuild check vs 20260317 summary: max |CS diff| =', max(abs(chk$CS.x - chk$CS.y)),
    '| max |CVS diff| =', max(abs(chk$CVS.x - chk$CVS.y)), '\n')

## 1. Table S5 CS/CVS vs each run
################################################################################
rel = function(a, b) abs(a - b) / abs(a)
cmp = bind_rows(lapply(names(cands), function(nm) {
  m = inner_join(s5 %>% filter(!is.na(CS)), cands[[nm]], by = 'RBP', suffix = c('_pub', '_run'))
  m28 = m %>% filter(!is.na(IS))
  data.frame(run = nm,
             n28 = nrow(m28),
             CS_exact_28 = sum(rel(m28$CS_pub, m28$CS_run) < 1e-6),
             CS_median_pctdiff_28 = round(100 * median(rel(m28$CS_pub, m28$CS_run)), 2),
             CS_max_pctdiff_28 = round(100 * max(rel(m28$CS_pub, m28$CS_run)), 2),
             CVS_exact_28 = sum(rel(m28$CVS_pub, m28$CVS_run) < 1e-6),
             n_all = nrow(m),
             CS_exact_all = sum(rel(m$CS_pub, m$CS_run) < 1e-6))
}))
cat('\n1. Table S5 CS/CVS vs stored runs (28 = RBPs with RBNS and eCLIP)\n')
print(cmp)

## 2. Values quoted in the text
################################################################################
cat('\n2. CS quoted in text: RBFOX2 3.39, EIF4G2 16.5\n')
txt = s5 %>% filter(RBP %in% c('RBFOX2', 'EIF4G2', 'HNRNPC', 'PCBP2')) %>% select(RBP, TableS5 = CS)
for (nm in names(cands)) {
  txt = left_join(txt, cands[[nm]] %>% select(RBP, CS) %>% rename(!!nm := CS), by = 'RBP')
}
print(txt %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

## 3. Fig. 2C/2D correlations (text: H IS-CS r = -0.32; L IS-CS r = 0.94; L VS-CVS r = 0.84)
################################################################################
fig2cd = function(nm, cand) {
  m = inner_join(s5 %>% filter(!is.na(IS)) %>% select(RBP, IS, VS),
                 cand %>% select(RBP, CS, CVS), by = 'RBP') %>%
    inner_join(structure, by = 'RBP')
  H = m %>% filter(structure_context == 'H')
  L = m %>% filter(structure_context == 'L')
  data.frame(run = nm, nH = nrow(H), r_IS_CS_H = cor(H$IS, H$CS),
             nL = nrow(L), r_IS_CS_L = cor(L$IS, L$CS), r_VS_CVS_L = cor(L$VS, L$CVS))
}
cat('\n3. Fig. 2C/2D Pearson r recomputed from each run (IS/VS from Table S5)\n')
print(bind_rows(fig2cd('Table S5 (published)', s5 %>% select(RBP, CS, CVS)),
                bind_rows(lapply(names(cands), function(nm) fig2cd(nm, cands[[nm]])))) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))))
