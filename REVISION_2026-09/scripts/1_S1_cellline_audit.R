################################################################################
## Task B: Table S1 cell-line column (Reviewer 2, comment 2)
## Genome Biology revision, 2026-09
##
## Table S1 = the ENCODE_ACCESSION sheet (Dataset/Table_S1.xlsx; byte-identical to
## the Genome Biology submission media-2.xlsx and the CellReports resubmission copy).
## It has separate "K562 eCLIP" and "HepG2 eCLIP" accession columns, so cell line is
## implicit, and merged RBPs are not marked. This builds a per-RBP table with an
## explicit cell_line column and cross-checks it against two independent sources.
##
## Merge rule follows the published analysis (2026-03-16 run = Table S5): an RBP with
## eCLIP in both cell lines is 'K562+HepG2 merged'. RBNS-only RBPs get 'none (RBNS only)'.
################################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
})
options(width = 200)

repo   = '/Users/soonyi/Repos/BITS_Specificity/'
dataDir = paste0(repo, 'Dataset/Analysis/')
outDir = paste0(repo, 'REVISION_2026-09/')

## 1. Parse Table S1 (offset header, merged RBP cells)
################################################################################
raw = read_excel(paste0(repo, 'Dataset/Table_S1.xlsx'), sheet = 'ENCODE_ACCESSION',
                 col_names = FALSE, .name_repair = 'minimal')
raw = as.data.frame(raw)[, 1:5]
colnames(raw) = c('RBP', 'Kmer', 'RBNS', 'K562_eCLIP', 'HepG2_eCLIP')
raw = raw[-c(1, 2), ]                       # drop title + header rows
raw$RBP_ff = raw$RBP
for (i in seq_len(nrow(raw))) if (is.na(raw$RBP_ff[i]) && i > 1) raw$RBP_ff[i] = raw$RBP_ff[i - 1]

s1 = raw %>%
  group_by(RBP = RBP_ff) %>%
  summarise(
    RBNS_kmers   = paste(Kmer[!is.na(RBNS)], collapse = ','),
    K562_eCLIP   = paste(unique(K562_eCLIP[!is.na(K562_eCLIP)]), collapse = '|'),
    HepG2_eCLIP  = paste(unique(HepG2_eCLIP[!is.na(HepG2_eCLIP)]), collapse = '|'),
    .groups = 'drop'
  ) %>%
  mutate(across(c(K562_eCLIP, HepG2_eCLIP, RBNS_kmers), ~ ifelse(.x == '', NA, .x))) %>%
  filter(!is.na(RBP))

# cell_line from Table S1's own accession columns
s1 = s1 %>% mutate(
  has_K562  = !is.na(K562_eCLIP),
  has_HepG2 = !is.na(HepG2_eCLIP),
  cell_line = case_when(
    has_K562 &  has_HepG2 ~ 'K562+HepG2 merged',
    has_K562 & !has_HepG2 ~ 'K562',
    !has_K562 & has_HepG2 ~ 'HepG2',
    TRUE                  ~ 'none (RBNS only)'
  )
)

## 2. Cross-check source A: Table S5 K562/HepG2 flags + published CS
################################################################################
s5 = read_excel(paste0(repo, 'Dataset/Table_S5.xlsx'))
s5 = data.frame(RBP = s5$RBP, S5_eCLIP = s5$eCLIP, S5_K562 = s5$K562, S5_HepG2 = s5$HepG2,
                S5_CS = suppressWarnings(as.numeric(s5$CS)))
s5$S5_cell_line = with(s5, case_when(
  S5_eCLIP != 'Yes'          ~ 'none (RBNS only)',
  S5_K562 == 'Yes' & S5_HepG2 == 'Yes' ~ 'K562+HepG2 merged',
  S5_K562 == 'Yes'           ~ 'K562',
  S5_HepG2 == 'Yes'          ~ 'HepG2',
  TRUE                       ~ NA_character_
))

## 3. Cross-check source B: eCLIP_list.csv (the 2026-03-16 analysis input)
################################################################################
el = read.csv(paste0(dataDir, 'eCLIP_all/eCLIP_list.csv'), fileEncoding = 'UTF-8-BOM')
el_cell = el %>% group_by(RBP) %>%
  summarise(
    input_K562  = 'K562'  %in% CELL,
    input_HepG2 = 'HepG2' %in% CELL,
    .groups = 'drop'
  ) %>%
  mutate(input_cell_line = case_when(
    input_K562 & input_HepG2 ~ 'K562+HepG2 merged',
    input_K562               ~ 'K562',
    input_HepG2              ~ 'HepG2',
    TRUE                     ~ NA_character_
  ))

## 4. Reconcile
##    cell_line (authoritative) = what was actually analyzed. For eCLIP RBPs that is the
##    analysis input list (eCLIP_list.csv), verified against the narrowPeak files and
##    Table S5. For RBPs with no eCLIP it is 'none (RBNS only)'. The value implied by
##    Table S1's own accession columns is kept alongside, with a correction flag.
################################################################################
out = s1 %>%
  rename(cell_line_as_listed_in_TableS1 = cell_line) %>%
  left_join(s5 %>% select(RBP, S5_cell_line, S5_CS, S5_eCLIP), by = 'RBP') %>%
  left_join(el_cell %>% select(RBP, input_cell_line), by = 'RBP') %>%
  mutate(
    cell_line = case_when(
      !is.na(input_cell_line) ~ input_cell_line,                   # analyzed eCLIP RBP
      TRUE                    ~ 'none (RBNS only)'                 # no eCLIP
    ),
    in_published_CS = !is.na(S5_CS),
    TableS1_correction_needed = cell_line != cell_line_as_listed_in_TableS1,
    correction_note = case_when(
      RBP == 'GRSF1' ~ 'Table S1 lists ENCFF929AWR under K562; it is a HepG2 experiment (analyzed as HepG2).',
      RBP == 'MBNL1' ~ 'Table S1 has no eCLIP accession; MBNL1 was analyzed in K562 (ENCFF603WDI, missing from Table S1).',
      TRUE           ~ NA_character_
    ),
    agrees_TableS5 = is.na(S5_cell_line) | cell_line == S5_cell_line
  ) %>%
  select(RBP, cell_line, cell_line_as_listed_in_TableS1, TableS1_correction_needed, correction_note,
         K562_eCLIP, HepG2_eCLIP, RBNS_kmers,
         in_published_CS, published_CS = S5_CS,
         cell_line_TableS5 = S5_cell_line, agrees_TableS5) %>%
  arrange(RBP)

write_csv(out, paste0(outDir, '1_S1_cellline_audit.csv'), na = '')

## 5. Report
################################################################################
cat('Table S1 rows (RBPs):', nrow(out), '\n')
cat('\nAuthoritative cell-line breakdown (what was actually analyzed):\n')
print(table(out$cell_line))
cat('\nMerged (both lines):', sum(out$cell_line == 'K562+HepG2 merged'),
    ' K562-only:', sum(out$cell_line == 'K562'),
    ' HepG2-only:', sum(out$cell_line == 'HepG2'),
    ' no eCLIP:', sum(out$cell_line == 'none (RBNS only)'), '\n')

cat('\nRows where Table S1 needs correction:\n')
print(as.data.frame(out %>% filter(TableS1_correction_needed) %>%
  select(RBP, cell_line, cell_line_as_listed_in_TableS1, K562_eCLIP, HepG2_eCLIP, correction_note)))

cat('\nAuthoritative cell_line vs Table S5 flags — disagreements:\n')
d5 = out %>% filter(!agrees_TableS5)
if (nrow(d5) == 0) cat('  none (all agree)\n') else print(as.data.frame(d5 %>% select(RBP, cell_line, cell_line_TableS5)))

cat('\nHead of output:\n')
print(as.data.frame(head(out, 12)))
