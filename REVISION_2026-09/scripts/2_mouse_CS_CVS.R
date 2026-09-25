################################################################################
## Task C step 5: provisional CS/CVS for the mouse eIF4G2 HITS-CLIP peaks
## (Reviewer 2, comment 1)
##
## PROVISIONAL: the mouse set is thin (1237 COMBINED resting peaks, ~329k uniquely
## mapped reads, short reads). Treat the values as qualitative.
##
## Exact-procedure replication: RBPSpecificity 0.99.0's motifEnrichment() cannot map
## "mm39" via selectGenome(), so its body is replicated here by calling the same internal
## functions with a manually-loaded mm39 BSgenome injected in place of selectGenome().
## No package source is modified. Parameters match the ENCODE eCLIP procedure used for
## Table S5: K=5, extension c(25,0), method "subtract", bkg_iter 100, dist 500-1000,
## scramble FALSE, then min-max to [1,e] and natural log. BED coords are passed to
## peakParse as-is (0-based start treated as 1-based), exactly as Figure_2.R fed eCLIP peaks.
################################################################################

suppressPackageStartupMessages({
  library(RBPSpecificity); library(BSgenome.Mmusculus.UCSC.mm39)
  library(data.table); library(dplyr); library(readr)
})
options(width = 150); set.seed(1)

dataDir = '/Users/soonyi/Repos/BITS_Specificity/Dataset/Analysis/'
outDir  = '/Users/soonyi/Repos/BITS_Specificity/REVISION_2026-09/'
peakBed = '/Users/soonyi/Desktop/Specificity_Genome_Biology_Resubmission/EIF4G2_HITS_CLIP/mouse_EIF4G2_HITS_CLIP_fastq_output/04_PEAKS/COMBINED_PEAKS/FINAL_COMBINED_PEAKS.bed'
genome  = BSgenome.Mmusculus.UCSC.mm39

gc_frac = function(m) nchar(gsub('[^GC]', '', m)) / nchar(m)

## --- replicate motifEnrichment() body with mm39 genome injected ---------------
motifEnrichment_injected = function(coordinates, genome_obj, K = 5, extension = c(25, 0),
                                    bkg_iter = 100, bkg_min_dist = 500, bkg_max_dist = 1000,
                                    scramble_bkg = FALSE, nucleic_acid_type = 'DNA') {
  ns = asNamespace('RBPSpecificity')
  peak_gr = ns$peakParse(input = coordinates)
  peak_seq = ns$getSequence(granges_obj = peak_gr, genome_obj = genome_obj,
                            extension = as.integer(extension), min_length = K)
  peak_counts = ns$countKmers(sequences = peak_seq, K = K, type = nucleic_acid_type)
  bkg = ns$countKmersBkg(original_peak_gr = peak_gr, K = K, type = nucleic_acid_type,
                         genome_obj = genome_obj, bkg_min_dist = bkg_min_dist,
                         bkg_iter = bkg_iter, bkg_max_dist = bkg_max_dist, scramble_bkg = scramble_bkg)
  enr = ns$calEnrichment(peak_kmer_counts_df = peak_counts, avg_bkg_counts_df = bkg, method = 'subtract')
  scaled = ns$normalizeScores(enr$EnrichmentScore, method = 'min_max', a = 1, b = exp(1))
  final = log(scaled)
  res = data.frame(MOTIF = enr$MOTIF, Score = final, stringsAsFactors = FALSE)
  res$Score[is.nan(res$Score)] = 0
  res
}

## --- run on mouse COMBINED peaks ----------------------------------------------
Peak = fread(peakBed, select = 1:6)
setnames(Peak, c('chr', 'start', 'end', 'name', 'score', 'strand'))
Peak = Peak[chr %in% paste0('chr', c(1:19, 'X', 'Y', 'M'))]
Peak = as.data.frame(Peak)
enr = motifEnrichment_injected(Peak, genome, K = 5, extension = c(25, 0), scramble_bkg = FALSE)

top = enr$MOTIF[which.max(enr$Score)]
CS  = suppressMessages(returnIS(enr, top))
CVS = suppressMessages(returnMS(enr, top, output_type = 'number'))

## --- top 20 5-mers + GC -------------------------------------------------------
enrU = enr %>% mutate(MOTIF = gsub('T', 'U', MOTIF)) %>% arrange(desc(Score))
top20 = enrU %>% head(20) %>% mutate(rank = row_number(), GC = sapply(MOTIF, gc_frac))

## --- correlations vs ENCODE eCLIP and RBNS ------------------------------------
readvec = function(path, col) {
  x = read_csv(path, show_col_types = FALSE)
  if (is.null(col)) { d = data.frame(MOTIF = x$MOTIF, S = x$Score) }
  else { d = data.frame(MOTIF = x$Motif, S = x[[col]]) }
  d$MOTIF = gsub('T', 'U', d$MOTIF); d %>% arrange(MOTIF)
}
mouse_v = enrU %>% arrange(MOTIF)
enc_v = readvec(paste0(dataDir, 'eCLIP_all/output/20260316/K562_eCLIP_Enrichment_EIF4G2_5mer.csv'), NULL)
rbns_v = readvec(paste0(dataDir, 'RBNS/RBNS_normalized_5mer.csv'), 'EIF4G2')
stopifnot(identical(mouse_v$MOTIF, enc_v$MOTIF), identical(mouse_v$MOTIF, rbns_v$MOTIF))

cor_row = function(lab, a, b) data.frame(record_type = 'correlation', item = lab,
  pearson = round(cor(a, b, method = 'pearson'), 4),
  spearman = round(cor(a, b, method = 'spearman'), 4), n = length(a))
cors = bind_rows(
  cor_row('mouse_vs_ENCODE_EIF4G2_eCLIP', mouse_v$Score, enc_v$S),
  cor_row('mouse_vs_EIF4G2_RBNS',         mouse_v$Score, rbns_v$S)
)

## --- write CSV ----------------------------------------------------------------
summary_rows = data.frame(record_type = 'summary',
  item = c('CS', 'CVS', 'top_motif', 'meanGC_top20', 'meanGC_all1024', 'n_peaks', 'n_top20_GCge0.6'),
  value = c(round(CS, 3), round(CVS, 4), top, round(mean(top20$GC), 3),
            round(mean(sapply(enrU$MOTIF, gc_frac)), 3), nrow(Peak), sum(top20$GC >= 0.6)))
top20_out = top20 %>% transmute(record_type = 'top5mer', rank, motif = MOTIF,
  score = round(Score, 4), GC_fraction = GC)
res_out = bind_rows(summary_rows, cors, top20_out)
write_csv(res_out, paste0(outDir, '2_mouse_CS_CVS.csv'), na = '')

## --- report -------------------------------------------------------------------
cat('=== PROVISIONAL mouse eIF4G2 HITS-CLIP CS/CVS (1237 resting peaks) ===\n')
cat('CS =', round(CS, 3), ' CVS =', round(CVS, 4), ' top 5-mer =', top,
    ' mean GC(top20) =', round(mean(top20$GC), 3), '\n\n')
cat('Top 20 5-mers by enrichment:\n'); print(as.data.frame(top20 %>% select(rank, MOTIF, Score, GC)))
cat('\nCorrelation of mouse 5-mer vector vs:\n'); print(cors)
cat('\nFor reference: ENCODE EIF4G2 eCLIP top = GUGUG (CS 16.4); RBNS EIF4G2 top = GUUGC (IS 2.04).\n')
