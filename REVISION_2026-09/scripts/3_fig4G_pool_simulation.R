################################################################################
## Task M/N extended to the FULL Sutandy in-vitro pool (all transcripts with iCLIP
## coverage), so the reviewer's "does it hold across many regions" is answered across
## the whole assay, not just PTBP2.
##
## Coordinates are hg19 (verified: the provided PTBP2.fa == hg19 chr1 + strand at its
## header coords). PTBP2 uses its provided FASTA (anchor); every other transcript's
## fragment is fetched from hg19 + strand over [site_min-FLANK, site_max+FLANK] (PTBP2's
## fragment sits ~35 nt outside its sites, so FLANK=50). Same convention as Figure_3_G.R:
## + strand genomic sequence + + strand bedgraph coords, tx_pos = genomic - start + 1.
## Sim params identical to the Fig 4G baseline (see 3_fig4G_cotarget_violin).
##
## CAVEATS: (1) MIRLET7A has no bedgraph coverage in its locus -> dropped. (2) fragments
## fetched as + strand genomic to match the PTBP2 convention; genes whose sense is the
## - strand would need the reverse complement (flagged, not corrected here). Validation:
## PTBP2 must reproduce the ~0.40 -> 0.50 baseline.
##
## Output: REVISION_2026-09/3_fig4G_pool_per_transcript.csv (per-transcript + pooled) ; objects
## pool_violin, pool_bytx.
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(Biostrings); library(BSgenome.Hsapiens.UCSC.hg19); library(RBPEqBind)
})
WIN = 200; MINPTS = 20; FLANK = 50
revDir = '/Users/soonyi/Repos/BITS_Specificity/REVISION_2026-09/'
simDir = path.expand('~/Repos/BITS_Specificity/Dataset/Analysis/RBPEqBind_Simulation')
modelFile = file.path(simDir, 'DATA_PROCESSED', 'rnacompete_affinity_scores.csv')
bgBedDir  = file.path(simDir, 'RESULTS_ANALYSIS', 'BEDGRAPH_PROCESSED')
sitesFile = file.path(simDir, 'DATA_PROCESSED', 'binding_sites.csv')
ptbp2_fa  = file.path(simDir, 'DATA_PROCESSED', 'FASTA', 'PTBP2.fa')
genome = BSgenome.Hsapiens.UCSC.hg19
u2af2_conc=500; hnrnpc_conc=200; ptbp1_conc=200; rna_conc=6.75; k=7
scratch = tempdir()

## ---- model + sim wrapper -----------------------------------------------------
model_raw = read_csv(modelFile, show_col_types = FALSE); colnames(model_raw)[1] = 'Motif'
tmpm = file.path(scratch, 'model.csv'); write_csv(model_raw, tmpm); raw_models = loadModel(tmpm)
sim_one = function(fa, rbps, concs) {
  rm = setModel(raw_models, max_affinity = setNames(1/rep(50, length(rbps)), rbps), min_affinity = 1e-5)
  simulateBindingF(fasta_file = fa, rbp_models = rm[rbps], protein_concs = concs, rna_conc = rna_conc, k = k)
}
bg_all = lapply(c(U='bg_u2af2_500nM.bedgraph', H='bg_u2af2_500nM_hnrnpc_200nM.bedgraph',
                  P='bg_u2af2_500nM_ptbp1_200nM.bedgraph'),
                function(f) as.data.table(read_tsv(file.path(bgBedDir, f),
                     col_names = c('chr','start','end','val'), show_col_types = FALSE)))
bg_track = function(cond, chrom, gstart, gend) {        # clean per-base, tx_pos = genomic - gstart + 1
  d = bg_all[[cond]][chr == chrom & start >= gstart - 1 & end <= gend + 1]
  if (!nrow(d)) return(data.table(pos = integer(), signal = numeric()))
  pos = unlist(lapply(seq_len(nrow(d)), function(i) seq.int(d$start[i], d$end[i]-1L)))
  data.table(pos = pos - gstart + 1L, signal = rep(d$val, d$end - d$start))[!duplicated(pos)]
}

## ---- transcript loci from binding_sites --------------------------------------
sites = read_csv(sitesFile, show_col_types = FALSE) %>%
  mutate(chr = ifelse(grepl('^chr', chromosome), chromosome, paste0('chr', tolower(sub('Chr','',chromosome)))),
         start = as.numeric(start), end = as.numeric(end))
loci = sites %>% group_by(transcript) %>%
  summarise(chr = chr[1], site_min = min(start), site_max = max(end), n_sites = n(), .groups='drop')

## ---- per-transcript pipeline -------------------------------------------------
corr_win = function(a, lo, hi) { s = a[pos >= lo & pos < hi]
  if (nrow(s) < MINPTS || sd(s$density)==0 || sd(s$signal)==0) return(NA)
  cor(s$density, s$signal, method = 'pearson') }

run_tx = function(tx) {
  L = loci[loci$transcript == tx, ]
  chr = L$chr
  if (tx == 'PTBP2') { gstart = 97269727; gend = 97272451; fa = ptbp2_fa
  } else {
    gstart = round(L$site_min - FLANK); gend = round(L$site_max + FLANK)
    dss = DNAStringSet(getSeq(genome, chr, gstart, gend))          # + strand hg19
    fa = file.path(scratch, paste0(tx, '.fa'))
    writeXStringSet(setNames(dss,
      sprintf('%s|%s:%d-%d|+|len=%d', tx, chr, gstart, gend, width(dss))), fa)
  }
  tx_name = sub('^>', '', readLines(fa, n = 1))
  ali = lapply(c('U','H','P'), function(cond) {
    conc = switch(cond, U=c(U2AF2=u2af2_conc), H=c(U2AF2=u2af2_conc,HNRNPC=hnrnpc_conc), P=c(U2AF2=u2af2_conc,PTBP1=ptbp1_conc))
    rbps = switch(cond, U='U2AF2', H=c('U2AF2','HNRNPC'), P=c('U2AF2','PTBP1'))
    sd = sim_one(fa, rbps, conc)[transcript == tx_name, .(pos, density = U2AF2_density)]
    merge(sd, bg_track(cond, chr, gstart, gend), by = 'pos') })
  names(ali) = c('U','H','P')
  if (any(sapply(ali, nrow) < MINPTS)) return(NULL)              # no coverage -> drop (e.g. MIRLET7A)
  # binding-site-centered windows
  ctr = round((sites$start[sites$transcript==tx] + sites$end[sites$transcript==tx]) / 2) - gstart + 1
  cov = range(unlist(lapply(ali, function(a) a$pos)))
  wins = data.frame(lo = ctr - WIN/2, hi = ctr + WIN/2) %>% filter(lo >= cov[1], hi <= cov[2])
  w = bind_rows(lapply(seq_len(nrow(wins)), function(i)
        data.frame(U = corr_win(ali$U, wins$lo[i], wins$hi[i]),
                   H = corr_win(ali$H, wins$lo[i], wins$hi[i]),
                   P = corr_win(ali$P, wins$lo[i], wins$hi[i])))) %>%
      filter(!is.na(U), !is.na(H), !is.na(P))
  if (!nrow(w)) return(NULL)
  w$transcript = tx
  attr(w, 'roi_diag') = if (tx == 'PTBP2') sapply(ali, function(a){ s=a[pos>=1854 & pos<2054]; cor(s$density,s$signal) }) else NULL
  w
}

res = lapply(loci$transcript, function(t) tryCatch(run_tx(t), error = function(e) { message(t,': ',conditionMessage(e)); NULL }))
res = Filter(Negate(is.null), res)
pool = bind_rows(res)
write_csv(pool %>% mutate(across(where(is.numeric), ~round(.x, 5))),   # per-window (U/H/P) for downstream splits
          paste0(revDir, '3_fig4G_pool_per_window.csv'), na = '')

## ---- PTBP2 validation --------------------------------------------------------
p2 = res[[which(sapply(res, function(w) w$transcript[1] == 'PTBP2'))]]
cat('PTBP2 ROI diagnostic (expect ~0.38/0.38/0.46):',
    paste(sprintf('%.3f', attr(p2, 'roi_diag')), collapse = ' / '), '\n\n')

## ---- per-transcript + pooled summary -----------------------------------------
wp = function(a,b) if (length(a) < 3) NA else suppressWarnings(wilcox.test(a, b, paired = TRUE)$p.value)
bytx = pool %>% group_by(transcript) %>% summarise(n_win = n(),
         medP_U2AF2 = median(U), medP_PTBP1 = median(P), medP_HNRNPC = median(H),
         dP_PTBP1 = median(P - U), p_PTBP1 = wp(P, U),
         dP_HNRNPC = median(H - U), p_HNRNPC = wp(H, U), .groups='drop') %>% arrange(desc(n_win))
pooled = data.frame(transcript = 'ALL POOLED', n_win = nrow(pool),
         medP_U2AF2 = median(pool$U), medP_PTBP1 = median(pool$P), medP_HNRNPC = median(pool$H),
         dP_PTBP1 = median(pool$P - pool$U), p_PTBP1 = wp(pool$P, pool$U),
         dP_HNRNPC = median(pool$H - pool$U), p_HNRNPC = wp(pool$H, pool$U))
out = bind_rows(pooled, bytx) %>% mutate(across(where(is.numeric), ~round(.x, 4)))
write_csv(out, paste0(revDir, '3_fig4G_pool_per_transcript.csv'), na = '')
cat('transcripts kept:', nrow(bytx), 'of', nrow(loci), ' | pooled binding-site windows:', nrow(pool), '\n\n')
print(as.data.frame(out %>% select(transcript, n_win, medP_U2AF2, medP_PTBP1, dP_PTBP1, p_PTBP1, dP_HNRNPC, p_HNRNPC)))

## ---- plots -------------------------------------------------------------------
cols = c('U2AF2 only'='#34495E','U2AF2 + HNRNPC'='#E91E63','U2AF2 + PTBP1'='#F39C12')
long = pool %>% transmute(transcript, `U2AF2 only`=U, `U2AF2 + HNRNPC`=H, `U2AF2 + PTBP1`=P) %>%
  pivot_longer(-transcript, names_to='Condition', values_to='Pearson') %>%
  mutate(Condition = factor(Condition, levels = names(cols)))
pool_violin = ggplot(long, aes(Condition, Pearson, fill = Condition)) +
  geom_violin(alpha=0.5, colour=NA) + geom_boxplot(width=0.15, outlier.shape=NA) +
  geom_jitter(width=0.08, size=0.4, alpha=0.3) + scale_fill_manual(values=cols) +
  labs(title = sprintf('Sutandy pool: per-binding-site correlation across %d transcripts (%d windows)',
                       nrow(bytx), nrow(pool)), x=NULL, y='Per-window Pearson r') +
  theme_bw() + theme(legend.position='none', plot.title=element_text(hjust=0.5, size=12),
                     axis.text=element_text(size=12), axis.title=element_text(size=12, face='bold'))
pool_bytx = ggplot(bytx, aes(x = reorder(transcript, dP_PTBP1), y = dP_PTBP1, fill = p_PTBP1 < 0.05)) +
  geom_col() + geom_hline(yintercept = 0, linetype='dashed') + coord_flip() +
  scale_fill_manual(values = c('FALSE'='grey70','TRUE'='#F39C12'), name='p<0.05') +
  labs(title = 'Per-transcript +PTBP1 improvement', x=NULL, y='median dPearson (+PTBP1 - U2AF2 only)') +
  theme_bw() + theme(plot.title=element_text(hjust=0.5, size=12),
                     axis.text=element_text(size=11), axis.title=element_text(size=12, face='bold'))
