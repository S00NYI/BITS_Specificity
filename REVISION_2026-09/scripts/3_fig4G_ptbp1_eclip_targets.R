################################################################################
## Define PTBP1 targets from REAL PTBP1 eCLIP (ENCODE, hg38) instead of the sim/motif
## proxies, and ask whether PTBP1 targeting explains the +PTBP1 helps/hurts split.
##
## Pool binding sites are hg19; PTBP1 eCLIP peaks are hg38 (manuscript ran
## motifEnrichment(...,'hg38')). We match at the GENE level (avoids liftover): count
## PTBP1 eCLIP peaks (K562 + HepG2) overlapping each pool gene's hg38 span, and total
## peak signal, then compare to the per-transcript dP_PTBP1.
################################################################################

suppressPackageStartupMessages({ library(data.table); library(readr); library(dplyr); library(GenomicRanges) })
revDir = '/Users/soonyi/Repos/BITS_Specificity/REVISION_2026-09/'
npDir  = '/Users/soonyi/Repos/BITS_Specificity/Dataset/Analysis/eCLIP_all/narrowPeaks/'
genesF = '/private/tmp/claude-501/-Users-soonyi-Repos-BITS-Specificity/ea84054a-8c3b-4e31-9a1c-8b844676e48a/scratchpad/gencode_genes.tsv'

## --- pool gene hg38 spans (TENT2 == PAPD4) ------------------------------------
gt = fread(genesF, header = FALSE, sep = '\t')
gt[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', V9)]
pool_syms = c('MYC','PTBP2','NF1','TENT2','MAT2A','MYL6','C4BPB','PCBP2','MALAT1')
gg = gt[gene_name %in% pool_syms, .(gene_name, chr = V1, start = V4, end = V5, strand = V7)]
gg[, transcript := ifelse(gene_name == 'TENT2', 'PAPD4', gene_name)]
genes_gr = GRanges(gg$chr, IRanges(gg$start, gg$end), transcript = gg$transcript)

## --- PTBP1 eCLIP peaks (K562 + HepG2), hg38 -----------------------------------
read_np = function(f, cell) { d = fread(f, header = FALSE)
  data.table(chr = d$V1, start = d$V2, end = d$V3, signal = d$V7, cell = cell) }
pk = rbind(read_np(paste0(npDir, 'PTBP1_K562_narrowPeak.bed'), 'K562'),
           read_np(paste0(npDir, 'PTBP1_HepG2_narrowPeak.bed'), 'HepG2'))
pk_gr = GRanges(pk$chr, IRanges(pk$start, pk$end), signal = pk$signal, cell = pk$cell)

## --- count PTBP1 peaks / signal per pool gene ---------------------------------
ov = findOverlaps(pk_gr, genes_gr)
tab = data.table(transcript = genes_gr$transcript[subjectHits(ov)],
                 cell = pk$cell[queryHits(ov)], signal = pk$signal[queryHits(ov)])
per_gene = tab[, .(ptbp1_peaks = .N,
                   ptbp1_peaks_K562 = sum(cell == 'K562'), ptbp1_peaks_HepG2 = sum(cell == 'HepG2'),
                   ptbp1_total_signal = sum(signal), ptbp1_max_signal = max(signal)), by = transcript]
# genes with zero PTBP1 peaks
per_gene = merge(data.table(transcript = gg$transcript), per_gene, by = 'transcript', all.x = TRUE)
per_gene[is.na(per_gene)] = 0
per_gene[, gene_kb := round((gg$end[match(transcript, gg$transcript)] - gg$start[match(transcript, gg$transcript)]) / 1000)]
per_gene[, ptbp1_peaks_per_kb := round(ptbp1_peaks / gene_kb, 3)]

## --- join to the helps/hurts result -------------------------------------------
pool = read_csv(paste0(revDir, '3_fig4G_pool_per_transcript.csv'), show_col_types = FALSE) %>%
  filter(transcript != 'ALL POOLED') %>% select(transcript, dP_PTBP1, p_PTBP1)
res = pool %>% left_join(as.data.frame(per_gene), by = 'transcript') %>%
  mutate(group = ifelse(dP_PTBP1 > 0, 'helps', 'hurts')) %>% arrange(desc(dP_PTBP1))
write_csv(res, paste0(revDir, '3_fig4G_ptbp1_eclip_targets.csv'), na = '')
print(as.data.frame(res %>% select(transcript, dP_PTBP1, group, gene_kb, ptbp1_peaks,
                                   ptbp1_peaks_per_kb, ptbp1_total_signal)), row.names = FALSE)

cat('\n--- group medians ---\n')
print(res %>% group_by(group) %>% summarise(n = n(), ptbp1_peaks = median(ptbp1_peaks),
        ptbp1_peaks_per_kb = median(ptbp1_peaks_per_kb), ptbp1_total_signal = median(ptbp1_total_signal)))
cat('\nSpearman(dP_PTBP1, PTBP1 eCLIP peaks):', round(cor(res$dP_PTBP1, res$ptbp1_peaks, method='spearman'),3),
    '| per-kb:', round(cor(res$dP_PTBP1, res$ptbp1_peaks_per_kb, method='spearman'),3),
    '| total signal:', round(cor(res$dP_PTBP1, res$ptbp1_total_signal, method='spearman'),3), '\n')
