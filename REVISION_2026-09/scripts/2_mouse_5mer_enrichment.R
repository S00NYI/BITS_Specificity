################################################################################
## Task C byproduct: dump the FULL 1024-mer mouse eIF4G2 HITS-CLIP enrichment
## vector (not just the top 20) so the affinity-distribution figure (2_mouse_affinity_distribution) can be
## drawn without re-running the genome compute each time.
##
## Identical procedure to 2_mouse_CS_CVS.R (RBPSpecificity 0.99.0 internals with
## mm39 BSgenome injected in place of selectGenome(); no package source modified).
## Writes: REVISION_2026-09/2_mouse_5mer_enrichment.csv  (MOTIF in U alphabet, Score)
## Prints CS / top motif as a reproduction check against 2_mouse_CS_CVS.csv.
################################################################################

suppressPackageStartupMessages({
  library(RBPSpecificity); library(BSgenome.Mmusculus.UCSC.mm39)
  library(data.table); library(dplyr); library(readr)
})
options(width = 150); set.seed(1)

outDir  = '/Users/soonyi/Repos/BITS_Specificity/REVISION_2026-09/'
peakBed = '/Users/soonyi/Desktop/Specificity_Genome_Biology_Resubmission/EIF4G2_HITS_CLIP/mouse_EIF4G2_HITS_CLIP_fastq_output/04_PEAKS/COMBINED_PEAKS/FINAL_COMBINED_PEAKS.bed'
genome  = BSgenome.Mmusculus.UCSC.mm39

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

Peak = fread(peakBed, select = 1:6)
setnames(Peak, c('chr', 'start', 'end', 'name', 'score', 'strand'))
Peak = Peak[chr %in% paste0('chr', c(1:19, 'X', 'Y', 'M'))]
Peak = as.data.frame(Peak)

enr = motifEnrichment_injected(Peak, genome, K = 5, extension = c(25, 0), scramble_bkg = FALSE)

## report U-alphabet motifs, ranked high->low, to match the figure convention
out = enr %>%
  mutate(MOTIF = gsub('T', 'U', MOTIF)) %>%
  arrange(desc(Score)) %>%
  transmute(MOTIF, Score = round(Score, 6))
write_csv(out, paste0(outDir, '2_mouse_5mer_enrichment.csv'))

## reproduction check
top = enr$MOTIF[which.max(enr$Score)]
CS  = suppressMessages(returnIS(enr, top))
cat('rows written :', nrow(out), '(expect 1024)\n')
cat('top 5-mer    :', gsub('T', 'U', top), '(expect GAGGA)\n')
cat('CS           :', round(CS, 3), '(expect 11.588)\n')
