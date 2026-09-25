################################################################################
## Task E: genomic region distribution of EIF4G2 CLIP peaks (Reviewer 2, comment 1)
##
## PART 1: human ENCODE EIF4G2 K562 eCLIP IDR peaks (ENCFF207BXK, both reps merged).
## PART 2: mouse eIF4G2 HITS-CLIP resting peaks (Hacisuleyman 2024, GSE213082),
##         processed on the user's Ubuntu box with CLIPittyClip (single-end, no UMI,
##         --no-dedup) against the mm39 STAR index -> COMBINED peak set over 6 resting
##         replicates (FINAL_COMBINED_PEAKS.bed, 1237 peaks).
##
## Classes match Ext. Data Fig. 7e of Hacisuleyman 2024: 5'UTR, CDS, 3'UTR, intron,
## upstream 10kb, downstream 10kb. Annotation = the manuscript's own method (GENCODE
## TxDb feature overlap, minoverlap=1) extended with strand-aware 10kb flanks;
## hierarchical priority 5UTR > CDS > 3UTR > intron > upstream_10kb > downstream_10kb >
## intergenic; peaks assigned by full range. Human uses GENCODE v49 (hg38), mouse uses
## GENCODE vM38 (mm39) -- each the build its peaks were called on.
################################################################################

suppressPackageStartupMessages({
  library(GenomicFeatures); library(AnnotationDbi); library(GenomicRanges)
  library(data.table); library(dplyr); library(readr)
})
options(width = 150)

outDir  = '/Users/soonyi/Repos/BITS_Specificity/REVISION_2026-09/'
PRIORITY = c('5UTR', 'CDS', '3UTR', 'intron', 'upstream_10kb', 'downstream_10kb')

build_features = function(txdb, std) {
  seqlevels(txdb, pruning.mode = 'coarse') = intersect(seqlevels(txdb), std)
  genes = genes(txdb, single.strand.genes.only = TRUE)
  list(
    '5UTR'            = unique(unlist(fiveUTRsByTranscript(txdb, use.names = TRUE))),
    'CDS'             = unique(unlist(cdsBy(txdb, 'tx', use.names = TRUE))),
    '3UTR'            = unique(unlist(threeUTRsByTranscript(txdb, use.names = TRUE))),
    'intron'          = unique(unlist(intronsByTranscript(txdb, use.names = TRUE))),
    'upstream_10kb'   = trim(flank(genes, width = 10000, start = TRUE)),
    'downstream_10kb' = trim(flank(genes, width = 10000, start = FALSE))
  )
}

annotate_peaks = function(bed, label, feat, std) {
  P = fread(bed, select = 1:6)
  setnames(P, c('chr', 'start', 'end', 'name', 'score', 'strand'))
  P = P[chr %in% std]
  gr = GRanges(P$chr, IRanges(P$start + 1, P$end), strand = P$strand)  # BED 0-based -> 1-based
  hit = as.data.frame(lapply(feat, function(f)
    !is.na(GenomicRanges::findOverlaps(gr, f, minoverlap = 1L, select = 'first', ignore.strand = FALSE))))
  colnames(hit) = names(feat)
  cls = apply(hit, 1, function(r) { w = PRIORITY[PRIORITY %in% names(which(r))]; if (length(w)) w[1] else 'intergenic' })
  tab = table(factor(cls, levels = c(PRIORITY, 'intergenic')))
  data.frame(dataset = label, region = names(tab), n = as.integer(tab),
             fraction = round(as.integer(tab) / length(gr), 4), n_peaks_total = length(gr))
}

## PART 1: human
################################################################################
std_h = paste0('chr', c(1:22, 'X', 'Y', 'M'))
txdb_h = loadDb('/Users/soonyi/Desktop/Genomics/Annotations/Human/gencode.v49.txdb.sqlite')
feat_h = build_features(txdb_h, std_h)
human = annotate_peaks(
  '/Users/soonyi/Repos/BITS_Specificity/Dataset/Analysis/eCLIP_peak/raw/ENCFF207BXK.bed',
  'human_EIF4G2_K562_eCLIP', feat_h, std_h)

## PART 2: mouse
################################################################################
std_m = paste0('chr', c(1:19, 'X', 'Y', 'M'))
txdb_m = loadDb('/Volumes/1TB_Data/Annotations/Mouse/gencode.vM38.txdb.sqlite')
feat_m = build_features(txdb_m, std_m)
mouse = annotate_peaks(
  '/Users/soonyi/Desktop/Specificity_Genome_Biology_Resubmission/EIF4G2_HITS_CLIP/mouse_EIF4G2_HITS_CLIP_fastq_output/04_PEAKS/COMBINED_PEAKS/FINAL_COMBINED_PEAKS.bed',
  'mouse_eIF4G2_HITSCLIP_resting', feat_m, std_m)

## PART 2b: mouse TAG-level distribution (fraction of aligned tags per region).
## The paper's "~50% in 5'UTRs" and its deposited "CLIP_tags_on_5UTRs" files are tag-level,
## not peak-level, so this is the directly comparable basis.
################################################################################
tagbeds = Sys.glob('/Users/soonyi/Desktop/Specificity_Genome_Biology_Resubmission/EIF4G2_HITS_CLIP/mouse_EIF4G2_HITS_CLIP_fastq_output/02_COLLAPSED_BED/*.bed')
pooled = rbindlist(lapply(tagbeds, function(f) fread(f, select = 1:6)))
tmp = tempfile(fileext = '.bed'); fwrite(pooled, tmp, sep = '\t', col.names = FALSE)
mouse_tags = annotate_peaks(tmp, 'mouse_eIF4G2_HITSCLIP_resting_TAGS', feat_m, std_m)

## Combined table
################################################################################
res = bind_rows(human, mouse, mouse_tags)
write_csv(res, paste0(outDir, '2_mouse_region_distribution.csv'), na = '')

fmt = function(d) { w = reshape(d[, c('region','fraction')], idvar='region', timevar=NULL, direction='wide'); d }
cat('=== EIF4G2 CLIP region distribution: human ENCODE eCLIP vs mouse HITS-CLIP ===\n')
cat('priority 5UTR>CDS>3UTR>intron>upstream10kb>downstream10kb>intergenic; peaks by full range\n\n')
wide = merge(human[,c('region','fraction')], mouse[,c('region','fraction')], by='region', suffixes=c('_human','_mouse'))
wide = wide[match(c(PRIORITY,'intergenic'), wide$region), ]
wide$human_pct = sprintf('%.1f%%', 100*wide$fraction_human)
wide$mouse_pct = sprintf('%.1f%%', 100*wide$fraction_mouse)
print(wide[, c('region','human_pct','mouse_pct')], row.names = FALSE)
cat('\nhuman peaks:', human$n_peaks_total[1], ' (ENCFF207BXK, K562 IDR, both reps)\n')
cat('mouse peaks:', mouse$n_peaks_total[1], ' (COMBINED over 6 resting replicates)\n')
cat('Paper (Hacisuleyman 2024) reports ~50% of mouse eIF4G2 peaks in 5UTRs.\n')
