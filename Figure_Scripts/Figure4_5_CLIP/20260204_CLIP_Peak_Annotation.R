################################################################################
## CLIP Peak Normalization, Biological Complexity, and Annotation
## Author: Soon Yi
## Date: February 2026
################################################################################

# === Load Libraries ===
library(data.table)
library(tibble)
library(tidyverse)
library(tidyr)
library(readr)
library(stringr)
library(biomaRt)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(txdbmaker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## 1. Basic Setup:
################################################################################
# BASE_DIR = "F:/Specificity/CLIP/Analysis/"
# INPUT_DIR = "F:/Specificity/CLIP/ALL_POOL_2/"
# OUTPUT_DIR = paste0(BASE_DIR, "output/")
# ANNO_DIR = "F:/Annotations/"

BASE_DIR = "/Volumes/1TB_Data/Specificity/CLIP/Analysis/"
INPUT_DIR = "/Volumes/1TB_Data/Specificity/CLIP/ALL_POOL_2/"
OUTPUT_DIR = paste0(BASE_DIR, "output/")
ANNO_DIR = "/Volumes/1TB_Data/Annotations/"

PEAK_MATRIX_FILE = paste0(INPUT_DIR, "meta_peakMatrix.txt")
DEPTH_FILE = paste0(INPUT_DIR, "meta_peakDepth.txt")
OUTPUT_FILE = paste0(OUTPUT_DIR, "meta_peakMatrix_normalized_annotated.txt")

GTF_FILE = paste0(ANNO_DIR, "gencode.v49.primary_assembly.basic.annotation.gtf")
RMSK_FILE = paste0(ANNO_DIR, "Rmsk_hg38_2026Jan.txt")

rmchr = function(gr) {
  seqlevels(gr) = gsub("^chr", "", seqlevels(gr))
  return(gr)
}
################################################################################

## 2. Load and Normalize Data:
################################################################################
peaksMatrix = fread(PEAK_MATRIX_FILE)
depthsTable = fread(DEPTH_FILE)

# Identify Sample Columns (those starting with TC_)
tcCols = grep("^TC_", names(peaksMatrix), value = TRUE)
sampleIds = gsub("^TC_", "", tcCols)

# --- RPM Normalization ---
normPeaksMatrix = copy(peaksMatrix)

for (i in seq_along(tcCols)) {
  colName = tcCols[i]
  sampleId = sampleIds[i]
  
  # Find corresponding depth
  d = depthsTable[sample == sampleId, depth]
  
  if (length(d) == 1) {
    # Calculate RPM: (TC / Depth) * 1e6
    normPeaksMatrix[[colName]] = (peaksMatrix[[colName]] / d) * 1e6
    
    # Rename column to indicate it is RPM
    setnames(normPeaksMatrix, colName, paste0("nTC_", sampleId))
  } else {
    warning(paste("Depth not found for sample:", sampleId))
  }
}

# --- Biological Complexity (BC Counts) ---
groups = unique(depthsTable$group)

for (grp in groups) {
  # Get samples belonging to this group
  grpSamples = depthsTable[group == grp, sample]
  grpTcCols = paste0("TC_", grpSamples)
  grpNtcCols = paste0("nTC_", grpSamples)
  
  # Filter to columns that actually exist in the matrix
  grpTcCols = intersect(grpTcCols, names(peaksMatrix))
  grpNtcCols = intersect(grpNtcCols, names(normPeaksMatrix))
  
  if (length(grpTcCols) > 0) {
    # BC = count of replicates with TC > 0
    normPeaksMatrix[[paste0("BC_", grp)]] = rowSums(peaksMatrix[, ..grpTcCols] > 0)
  }
  
  if (length(grpNtcCols) > 0) {
    # Calculate average nTC for the group
    normPeaksMatrix[[paste0("nTC_", grp)]] = rowMeans(normPeaksMatrix[, ..grpNtcCols], na.rm = TRUE)
  }
}
################################################################################

## 3. Annotation - Feature Setup:
################################################################################
gtf = import(GTF_FILE)
gtf_df = as.data.frame(gtf)

## Subset the data_frame by gene/transcript types:
protCoding = subset(gtf_df, gene_type == "protein_coding")
miR = subset(gtf_df, transcript_type == "miRNA" & type == "transcript")
lncRNA = subset(gtf_df, transcript_type == "lncRNA" & type == "transcript")
rRNA = subset(gtf_df, transcript_type == "rRNA" & type == "transcript")
snoRNA = subset(gtf_df, transcript_type == "snoRNA" & type == "transcript")
scaRNA = subset(gtf_df, transcript_type == "scaRNA" & type == "transcript")
snRNA = subset(gtf_df, transcript_type == "snRNA" & type == "transcript")
miscRNA = subset(gtf_df, transcript_type == "misc_RNA" & type == "transcript")
Prot_retained_int = subset(gtf_df, transcript_type == "retained_intron" & gene_type == "protein_coding" & type == "transcript")
nc_retained_int = subset(gtf_df, transcript_type == "retained_intron" & gene_type != "protein_coding" & type == "transcript")

## Make Granges for GTF derived categories above:
miR.gr = GRanges(seqnames=miR$seqnames, ranges=IRanges(start=miR$start, end=miR$end, names=miR$gene_name), strand=miR$strand)
lncRNA.gr = GRanges(seqnames=lncRNA$seqnames, ranges=IRanges(start=lncRNA$start, end=lncRNA$end, names=lncRNA$gene_name), strand=lncRNA$strand)
rRNA.gr = GRanges(seqnames=rRNA$seqnames, ranges=IRanges(start=rRNA$start, end=rRNA$end, names=rRNA$gene_name), strand=rRNA$strand)
snoRNA.gr = GRanges(seqnames=snoRNA$seqnames, ranges=IRanges(start=snoRNA$start, end=snoRNA$end, names=snoRNA$gene_name), strand=snoRNA$strand)
scaRNA.gr = GRanges(seqnames=scaRNA$seqnames, ranges=IRanges(start=scaRNA$start, end=scaRNA$end, names=scaRNA$gene_name), strand=scaRNA$strand)
snRNA.gr = GRanges(seqnames=snRNA$seqnames, ranges=IRanges(start=snRNA$start, end=snRNA$end, names=snRNA$gene_name), strand=snRNA$strand)
miscRNA.gr = GRanges(seqnames=miscRNA$seqnames, ranges=IRanges(start=miscRNA$start, end=miscRNA$end, names=miscRNA$gene_name), strand=miscRNA$strand)
Prot_retained_int.gr = GRanges(seqnames=Prot_retained_int$seqnames, ranges=IRanges(start=Prot_retained_int$start, end=Prot_retained_int$end, names=Prot_retained_int$gene_name), strand=Prot_retained_int$strand)
nc_retained_int.gr = GRanges(seqnames=nc_retained_int$seqnames, ranges=IRanges(start=nc_retained_int$start, end=nc_retained_int$end, names=nc_retained_int$gene_name), strand=nc_retained_int$strand)

## Use the same GTF to make txdb
txdb = makeTxDbFromGFF(GTF_FILE, format="gtf")
hg38.tx = transcriptsBy(txdb, "gene")

## Extract features of interest from hg38 TxDb:
exons = unique(unlist(exonsBy(txdb, "tx", use.names=T)))    
CDS = unique(unlist(cdsBy(txdb, "tx", use.names=T)))
introns = unique(unlist(intronsByTranscript(txdb, use.names=T)))
fiveUTRs = unique(unlist(fiveUTRsByTranscript(txdb,use.names=T)))
threeUTRs = unique(unlist(threeUTRsByTranscript(txdb,use.names=T)))

## RepeatMasker features:
repMask_Main = read.delim(RMSK_FILE)
repMask_Main = repMask_Main[, c('genoName', 'genoStart', 'genoEnd', 'strand', 'repName', 'repClass', 'repFamily')]
colnames(repMask_Main) = c('chr', 'start', 'end', 'strand', 'name', 'repClass', 'repFamily')

tRNA.bed = subset(repMask_Main, repClass == 'tRNA')
LINE.bed = subset(repMask_Main, repClass == 'LINE')
LC_SR.bed = subset(repMask_Main, (repClass == 'Low_complexity' | repClass == 'Simple_repeat'))
LTR.bed = subset(repMask_Main, repClass == 'LTR')
Satellite.bed = subset(repMask_Main, repClass == 'Satellite')
SINE.bed = subset(repMask_Main, repClass == 'SINE')

tRNA.gr = GRanges(seqnames=tRNA.bed$chr, ranges=IRanges(start=tRNA.bed$start, end=tRNA.bed$end, names=tRNA.bed$name), strand=tRNA.bed$strand)
LINE.gr = GRanges(seqnames=LINE.bed$chr, ranges=IRanges(start=LINE.bed$start, end=LINE.bed$end, names=LINE.bed$name), strand=LINE.bed$strand)
LC_SR.gr = GRanges(seqnames=LC_SR.bed$chr, ranges=IRanges(start=LC_SR.bed$start, end=LC_SR.bed$end, names=LC_SR.bed$name), strand=LC_SR.bed$strand)
LTR.gr = GRanges(seqnames=LTR.bed$chr, ranges=IRanges(start=LTR.bed$start, end=LTR.bed$end, names=LTR.bed$name), strand=LTR.bed$strand)
Satellite.gr = GRanges(seqnames=Satellite.bed$chr, ranges=IRanges(start=Satellite.bed$start, end=Satellite.bed$end, names=Satellite.bed$name), strand=Satellite.bed$strand)
SINE.gr = GRanges(seqnames=SINE.bed$chr, ranges=IRanges(start=SINE.bed$start, end=SINE.bed$end, names=SINE.bed$name), strand=SINE.bed$strand)

## Extract gene regions for annotation of genes:
genes = genes(txdb)
genes.bed = data.frame(chr=seqnames(genes), start=start(genes)-1, end=end(genes), name=unlist(genes$gene_id), score=0, strand=strand(genes))
################################################################################

## 4. Annotation - Overlap Detection:
################################################################################
## Create GRanges object for peaks
peaksGR = GRanges(seqnames=normPeaksMatrix$chr, ranges=IRanges(start=normPeaksMatrix$start, end=normPeaksMatrix$end, names=normPeaksMatrix$name), strand=normPeaksMatrix$strand)

## Overlap with features:
fiveUTRs.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=fiveUTRs,  minoverlap=1, select="first"))
threeUTRs.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=threeUTRs,  minoverlap=1, select="first"))
CDS.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=CDS,  minoverlap=1, select="first"))
introns.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=introns,  minoverlap=1, select="first"))
genes.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=genes,  minoverlap=1, select="first"))
tRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=tRNA.gr,  minoverlap=1, select="first"))
LINE.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=LINE.gr,  minoverlap=1, select="first"))
LC_SR.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=LC_SR.gr,  minoverlap=1, select="first"))
LTR.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=LTR.gr,  minoverlap=1, select="first"))
Satellite.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=Satellite.gr,  minoverlap=1, select="first"))
SINE.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=SINE.gr,  minoverlap=1, select="first"))
miR.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=miR.gr,  minoverlap=1, select="first"))
lncRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=lncRNA.gr,  minoverlap=1, select="first"))
rRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=rRNA.gr,  minoverlap=1, select="first"))
snoRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=snoRNA.gr,  minoverlap=1, select="first"))
scaRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=scaRNA.gr,  minoverlap=1, select="first"))
snRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=snRNA.gr,  minoverlap=1, select="first"))
miscRNA.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=miscRNA.gr,  minoverlap=1, select="first"))
Prot_retained_int.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=Prot_retained_int.gr,  minoverlap=1, select="first"))
nc_retained_int.peaks = as.data.frame(findOverlaps(query=peaksGR, subject=nc_retained_int.gr,  minoverlap=1, select="first"))

## Assign feature data to matrix:

normPeaksMatrix$fiveUTRs = ifelse(is.na(fiveUTRs.peaks[,1]), NA, "5'UTR")
normPeaksMatrix$threeUTRs = ifelse(is.na(threeUTRs.peaks[,1]), NA, "3'UTR")
normPeaksMatrix$CDS = ifelse(is.na(CDS.peaks[,1]), NA, "CDS")
normPeaksMatrix$introns = ifelse(is.na(introns.peaks[,1]), NA, "intron")
normPeaksMatrix$gene = genes.bed[genes.peaks[,1], "name"]

normPeaksMatrix$tRNA1 = ifelse(is.na(tRNA.peaks[,1]), NA, "tRNA")
normPeaksMatrix$LINE1 = ifelse(is.na(LINE.peaks[,1]), NA, "TE")
normPeaksMatrix$LTR1 = ifelse(is.na(LTR.peaks[,1]), NA, "Other")
normPeaksMatrix$LC_SR1 = ifelse(is.na(LC_SR.peaks[,1]), NA, "Other")
normPeaksMatrix$Satellite1 = ifelse(is.na(Satellite.peaks[,1]), NA, "Other")
normPeaksMatrix$SINE1 = ifelse(is.na(SINE.peaks[,1]), NA, "TE")
normPeaksMatrix$miR1 = ifelse(is.na(miR.peaks[,1]), NA, "miRNA")
normPeaksMatrix$lncRNA1 = ifelse(is.na(lncRNA.peaks[,1]), NA, "lncRNA")
normPeaksMatrix$rRNA1 = ifelse(is.na(rRNA.peaks[,1]), NA, "rRNA")
normPeaksMatrix$snoRNA1 = ifelse(is.na(snoRNA.peaks[,1]), NA, "snoRNA")
normPeaksMatrix$scaRNA1 = ifelse(is.na(scaRNA.peaks[,1]), NA, "scaRNA")
normPeaksMatrix$snRNA1 = ifelse(is.na(snRNA.peaks[,1]), NA, "snRNA")
normPeaksMatrix$miscRNA1 = ifelse(is.na(miscRNA.peaks[,1]), NA, "Other")
normPeaksMatrix$Prot_retained_int1 = ifelse(is.na(Prot_retained_int.peaks[,1]), NA, "CDS_Retained_intron")
normPeaksMatrix$nc_retained_int1 = ifelse(is.na(nc_retained_int.peaks[,1]), NA, "ncRNA_Retained_intron")
################################################################################

## 5. Annotation - Consolidation:
################################################################################
# Columns to collapse
annoCols = c("fiveUTRs", "threeUTRs", "CDS", "introns", "tRNA1", "LTR1", 
            "LINE1", "SINE1", "Satellite1", "LC_SR1", "miR1", "lncRNA1", 
            "rRNA1", "snoRNA1", "scaRNA1", "snRNA1", "miscRNA1", 
            "Prot_retained_int1", "nc_retained_int1")

# Collapse all overlapping features into a single column using vectorized unite
normPeaksMatrix = normPeaksMatrix %>%
  unite("annotation", all_of(annoCols), sep = "|", na.rm = TRUE, remove = FALSE)

# Priority list
priority_list = c("5'UTR", "CDS", "3'UTR", "intron",
                  "miRNA", "lncRNA",  "rRNA", "snoRNA",  "scaRNA",  "snRNA",  "miscRNA",  "tRNA",  "TE",  "Other",  
                  "CDS_Retained_intron", "ncRNA_Retained_intron", 
                  "unannotated")

# Extract the highest priority term for finalized_annotation
normPeaksMatrix = normPeaksMatrix %>%
  mutate(finalized_annotation = sapply(strsplit(annotation, "\\|"), function(terms) {
    for (term in priority_list) {
      if (term %in% terms) {
        return(term)
      }
    }
    return(NA)
  }))

# Map terms to grouped_annotation
normPeaksMatrix = normPeaksMatrix %>%
  mutate(grouped_annotation = ifelse(finalized_annotation %in% c("5'UTR", "CDS", "3'UTR", "intron", "Other", "unannotated"), finalized_annotation,
                                     ifelse(finalized_annotation %in% c("miRNA", "lncRNA", "rRNA", "scaRNA", "snRNA", "miscRNA", "tRNA", "snoRNA"), "ncRNA",
                                            ifelse(finalized_annotation %in% c("CDS_Retained_intron", "ncRNA_Retained_intron"), "retained_intron", "unannotated"))),
         annotation_count = sapply(strsplit(annotation, "\\|"), length))

normPeaksMatrix[normPeaksMatrix == '' | normPeaksMatrix == 'NA'] = NA

# Cleanup intermediate annotation columns
drops = c('fiveUTRs', 'threeUTRs', 'CDS', 'introns',
          'tRNA1', 'LINE1', 'LTR1', 'LC_SR1', 'Satellite1', 'SINE1',
          'miR1', 'lncRNA1', 'rRNA1', 'snoRNA1', 'scaRNA1', 'snRNA1', 'miscRNA1',
          'Prot_retained_int1', 'nc_retained_int1')
setDT(normPeaksMatrix); normPeaksMatrix[, (drops) := NULL]

## Use biomaRt for gene mapping:
mart.hs = useMart("ensembl", host = "https://useast.ensembl.org", dataset="hsapiens_gene_ensembl")
gene_names = getBM(attributes = c("ensembl_gene_id_version", "external_gene_name"),
                   filters = "ensembl_gene_id_version",
                   values = normPeaksMatrix$gene,
                   mart = mart.hs)

normPeaksMatrix = normPeaksMatrix %>% left_join(gene_names, by = c("gene" = "ensembl_gene_id_version"), relationship = "many-to-many")
normPeaksMatrix = normPeaksMatrix %>% mutate(finalized_annotation = ifelse(is.na(finalized_annotation), 'unannotated', finalized_annotation))
################################################################################

## 6. Save Results:
################################################################################
fwrite(normPeaksMatrix, OUTPUT_FILE, sep = "\t")
################################################################################
