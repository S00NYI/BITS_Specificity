################################################################################
## CLIP Metagene Analysis
## Author: Soon Yi with Antigravity
## Date: February 2026
################################################################################

# === Load Libraries ===
library(dplyr)
library(tidyr)
library(data.table)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(stringr)

## 1. Basic Setup:
################################################################################
# BASE_DIR = "F:/Specificity/CLIP/Analysis/"
# INPUT_DIR = paste0(BASE_DIR, "output/")
# OUTPUT_DIR = paste0(BASE_DIR, "output/metagene/")
# ANNOTATION_DIR = "F:/Annotations/"

BASE_DIR = "/Volumes/1TB_Data/Specificity/CLIP/Analysis/"
INPUT_DIR = paste0(BASE_DIR, "output/")
OUTPUT_DIR = paste0(BASE_DIR, "output/metagene/")
ANNOTATION_DIR = "/Volumes/1TB_Data/Annotations/"

GTF_FILE = paste0(ANNOTATION_DIR, "Ensembl_Homo_sapiens.GRCh38.115.gtf")

genome = BSgenome.Hsapiens.UCSC.hg38

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

CLIP_List = c('hnRNPC_WT', 'hnRNPC_Mut', 'hnRNPC_WT_inRBM25_WT', 
              'hnRNPC_WT_inRBM25_Mut', 'hnRNPC_WT_inKD', 'hnRNPC_Mut_inKD',
              'RBM25_WT', 'RBM25_Mut')
################################################################################

## 2. Custom Functions:
################################################################################
# Normalize densities by base-correcting.
normalize_density = function(density_data) {
  density_data$density = density_data$density - min(density_data$density)
  return(density_data)
}

# Function to calculate peak densities for metagene:
metaDensity = function(feature_df, peaks_df, feature_center, window_width = 20, window_definition = 500) {
  feature_df_sense = feature_df %>% filter(strand == '+')
  feature_df_antisense = feature_df %>% filter(strand == '-')
  
  peaks_df_sense = peaks_df %>% filter(strand == '+')
  peaks_df_antisense = peaks_df %>% filter(strand == '-')
  
  if (feature_center == 'start') {
    feature_center_sense = 'start'
    feature_center_antisense = 'end'
  } else if (feature_center == 'end') {
    feature_center_sense = 'end'
    feature_center_antisense = 'start'
  }
  
  ## Make bins for density calculation:  
  bins = list()
  bin_counts = floor(window_definition / window_width)
  for (i in -(bin_counts):(bin_counts-1)) {
    start = 0 + (i * window_width)
    end = start + window_width - 1
    bins[[i + bin_counts + 1]] = c(start, end)
  }
  
  ## For Sense Strands:
  counts_sense = matrix(0, nrow = 1, ncol = length(bins))
  
  if (nrow(feature_df_sense) > 0) {
    chroms = feature_df_sense$seqid
    centers = unlist(feature_df_sense[, feature_center_sense])
    strands = feature_df_sense$strand
    
    peaksGR = GRanges(seqnames = peaks_df_sense$chr, 
                      ranges = IRanges(start = peaks_df_sense$start, 
                                       end = peaks_df_sense$end), 
                      strand = peaks_df_sense$strand)
    
    for (bin_id in 1:length(bins)) {
      bin = bins[bin_id]
      bin_Ranges = matrix(0, nrow = nrow(feature_df_sense), ncol = 4)
      for (center_id in 1:length(centers)) {
        bin_Ranges[center_id, 1] = chroms[center_id]
        bin_Ranges[center_id, 2] = centers[center_id] + bin[[1]][1]
        bin_Ranges[center_id, 3] = centers[center_id] + bin[[1]][2]
        bin_Ranges[center_id, 4] = strands[center_id]
      }
      bin_Ranges = data.frame(bin_Ranges)
      colnames(bin_Ranges) = c('chrom', 'start', 'end', 'strand')
      bin_Ranges_GR = GRanges(bin_Ranges)
      
      counts_sense[1, bin_id] = sum(countOverlaps(peaksGR, bin_Ranges_GR, minoverlap = 10))
    }
  }
  
  
  ## For AntiSense Strands:
  counts_antisense = matrix(0, nrow = 1, ncol = length(bins))
  
  if (nrow(feature_df_antisense) > 0) {
    chroms = feature_df_antisense$seqid
    centers = unlist(feature_df_antisense[, feature_center_antisense])
    strands = feature_df_antisense$strand
    
    peaksGR = GRanges(seqnames = peaks_df_antisense$chr, 
                      ranges = IRanges(start = peaks_df_antisense$start, 
                                       end = peaks_df_antisense$end), 
                      strand = peaks_df_antisense$strand)
    
    bins_antisense = rev(lapply(bins, function(bin) -(bin + 1)))
    
    for (bin_id in 1:length(bins)) {
      bin = bins_antisense[bin_id]
      bin_Ranges = matrix(0, nrow = nrow(feature_df_antisense), ncol = 4)
      for (center_id in 1:length(centers)) {
        bin_Ranges[center_id, 1] = chroms[center_id]
        bin_Ranges[center_id, 2] = centers[center_id] - bin[[1]][1]
        bin_Ranges[center_id, 3] = centers[center_id] - bin[[1]][2]
        bin_Ranges[center_id, 4] = strands[center_id]
      }
      bin_Ranges = data.frame(bin_Ranges)
      colnames(bin_Ranges) = c('chrom', 'start', 'end', 'strand')
      bin_Ranges_GR = GRanges(bin_Ranges)
      
      counts_antisense[1, bin_id] = sum(countOverlaps(peaksGR, bin_Ranges_GR, minoverlap = 10))
    }
  }
  
  counts = data.frame(sense = t(data.frame(counts_sense)), 
                      antisense = t(data.frame(counts_antisense)))
  
  counts$both = counts$sense + counts$antisense
  
  counts$density = counts$both
  
  counts$midpoint = unlist(lapply(bins, function(bin) mean(bin)))
  
  counts = counts[, c('midpoint', 'sense', 'antisense', 'both', 'density')]
  
  return(counts)
}

# Function to count how many peaks were in search space:
countPeaks = function(feature_df, peaks_df, feature_center, window_definition = 500) {
  
  if (feature_center == 'start') {
    feature_GR = GRanges(seqnames = feature_df$seqid, 
                         ranges = IRanges(start = feature_df$start - window_definition, 
                                          end = feature_df$start + window_definition), 
                         strand = feature_df$strand)
  } else if (feature_center == 'end') {
    feature_GR = GRanges(seqnames = feature_df$seqid, 
                         ranges = IRanges(start = feature_df$end - window_definition, 
                                          end = feature_df$end + window_definition), 
                         strand = feature_df$strand)
  }
  
  peaksGR = GRanges(seqnames = peaks_df$chr, 
                    ranges = IRanges(start = peaks_df$start, 
                                     end = peaks_df$end), 
                    strand = peaks_df$strand)
  
  counts = sum(countOverlaps(feature_GR, peaksGR))
  
  return(counts)
}

# Function to compute nucleotide content around genomic features (vectorized):
nucleotideContentAroundFeature = function(feature_df, genome, target_nucleotide = "T", 
                                          feature_center = 'start',
                                          window_width = 20, window_definition = 500) {
  
  feature_df_sense = feature_df %>% filter(strand == '+')
  feature_df_antisense = feature_df %>% filter(strand == '-')
  
  ## Make bins:
  bins = list()
  bin_counts = floor(window_definition / window_width)
  for (i in -(bin_counts):(bin_counts-1)) {
    start = 0 + (i * window_width)
    end = start + window_width - 1
    bins[[i + bin_counts + 1]] = c(start, end)
  }
  
  midpoints = unlist(lapply(bins, function(bin) mean(bin)))
  n_bins = length(bins)
  bin_starts = sapply(bins, function(b) b[1])
  bin_ends = sapply(bins, function(b) b[2])
  fraction_totals = rep(0, n_bins)
  feature_count = 0
  
  ## Complement mapping for antisense (counting complement = counting target in revcomp)
  comp_map = c("A" = "T", "T" = "A", "G" = "C", "C" = "G")
  complement_nucleotide = comp_map[target_nucleotide]
  
  ## Determine which column to use as center based on feature_center and strand:
  if (feature_center == 'start') {
    sense_center_col = 'start'
    antisense_center_col = 'end'
  } else if (feature_center == 'end') {
    sense_center_col = 'end'
    antisense_center_col = 'start'
  }
  
  ## Sense strand features: vectorized per-bin
  if (nrow(feature_df_sense) > 0) {
    chroms = as.character(feature_df_sense$seqid)
    centers = feature_df_sense[[sense_center_col]]
    
    ## Filter to features on chromosomes present in genome
    valid_chr = chroms %in% seqnames(genome)
    chroms = chroms[valid_chr]
    centers = centers[valid_chr]
    chr_lens = seqlengths(genome)[chroms]
    
    if (length(centers) > 0) {
      for (bin_id in 1:n_bins) {
        all_starts = centers + bin_starts[bin_id]
        all_ends = centers + bin_ends[bin_id]
        
        valid = all_starts >= 1 & all_ends <= chr_lens
        
        if (sum(valid) > 0) {
          gr = GRanges(seqnames = chroms[valid],
                       ranges = IRanges(start = all_starts[valid], end = all_ends[valid]))
          seqs = getSeq(genome, gr)
          freqs = letterFrequency(seqs, letters = target_nucleotide, as.prob = TRUE)[, 1]
          fraction_totals[bin_id] = fraction_totals[bin_id] + sum(freqs)
        }
      }
      feature_count = feature_count + length(centers)
    }
  }
  
  ## Antisense strand features: vectorized per-bin
  if (nrow(feature_df_antisense) > 0) {
    chroms = as.character(feature_df_antisense$seqid)
    centers = feature_df_antisense[[antisense_center_col]]
    
    valid_chr = chroms %in% seqnames(genome)
    chroms = chroms[valid_chr]
    centers = centers[valid_chr]
    chr_lens = seqlengths(genome)[chroms]
    
    ## Reverse bins for antisense
    bins_antisense = rev(lapply(bins, function(bin) -(bin + 1)))
    bin_starts_as = sapply(bins_antisense, function(b) b[1])
    bin_ends_as = sapply(bins_antisense, function(b) b[2])
    
    if (length(centers) > 0) {
      for (bin_id in 1:n_bins) {
        all_raw_starts = centers - bin_starts_as[bin_id]
        all_raw_ends = centers - bin_ends_as[bin_id]
        all_starts = pmin(all_raw_starts, all_raw_ends)
        all_ends = pmax(all_raw_starts, all_raw_ends)
        
        valid = all_starts >= 1 & all_ends <= chr_lens
        
        if (sum(valid) > 0) {
          gr = GRanges(seqnames = chroms[valid],
                       ranges = IRanges(start = all_starts[valid], end = all_ends[valid]))
          seqs = getSeq(genome, gr)
          ## Count complement nucleotide (equivalent to target in reverse complement)
          freqs = letterFrequency(seqs, letters = complement_nucleotide, as.prob = TRUE)[, 1]
          fraction_totals[bin_id] = fraction_totals[bin_id] + sum(freqs)
        }
      }
      feature_count = feature_count + length(centers)
    }
  }
  
  ## Average across all features:
  result = data.frame(
    midpoint = midpoints,
    fraction = fraction_totals / feature_count
  )
  
  return(result)
}

# Function to compute nucleotide content around peak centers (vectorized):
nucleotideContentAroundPeaks = function(peaks_df, genome, target_nucleotide = "T", 
                                         window_width = 20, window_definition = 500) {
  
  peaks_df_sense = peaks_df %>% filter(strand == '+')
  peaks_df_antisense = peaks_df %>% filter(strand == '-')
  
  ## Make bins:
  bins = list()
  bin_counts = floor(window_definition / window_width)
  for (i in -(bin_counts):(bin_counts-1)) {
    start = 0 + (i * window_width)
    end = start + window_width - 1
    bins[[i + bin_counts + 1]] = c(start, end)
  }
  
  midpoints = unlist(lapply(bins, function(bin) mean(bin)))
  n_bins = length(bins)
  bin_starts = sapply(bins, function(b) b[1])
  bin_ends = sapply(bins, function(b) b[2])
  fraction_totals = rep(0, n_bins)
  peak_count = 0
  
  ## Complement mapping for antisense
  comp_map = c("A" = "T", "T" = "A", "G" = "C", "C" = "G")
  complement_nucleotide = comp_map[target_nucleotide]
  
  ## Sense strand peaks: vectorized per-bin
  if (nrow(peaks_df_sense) > 0) {
    chroms = as.character(peaks_df_sense$chr)
    centers = floor((peaks_df_sense$start + peaks_df_sense$end) / 2)
    
    valid_chr = chroms %in% seqnames(genome)
    chroms = chroms[valid_chr]
    centers = centers[valid_chr]
    chr_lens = seqlengths(genome)[chroms]
    
    if (length(centers) > 0) {
      for (bin_id in 1:n_bins) {
        all_starts = centers + bin_starts[bin_id]
        all_ends = centers + bin_ends[bin_id]
        
        valid = all_starts >= 1 & all_ends <= chr_lens
        
        if (sum(valid) > 0) {
          gr = GRanges(seqnames = chroms[valid],
                       ranges = IRanges(start = all_starts[valid], end = all_ends[valid]))
          seqs = getSeq(genome, gr)
          freqs = letterFrequency(seqs, letters = target_nucleotide, as.prob = TRUE)[, 1]
          fraction_totals[bin_id] = fraction_totals[bin_id] + sum(freqs)
        }
      }
      peak_count = peak_count + length(centers)
    }
  }
  
  ## Antisense strand peaks: vectorized per-bin
  if (nrow(peaks_df_antisense) > 0) {
    chroms = as.character(peaks_df_antisense$chr)
    centers = floor((peaks_df_antisense$start + peaks_df_antisense$end) / 2)
    
    valid_chr = chroms %in% seqnames(genome)
    chroms = chroms[valid_chr]
    centers = centers[valid_chr]
    chr_lens = seqlengths(genome)[chroms]
    
    bins_antisense = rev(lapply(bins, function(bin) -(bin + 1)))
    bin_starts_as = sapply(bins_antisense, function(b) b[1])
    bin_ends_as = sapply(bins_antisense, function(b) b[2])
    
    if (length(centers) > 0) {
      for (bin_id in 1:n_bins) {
        all_raw_starts = centers - bin_starts_as[bin_id]
        all_raw_ends = centers - bin_ends_as[bin_id]
        all_starts = pmin(all_raw_starts, all_raw_ends)
        all_ends = pmax(all_raw_starts, all_raw_ends)
        
        valid = all_starts >= 1 & all_ends <= chr_lens
        
        if (sum(valid) > 0) {
          gr = GRanges(seqnames = chroms[valid],
                       ranges = IRanges(start = all_starts[valid], end = all_ends[valid]))
          seqs = getSeq(genome, gr)
          freqs = letterFrequency(seqs, letters = complement_nucleotide, as.prob = TRUE)[, 1]
          fraction_totals[bin_id] = fraction_totals[bin_id] + sum(freqs)
        }
      }
      peak_count = peak_count + length(centers)
    }
  }
  
  ## Average across all peaks:
  result = data.frame(
    midpoint = midpoints,
    fraction = fraction_totals / peak_count
  )
  
  return(result)
}

## Plot density plot:
plot_Density = function(density_data, columns_list, xaxis_lims = NULL, yaxis_lims = NULL, custom_colors = NULL, densityType = NULL, sampleName = NULL, featureName = NULL, smoothing = NULL) {
  # Select the columns based on the columns_list
  plot_data = density_data[, c('position', columns_list)]
  
  # Create a long-format dataframe for better legend handling
  plot_data_long = plot_data %>%
    pivot_longer(cols = {{columns_list}}, names_to = "Data", values_to = "Density")
  
  # Define the order of levels for the Data factor
  plot_data_long$Data = factor(plot_data_long$Data, levels = columns_list)
  
  # Create the ggplot object
  plot = ggplot(plot_data_long, aes(x = position, y = Density, color = Data)) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = 'bold')) +
    labs(y = "Peak Density")
  
  # Add smoothed lines
  if (!is.null(smoothing)) {
    plot = plot + geom_smooth(span = smoothing, se = F, level = 0.9)
  } else {
    plot = plot + geom_line(linewidth = 1)
  }

  if (!is.null(custom_colors)) {
    plot = plot + scale_color_manual(values = custom_colors)
  }
  
  if (!is.null(xaxis_lims)) {
    plot = plot + xlim(xaxis_lims)
  }
  
  if (!is.null(yaxis_lims)) {
    plot = plot + ylim(yaxis_lims)
  }
  
  if (is.null(densityType)) {
    if (!is.null(sampleName)) {
      plot = plot + labs(title = paste0('Density for ', sampleName),
                         x = "Distance from center (nucleotides)")
      
    } else {
      plot = plot + labs(title = 'Density Plot',
                         x = "Distance from center (nucleotides)")
    }
    
  } else if (densityType == 'motif_density') {
    if (!is.null(sampleName)) {
      plot = plot + labs(title = paste0('Motif Density around Peaks for ', sampleName),
                         x = "Distance from peak center (nucleotides)")
    } else {
      plot = plot + labs(title = 'Motif Density around Peaks',
                         x = "Distance from peak center (nucleotides)")
    }
  } else if (densityType == 'feature_metagene') {
    if (!is.null(sampleName)) {
      if (!is.null(featureName)) {
        plot = plot + labs(title = paste0('Peaks density around ', featureName, ' for ', sampleName),
                           x = paste0("Distance from ", featureName, " (nucleotides)"))
      } else {
        plot = plot + labs(title = 'Peaks density around feature',
                           x = "Distance from feature (nucleotides)")
      }
    } else {
      if (!is.null(featureName)) {
        plot = plot + labs(title = paste0('Peaks density around ', featureName),
                           x = paste0("Distance from ", featureName, " (nucleotides)"))
      } else {
        plot = plot + labs(title = 'Peaks density around feature',
                           x = "Distance from feature (nucleotides)")
      }
    }
  } else if (densityType == 'nucleotide_content') {
    if (!is.null(sampleName)) {
      if (!is.null(featureName)) {
        plot = plot + labs(title = paste0(featureName, ' for ', sampleName),
                           x = "Distance from center (nucleotides)",
                           y = "Nucleotide Fraction")
      } else {
        plot = plot + labs(title = paste0('Nucleotide Content for ', sampleName),
                           x = "Distance from center (nucleotides)",
                           y = "Nucleotide Fraction")
      }
    } else {
      if (!is.null(featureName)) {
        plot = plot + labs(title = featureName,
                           x = "Distance from center (nucleotides)",
                           y = "Nucleotide Fraction")
      } else {
        plot = plot + labs(title = 'Nucleotide Content',
                           x = "Distance from center (nucleotides)",
                           y = "Nucleotide Fraction")
      }
    }
  }
  return(plot)
}
################################################################################

## 3. Load Annotations & Define RNA Features:
################################################################################
# Read in GTF file (Ensembl GRCh38.115, needs chr prefix added)
gtf_raw = readGFF(GTF_FILE)
gtf_raw$seqid = paste0("chr", gtf_raw$seqid)

# Filter out gene-level entries (keep transcripts and sub-features), focus on protein coding
gtf_raw = gtf_raw[gtf_raw$type != "gene", ]
gtf = subset(gtf_raw, transcript_biotype == "protein_coding")

## Filter to longest transcript with MANE_Select tag:
gtf_transcripts = subset(gtf, type == 'transcript')

## Gencode uses 'tag' column which can contain multiple tags; filter for MANE_Select
# gtf_transcripts = gtf_transcripts[grep("MANE_Select", gtf_transcripts$tag), ]
gtf_transcripts = gtf_transcripts %>% filter(gene_source == "ensembl_havana")
gtf_transcripts$length = gtf_transcripts$end - gtf_transcripts$start
gtf_transcripts = gtf_transcripts %>% group_by(gene_id) %>% filter(length == max(length)) %>% slice_head(n = 1)

nrow(gtf_transcripts) == length(unique(gtf_transcripts$gene_id))

gtf = gtf[gtf$transcript_id %in% unique(gtf_transcripts$transcript_id), ]
################################################################################

## Define RNA Features:
## TSS (Transcription Start Sites)
## CDS (Coding Start sites; a.k.a., translation start sites)
## TLS (Translation Stop sites)
## TTS (Transcription Termination Sites)
################################################################################
# Transcription Start Sites (TSS) — 5'UTR start  
TSS = subset(gtf, type == "five_prime_utr")

# Translation Start Sites (aka CDS starts)
CDS = subset(gtf, type == "start_codon")

# Translation Stop Sites (aka 3'UTR starts) 
# Transcription Termination Sites (aka 3'UTR ends)
UTR3 = subset(gtf, type == "three_prime_utr")
################################################################################

## Splice sites:
################################################################################
# Subset exons
exons = subset(gtf, type == "exon")

# Group by transcript and arrange by start position
exons$exon_number = as.numeric(exons$exon_number)
exons = exons %>% arrange(transcript_id, exon_number) %>% group_by(transcript_id)

# Calculate intron boundaries on sense strand
exons_sense = exons %>% filter(strand == '+')
introns_sense = exons_sense %>%
  mutate(next_start = lead(start), 
         next_end = lead(end)) %>%
  rowwise() %>%
  filter(!is.na(next_start)) %>%
  summarize(intron_start = end + 1,
            intron_end = next_start - 1,
            seqid = seqid[1],
            strand = strand[1],
            intron_number = exon_number) %>%
  ungroup()

# Calculate intron boundaries on anti-sense strand
exons_antisense = exons %>% filter(strand == '-')
introns_antisense = exons_antisense %>%
  mutate(next_start = lead(start), 
         next_end = lead(end)) %>%
  rowwise() %>%
  filter(!is.na(next_start)) %>%
  summarize(intron_start = next_end + 1,
            intron_end = start - 1,
            seqid = seqid[1],
            strand = strand[1],
            intron_number = exon_number) %>%
  ungroup()

INTRONs = rbind(introns_sense, introns_antisense)
INTRONs = INTRONs %>% arrange(transcript_id, intron_number)
colnames(INTRONs) = c('transcript_id', 'start', 'end', 'seqid', 'strand', 'intron_number')
################################################################################

## 4. Load Pre-filtered Peaks:
################################################################################
peaks_hnRNPC_WT = fread(paste0(INPUT_DIR, "peaks_hnRNPC_WT.txt"))
peaks_hnRNPC_Mut = fread(paste0(INPUT_DIR, "peaks_hnRNPC_Mut.txt"))
peaks_hnRNPC_WT_inRBM25_WT = fread(paste0(INPUT_DIR, "peaks_hnRNPC_WT_inRBM25_WT.txt"))
peaks_hnRNPC_WT_inRBM25_Mut = fread(paste0(INPUT_DIR, "peaks_hnRNPC_WT_inRBM25_Mut.txt"))
peaks_hnRNPC_WT_inKD = fread(paste0(INPUT_DIR, "peaks_hnRNPC_WT_inKD.txt"))
peaks_hnRNPC_Mut_inKD = fread(paste0(INPUT_DIR, "peaks_hnRNPC_Mut_inKD.txt"))
peaks_RBM25_WT = fread(paste0(INPUT_DIR, "peaks_RBM25_WT.txt"))
peaks_RBM25_Mut = fread(paste0(INPUT_DIR, "peaks_RBM25_Mut.txt"))
################################################################################

## 5. Peak Density Metagene:
################################################################################
search_space = 1000

## Transcription Start Site (TSS):
#########################################################################################################################
hnRNPC_WT_metaDensity = metaDensity(feature_df = TSS, peaks_df = peaks_hnRNPC_WT, feature_center = 'start', window_definition = search_space)
hnRNPC_Mut_metaDensity = metaDensity(feature_df = TSS, peaks_df = peaks_hnRNPC_Mut, feature_center = 'start', window_definition = search_space)
hnRNPC_WT_inRBM25_WT_metaDensity = metaDensity(feature_df = TSS, peaks_df = peaks_hnRNPC_WT_inRBM25_WT, feature_center = 'start', window_definition = search_space)
hnRNPC_WT_inRBM25_Mut_metaDensity = metaDensity(feature_df = TSS, peaks_df = peaks_hnRNPC_WT_inRBM25_Mut, feature_center = 'start', window_definition = search_space)
hnRNPC_WT_inKD_metaDensity = metaDensity(feature_df = TSS, peaks_df = peaks_hnRNPC_WT_inKD, feature_center = 'start', window_definition = search_space)
hnRNPC_Mut_inKD_metaDensity = metaDensity(feature_df = TSS, peaks_df = peaks_hnRNPC_Mut_inKD, feature_center = 'start', window_definition = search_space)
RBM25_WT_metaDensity = metaDensity(feature_df = TSS, peaks_df = peaks_RBM25_WT, feature_center = 'start', window_definition = search_space)
RBM25_Mut_metaDensity = metaDensity(feature_df = TSS, peaks_df = peaks_RBM25_Mut, feature_center = 'start', window_definition = search_space)

TSS_densities = data.frame(position = hnRNPC_WT_metaDensity$midpoint,
                           hnRNPC_WT = hnRNPC_WT_metaDensity$density,
                           hnRNPC_Mut = hnRNPC_Mut_metaDensity$density,
                           hnRNPC_WT_inRBM25_WT = hnRNPC_WT_inRBM25_WT_metaDensity$density,
                           hnRNPC_WT_inRBM25_Mut = hnRNPC_WT_inRBM25_Mut_metaDensity$density,
                           hnRNPC_WT_inKD = hnRNPC_WT_inKD_metaDensity$density,
                           hnRNPC_Mut_inKD = hnRNPC_Mut_inKD_metaDensity$density,
                           RBM25_WT = RBM25_WT_metaDensity$density,
                           RBM25_Mut = RBM25_Mut_metaDensity$density)
#########################################################################################################################

## Translation Start (Coding Start) Site (CDS):
#########################################################################################################################
hnRNPC_WT_metaDensity = metaDensity(feature_df = CDS, peaks_df = peaks_hnRNPC_WT, feature_center = 'start', window_definition = search_space)
hnRNPC_Mut_metaDensity = metaDensity(feature_df = CDS, peaks_df = peaks_hnRNPC_Mut, feature_center = 'start', window_definition = search_space)
hnRNPC_WT_inRBM25_WT_metaDensity = metaDensity(feature_df = CDS, peaks_df = peaks_hnRNPC_WT_inRBM25_WT, feature_center = 'start', window_definition = search_space)
hnRNPC_WT_inRBM25_Mut_metaDensity = metaDensity(feature_df = CDS, peaks_df = peaks_hnRNPC_WT_inRBM25_Mut, feature_center = 'start', window_definition = search_space)
hnRNPC_WT_inKD_metaDensity = metaDensity(feature_df = CDS, peaks_df = peaks_hnRNPC_WT_inKD, feature_center = 'start', window_definition = search_space)
hnRNPC_Mut_inKD_metaDensity = metaDensity(feature_df = CDS, peaks_df = peaks_hnRNPC_Mut_inKD, feature_center = 'start', window_definition = search_space)
RBM25_WT_metaDensity = metaDensity(feature_df = CDS, peaks_df = peaks_RBM25_WT, feature_center = 'start', window_definition = search_space)
RBM25_Mut_metaDensity = metaDensity(feature_df = CDS, peaks_df = peaks_RBM25_Mut, feature_center = 'start', window_definition = search_space)

CDS_densities = data.frame(position = hnRNPC_WT_metaDensity$midpoint,
                           hnRNPC_WT = hnRNPC_WT_metaDensity$density,
                           hnRNPC_Mut = hnRNPC_Mut_metaDensity$density,
                           hnRNPC_WT_inRBM25_WT = hnRNPC_WT_inRBM25_WT_metaDensity$density,
                           hnRNPC_WT_inRBM25_Mut = hnRNPC_WT_inRBM25_Mut_metaDensity$density,
                           hnRNPC_WT_inKD = hnRNPC_WT_inKD_metaDensity$density,
                           hnRNPC_Mut_inKD = hnRNPC_Mut_inKD_metaDensity$density,
                           RBM25_WT = RBM25_WT_metaDensity$density,
                           RBM25_Mut = RBM25_Mut_metaDensity$density)
#########################################################################################################################

## 5' Splice Site (SS5):
#########################################################################################################################
hnRNPC_WT_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = peaks_hnRNPC_WT, feature_center = 'start', window_definition = search_space)
hnRNPC_Mut_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = peaks_hnRNPC_Mut, feature_center = 'start', window_definition = search_space)
hnRNPC_WT_inRBM25_WT_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = peaks_hnRNPC_WT_inRBM25_WT, feature_center = 'start', window_definition = search_space)
hnRNPC_WT_inRBM25_Mut_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = peaks_hnRNPC_WT_inRBM25_Mut, feature_center = 'start', window_definition = search_space)
hnRNPC_WT_inKD_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = peaks_hnRNPC_WT_inKD, feature_center = 'start', window_definition = search_space)
hnRNPC_Mut_inKD_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = peaks_hnRNPC_Mut_inKD, feature_center = 'start', window_definition = search_space)
RBM25_WT_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = peaks_RBM25_WT, feature_center = 'start', window_definition = search_space)
RBM25_Mut_metaDensity = metaDensity(feature_df = INTRONs, peaks_df = peaks_RBM25_Mut, feature_center = 'start', window_definition = search_space)

SS5_densities = data.frame(position = hnRNPC_WT_metaDensity$midpoint,
                           hnRNPC_WT = hnRNPC_WT_metaDensity$density,
                           hnRNPC_Mut = hnRNPC_Mut_metaDensity$density,
                           hnRNPC_WT_inRBM25_WT = hnRNPC_WT_inRBM25_WT_metaDensity$density,
                           hnRNPC_WT_inRBM25_Mut = hnRNPC_WT_inRBM25_Mut_metaDensity$density,
                           hnRNPC_WT_inKD = hnRNPC_WT_inKD_metaDensity$density,
                           hnRNPC_Mut_inKD = hnRNPC_Mut_inKD_metaDensity$density,
                           RBM25_WT = RBM25_WT_metaDensity$density,
                           RBM25_Mut = RBM25_Mut_metaDensity$density)
#########################################################################################################################

## 3' Splice Site (SS3):
#########################################################################################################################
INTRONs_noLast = INTRONs %>% group_by(transcript_id) %>% filter(intron_number != max(intron_number)) %>% ungroup()

hnRNPC_WT_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = peaks_hnRNPC_WT, feature_center = 'end', window_definition = search_space)
hnRNPC_Mut_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = peaks_hnRNPC_Mut, feature_center = 'end', window_definition = search_space)
hnRNPC_WT_inRBM25_WT_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = peaks_hnRNPC_WT_inRBM25_WT, feature_center = 'end', window_definition = search_space)
hnRNPC_WT_inRBM25_Mut_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = peaks_hnRNPC_WT_inRBM25_Mut, feature_center = 'end', window_definition = search_space)
hnRNPC_WT_inKD_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = peaks_hnRNPC_WT_inKD, feature_center = 'end', window_definition = search_space)
hnRNPC_Mut_inKD_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = peaks_hnRNPC_Mut_inKD, feature_center = 'end', window_definition = search_space)
RBM25_WT_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = peaks_RBM25_WT, feature_center = 'end', window_definition = search_space)
RBM25_Mut_metaDensity = metaDensity(feature_df = INTRONs_noLast, peaks_df = peaks_RBM25_Mut, feature_center = 'end', window_definition = search_space)

SS3_densities = data.frame(position = hnRNPC_WT_metaDensity$midpoint,
                           hnRNPC_WT = hnRNPC_WT_metaDensity$density,
                           hnRNPC_Mut = hnRNPC_Mut_metaDensity$density,
                           hnRNPC_WT_inRBM25_WT = hnRNPC_WT_inRBM25_WT_metaDensity$density,
                           hnRNPC_WT_inRBM25_Mut = hnRNPC_WT_inRBM25_Mut_metaDensity$density,
                           hnRNPC_WT_inKD = hnRNPC_WT_inKD_metaDensity$density,
                           hnRNPC_Mut_inKD = hnRNPC_Mut_inKD_metaDensity$density,
                           RBM25_WT = RBM25_WT_metaDensity$density,
                           RBM25_Mut = RBM25_Mut_metaDensity$density)
#########################################################################################################################

## Translation Stop Site (TLS):
#########################################################################################################################
hnRNPC_WT_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_hnRNPC_WT, feature_center = 'start', window_definition = search_space)
hnRNPC_Mut_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_hnRNPC_Mut, feature_center = 'start', window_definition = search_space)
hnRNPC_WT_inRBM25_WT_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_hnRNPC_WT_inRBM25_WT, feature_center = 'start', window_definition = search_space)
hnRNPC_WT_inRBM25_Mut_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_hnRNPC_WT_inRBM25_Mut, feature_center = 'start', window_definition = search_space)
hnRNPC_WT_inKD_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_hnRNPC_WT_inKD, feature_center = 'start', window_definition = search_space)
hnRNPC_Mut_inKD_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_hnRNPC_Mut_inKD, feature_center = 'start', window_definition = search_space)
RBM25_WT_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_RBM25_WT, feature_center = 'start', window_definition = search_space)
RBM25_Mut_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_RBM25_Mut, feature_center = 'start', window_definition = search_space)

TLS_densities = data.frame(position = hnRNPC_WT_metaDensity$midpoint,
                           hnRNPC_WT = hnRNPC_WT_metaDensity$density,
                           hnRNPC_Mut = hnRNPC_Mut_metaDensity$density,
                           hnRNPC_WT_inRBM25_WT = hnRNPC_WT_inRBM25_WT_metaDensity$density,
                           hnRNPC_WT_inRBM25_Mut = hnRNPC_WT_inRBM25_Mut_metaDensity$density,
                           hnRNPC_WT_inKD = hnRNPC_WT_inKD_metaDensity$density,
                           hnRNPC_Mut_inKD = hnRNPC_Mut_inKD_metaDensity$density,
                           RBM25_WT = RBM25_WT_metaDensity$density,
                           RBM25_Mut = RBM25_Mut_metaDensity$density)
#########################################################################################################################

## Transcription Termination Site (TTS):
#########################################################################################################################
hnRNPC_WT_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_hnRNPC_WT, feature_center = 'end', window_definition = search_space)
hnRNPC_Mut_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_hnRNPC_Mut, feature_center = 'end', window_definition = search_space)
hnRNPC_WT_inRBM25_WT_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_hnRNPC_WT_inRBM25_WT, feature_center = 'end', window_definition = search_space)
hnRNPC_WT_inRBM25_Mut_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_hnRNPC_WT_inRBM25_Mut, feature_center = 'end', window_definition = search_space)
hnRNPC_WT_inKD_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_hnRNPC_WT_inKD, feature_center = 'end', window_definition = search_space)
hnRNPC_Mut_inKD_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_hnRNPC_Mut_inKD, feature_center = 'end', window_definition = search_space)
RBM25_WT_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_RBM25_WT, feature_center = 'end', window_definition = search_space)
RBM25_Mut_metaDensity = metaDensity(feature_df = UTR3, peaks_df = peaks_RBM25_Mut, feature_center = 'end', window_definition = search_space)

TTS_densities = data.frame(position = hnRNPC_WT_metaDensity$midpoint,
                           hnRNPC_WT = hnRNPC_WT_metaDensity$density,
                           hnRNPC_Mut = hnRNPC_Mut_metaDensity$density,
                           hnRNPC_WT_inRBM25_WT = hnRNPC_WT_inRBM25_WT_metaDensity$density,
                           hnRNPC_WT_inRBM25_Mut = hnRNPC_WT_inRBM25_Mut_metaDensity$density,
                           hnRNPC_WT_inKD = hnRNPC_WT_inKD_metaDensity$density,
                           hnRNPC_Mut_inKD = hnRNPC_Mut_inKD_metaDensity$density,
                           RBM25_WT = RBM25_WT_metaDensity$density,
                           RBM25_Mut = RBM25_Mut_metaDensity$density)
#########################################################################################################################

## 6. Normalize Peak Densities:
#########################################################################################################################
TSS_densities_norm = TSS_densities
CDS_densities_norm = CDS_densities
SS5_densities_norm = SS5_densities
SS3_densities_norm = SS3_densities
TLS_densities_norm = TLS_densities
TTS_densities_norm = TTS_densities

for (set in CLIP_List) {
  densitySum = sum(TSS_densities_norm[, set]) +
    sum(CDS_densities_norm[, set])+
    sum(SS5_densities_norm[, set])+
    sum(SS3_densities_norm[, set])+
    sum(TLS_densities_norm[, set])+
    sum(TTS_densities_norm[, set])
  
  TSS_densities_norm[, set] = TSS_densities_norm[, set] / densitySum
  CDS_densities_norm[, set] = CDS_densities_norm[, set] / densitySum
  SS5_densities_norm[, set] = SS5_densities_norm[, set] / densitySum
  SS3_densities_norm[, set] = SS3_densities_norm[, set] / densitySum
  TLS_densities_norm[, set] = TLS_densities_norm[, set] / densitySum
  TTS_densities_norm[, set] = TTS_densities_norm[, set] / densitySum
}
#########################################################################################################################

## 7. Peak Density Plots (Normalized):
#########################################################################################################################
smoothing = 0.5
y_lim = c(-0.001, 0.05)
x_lim = c(-500, 500)

## hnRNPC WT vs Mut:
hnRNPC_colors = c('dodgerblue3', 'firebrick3')
hnRNPC_list = c('hnRNPC_WT', 'hnRNPC_Mut')

feature_names = c("Transcription Start Sites", "Translation Start Sites", "5' Splice Site", 
                  "3' Splice Site", "Translation Stop Sites", "Transcription Termination Sites")
density_dfs = list(TSS_densities_norm, CDS_densities_norm, SS5_densities_norm, 
                   SS3_densities_norm, TLS_densities_norm, TTS_densities_norm)

for (i in 1:length(feature_names)) {
  print(plot_Density(density_data = density_dfs[[i]],
               columns_list = hnRNPC_list,
               custom_colors = hnRNPC_colors,
               densityType = 'feature_metagene',
               featureName = feature_names[i],
               yaxis_lims = y_lim,
               xaxis_lims = x_lim,
               smoothing = NULL,
               sampleName = 'hnRNPC WT vs Mut'))
}

## hnRNPC in RBM25 WT vs Mut:
hnRNPC_inRBM25_colors = c('dodgerblue3', 'firebrick3')
hnRNPC_inRBM25_list = c('hnRNPC_WT_inRBM25_WT', 'hnRNPC_WT_inRBM25_Mut')

for (i in 1:length(feature_names)) {
  print(plot_Density(density_data = density_dfs[[i]],
               columns_list = hnRNPC_inRBM25_list,
               custom_colors = hnRNPC_inRBM25_colors,
               densityType = 'feature_metagene',
               featureName = feature_names[i],
               yaxis_lims = y_lim,
               xaxis_lims = x_lim,
               smoothing = smoothing,
               sampleName = 'hnRNPC in RBM25OE WT vs Mut'))
}

## RBM25 WT vs Mut:
RBM25_colors = c('dodgerblue3', 'firebrick3')
RBM25_list = c('RBM25_WT', 'RBM25_Mut')

for (i in 1:length(feature_names)) {
  print(plot_Density(density_data = density_dfs[[i]],
               columns_list = RBM25_list,
               custom_colors = RBM25_colors,
               densityType = 'feature_metagene',
               featureName = feature_names[i],
               yaxis_lims = c(-0.01, 0.025),
               xaxis_lims = NULL,
               smoothing = 0.5,
               sampleName = 'RBM25 WT vs Mut'))
}

## hnRNPC WT vs Mut in KD:
hnRNPC_inKD_colors = c('dodgerblue3', 'firebrick3')
hnRNPC_inKD_list = c('hnRNPC_WT_inKD', 'hnRNPC_Mut_inKD')

for (i in 1:length(feature_names)) {
  print(plot_Density(density_data = density_dfs[[i]],
               columns_list = hnRNPC_inKD_list,
               custom_colors = hnRNPC_inKD_colors,
               densityType = 'feature_metagene',
               featureName = feature_names[i],
               yaxis_lims = y_lim,
               xaxis_lims = x_lim,
               smoothing = smoothing,
               sampleName = 'hnRNPC WT vs Mut in KD'))
}


## All hnRNPC samples together:
all_hnRNPC_colors = c('dodgerblue3', 'firebrick3', 'skyblue', 'salmon', 'mediumpurple3', 'orchid3')
all_hnRNPC_list = c('hnRNPC_WT', 'hnRNPC_Mut', 'hnRNPC_WT_inRBM25_WT', 'hnRNPC_WT_inRBM25_Mut', 'hnRNPC_WT_inKD', 'hnRNPC_Mut_inKD')

y_lim = c(-0.001, 0.015)
y_lim = NULL
x_lim = c(-1000, 1000)
x_lim = NULL
smoothing = 0.5

for (i in 1:length(feature_names)) {
  print(plot_Density(density_data = density_dfs[[i]],
               columns_list = all_hnRNPC_list,
               custom_colors = all_hnRNPC_colors,
               densityType = 'feature_metagene',
               featureName = feature_names[i],
               yaxis_lims = y_lim,
               xaxis_lims = x_lim,
               smoothing = smoothing,
               sampleName = 'All hnRNPC'))
}

## All samples:
all_colors = c('dodgerblue3', 'firebrick3', 'skyblue', 'salmon', 'mediumpurple3', 'orchid3', 'darkseagreen3', 'darkseagreen4')
all_list = CLIP_List

y_lim = c(-0.001, 0.01)
y_lim = NULL
x_lim = c(-100, 100)
x_lim = NULL
smoothing = 0.5

for (i in 1:length(feature_names)) {
  print(plot_Density(density_data = density_dfs[[i]],
               columns_list = all_list,
               custom_colors = all_colors,
               densityType = 'feature_metagene',
               featureName = feature_names[i],
               yaxis_lims = y_lim,
               xaxis_lims = x_lim,
               smoothing = smoothing,
               sampleName = 'All Samples'))
}
#########################################################################################################################

## 8. Peak Density Plots (Raw Counts):
#########################################################################################################################
y_lim_raw = NULL
smoothing_raw = 0.5

for (i in 1:length(feature_names)) {
  print(plot_Density(density_data = list(TSS_densities, CDS_densities, SS5_densities, 
                                         SS3_densities, TLS_densities, TTS_densities)[[i]],
               columns_list = all_list,
               custom_colors = all_colors,
               densityType = 'feature_metagene',
               featureName = feature_names[i],
               yaxis_lims = c(-5, 200),
               xaxis_lims = x_lim,
               smoothing = smoothing_raw,
               sampleName = 'All Samples (Raw)'))
}
#########################################################################################################################

## 9. Peak Counts Summary:
#########################################################################################################################
peakCounts = data.frame(locus = c('TSS', 'CDS', 'SS5', 'SS3', 'TLS', 'TTS'))
for (sample in CLIP_List) {
  peakCounts[[sample]] = 0
}

samples = CLIP_List
locus = c('TSS', 'CDS', 'SS5', 'SS3', 'TLS', 'TTS')

for (sample in samples) {
  for (loci in locus) {
    peakData = get(paste0('peaks_', sample))
    if (loci == 'SS5') {
      lociData = INTRONs
      peakCounts[which(peakCounts$locus == loci), sample] = countPeaks(lociData, peakData, 'start', 500)
    } else if (loci == 'SS3') {
      lociData = INTRONs_noLast
      peakCounts[which(peakCounts$locus == loci), sample] = countPeaks(lociData, peakData, 'end', 500)
    } else if (loci == 'TLS') {
      lociData = UTR3
      peakCounts[which(peakCounts$locus == loci), sample] = countPeaks(lociData, peakData, 'start', 500)
    } else if (loci == 'TTS') {
      lociData = UTR3
      peakCounts[which(peakCounts$locus == loci), sample] = countPeaks(lociData, peakData, 'end', 500)
    } else {
      lociData = get(loci)
      peakCounts[which(peakCounts$locus == loci), sample] = countPeaks(lociData, peakData, 'start', 500)
    }
  }
}

print(peakCounts)
#########################################################################################################################

## 10. Nucleotide Content Around Genomic Features:
#########################################################################################################################
## These are sequence-intrinsic properties (not sample-dependent).
## Compute U-content (T on DNA) and G-content around each feature.

nc_search_space = 1000
nc_window = 10

## TSS:
TSS_Ucontent = nucleotideContentAroundFeature(TSS, genome, "T", feature_center = 'start', nc_window, nc_search_space)
TSS_Gcontent = nucleotideContentAroundFeature(TSS, genome, "G", feature_center = 'start', nc_window, nc_search_space)

## CDS:
CDS_Ucontent = nucleotideContentAroundFeature(CDS, genome, "T", feature_center = 'start', nc_window, nc_search_space)
CDS_Gcontent = nucleotideContentAroundFeature(CDS, genome, "G", feature_center = 'start', nc_window, nc_search_space)

## SS5 (5' Splice Site — intron start):
SS5_Ucontent = nucleotideContentAroundFeature(INTRONs, genome, "T", feature_center = 'start', nc_window, nc_search_space)
SS5_Gcontent = nucleotideContentAroundFeature(INTRONs, genome, "G", feature_center = 'start', nc_window, nc_search_space)

## SS3 (3' Splice Site — intron end):
SS3_Ucontent = nucleotideContentAroundFeature(INTRONs_noLast, genome, "T", feature_center = 'end', nc_window, nc_search_space)
SS3_Gcontent = nucleotideContentAroundFeature(INTRONs_noLast, genome, "G", feature_center = 'end', nc_window, nc_search_space)

## TLS (Translation Stop — 3'UTR start):
TLS_Ucontent = nucleotideContentAroundFeature(UTR3, genome, "T", feature_center = 'start', nc_window, nc_search_space)
TLS_Gcontent = nucleotideContentAroundFeature(UTR3, genome, "G", feature_center = 'start', nc_window, nc_search_space)

## TTS (Transcription Termination — 3'UTR end):
TTS_Ucontent = nucleotideContentAroundFeature(UTR3, genome, "T", feature_center = 'end', nc_window, nc_search_space)
TTS_Gcontent = nucleotideContentAroundFeature(UTR3, genome, "G", feature_center = 'end', nc_window, nc_search_space)

## Assemble into combined data frames for plotting:
nc_feature_names = c('TSS', 'CDS', 'SS5', 'SS3', 'TLS', 'TTS')
nc_feature_labels = c("Transcription Start Sites", "Translation Start Sites", "5' Splice Site", 
                      "3' Splice Site", "Translation Stop Sites", "Transcription Termination Sites")

Ucontent_features = data.frame(position = TSS_Ucontent$midpoint,
                               TSS = TSS_Ucontent$fraction,
                               CDS = CDS_Ucontent$fraction,
                               SS5 = SS5_Ucontent$fraction,
                               SS3 = SS3_Ucontent$fraction,
                               TLS = TLS_Ucontent$fraction,
                               TTS = TTS_Ucontent$fraction)

Gcontent_features = data.frame(position = TSS_Gcontent$midpoint,
                               TSS = TSS_Gcontent$fraction,
                               CDS = CDS_Gcontent$fraction,
                               SS5 = SS5_Gcontent$fraction,
                               SS3 = SS3_Gcontent$fraction,
                               TLS = TLS_Gcontent$fraction,
                               TTS = TTS_Gcontent$fraction)

#########################################################################################################################

## 10-A. Plot Nucleotide Content Around Genomic Features:
#########################################################################################################################
## Plot U-content per feature:
for (i in 1:length(nc_feature_names)) {
  uc_df = data.frame(position = Ucontent_features$position,
                     U_content = Ucontent_features[[nc_feature_names[i]]])
  
  print(ggplot(uc_df, aes(x = position, y = U_content)) +
    geom_vline(xintercept = 0.0, color = "red", linetype = "dashed") +
    geom_smooth(span = 0.2, se = FALSE, color = "darkblue", linewidth = 1) +
    ylim(c(0.1, 0.5)) +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = 'bold')) +
    labs(title = paste0('U-content around ', nc_feature_labels[i]),
         x = paste0("Distance from ", nc_feature_labels[i], " (nucleotides)"),
         y = "U Fraction"))
}

## Plot G-content per feature:
for (i in 1:length(nc_feature_names)) {
  gc_df = data.frame(position = Gcontent_features$position,
                     G_content = Gcontent_features[[nc_feature_names[i]]])
  
  print(ggplot(gc_df, aes(x = position, y = G_content)) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(span = 0.2, se = FALSE, color = "darkgreen", linewidth = 1) +
    ylim(c(0.1, 0.5)) +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = 'bold')) +
    labs(title = paste0('G-content around ', nc_feature_labels[i]),
         x = paste0("Distance from ", nc_feature_labels[i], " (nucleotides)"),
         y = "G Fraction"))
}

## Combined U and G overlay per feature:
for (i in 1:length(nc_feature_names)) {
  combined_df = data.frame(position = Ucontent_features$position,
                           U = Ucontent_features[[nc_feature_names[i]]],
                           G = Gcontent_features[[nc_feature_names[i]]])
  combined_long = combined_df %>% pivot_longer(cols = c(U, G), names_to = "Nucleotide", values_to = "Fraction")
  
  print(ggplot(combined_long, aes(x = position, y = Fraction, color = Nucleotide)) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(span = 0.2, se = FALSE, linewidth = 1) +
    scale_color_manual(values = c("U" = "darkblue", "G" = "darkgreen")) +
    ylim(c(0.05, 0.5)) +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = 'bold')) +
    labs(title = paste0('U & G content around ', nc_feature_labels[i]),
         x = paste0("Distance from ", nc_feature_labels[i], " (nucleotides)"),
         y = "Nucleotide Fraction"))
}
#########################################################################################################################

## 11. Nucleotide Content Around Peak Centers:
#########################################################################################################################
## U-content around peak centers for each sample:
peak_nc_search_space = 1000
peak_nc_window = 10

Ucontent_peaks = data.frame(position = seq(from = -(peak_nc_search_space - peak_nc_window/2), 
                                            to = (peak_nc_search_space - peak_nc_window/2 - 1), 
                                            length.out = 2 * floor(peak_nc_search_space / peak_nc_window)))

for (sample in CLIP_List) {
  peakData = get(paste0('peaks_', sample))
  uc = nucleotideContentAroundPeaks(peakData, genome, "T", peak_nc_window, peak_nc_search_space)
  Ucontent_peaks[[sample]] = uc$fraction
  Ucontent_peaks$position = uc$midpoint
}

## G-content around peak centers for each sample:
Gcontent_peaks = data.frame(position = Ucontent_peaks$position)

for (sample in CLIP_List) {
  peakData = get(paste0('peaks_', sample))
  gc_result = nucleotideContentAroundPeaks(peakData, genome, "G", peak_nc_window, peak_nc_search_space)
  Gcontent_peaks[[sample]] = gc_result$fraction
}
#########################################################################################################################

## 11-A. Plot Nucleotide Content Around Peak Centers:
#########################################################################################################################
smoothing = 0.1

## Plot U-content around peak centers:
print(plot_Density(density_data = Ucontent_peaks,
             columns_list = CLIP_List,
             custom_colors = all_colors,
             densityType = 'nucleotide_content',
             featureName = "U-content around Peak Centers",
             smoothing = smoothing,
             yaxis_lims = c(0, 0.6),
             xaxis_lims = c(-500, 500),
             sampleName = 'All Samples'))

## Plot G-content around peak centers:
print(plot_Density(density_data = Gcontent_peaks,
                   columns_list = CLIP_List,
                   custom_colors = all_colors,
                   densityType = 'nucleotide_content',
                   featureName = "G-content around Peak Centers",
                   smoothing = smoothing,
                   yaxis_lims = c(0, 0.6),
                   xaxis_lims = c(-500, 500),
                   sampleName = 'All Samples'))



## hnRNPC WT vs Mut:
print(plot_Density(density_data = Ucontent_peaks,
             columns_list = hnRNPC_list,
             custom_colors = hnRNPC_colors,
             densityType = 'nucleotide_content',
             featureName = "U-content around Peak Centers",
             smoothing = smoothing,
             sampleName = 'hnRNPC WT vs Mut'))

## hnRNPC in RBM25 WT vs Mut:
print(plot_Density(density_data = Ucontent_peaks,
             columns_list = hnRNPC_inRBM25_list,
             custom_colors = hnRNPC_inRBM25_colors,
             densityType = 'nucleotide_content',
             featureName = "U-content around Peak Centers",
             smoothing = smoothing,
             sampleName = 'hnRNPC in RBM25OE WT vs Mut'))

## hnRNPC WT vs Mut in KD:
print(plot_Density(density_data = Ucontent_peaks,
             columns_list = hnRNPC_inKD_list,
             custom_colors = hnRNPC_inKD_colors,
             densityType = 'nucleotide_content',
             featureName = "U-content around Peak Centers",
             smoothing = smoothing,
             sampleName = 'hnRNPC WT vs Mut in KD'))

## RBM25 WT vs Mut:
print(plot_Density(density_data = Ucontent_peaks,
             columns_list = RBM25_list,
             custom_colors = RBM25_colors,
             densityType = 'nucleotide_content',
             featureName = "U-content around Peak Centers",
             smoothing = smoothing,
             sampleName = 'RBM25 WT vs Mut'))



## hnRNPC WT vs Mut:
print(plot_Density(density_data = Gcontent_peaks,
             columns_list = hnRNPC_list,
             custom_colors = hnRNPC_colors,
             densityType = 'nucleotide_content',
             featureName = "G-content around Peak Centers",
             smoothing = smoothing,
             sampleName = 'hnRNPC WT vs Mut'))

## hnRNPC in RBM25 WT vs Mut:
print(plot_Density(density_data = Gcontent_peaks,
             columns_list = hnRNPC_inRBM25_list,
             custom_colors = hnRNPC_inRBM25_colors,
             densityType = 'nucleotide_content',
             featureName = "G-content around Peak Centers",
             smoothing = smoothing,
             sampleName = 'hnRNPC in RBM25OE WT vs Mut'))

## hnRNPC WT vs Mut in KD:
print(plot_Density(density_data = Gcontent_peaks,
             columns_list = hnRNPC_inKD_list,
             custom_colors = hnRNPC_inKD_colors,
             densityType = 'nucleotide_content',
             featureName = "G-content around Peak Centers",
             smoothing = smoothing,
             sampleName = 'hnRNPC WT vs Mut in KD'))

## RBM25 WT vs Mut:
print(plot_Density(density_data = Gcontent_peaks,
             columns_list = RBM25_list,
             custom_colors = RBM25_colors,
             densityType = 'nucleotide_content',
             featureName = "G-content around Peak Centers",
             smoothing = smoothing,
             sampleName = 'RBM25 WT vs Mut'))
#########################################################################################################################


