################################################################################
## CLIP Peak Filtering
## Author: Soon Yi with Antigravity
## Date: February 2026
################################################################################

# === Load Libraries ===
library(data.table)
library(tidyverse)
library(tibble)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggsignif)
library(eulerr)

## 1. Basic Setup:
################################################################################
BASE_DIR = "F:/Specificity/CLIP/Analysis/"
# INPUT_DIR = paste0(BASE_DIR, "output/")
# OUTPUT_DIR = paste0(BASE_DIR, "output/")

# BASE_DIR = "/Volumes/1TB_Data/Specificity/CLIP/Analysis/"
INPUT_DIR = paste0(BASE_DIR, "output/")
OUTPUT_DIR = paste0(BASE_DIR, "output/")

PEAK_MATRIX_FILE = paste0(INPUT_DIR, "meta_peakMatrix_normalized_annotated.txt")

################################################################################

## 2. Load and filter annotated peak matrix:
################################################################################
peaksMatrix = fread(PEAK_MATRIX_FILE)
# peaksMatrix = peaksMatrix %>% filter(grouped_annotation != "unannotated")

peaks_hnRNPC_WT = peaksMatrix %>% filter(BC_hnRNPC_WT >= 4)
peaks_hnRNPC_WT = peaks_hnRNPC_WT %>% filter(nTC_hnRNPC_WT >= quantile(peaks_hnRNPC_WT$nTC_hnRNPC_WT)[3])

peaks_hnRNPC_Mut = peaksMatrix %>% filter(BC_hnRNPC_Mut >= 4)
peaks_hnRNPC_Mut = peaks_hnRNPC_Mut %>% filter(nTC_hnRNPC_Mut >= quantile(peaks_hnRNPC_Mut$nTC_hnRNPC_Mut)[3])

peaks_hnRNPC_WT_inRBM25_WT = peaksMatrix %>% filter(BC_hnRNPC_WT_inRBM25_WT >= 3)
peaks_hnRNPC_WT_inRBM25_WT = peaks_hnRNPC_WT_inRBM25_WT %>% filter(nTC_hnRNPC_WT_inRBM25_WT >= 3*quantile(peaks_hnRNPC_WT_inRBM25_WT$nTC_hnRNPC_WT_inRBM25_WT)[4])

peaks_hnRNPC_WT_inRBM25_Mut = peaksMatrix %>% filter(BC_hnRNPC_WT_inRBM25_Mut >= 3)
peaks_hnRNPC_WT_inRBM25_Mut = peaks_hnRNPC_WT_inRBM25_Mut %>% filter(nTC_hnRNPC_WT_inRBM25_Mut >= 3*quantile(peaks_hnRNPC_WT_inRBM25_Mut$nTC_hnRNPC_WT_inRBM25_Mut)[4])

peaks_hnRNPC_WT_inKD = peaksMatrix %>% filter(BC_hnRNPC_WT_inKD >= 3)
peaks_hnRNPC_WT_inKD = peaks_hnRNPC_WT_inKD %>% filter(nTC_hnRNPC_WT_inKD >= quantile(peaks_hnRNPC_WT_inKD$nTC_hnRNPC_WT_inKD)[3])

peaks_hnRNPC_Mut_inKD = peaksMatrix %>% filter(BC_hnRNPC_Mut_inKD >= 3)
peaks_hnRNPC_Mut_inKD = peaks_hnRNPC_Mut_inKD %>% filter(nTC_hnRNPC_Mut_inKD >= quantile(peaks_hnRNPC_Mut_inKD$nTC_hnRNPC_Mut_inKD)[3])

peaks_RBM25_WT = peaksMatrix %>% filter(BC_RBM25_WT >= 3)
peaks_RBM25_WT = peaks_RBM25_WT %>% filter(nTC_RBM25_WT >= quantile(peaks_RBM25_WT$nTC_RBM25_WT)[3])

peaks_RBM25_Mut = peaksMatrix %>% filter(BC_RBM25_Mut >= 3)
peaks_RBM25_Mut = peaks_RBM25_Mut %>% filter(nTC_RBM25_Mut >= quantile(peaks_RBM25_Mut$nTC_RBM25_Mut)[3])


peaks_hnRNPC_WT = peaks_hnRNPC_WT %>% distinct()
peaks_hnRNPC_Mut = peaks_hnRNPC_Mut %>% distinct()
peaks_hnRNPC_WT_inRBM25_WT = peaks_hnRNPC_WT_inRBM25_WT %>% distinct()
peaks_hnRNPC_WT_inRBM25_Mut = peaks_hnRNPC_WT_inRBM25_Mut %>% distinct()
peaks_hnRNPC_WT_inKD = peaks_hnRNPC_WT_inKD %>% distinct()
peaks_hnRNPC_Mut_inKD = peaks_hnRNPC_Mut_inKD %>% distinct()
peaks_RBM25_WT = peaks_RBM25_WT %>% distinct()
peaks_RBM25_Mut = peaks_RBM25_Mut %>% distinct()

# fwrite(peaks_hnRNPC_WT, paste0(OUTPUT_DIR, "peaks_hnRNPC_WT", ".txt"), sep = "\t")
# fwrite(peaks_hnRNPC_Mut, paste0(OUTPUT_DIR, "peaks_hnRNPC_Mut", ".txt"), sep = "\t")
# fwrite(peaks_hnRNPC_WT_inRBM25_WT, paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inRBM25_WT", ".txt"), sep = "\t")
# fwrite(peaks_hnRNPC_WT_inRBM25_Mut, paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inRBM25_Mut", ".txt"), sep = "\t")
# fwrite(peaks_hnRNPC_WT_inKD, paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inKD", ".txt"), sep = "\t")
# fwrite(peaks_hnRNPC_Mut_inKD, paste0(OUTPUT_DIR, "peaks_hnRNPC_Mut_inKD", ".txt"), sep = "\t")
# fwrite(peaks_RBM25_WT, paste0(OUTPUT_DIR, "peaks_RBM25_WT", ".txt"), sep = "\t")
# fwrite(peaks_RBM25_Mut, paste0(OUTPUT_DIR, "peaks_RBM25_Mut", ".txt"), sep = "\t")

################################################################################

## 3. Load peak matrix:
################################################################################
peaks_hnRNPC_WT = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_WT.txt"))
peaks_hnRNPC_Mut = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_Mut.txt"))
peaks_hnRNPC_WT_inRBM25_WT = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inRBM25_WT.txt"))
peaks_hnRNPC_WT_inRBM25_Mut = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inRBM25_Mut.txt"))
peaks_hnRNPC_WT_inKD = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inKD.txt"))
peaks_hnRNPC_Mut_inKD = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_Mut_inKD.txt"))
peaks_RBM25_WT = fread(paste0(OUTPUT_DIR, "peaks_RBM25_WT.txt"))
peaks_RBM25_Mut = fread(paste0(OUTPUT_DIR, "peaks_RBM25_Mut.txt"))
################################################################################

## 4. Visualization - Stacked Bar Graphs:
################################################################################
## Custom Functions 
# Get Annotation counts:
countAnnotation = function(peak_matrix, annotation_column, new_column_name = NULL, annotation_to_skip = NULL, fraction = NULL) {
  temp = data.frame(table(peak_matrix[[annotation_column]]), row.names = 1)
  if(!is.null(new_column_name)) {
    colnames(temp) = new_column_name
  }
  
  if(!is.null(annotation_to_skip)) {
    temp = temp[rownames(temp) != annotation_to_skip, , drop = FALSE]
  }
  
  if(!is.null(fraction)) {
    temp = temp/sum(temp)
  }
  
  return(temp)
}

fillAnnotation = function(annotation_counts, annotation_list) {
  temp = data.frame(Sample = numeric(length(annotation_list)))
  rownames(temp) = annotation_list
  temp2 = merge(temp, annotation_counts, by = "row.names", all = TRUE)
  temp2[is.na(temp2)] = 0 
  temp2$Sample = NULL
  rownames(temp2) = temp2$Row.names
  temp2 = temp2[annotation_list, -1, drop = FALSE]
  return(temp2)
}

## Plot Stacked bar:
plotStackedBar = function(annotation_counts, sample_list, sample_label, title, y_lim = NULL, y_tick = NULL) {
  plot = ggplot(annotation_counts %>% filter(Source %in% sample_list), aes(fill = Annotation, y=Freq, x=Source)) + 
    geom_bar(position='stack', stat='identity') +
    scale_x_discrete(labels = sample_label) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Set3") +
    theme_bw() + 
    theme(axis.text = element_text(size=14), 
          axis.title = element_text(size=14, face = 'bold'), 
          legend.text = element_text(size=14))
  
  if(!is.null(y_lim)) {
    plot = plot + ylim(y_lim) 
  }
  
  if (!is.null(y_tick)) {
    plot = plot + scale_y_continuous(breaks = seq(0, y_lim[2], by=y_tick), limits=c(0, y_lim[2]))
  }
  
  return(plot)
}

## Global Variables:
CLIP_List = c('hnRNPC_WT', 'hnRNPC_Mut', 'hnRNPC_WT_inRBM25_WT', 'hnRNPC_WT_inRBM25_Mut', 'hnRNPC_WT_inKD', 'hnRNPC_Mut_inKD', 'RBM25_WT', 'RBM25_Mut')
All_Annotation_List = c("5'UTR", "CDS", "3'UTR", "intron", "snoRNA", 'ncRNA', "TE", "Other", "retained_intron", 'unannotated')

## Peak Counts Distribution Stacked Bar Graph
PeakDistribution_hnRNPC_WT = countAnnotation(peaks_hnRNPC_WT, 'grouped_annotation', 'hnRNPC_WT', 'unannotated')
PeakDistribution_hnRNPC_Mut = countAnnotation(peaks_hnRNPC_Mut, 'grouped_annotation', 'hnRNPC_Mut', 'unannotated')
PeakDistribution_RBM25_WT = countAnnotation(peaks_RBM25_WT, 'grouped_annotation', 'RBM25_WT', 'unannotated')
PeakDistribution_RBM25_Mut = countAnnotation(peaks_RBM25_Mut, 'grouped_annotation', 'RBM25_Mut', 'unannotated')
PeakDistribution_hnRNPC_WT_inRBM25_WT = countAnnotation(peaks_hnRNPC_WT_inRBM25_WT, 'grouped_annotation', 'hnRNPC_WT_inRBM25_WT', 'unannotated')
PeakDistribution_hnRNPC_WT_inRBM25_Mut = countAnnotation(peaks_hnRNPC_WT_inRBM25_Mut, 'grouped_annotation', 'hnRNPC_WT_inRBM25_Mut', 'unannotated')
PeakDistribution_hnRNPC_WT_inKD = countAnnotation(peaks_hnRNPC_WT_inKD, 'grouped_annotation', 'hnRNPC_WT_inKD', 'unannotated')
PeakDistribution_hnRNPC_Mut_inKD = countAnnotation(peaks_hnRNPC_Mut_inKD, 'grouped_annotation', 'hnRNPC_Mut_inKD', 'unannotated')

PeakDistribution_hnRNPC_Mut_inKD = fillAnnotation(PeakDistribution_hnRNPC_Mut_inKD, rownames(PeakDistribution_hnRNPC_WT))
PeakDistribution_RBM25_WT = fillAnnotation(PeakDistribution_RBM25_WT, rownames(PeakDistribution_hnRNPC_WT))
PeakDistribution_RBM25_Mut = fillAnnotation(PeakDistribution_RBM25_Mut, rownames(PeakDistribution_hnRNPC_WT))

PeakDistribution_combined = cbind(PeakDistribution_hnRNPC_WT, PeakDistribution_hnRNPC_Mut, 
                                  PeakDistribution_hnRNPC_WT_inRBM25_WT, PeakDistribution_hnRNPC_WT_inRBM25_Mut,
                                  PeakDistribution_hnRNPC_WT_inKD, PeakDistribution_hnRNPC_Mut_inKD,
                                  PeakDistribution_RBM25_WT, PeakDistribution_RBM25_Mut)
PeakDistribution_combined$Annotation = rownames(PeakDistribution_combined)

PeakDistribution_combined = PeakDistribution_combined %>% gather(key = "Source", value = "Freq", all_of(CLIP_List)) %>% dplyr::select(Source, Freq, Annotation)
PeakDistribution_combined$Source = factor(PeakDistribution_combined$Source, levels = CLIP_List)
PeakDistribution_combined$Annotation = factor(PeakDistribution_combined$Annotation, levels = All_Annotation_List)

peakDistribution_hnRNPC = plotStackedBar(PeakDistribution_combined, c('hnRNPC_WT', 'hnRNPC_Mut'), c('HNRNPC WT', 'HNRNPC Mut'), 'HNRNPC WT vs Mut Peak Distribution') 
peakDistribution_hnRNPC_inRBM25 = plotStackedBar(PeakDistribution_combined, c('hnRNPC_WT_inRBM25_WT', 'hnRNPC_WT_inRBM25_Mut'), c('HNRNPC in RBM25 WT', 'HNRNPC in RBM25 Mut'), 'HNRNPC in RBM25 WT vs Mut Peak Distribution') 
peakDistribution_hnRNPC_inKD = plotStackedBar(PeakDistribution_combined, c('hnRNPC_WT_inKD', 'hnRNPC_Mut_inKD'), c('HNRNPC WT in KD', 'HNRNPC Mut in KD'), 'HNRNPC WT vs Mut in KD Peak Distribution')
peakDistribution_RBM25 = plotStackedBar(PeakDistribution_combined, c('RBM25_WT', 'RBM25_Mut'), c('RBM25 WT', 'RBM25 Mut'), 'RBM25 WT vs Mut Peak Distribution')

## Add proportional plot for all samples:
PeakDistribution_prop = PeakDistribution_combined %>% 
  group_by(Source) %>% 
  mutate(Freq = Freq / sum(Freq)) %>% 
  ungroup()

peakDistribution_group = plotStackedBar(PeakDistribution_prop, CLIP_List, 
                                        c('HNRNPC WT', 'HNRNPC Mut', 'HNRNPC WT in RBM25 WT', 'HNRNPC WT in RBM25 Mut', 'HNRNPC WT in KD', 'HNRNPC Mut in KD', 'RBM25 WT', 'RBM25 Mut'), 
                                        'All Samples Peak Distribution (Proportion)')

# Display plots in RStudio
print(peakDistribution_hnRNPC)
print(peakDistribution_hnRNPC_inRBM25)
print(peakDistribution_hnRNPC_inKD)
print(peakDistribution_RBM25)
print(peakDistribution_group)

################################################################################

## 5. Peak Overlap Venn Diagram:
################################################################################
library(VennDiagram)
library(RColorBrewer)

# Standardize Venn plotting to RStudio Plots pane
plotOverlapVenn = function(list_of_sets, titles, filename) {
  # Create the Venn diagram
  venn = venn.diagram(
    x = list_of_sets,
    category.names = titles,
    filename = NULL,
    output = TRUE,
    
    # Aesthetics
    fill = brewer.pal(max(3, length(list_of_sets)), "Pastel1")[1:length(list_of_sets)], # nolint
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1, # Slightly smaller to fit long names
    main = gsub("_", " ", filename),
    main.cex = 1.5,
    margin = 0.1
  )
  
  # Remove the log file generated by VennDiagram
  if (file.exists(list.files(pattern = "VennDiagram.*\\.log"))) {
    file.remove(list.files(pattern = "VennDiagram.*\\.log"))
  }
  
  # Draw to the current device (RStudio Plots pane)
  grid.newpage()
  grid.draw(venn)
}

# 1. HNRNPC WT vs MUT
plotOverlapVenn(
  list(peaks_hnRNPC_WT$name, peaks_hnRNPC_Mut$name),
  c("HNRNPC WT", "HNRNPC Mut"),
  "hnRNPC WT vs Mut"
)

plotOverlapVenn(
  list(peaks_hnRNPC_WT$name, peaks_hnRNPC_WT_inRBM25_WT$name),
  c("HNRNPC WT", "HNRNPC WT in RBM25 WT"),
  "hnRNPC WT vs hnRNPC WT in RBM25 WT"
)

plotOverlapVenn(
  list(peaks_hnRNPC_WT$name, peaks_hnRNPC_WT_inRBM25_Mut$name),
  c("HNRNPC WT", "HNRNPC WT in RBM25 Mut"),
  "hnRNPC WT vs hnRNPC WT in RBM25 Mut"
)

plotOverlapVenn(
  list(peaks_hnRNPC_WT$name, peaks_hnRNPC_WT_inKD$name),
  c("HNRNPC WT", "HNRNPC WT in KD"),
  "hnRNPC WT vs hnRNPC WT in KD"
)

plotOverlapVenn(
  list(peaks_hnRNPC_WT$name, peaks_hnRNPC_Mut_inKD$name),
  c("HNRNPC WT", "HNRNPC Mut in KD"),
  "hnRNPC WT vs hnRNPC Mut in KD"
)

# 2. RBM25 WT vs MUT
plotOverlapVenn(
  list(peaks_RBM25_WT$name, peaks_RBM25_Mut$name),
  c("RBM25 WT", "RBM25 Mut"),
  "Overlap_RBM25_WT_vs_Mut"
)

# 3. HNRNPC WT vs HNRNPC in RBM25 WT vs HNRNPC in RBM25 Mut
plotOverlapVenn(
  list(peaks_hnRNPC_WT$name, peaks_hnRNPC_WT_inRBM25_WT$name, peaks_hnRNPC_WT_inRBM25_Mut$name),
  c("HNRNPC WT", "HNRNPC in RBM25 WT", "HNRNPC in RBM25 Mut"),
  "Overlap_HNRNPC_3way"
)

# 4. HNRNPC WT in RBM25 WT vs RBM25 WT
plotOverlapVenn(
  list(peaks_hnRNPC_WT_inRBM25_WT$name, peaks_RBM25_WT$name),
  c("HNRNPC in RBM25 WT", "RBM25 WT"),
  "Overlap_HNRNPC_inRBM25WT_vs_RBM25WT"
)

# 5. HNRNPC WT in RBM25 Mut vs RBM25 Mut
plotOverlapVenn(
  list(peaks_hnRNPC_WT_inRBM25_Mut$name, peaks_RBM25_Mut$name),
  c("HNRNPC in RBM25 Mut", "RBM25 Mut"),
  "Overlap_HNRNPC_inRBM25Mut_vs_RBM25Mut"
)


# 6. HNRNPC WT vs RBM25 WT
plotOverlapVenn(
  list(peaks_hnRNPC_WT$name, peaks_RBM25_WT$name),
  c("HNRNPC WT", "RBM25 WT"),
  "Overlap_HNRNPC_vs_RBM25WT"
)

# 7. HNRNPC WT vs RBM25 WT
plotOverlapVenn(
  list(peaks_hnRNPC_WT$name, peaks_RBM25_Mut$name),
  c("HNRNPC WT", "RBM25 Mut"),
  "Overlap_HNRNPC_vs_RBM25Mut"
)

################################################################################

################################################################################
peak_colors = c(
  "5'UTR"       = "#332f85",
  "CDS"         = "#89ccee",
  "3'UTR"       = "#44aa99",
  "intron"      = "#cc6677",
  "ncRNA"       = "#ddcc76",
  "unannotated" = "#dddddd"
)


# Join all dataframes
peaks_union_df = peaks_hnRNPC_WT %>%
  select(name, external_gene_name, grouped_annotation, nTC_hnRNPC_WT) %>%
  full_join(
    peaks_hnRNPC_WT_inRBM25_WT %>% select(name, nTC_hnRNPC_WT_inRBM25_WT),
    by = "name"
  ) %>%
  full_join(
    peaks_hnRNPC_WT_inRBM25_Mut %>% select(name, nTC_hnRNPC_WT_inRBM25_Mut),
    by = "name"
  )

peaks_union_df = peaks_union_df %>% filter(external_gene_name != "")

# Targeted NA replacement for numeric columns only
peaks_union_df = peaks_union_df %>%
  mutate(across(starts_with("nTC"), ~replace_na(., 0)))



#     Null    RBM25_WT_OE   RBM25_MUT_OE
# A    +           +             +
# B    +           +             -
# C    +           -             +
# D    +           -             -

peaks_A = peaks_union_df %>% filter(nTC_hnRNPC_WT > 0 & 
                                    nTC_hnRNPC_WT_inRBM25_WT > 0 & 
                                    nTC_hnRNPC_WT_inRBM25_Mut > 0)

peaks_A$L2FC_WT = log2(peaks_A$nTC_hnRNPC_WT_inRBM25_WT/peaks_A$nTC_hnRNPC_WT)
peaks_A$L2FC_Mut = log2(peaks_A$nTC_hnRNPC_WT_inRBM25_Mut/peaks_A$nTC_hnRNPC_WT)

peaks_B = peaks_union_df %>% filter(nTC_hnRNPC_WT > 0 & 
                                      nTC_hnRNPC_WT_inRBM25_WT > 0 & 
                                      nTC_hnRNPC_WT_inRBM25_Mut == 0)

peaks_B$L2FC_WT = log2(peaks_B$nTC_hnRNPC_WT_inRBM25_WT/peaks_B$nTC_hnRNPC_WT)

peaks_C = peaks_union_df %>% filter(nTC_hnRNPC_WT > 0 & 
                                      nTC_hnRNPC_WT_inRBM25_WT == 0 & 
                                      nTC_hnRNPC_WT_inRBM25_Mut > 0)

peaks_C$L2FC_Mut = log2(peaks_C$nTC_hnRNPC_WT_inRBM25_Mut/peaks_C$nTC_hnRNPC_WT)

peaks_D = peaks_union_df %>% filter(nTC_hnRNPC_WT > 0 & 
                                      nTC_hnRNPC_WT_inRBM25_WT == 0 & 
                                      nTC_hnRNPC_WT_inRBM25_Mut == 0)


peaks_A %>% filter(name %in% peaks_hnRNPC_WT$name)
peaks_A %>% filter(name %in% peaks_hnRNPC_WT_inRBM25_WT$name)
peaks_A %>% filter(name %in% peaks_hnRNPC_WT_inRBM25_Mut$name)

peaks_B %>% filter(name %in% peaks_hnRNPC_WT$name)
peaks_B %>% filter(name %in% peaks_hnRNPC_WT_inRBM25_WT$name)
peaks_B %>% filter(name %in% peaks_hnRNPC_WT_inRBM25_Mut$name)

peaks_C %>% filter(name %in% peaks_hnRNPC_WT$name)
peaks_C %>% filter(name %in% peaks_hnRNPC_WT_inRBM25_WT$name)
peaks_C %>% filter(name %in% peaks_hnRNPC_WT_inRBM25_Mut$name)

peaks_D %>% filter(name %in% peaks_hnRNPC_WT$name)
peaks_D %>% filter(name %in% peaks_hnRNPC_WT_inRBM25_WT$name)
peaks_D %>% filter(name %in% peaks_hnRNPC_WT_inRBM25_Mut$name)

## peaks_A plotting
# Ensure Annotation is a factor with the desired order for the grid
peaks_A$grouped_annotation = factor(
  peaks_A$grouped_annotation, 
  levels = c("5'UTR", "CDS", "3'UTR", "intron", "ncRNA", "unannotated")
)

# Figure 1: Overlaid Scatterplot
ggplot(peaks_A, aes(x = L2FC_WT, y = L2FC_Mut, color = grouped_annotation)) +
  geom_hline(yintercept = c(1, -1), color = "red", linetype = 'dotted') +
  geom_vline(xintercept = c(1, -1), color = "red", linetype = 'dotted') +
  geom_point(pch = 16, size = 3, alpha = 0.6) +
  scale_color_manual(values = peak_colors) +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 2.5)) +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 2.5)) +
  labs(title = "Overlaid Regions",
       x = "Log2(HNRNPC WT in RBM25 WT OE / HNRNPC WT)", 
       y = "Log2(HNRNPC WT in RBM25 MUT OE / HNRNPC WT)") +
  theme_bw() + 
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14, face = 'bold'))

# Figure 2: 2x3 Faceted Panels
ggplot(peaks_A, aes(x = L2FC_WT, y = L2FC_Mut, color = grouped_annotation)) +
  geom_hline(yintercept = c(1, -1), color = "red", linetype = 'dotted') +
  geom_vline(xintercept = c(1, -1), color = "red", linetype = 'dotted') +
  geom_point(pch = 16, size = 2, alpha = 0.6) +
  # Facet into 2 rows and 3 columns
  facet_wrap(~grouped_annotation, nrow = 2, ncol = 3) +
  scale_color_manual(values = peak_colors) +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 2.5)) +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 2.5)) +
  labs(x = "Log2(HNRNPC WT in RBM25 WT OE / HNRNPC WT)", 
       y = "Log2(HNRNPC WT in RBM25 MUT OE / HNRNPC WT)") +
  theme_bw() + 
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=14, face = 'bold'),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")



## genes that has shared peaks with greater than 2 fold increase/reduction 
peaks_A_focused = peaks_A %>% filter(L2FC_WT >= 1 & L2FC_Mut >= 1)
focused_gene_summary_table = peaks_A_focused %>%
  count(external_gene_name, name = "peak_count") %>%
  arrange(desc(peak_count))

peaks_A_masked = peaks_A %>% filter(L2FC_WT <= -1 & L2FC_Mut <= -1)
masked_gene_summary_table = peaks_A_masked %>%
  count(external_gene_name, name = "peak_count") %>%
  arrange(desc(peak_count))


## Master gene count table
# 1. Calculate peak counts per gene (dropping empty names)
counts_Null = peaks_hnRNPC_WT %>% 
  filter(external_gene_name != "") %>% 
  count(external_gene_name, name = "Null")

counts_WT_OE = peaks_hnRNPC_WT_inRBM25_WT %>% 
  filter(external_gene_name != "") %>% 
  count(external_gene_name, name = "RBM25_WT_OE")

counts_Mut_OE = peaks_hnRNPC_WT_inRBM25_Mut %>% 
  filter(external_gene_name != "") %>% 
  count(external_gene_name, name = "RBM25_Mut_OE")

# 2. Join into a master reference table
master_gene_counts = counts_Null %>%
  full_join(counts_WT_OE, by = "external_gene_name") %>%
  full_join(counts_Mut_OE, by = "external_gene_name") %>%
  # Replace NAs with 0 for genes not found in a specific condition
  mutate(across(c(Null, RBM25_WT_OE, RBM25_Mut_OE), ~replace_na(., 0))) %>%
  rename(genes = external_gene_name)

# 3. Create peak_focusing and peak_masking tables
focused_genes = unique(peaks_A_focused$external_gene_name)
masked_genes  = unique(peaks_A_masked$external_gene_name)

peak_focusing = master_gene_counts %>%
  filter(genes %in% focused_genes) %>%
  arrange(desc(RBM25_WT_OE))

peak_masking = master_gene_counts %>%
  filter(genes %in% masked_genes) %>%
  arrange(desc(Null))


write.csv(focused_gene_summary_table, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/GroupA_focused_peaks_count_per_gene.csv", row.names = F)
write.csv(masked_gene_summary_table, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/GroupA_masked_peaks_count_per_gene.csv", row.names = F)

write.csv(master_gene_counts, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/Union_peaks_total_count_per_gene.csv", row.names = F)
write.csv(peak_focusing, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/GroupA_peaks_total_count_per_gene_focused.csv", row.names = F)
write.csv(peak_masking, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/GroupA_peaks_total_count_per_gene_masked.csv", row.names = F)

################################################################################

################################################################################

peaks_B_focused = peaks_B %>% filter(L2FC_WT >= 1)
focused_gene_summary_table = peaks_B_focused %>%
  count(external_gene_name, name = "peak_count") %>%
  arrange(desc(peak_count))

peaks_B_masked = peaks_B %>% filter(L2FC_WT <= -1)
masked_gene_summary_table = peaks_B_masked %>%
  count(external_gene_name, name = "peak_count") %>%
  arrange(desc(peak_count))

focused_genes = unique(peaks_B_focused$external_gene_name)
masked_genes  = unique(peaks_B_masked$external_gene_name)

peak_focusing = master_gene_counts %>%
  filter(genes %in% focused_genes) %>%
  arrange(desc(RBM25_WT_OE))

peak_masking = master_gene_counts %>%
  filter(genes %in% masked_genes) %>%
  arrange(desc(Null))

write.csv(focused_gene_summary_table, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/GroupB_focused_peaks_count_per_gene.csv", row.names = F)
write.csv(masked_gene_summary_table, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/GroupB_masked_peaks_count_per_gene.csv", row.names = F)
write.csv(peak_focusing, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/GroupB_peaks_total_count_per_gene_focused.csv", row.names = F)
write.csv(peak_masking, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/GroupB_peaks_total_count_per_gene_masked.csv", row.names = F)

################################################################################

################################################################################

peaks_C_focused = peaks_C %>% filter(L2FC_Mut >= 1)
focused_gene_summary_table = peaks_C_focused %>%
  count(external_gene_name, name = "peak_count") %>%
  arrange(desc(peak_count))

peaks_C_masked = peaks_C %>% filter(L2FC_Mut <= -1)
masked_gene_summary_table = peaks_C_masked %>%
  count(external_gene_name, name = "peak_count") %>%
  arrange(desc(peak_count))

focused_genes = unique(peaks_C_focused$external_gene_name)
masked_genes  = unique(peaks_C_masked$external_gene_name)

peak_focusing = master_gene_counts %>%
  filter(genes %in% focused_genes) %>%
  arrange(desc(RBM25_WT_OE))

peak_masking = master_gene_counts %>%
  filter(genes %in% masked_genes) %>%
  arrange(desc(Null))

write.csv(focused_gene_summary_table, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/GroupC_focused_peaks_count_per_gene.csv", row.names = F)
write.csv(masked_gene_summary_table, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/GroupC_masked_peaks_count_per_gene.csv", row.names = F)
write.csv(peak_focusing, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/GroupC_peaks_total_count_per_gene_focused.csv", row.names = F)
write.csv(peak_masking, "~/Repos/BITS_Specificity/Dataset/Analysis/CompCLIP/GroupC_peaks_total_count_per_gene_masked.csv", row.names = F)

################################################################################

## Focusing/Masking Analysis
################################################################################
## 1. Classification & Data Prep
# Defined standard group palette to keep it consistent across all plots
group_palette = c("A" = "#999999", "B" = "#E69F00", "C" = "#56B4E9", "D" = "#009E73")

peaks_union_df = peaks_union_df %>%
  mutate(Group = case_when(
    nTC_hnRNPC_WT > 0 & nTC_hnRNPC_WT_inRBM25_WT > 0 & nTC_hnRNPC_WT_inRBM25_Mut > 0 ~ "A",
    nTC_hnRNPC_WT > 0 & nTC_hnRNPC_WT_inRBM25_WT > 0 & nTC_hnRNPC_WT_inRBM25_Mut == 0 ~ "B",
    nTC_hnRNPC_WT > 0 & nTC_hnRNPC_WT_inRBM25_WT == 0 & nTC_hnRNPC_WT_inRBM25_Mut > 0 ~ "C",
    nTC_hnRNPC_WT > 0 & nTC_hnRNPC_WT_inRBM25_WT == 0 & nTC_hnRNPC_WT_inRBM25_Mut == 0 ~ "D",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group)) %>%
  mutate(Group = factor(Group, levels = c("A", "B", "C", "D")))

## 2. Peak Composition (Single Bar - Absolute & Proportional)
group_summary = peaks_union_df %>% count(Group)

# Base plot for composition
p_comp = ggplot(group_summary, aes(x = "All Peaks", y = n, fill = Group)) +
  scale_fill_manual(values = group_palette) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"))

# Absolute
p_comp_abs = p_comp + 
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 5, color = "white", fontface = "bold") +
  labs(title = "Composition of Overlap Groups (A-D)", x = "", y = "Number of Peaks")

# Proportional
p_comp_prop = p_comp + 
  geom_bar(stat = "identity", position = "fill", width = 0.5) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Proportion of Overlap Groups (A-D)", x = "", y = "Percentage")

## 3. 4-Way Gene Venn (Euler)
# Streamlined: Create list directly from the grouped dataframe
gene_overlap_list = peaks_union_df %>%
  split(.$Group) %>%
  map(~unique(.$external_gene_name))

plot(euler(gene_overlap_list), 
     quantities = TRUE, fill = group_palette, alpha = 0.6,
     main = "Gene Overlap across Peak Groups A-D")

## 4. Peak Annotation Distributions (Groups A-D)
# Aggregated cleaning and annotation prep
PeakDist_ABCD = peaks_union_df %>%
  count(Group, grouped_annotation) %>%
  rename(Source = Group, Annotation = grouped_annotation, Freq = n) %>%
  mutate(Annotation = factor(Annotation, levels = All_Annotation_List))

# Absolute Plot
plot_ABCD_abs = plotStackedBar(PeakDist_ABCD, c("A", "B", "C", "D"), 
                               paste("Group", c("A", "B", "C", "D")), 
                               "Annotation Distribution (Counts)") +
  scale_fill_manual(values = peak_colors)

# Proportional Plot
plot_ABCD_prop = PeakDist_ABCD %>%
  group_by(Source) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  plotStackedBar(c("A", "B", "C", "D"), 
                 paste("Group", c("A", "B", "C", "D")), 
                 "Annotation Distribution (Proportion)") +
  scale_fill_manual(values = peak_colors)

print(p_comp_abs); print(p_comp_prop); print(plot_ABCD_abs); print(plot_ABCD_prop)
################################################################################





## 6. Peak Overlap TC comparison for hnRNPC in RBM25 WT/Mut OE:
################################################################################

# Identify overlapping peaks
peaks_overlap = intersect(peaks_hnRNPC_WT$name, 
                          intersect(peaks_hnRNPC_WT_inRBM25_WT$name, 
                                    peaks_hnRNPC_WT_inRBM25_Mut$name))

# Filter and align dataframes (Keep extra columns)
peaks_hnRNPC_WT_overlap = peaks_hnRNPC_WT %>% 
  filter(name %in% peaks_overlap) %>% 
  arrange(name)

peaks_hnRNPC_WT_inRBM25_WT_overlap = peaks_hnRNPC_WT_inRBM25_WT %>% 
  filter(name %in% peaks_overlap) %>% 
  arrange(name)

peaks_hnRNPC_WT_inRBM25_Mut_overlap = peaks_hnRNPC_WT_inRBM25_Mut %>% 
  filter(name %in% peaks_overlap) %>% 
  arrange(name)

# Construct peaks_overlap_df with gene names
# Assuming external_gene_name is consistent across dataframes for the same peak name
peaks_overlap_df = data.frame(
  peak = peaks_hnRNPC_WT_overlap$name,
  gene = peaks_hnRNPC_WT_overlap$external_gene_name,
  annotation = peaks_hnRNPC_WT_overlap$grouped_annotation,
  inRBM25WT = peaks_hnRNPC_WT_inRBM25_WT_overlap$nTC_hnRNPC_WT_inRBM25_WT / peaks_hnRNPC_WT_overlap$nTC_hnRNPC_WT,
  inRBM25Mut = peaks_hnRNPC_WT_inRBM25_Mut_overlap$nTC_hnRNPC_WT_inRBM25_Mut / peaks_hnRNPC_WT_overlap$nTC_hnRNPC_WT
)

 # Legend is redundant with facet titles
################################################################################



## 7. Unique peak comparison for hnRNPC in RBM25 WT/Mut OE:
################################################################################

peaks_unique = setdiff(peaks_hnRNPC_WT_inRBM25_WT$name, intersect(peaks_hnRNPC_WT_inRBM25_WT$name, peaks_hnRNPC_WT$name))
peaks_unique = setdiff(peaks_unique, intersect(peaks_hnRNPC_WT_inRBM25_Mut$name, peaks_unique))
peaks_hnRNPC_WT_inRBM25_WT_unique = peaks_hnRNPC_WT_inRBM25_WT %>% filter(name %in% peaks_unique) %>% arrange(name)

peaks_unique = setdiff(peaks_hnRNPC_WT_inRBM25_Mut$name, intersect(peaks_hnRNPC_WT_inRBM25_Mut$name, peaks_hnRNPC_WT$name))
peaks_unique = setdiff(peaks_unique, intersect(peaks_hnRNPC_WT_inRBM25_WT$name, peaks_unique))
peaks_hnRNPC_WT_inRBM25_Mut_unique = peaks_hnRNPC_WT_inRBM25_Mut %>% filter(name %in% peaks_unique) %>% arrange(name)

plot_data = data.frame(
  values = c(peaks_hnRNPC_WT_inRBM25_WT_unique$nTC_hnRNPC_WT_inRBM25_WT, 
             peaks_hnRNPC_WT_inRBM25_Mut_unique$nTC_hnRNPC_WT_inRBM25_Mut),
  group = c(rep("WT", length(peaks_hnRNPC_WT_inRBM25_WT_unique$nTC_hnRNPC_WT_inRBM25_WT)), 
            rep("Mut", length(peaks_hnRNPC_WT_inRBM25_Mut_unique$nTC_hnRNPC_WT_inRBM25_Mut)))
)

ggplot(plot_data, aes(x = group, y = values, fill = group)) +
  geom_boxplot(outlier.size = 1, alpha = 0.7, notch = T, outliers = F) +
  scale_fill_brewer(palette = "Set3") +
  # scale_y_continuous(transform = 'log10', limits = c(1, 1200)) +
  labs(x = "Condition", y = "nTC") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.position = "none") + 
  geom_signif(comparisons = list(c("WT", "Mut")),
              test = c('t.test'),
              map_signif_level = TRUE,
              textsize = 5, y_position = 1.7, tip_length = 0.001)




################################################################################

## 8. Peak Overlap TC comparison for RBM25:
################################################################################

peaks_overlap = intersect(peaks_RBM25_WT$name, peaks_RBM25_Mut$name)

peaks_RBM25_WT_overlap = peaks_RBM25_WT %>% filter(name %in% peaks_overlap) %>% arrange(name)
peaks_RBM25_Mut_overlap = peaks_RBM25_Mut %>% filter(name %in% peaks_overlap) %>% arrange(name)

TC_RBM25 = data.frame(WT = peaks_RBM25_WT_overlap$nTC_RBM25_WT,
                      Mut = peaks_RBM25_Mut_overlap$nTC_RBM25_Mut)


ggplot(TC_RBM25, aes(x = WT, y = Mut)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid") + 
  geom_hline(yintercept = 1, color = "red", linetype = 'dotted') +
  geom_vline(xintercept = 1, color = "red", linetype = 'dotted') +
  geom_point(pch = 16, size = 3, alpha = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  scale_x_continuous(transform = 'log10', limits = c(1, 1200)) +
  scale_y_continuous(transform = 'log10', limits = c(1, 1200)) +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))


TC_RBM25_long = TC_RBM25 %>%
  pivot_longer(cols = c(WT, Mut), names_to = "Condition", values_to = "nTC")

ggplot(TC_RBM25_long, aes(x = Condition, y = nTC, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, notch = T) + 
  geom_jitter(width = 0.2, alpha = 0.1)+
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(transform = 'log10', limits = c(1, 1200)) +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_signif(comparisons = list(c("WT", "Mut")),
              test = c('t.test'),
              map_signif_level = TRUE,
              textsize = 5, y_position = 2.9, tip_length = 0.05)

################################################################################

## 9. Peak Overlap TC comparison for RBM25 peaks that overlaps with hnRNPC peaks:
################################################################################

peaks_overlap = intersect(peaks_RBM25_WT$name, peaks_hnRNPC_WT$name)
peaks_RBM25_WT_overlap = peaks_RBM25_WT %>% filter(name %in% peaks_overlap) %>% arrange(name)

peaks_overlap = intersect(peaks_RBM25_Mut$name, peaks_hnRNPC_WT$name)
peaks_RBM25_Mut_overlap = peaks_RBM25_Mut %>% filter(name %in% peaks_overlap) %>% arrange(name)

# peaks_overlap = setdiff(intersect(peaks_RBM25_WT$name, peaks_hnRNPC_WT$name), peaks_RBM25_Mut$name)
# peaks_RBM25_WT_overlap = peaks_RBM25_WT %>% filter(name %in% peaks_overlap) %>% arrange(name)
# 
# peaks_overlap = setdiff(intersect(peaks_RBM25_Mut$name, peaks_hnRNPC_WT$name), peaks_RBM25_WT$name)
# peaks_RBM25_Mut_overlap = peaks_RBM25_Mut %>% filter(name %in% peaks_overlap) %>% arrange(name)

plot_data = data.frame(
  values = c(peaks_RBM25_WT_overlap$nTC_RBM25_WT, 
             peaks_RBM25_Mut_overlap$nTC_RBM25_Mut),
  group = c(rep("WT", length(peaks_RBM25_WT_overlap$nTC_RBM25_WT)), 
            rep("Mut", length(peaks_RBM25_Mut_overlap$nTC_RBM25_Mut)))
)

ggplot(plot_data, aes(x = group, y = values, fill = group)) +
  geom_boxplot(outlier.size = 1, alpha = 0.7, notch = T) +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(transform = 'log10', limits = c(1, 1200)) +
  labs(x = "Condition", y = "nTC") +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.position = "none") + 
  geom_signif(comparisons = list(c("WT", "Mut")),
              test = c('t.test'),
              map_signif_level = TRUE,
              textsize = 5, y_position = 2.9, tip_length = 0.001)

################################################################################