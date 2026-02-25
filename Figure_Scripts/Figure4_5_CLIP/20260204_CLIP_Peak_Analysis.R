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

## 1. Basic Setup:
################################################################################
# BASE_DIR = "F:/Specificity/CLIP/Analysis/"
# INPUT_DIR = paste0(BASE_DIR, "output/")
# OUTPUT_DIR = paste0(BASE_DIR, "output/")

BASE_DIR = "/Volumes/1TB_Data/Specificity/CLIP/Analysis/"
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

fwrite(peaks_hnRNPC_WT, paste0(OUTPUT_DIR, "peaks_hnRNPC_WT", ".txt"), sep = "\t")
fwrite(peaks_hnRNPC_Mut, paste0(OUTPUT_DIR, "peaks_hnRNPC_Mut", ".txt"), sep = "\t")
fwrite(peaks_hnRNPC_WT_inRBM25_WT, paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inRBM25_WT", ".txt"), sep = "\t")
fwrite(peaks_hnRNPC_WT_inRBM25_Mut, paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inRBM25_Mut", ".txt"), sep = "\t")
fwrite(peaks_hnRNPC_WT_inKD, paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inKD", ".txt"), sep = "\t")
fwrite(peaks_hnRNPC_Mut_inKD, paste0(OUTPUT_DIR, "peaks_hnRNPC_Mut_inKD", ".txt"), sep = "\t")
fwrite(peaks_RBM25_WT, paste0(OUTPUT_DIR, "peaks_RBM25_WT", ".txt"), sep = "\t")
fwrite(peaks_RBM25_Mut, paste0(OUTPUT_DIR, "peaks_RBM25_Mut", ".txt"), sep = "\t")

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

## 6. Peak Overlap TC comparison for hnRNPC in RBM25 WT/Mut OE:
################################################################################

peaks_overlap = intersect(peaks_hnRNPC_WT$name, intersect(peaks_hnRNPC_WT_inRBM25_WT$name, peaks_hnRNPC_WT_inRBM25_Mut$name))

peaks_hnRNPC_WT_overlap = peaks_hnRNPC_WT %>% filter(name %in% peaks_overlap) %>% arrange(name)
peaks_hnRNPC_WT_inRBM25_WT_overlap = peaks_hnRNPC_WT_inRBM25_WT %>% filter(name %in% peaks_overlap) %>% arrange(name)
peaks_hnRNPC_WT_inRBM25_Mut_overlap = peaks_hnRNPC_WT_inRBM25_Mut %>% filter(name %in% peaks_overlap) %>% arrange(name)

L2FC_hnRNPC_inRBM25 = data.frame(WT = log2(peaks_hnRNPC_WT_inRBM25_WT_overlap$nTC_hnRNPC_WT_inRBM25_WT / peaks_hnRNPC_WT_overlap$nTC_hnRNPC_WT),
                                 Mut = log2(peaks_hnRNPC_WT_inRBM25_Mut_overlap$nTC_hnRNPC_WT_inRBM25_Mut / peaks_hnRNPC_WT_overlap$nTC_hnRNPC_WT))

ggplot(L2FC_hnRNPC_inRBM25, aes(x = WT, y = Mut)) +
  geom_hline(yintercept = 0, color = "red", linetype = 'dotted') +
  geom_vline(xintercept = 0, color = "red", linetype = 'dotted') +
  geom_point(pch = 16, size = 3, alpha = 0.1) +
  scale_fill_brewer(palette = "Set3") +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 2.5)) +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 2.5)) +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14))

L2FC_hnRNPC_inRBM25_long = L2FC_hnRNPC_inRBM25 %>%
  pivot_longer(cols = c(WT, Mut), names_to = "Condition", values_to = "nTC")

ggplot(L2FC_hnRNPC_inRBM25_long, aes(x = Condition, y = nTC, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, notch = T) + 
  # geom_jitter(width = 0.2, alpha = 0.1)+
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2.5, 2.5, by = 0.5)) +
  theme_bw() + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=14, face = 'bold'), 
        legend.text = element_text(size=14)) + 
  geom_signif(comparisons = list(c("WT", "Mut")),
              test = c('t.test'),
              map_signif_level = TRUE,
              textsize = 5, y_position = 2, tip_length = 0.05)



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