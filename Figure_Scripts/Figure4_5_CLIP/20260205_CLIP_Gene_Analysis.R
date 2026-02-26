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
library(purrr)
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

## 2. Load peak matrix:
################################################################################
peaks_hnRNPC_WT = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_WT.txt"))
peaks_RBM25_WT = fread(paste0(OUTPUT_DIR, "peaks_RBM25_WT.txt"))
peaks_RBM25_Mut = fread(paste0(OUTPUT_DIR, "peaks_RBM25_Mut.txt"))
peaks_hnRNPC_WT_inRBM25_WT = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inRBM25_WT.txt"))
peaks_hnRNPC_WT_inRBM25_Mut = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inRBM25_Mut.txt"))

# peaks_hnRNPC_Mut = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_Mut.txt"))
# peaks_hnRNPC_WT_inKD = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_WT_inKD.txt"))
# peaks_hnRNPC_Mut_inKD = fread(paste0(OUTPUT_DIR, "peaks_hnRNPC_Mut_inKD.txt"))
################################################################################

## 3. Unique target genes per sample:
################################################################################
genes_hnRNPC_WT = unique(peaks_hnRNPC_WT$external_gene_name); genes_hnRNPC_WT = genes_hnRNPC_WT[genes_hnRNPC_WT != ""]
genes_RBM25_WT = unique(peaks_RBM25_WT$external_gene_name); genes_RBM25_WT = genes_RBM25_WT[genes_RBM25_WT != ""]
genes_RBM25_Mut = unique(peaks_RBM25_Mut$external_gene_name); genes_RBM25_Mut = genes_RBM25_Mut[genes_RBM25_Mut != ""]
genes_hnRNPC_WT_inRBM25_WT = unique(peaks_hnRNPC_WT_inRBM25_WT$external_gene_name); genes_hnRNPC_WT_inRBM25_WT = genes_hnRNPC_WT_inRBM25_WT[genes_hnRNPC_WT_inRBM25_WT != ""]
genes_hnRNPC_WT_inRBM25_Mut = unique(peaks_hnRNPC_WT_inRBM25_Mut$external_gene_name); genes_hnRNPC_WT_inRBM25_Mut = genes_hnRNPC_WT_inRBM25_Mut[genes_hnRNPC_WT_inRBM25_Mut != ""]

# genes_hnRNPC_Mut = unique(peaks_hnRNPC_Mut$external_gene_name); genes_hnRNPC_Mut = genes_hnRNPC_Mut[genes_hnRNPC_Mut != ""]
# genes_hnRNPC_WT_inKD = unique(peaks_hnRNPC_WT_inKD$external_gene_name); genes_hnRNPC_WT_inKD = genes_hnRNPC_WT_inKD[genes_hnRNPC_WT_inKD != ""]
# genes_hnRNPC_Mut_inKD = unique(peaks_hnRNPC_Mut_inKD$external_gene_name); genes_hnRNPC_Mut_inKD = genes_hnRNPC_Mut_inKD[genes_hnRNPC_Mut_inKD != ""]
################################################################################

## 4. Target genes overlap:
################################################################################
genes_hnRNPC_WT_overlap = intersect(genes_hnRNPC_WT, intersect(genes_hnRNPC_WT_inRBM25_WT, genes_hnRNPC_WT_inRBM25_Mut))

genes_RBM25_WT_hnRNPC_WT_target = intersect(genes_RBM25_WT, genes_hnRNPC_WT_overlap)
genes_RBM25_Mut_hnRNPC_WT_target = intersect(genes_RBM25_Mut, genes_hnRNPC_WT_overlap)

genes_RBM25_WT_target_only = setdiff(genes_RBM25_WT_hnRNPC_WT_target, genes_RBM25_Mut_hnRNPC_WT_target)
genes_RBM25_Mut_target_only = setdiff(genes_RBM25_Mut_hnRNPC_WT_target, genes_RBM25_WT_hnRNPC_WT_target)
genes_RBM25_target_overlap = intersect(genes_RBM25_WT_hnRNPC_WT_target, genes_RBM25_Mut_hnRNPC_WT_target)


fwrite(as.data.frame(genes_RBM25_WT_target_only), paste0(OUTPUT_DIR, "genes_RBM25_WT_target_only", ".txt"), sep = "\t")
fwrite(as.data.frame(genes_RBM25_Mut_target_only), paste0(OUTPUT_DIR, "genes_RBM25_Mut_target_only", ".txt"), sep = "\t")
fwrite(as.data.frame(genes_RBM25_target_overlap), paste0(OUTPUT_DIR, "genes_RBM25_target_overlap", ".txt"), sep = "\t")


################################################################################
