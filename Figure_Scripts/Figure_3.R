################################################################################
## eCLIP Analysis Using RBPSpecificity Package:
## Written by Soon Yi
## Created: 2023-12-11
## Last Edited: 2026-02-21
## Figure 3
################################################################################

library(data.table)
library(tidyr)
library(RBPEqBind)

library(readr)
library(dplyr)

## Set up basic parameters:
################################################################################
baseDir = '~/Repos/BITS_Specificity/Dataset/Analysis/'
Model_RBPs = loadModel(paste0(baseDir, 'Model_RBP/Model_RBPs.csv'))

max_affinities = c("HH" = 100, "HL" = 100, "LH" = 10,  "LL" = 10)
min_affinities = c("HH" = 0.001, "HL" = 0.001, "LH" = 0.001, "LL" = 0.001)
rbp_Models = setModel(Model_RBPs, 
                     max_affinity = max_affinities, 
                     min_affinity = min_affinities)

viewModel(Model_RBPs, rbp = c("HH", "HL", "LH", "LL"), bins = 100, alpha = 0.4)


seq_Target = "UAGCGGUGCGAUUGGCCCGUGGACCGCGUUUUUGCACUCAUCGUUUCGCACUAAGUACAUAUAGUUGCGACAAAGCCGCUUAUGAGUUGGGGGUAUAUUC"
conc_Protein = c("HH" = 100, "HL" = 100, "LH" = 100, "LL" = 100)
conc_RNA = 10.0
################################################################################

## Run Simulation
################################################################################
simResults = simulateBinding(
  sequence = seq_Target,
  rbp_models = rbp_Models,
  protein_concs = conc_Protein,
  rna_conc = conc_RNA,
  k = 5
)

# Add transcript column for visualization
simResults$transcript = "CustomRNA"
################################################################################

## Visualization
################################################################################
xlim = NULL
xlim = c(20, 60)

plotBinding(simResults, 
            rbp = c("HH", "HL", "LH", "LL"), 
            metric = "density", 
            transcript = "CustomRNA", 
            xlim = xlim)

plotBinding(simResults, 
            rbp = c("HH", "HL", "LH", "LL"), 
            metric = "density", 
            transcript = "CustomRNA",
            window = 5,
            xlim = xlim)

plotHeatmap(simResults,
            metric = "occupancy_fc",
            transcript = "CustomRNA", 
            xlim = xlim, 
            xaxis_type = "both")
################################################################################

## Concentration Grid Sweep
################################################################################
grid_Protein = list(
  HH = c(0, 10, 25, 50, 100, 250, 500, 1000), 
  HL = c(0, 10, 25, 50, 100, 250, 500, 1000),
  LH = c(0, 10, 25, 50, 100, 250, 500, 1000),
  LL = c(0, 10, 25, 50, 100, 250, 500, 1000)
)
grid_RNA = c(10, 100, 500, 1000)

gridResults = simulateGrid(
  sequence = seq_Target,
  rbp_models = rbp_Models,
  protein_conc_grid = grid_Protein,
  rna_conc_grid = grid_RNA,
  k = 5,
  parallel = TRUE
)

# Increasing the non-specific protein concentration focuses the specific proteins
plotBinding(results = gridResults,
            rbp = c("HH", "HL", "LH", "LL"),
            rna_conc = 100,
            protein_conc = c(HH=100, HL=100, LH=0, LL=0),
            metric = "density",
            ylim = c(0, 0.05))

plotBinding(results = gridResults,
            rbp = c("HH", "HL", "LH", "LL"),
            rna_conc = 100,
            protein_conc = c(HH=100, HL=100, LH=10, LL=10),
            metric = "density",
            ylim = c(0, 0.05))

plotBinding(results = gridResults,
            rbp = c("HH", "HL", "LH", "LL"),
            rna_conc = 100,
            protein_conc = c(HH=100, HL=100, LH=50, LL=50),
            metric = "density",
            ylim = c(0, 0.05))

plotBinding(results = gridResults,
            rbp = c("HH", "HL", "LH", "LL"),
            rna_conc = 100,
            protein_conc = c(HH=100, HL=100, LH=100, LL=100),
            metric = "density",
            ylim = c(0, 0.05))




################################################################################

## Competition: HH vs LL on same sequence
################################################################################
# Scenario 1: HH is the stronger binder (HH: 0.001~100, LL: 0.001~10)
max_affinities_s1 = c("HH" = 100, "LL" = 10)
min_affinities_s1 = c("HH" = 0.001, "LL" = 0.001)
rbp_Models_s1 = setModel(Model_RBPs,
                         max_affinity = max_affinities_s1,
                         min_affinity = min_affinities_s1)
conc_LL_range = c(0, 10, 25, 50, 100, 250, 500, 1000)
xlim = c(20, 60)
xlim = NULL

for (conc_LL in conc_LL_range) {
  conc_Protein = c("HH" = 100, "LL" = conc_LL)
  simResults_s1 = simulateBinding(
    sequence = seq_Target,
    rbp_models = rbp_Models_s1,
    protein_concs = conc_Protein,
    rna_conc = 10.0,
    k = 5
  )
  simResults_s1$transcript = "CustomRNA"
  message("Scenario 1 | HH=100, LL=", conc_LL)
  print(plotBinding(simResults_s1,
                    rbp = c("HH", "LL"),
                    metric = "density",
                    transcript = "CustomRNA",
                    xlim = xlim,
                    ylim = c(0, 0.08)))
  print(plotHeatmap(simResults_s1,
                    metric = "density",
                    transcript = "CustomRNA",
                    xlim = xlim,
                    zlim = c(0, 0.08)))
}

# Scenario 2: LL is the stronger binder (HH: 0.001~10, LL: 0.001~100) — flipped
max_affinities_s2 = c("HH" = 10, "LL" = 100)
min_affinities_s2 = c("HH" = 0.001, "LL" = 0.001)
rbp_Models_s2 = setModel(Model_RBPs,
                         max_affinity = max_affinities_s2,
                         min_affinity = min_affinities_s2)

for (conc_LL in conc_LL_range) {
  conc_Protein = c("HH" = 100, "LL" = conc_LL)
  simResults_s2 = simulateBinding(
    sequence = seq_Target,
    rbp_models = rbp_Models_s2,
    protein_concs = conc_Protein,
    rna_conc = 10.0,
    k = 5
  )
  simResults_s2$transcript = "CustomRNA"
  message("Scenario 2 (flipped) | HH=100, LL=", conc_LL)
  print(plotBinding(simResults_s2,
                    rbp = c("HH", "LL"),
                    metric = "density",
                    transcript = "CustomRNA",
                    xlim = xlim,
                    ylim = c(0, 0.08)))
  print(plotHeatmap(simResults_s2,
                    metric = "density",
                    transcript = "CustomRNA",
                    xlim = xlim,
                    zlim = c(0, 0.08)))
}
################################################################################

