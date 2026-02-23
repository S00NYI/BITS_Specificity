## RBP Specificity Analysis Scripts
- Author: Soon Yi
- Last updated: 2026-02-22
---

## Repository Structure

### `PreProcessing_Scripts/`
Scripts for preprocessing raw RBNS data.

| Script | Description |
|--------|-------------|
| `0_1_RBNS_DataParsing.R` | RBNS data parsing based on the RBNS sample list |
| `0_2_RBNS_MetricCalculation.R` | IS and MS calculation across K-mers for RBNS data |

---

### `Figure_Scripts/`
Scripts for generating manuscript figures using the [RBPEqBind](https://github.com/S00NYI/RBPEqBind) and [RBPSpecificity](https://github.com/S00NYI/RBPSpecificity) packages.

| Script | Description |
|--------|-------------|
| `Figure_1.R` | Figure 1 panels |
| `Figure_2.R` | Figure 2 panels |
| `Figure_3.R` | Figure 3 panels — binding simulation with model RBPs (HH/HL/LH/LL), concentration sweeps, and competitive binding |
| `Figure_3_G.R` | Figure 3G — U2AF2 cofactor competition (HNRNPC, PTBP1) on PTBP2 transcript with iCLIP vs simulation comparison |
| `Figure_S1.R` | Supplemental Figure S1 panels |

---

### `Dataset/`
Input data and analysis outputs.

| Path | Description |
|------|-------------|
| `Table_S1.xlsx` | Supplemental Table S1 |
| `Table_S2.xlsx` | Supplemental Table S2 |
| `Table_S3.xlsx` | Supplemental Table S3 |
| `Table_S4.xlsx` | Supplemental Table S4 |
| `Analysis/sample_table_RBNS.csv` | RBNS sample metadata |
| `Analysis/sample_table_eCLIP.csv` | eCLIP sample metadata |
| `Analysis/Model_RBP/` | Model RBP affinity data |
| `Analysis/RBNS/` | Processed RBNS data |
| `Analysis/eCLIP_peak/` | eCLIP peak data |
| `Analysis/RBPEqBind_Simulation/` | U2AF2 simulation data (RNACompete model, FASTA, iCLIP bedgraphs, binding sites) |
| `AnalysisOutput/` | Analysis output tables and figures |

---

### `Deprecated/`
Legacy analysis scripts superseded by the current Figure/PreProcessing scripts.


---
For more information, please check our manuscript:
> Yi S, Singh SS, Ye X, Krishna R, Jankowsky E, Luna JM. (2025). 
> *Inherent Specificity and Mutational Sensitivity as Quantitative Metrics for RBP Binding.* 
> bioRxiv. https://www.biorxiv.org/content/10.1101/2025.03.28.646018v2

