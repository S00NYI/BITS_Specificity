# REVISION_2026-09

Analyses for the Genome Biology revision of the IS/VS manuscript. Scripts are in `scripts/`,
data outputs at the top level, figure PDFs in `figures/`. Files are grouped by the **four items
included in the reviewer response** (prefixes `1_`–`4_`); superseded / exploratory work is in
`_deprecated/`.

Run any script from the repo root, e.g.:

```bash
cd /Users/soonyi/Repos/BITS_Specificity && Rscript REVISION_2026-09/scripts/2_mouse_CS_CVS.R
```

## CS/CVS procedure used for every number here

The manuscript CS/CVS come from `Figure_Scripts/Figure_2.R` (Fig. 2 set) and a 20260316 run
read by `Figure_Scripts/Figure_6.R` (Fig. 6 set). The package is **RBPSpecificity 0.99.0**
as installed: built 2026-02-21 from commit 503e73f, function bodies identical to 8b121b3
(HEAD on 2026-03-17). Do not update the package; HEAD 0.99.5 uses a different algorithm.

- Peaks: ENCODE GRCh38 IDR-thresholded narrowPeak (reps 1+2); keep chr1-22,X,Y,M; no score filter.
- BED start read as 1-based (start offset by 1 nt). 25-nt 5' extension (strand-aware), none at 3'.
- Counts: all 1024 5-mers over all peak sequences.
- Background: 100 iterations, each peak shifted +/-500-1000 nt (drop shifts overlapping any peak),
  no scrambling, averaged.
- Enrichment: peak minus mean background; min-max to [1, e]; natural log -> [0, 1].
- CS = Score(top 5-mer) / median(Score). CVS = mean of (Score_top - Score_variant) over the 15
  single-nt variants (zeros excluded).
- Merge: average the K562 and HepG2 Score vectors, min-max to [0, 1], then compute CS/CVS
  (both sets; in the Fig. 2 set, XRCC6 uses K562 only).

---

## 1. Updated Table S5 (per-cell-line + merged CS/CVS) and Table S1 fixes

- **`1_published_CS_provenance.R`** — confirms Table S5 CS/CVS equal the 2026-03-16 RBPSpecificity
  0.99.0 run (168/168 RBPs exact match); the pre-package and 0.99.5 runs do not match.
- **`1_S1_cellline_audit.R`** -> `1_S1_cellline_audit.csv` — per-RBP cell-line audit of Table S1.
  Two accession fixes: **GRSF1** (ENCFF929AWR is a HepG2 experiment, listed under K562 -> move to
  HepG2); **MBNL1** (analyzed in K562 as ENCFF603WDI, missing from Table S1 -> add to K562).
- **`1_S5_percell_CS_CVS.R`** -> `1_S5_percell_CS_CVS.csv` — per-cell-line and merged CS/CVS for all
  168 eCLIP RBPs (76 merged / 63 K562-only / 29 HepG2-only), validated against
  `eCLIP_5mer_IS_VS_summary.csv` (0 mismatches).
- **`1_S5_paste_columns.R`** -> `1_S5_paste_columns.csv` — the four new columns (CS_K562, CVS_K562,
  CS_HepG2, CVS_HepG2) in Table S5 row order, ready to paste. *(reads `1_S5_percell_CS_CVS.csv`.)*

## 2. Mouse eIF4G2 HITS-CLIP

- **`2_mouse_CS_CVS.R`** -> `2_mouse_CS_CVS.csv` — provisional CS = 11.6, CVS = 0.64, top 5-mer
  **GAGGA** (mm39 BSgenome injected into RBPSpecificity 0.99.0 internals; no package edit).
- **`2_mouse_5mer_enrichment.R`** -> `2_mouse_5mer_enrichment.csv` — the full 1024-mer enrichment
  vector (input to the affinity-distribution figure).
- **`2_mouse_region_distribution.R`** -> `2_mouse_region_distribution.csv` — peak- and tag-level
  region distribution, mouse HITS-CLIP + human EIF4G2 K562 eCLIP.
- **`2_mouse_affinity_distribution.R`** -> `figures/2_mouse_affinity_distribution.pdf`.
- **`2_mouse_peak_distribution.R`** -> `figures/2_mouse_peak_distribution.pdf` (+ `2_mouse_tag_distribution.pdf`).
- `2_mouse_SRA_manifest.csv` — GEO/SRA provenance for the mouse data (GSE213082).

## 3. Fig 4G global analysis (PTBP1 co-target vs non-target)

- **`3_fig4G_pool_simulation.R`** -> `3_fig4G_pool_per_transcript.csv`, `3_fig4G_pool_per_window.csv`
  — RBPEqBind simulation across the Sutandy in-vitro-iCLIP pool (8 transcripts with coverage);
  per-window correlation of simulated U2AF2 occupancy vs measured iCLIP, three conditions.
- **`3_fig4G_ptbp1_eclip_targets.R`** -> `3_fig4G_ptbp1_eclip_targets.csv` — classifies transcripts
  as PTBP1 co-target / non-target using ENCODE PTBP1 eCLIP (K562+HepG2 IDR peaks, gene-level overlap).
- **`3_fig4G_cotarget_violin.R`** -> `3_fig4G_cotarget_zerofill.csv`,
  `figures/3_fig4G_cotarget_violin.pdf` — the reported figure. Zero-fill (uncovered positions = 0),
  200-nt windows on measured U2AF2 sites. Co-target 0.37 -> 0.42 (p = 7e-17), non-target
  0.40 -> 0.39 (p = 0.58); paired Wilcoxon. *(reads `3_fig4G_cotarget_classification.csv`.)*
- `3_fig4G_cotarget_classification.csv` — PTBP1 and HNRNPC eCLIP peak counts + labels per transcript
  (input to the violin).

## 4. Bootstrap / jackknife confidence intervals (Fig 1D, 2C, 2D, 6A, S6A)

- **`4_bootstrap_cellline_reproducibility.R`** -> `4_bootstrap_cellline_reproducibility.csv` —
  per-RBP CS/CVS, K562 vs HepG2 (supplies the Fig 6A/S6A correlations).
- **`4_bootstrap_fig2CD.R`** -> `4_bootstrap_fig2CD.csv` — bootstrap + jackknife for Fig 2C/2D.
- **`4_bootstrap_all_CI.R`** -> `4_bootstrap_all_CI.csv` — consolidated 95% percentile bootstrap CIs
  + jackknife for all five panels (exact enumeration for the n=9 groups, 10,000 Monte-Carlo else,
  fixed seed). *(reads `4_bootstrap_cellline_reproducibility.csv`.)*
- **`4_bootstrap_forest.R`** -> `figures/4_bootstrap_forest.pdf` — forest plot of the CIs.
  *(reads `4_bootstrap_all_CI.csv`.)*

---

## figures/
`2_mouse_affinity_distribution.pdf`, `2_mouse_peak_distribution.pdf`, `2_mouse_tag_distribution.pdf`,
`3_fig4G_cotarget_violin.pdf`, `4_bootstrap_forest.pdf`.

## _deprecated/
Exploratory or superseded work, and analyses for reviewer points not among the four figures:
human eCLIP motif-GC, off-diagonal leave-one-out, supp-table numbering audit, Mut->DSwap audit,
eIF4G2 cross-species scatter, the Fig 4G baseline reproduction, window-size sweeps, single-transcript
windowed analysis, sim/motif proxy, ROI per-site, and the combined + inner-join violins. Note:
`3_fig4G_cotarget_classification.csv` (kept, read by the zero-fill script) was originally produced by
`_deprecated/scripts/TaskN_pool_cotarget_violin.R`.
